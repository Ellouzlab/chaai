#pragma once

#include "db/types.h"
#include "aai/seeds.h"
#include "db/io.h"
#include "db/seedcache.h"
#include "aai/merge.h"
#include "aai/chain.h"
#include "aai/pairwise.h"
#include "aai/align.h"
#include "aai/pair.h"
#include "run/progress.h"
#include "run/interrupt.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <atomic>
#include <stdexcept>
#include <mutex>
#include <omp.h>

static int run_dist(
    const std::vector<std::string> &query_paths,
    const std::vector<std::string> &ref_paths,
    const std::string &output_path,
    int threads, bool is_self,
    const std::string &seed_shape, int seed_scaled,
    int max_seed_occ,
    int chain_band, int chain_max_qgap, int min_chain_seeds,
    float min_aai, int min_aligned_aa,
    int min_chain_aligned_aa,
    int min_matched_orfs,
    float min_orf_cov,
    bool aggressive_filter,
    int xdrop,
    double max_evalue,
    size_t max_pair_matches,
    int gene_scaled = 1,
    const std::string &temp_dir = "",
    bool no_triangle = false)
{
    InterruptScope guard;
    set_temp_dir(temp_dir);
    const PatSet SHAPE = parse_seedspec(seed_shape.c_str());
    if (SHAPE.n == 0) {
        throw std::invalid_argument(
            "seed_shape must be \"w:off,off,...\" (e.g. \"4:0,3,4,5\"); got \"" +
            seed_shape + "\"");
    }
    int total_pairs = 0;

    const int nq_total = static_cast<int>(query_paths.size());
    const int nr_total = static_cast<int>(ref_paths.size());

    DistParams p{max_seed_occ,
                 chain_band, chain_max_qgap, min_chain_seeds,
                 min_aai, min_aligned_aa,
                 min_chain_aligned_aa,
                 patset_weight(SHAPE), patset_span(SHAPE), min_matched_orfs,
                 min_orf_cov, aggressive_filter, xdrop,
                 max_evalue, max_pair_matches};
    if (gene_scaled < 1) throw std::invalid_argument("gene_scaled must be >= 1");
    p.gene_scaled = static_cast<uint32_t>(gene_scaled);

    std::cerr << "[dist] Seed params: shape=" << seed_shape
              << " scaled=" << seed_scaled
              << " xdrop=" << xdrop;
    if (gene_scaled > 1) std::cerr << " gene_scaled=" << gene_scaled;
    std::cerr << " [SW-align]\n";

    std::cerr << "[dist] " << nq_total << " queries x " << nr_total << " refs\n";

    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Cannot open: " + output_path);
    out << DIST_HEADER;

    std::mutex out_mutex;
    std::atomic<int> total_pass{0};
    std::atomic<uint64_t> total_filtered{0};
    int total_candidates = 0;

    constexpr size_t FLUSH_THRESHOLD = 1 << 16;

    std::cerr << "[dist] Loading " << nr_total << " refs\n";
    std::vector<GenomeDb> refs;
    load_and_seed(ref_paths, 0, nr_total, refs,
                  SHAPE, seed_scaled, max_seed_occ, threads, p.gene_scaled);

    std::vector<GenomeDb> queries_storage;
    if (!is_self) {
        load_and_seed(query_paths, 0, nq_total, queries_storage,
                      SHAPE, seed_scaled, max_seed_occ, threads, p.gene_scaled, 1);
    }
    auto &queries = is_self ? refs : queries_storage;
    const int nq = static_cast<int>(queries.size());
    const int nr = nr_total;
    const bool triangle = is_self && !no_triangle;

    uint64_t np64 = 0;
    for (int qi = 0; qi < nq; qi++) {
        if (queries[qi].seed_count() == 0) continue;
        const uint32_t first = triangle ? static_cast<uint32_t>(qi) + 1 : 0;
        if (first < static_cast<uint32_t>(nr)) np64 += nr - first - (is_self && !triangle ? 1 : 0);
    }
    const int np = static_cast<int>(std::min<uint64_t>(np64, INT32_MAX));
    total_candidates += np;
    if (np == 0) return 0;

    ProgressBar prog(np, "dist");
    std::atomic<int> progress{0};

#pragma omp parallel num_threads(threads)
    {
        PairBufs bufs;
        std::vector<QueryHashGroup> qgroups;
        std::string outbuf;
        outbuf.reserve(FLUSH_THRESHOLD * 2);
        char linebuf[1024];

        auto do_pair = [&](int qi, uint32_t rid) {
            DistResult res;
            if (!compare_pair(queries[qi], refs[rid], qgroups, p, bufs, res)) {
                total_filtered++;
            } else {
                outbuf.append(linebuf, format_dist_row(linebuf, sizeof(linebuf), res));
                total_pass++;
                if (outbuf.size() >= FLUSH_THRESHOLD) {
                    std::lock_guard<std::mutex> lock(out_mutex);
                    out.write(outbuf.data(), outbuf.size());
                    outbuf.clear();
                }
            }
            const int pg = ++progress;
            if ((pg & (PROGRESS_INTERVAL - 1)) == 0) prog.update(pg);
        };

        if (nq >= threads * 2) {
#pragma omp for schedule(dynamic, 1)
            for (int qi = 0; qi < nq; qi++) {
                if (interrupted()) continue;
                if (queries[qi].seed_count() == 0) continue;
                group_query_seeds(queries[qi].seed_span(), qgroups);
                const uint32_t first_rid =
                    triangle ? static_cast<uint32_t>(qi) + 1 : 0;
                for (uint32_t rid = first_rid; rid < static_cast<uint32_t>(nr); rid++) {
                    if (is_self && rid == static_cast<uint32_t>(qi)) continue;
                    do_pair(qi, rid);
                }
            }
        } else {
            for (int qi = 0; qi < nq; qi++) {
                if (interrupted()) continue;
                if (queries[qi].seed_count() == 0) continue;
                group_query_seeds(queries[qi].seed_span(), qgroups);
                const int64_t first_rid = triangle ? qi + 1 : 0;
#pragma omp for schedule(dynamic, 8)
                for (int64_t rid = first_rid; rid < static_cast<int64_t>(nr); rid++) {
                    if (is_self && rid == qi) continue;
                    do_pair(qi, static_cast<uint32_t>(rid));
                }
            }
        }

        if (!outbuf.empty()) {
            std::lock_guard<std::mutex> lock(out_mutex);
            out.write(outbuf.data(), outbuf.size());
        }
    }
    prog.finish();

    total_pairs = total_pass.load();
    std::cerr << "[dist] Pairs: " << total_candidates << "\n";
    std::cerr << "[dist] Filtered: " << total_filtered.load() << "\n";
    std::cerr << "[dist] Done. " << total_pairs
              << " -> " << output_path << "\n";

    if (interrupted()) throw Interrupted{};
    return total_pairs;
}
