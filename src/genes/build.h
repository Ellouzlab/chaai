#pragma once

#include <cmath>
#include <tuple>
#include <filesystem>
#include <iostream>
#include <csignal>
#include <omp.h>

#include "db/types.h"
#include "db/io.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include "genes/orf.h"
#include "genes/fasta.h"
#include "genes/gene_finder.h"
#include "genes/pipeline.h"
#include "run/progress.h"

namespace fs = std::filesystem;

static void db_sigint_handler(int) { g_interrupted = 1; }

static constexpr uint64_t WAVE_BYTES_LIMIT = 500'000'000ULL;

static void write_genome_db(
    const GenomeContigs &gc, GenomeCalls &calls,
    const std::string &outfolder, bool verbose, bool coding_filter)
{
    GenomeDb sk;
    sk.name = gc.name;
    sk.nuc_length = gc.nuc_length;
    sk.seg_headers = std::move(calls.headers);
    sk.compute_coding_nuc();

    std::vector<ContigEntry> ctgs;
    ctgs.reserve(gc.contig_names.size());
    for (size_t ci = 0; ci < gc.contig_names.size(); ci++) {
        ctgs.push_back({gc.contig_names[ci],
                        static_cast<uint64_t>(gc.offsets[ci]),
                        static_cast<uint64_t>(gc.lens[ci])});
    }

    write_db(sk.name, sk.nuc_length, sk.coding_nuc,
                 ctgs, sk.seg_headers, calls.packed,
                 outfolder + "/" + safe_filename(sk.name) + ".chaai");

    if (verbose) {
        const char *filt = coding_filter ? "applied" : "off";
        const double db_cod = sk.nuc_length > 0
            ? 100.0 * sk.coding_nuc / sk.nuc_length : 0;
#pragma omp critical
        fprintf(stderr,
                "\r[verbose] %-25s  %lubp  cod=%.1f%%  filter=%-7s  segs=%d\n",
                sk.name.c_str(), static_cast<unsigned long>(sk.nuc_length),
                db_cod, filt, static_cast<int>(sk.seg_headers.size()));
    }
}

static int build_from_folder(
    const std::string &infolder,
    const std::string &outfolder,
    int threads,
    bool coding_filter = false,
    int min_orf = MIN_SEG_LEN_DEFAULT,
    int map_window = 300,
    bool verbose = false)
{
    fs::create_directories(outfolder);
    g_interrupted = 0;
    auto prev_handler = std::signal(SIGINT, db_sigint_handler);

    struct FInfo { std::string path; uintmax_t size; };
    std::vector<FInfo> files;
    for (auto &e : fs::directory_iterator(infolder)) {
        if (e.is_regular_file() && is_fasta_ext(e.path().string())) {
            files.push_back({e.path().string(), e.file_size()});
        }
    }
    std::sort(files.begin(), files.end(),
              [](const FInfo &a, const FInfo &b) { return a.size > b.size; });

    const int n_files = static_cast<int>(files.size());
    int total_written = 0;
    omp_set_num_threads(threads);
    ProgressBar prog(n_files, "call");

    int fi = 0;
    while (fi < n_files && !g_interrupted) {
        std::vector<GenomeContigs> wave;
        uint64_t wave_bytes = 0;
        int wave_end = fi;

        while (wave_end < n_files && !g_interrupted) {
            const auto stem = fs::path(files[wave_end].path).stem().string();
            auto gc = read_genome_contigs(files[wave_end].path, stem, true);
            wave_bytes += gc.nuc_length;
            wave.push_back(std::move(gc));
            wave_end++;
            if (wave_bytes >= WAVE_BYTES_LIMIT) break;
            if (static_cast<int>(wave.size()) >= threads * 2) break;
        }
        fi = wave_end;

        const GeneFinderParams prm{min_orf, coding_filter, map_window, threads};
        std::vector<GenomeCalls> results(wave.size());
#pragma omp parallel for schedule(dynamic, 1) if (wave.size() > 1)
        for (int g = 0; g < static_cast<int>(wave.size()); g++) {
            if (g_interrupted) continue;
            try {
                results[g] = call_genome(wave[g], prm);
            } catch (const std::exception &e) {
#pragma omp critical
                std::cerr << "[call ERROR] " << wave[g].name << ": "
                          << e.what() << std::endl;
            }
            for (auto &s : wave[g].seqs) { std::string().swap(s); }
        }

        for (int g = 0; g < static_cast<int>(wave.size()) && !g_interrupted; g++) {
            write_genome_db(wave[g], results[g], outfolder,
                                verbose, coding_filter);
            total_written++;
            prog.update(total_written);
        }
    }
    prog.finish();
    std::signal(SIGINT, prev_handler);
    if (g_interrupted) throw std::runtime_error("Interrupted");

    return total_written;
}

static int build_from_fasta(
    const std::string &infasta,
    const std::string &outfolder,
    int threads,
    bool coding_filter = false,
    int min_orf = MIN_SEG_LEN_DEFAULT,
    int map_window = 300,
    bool verbose = false)
{
    fs::create_directories(outfolder);
    omp_set_num_threads(threads);
    g_interrupted = 0;
    auto prev_handler = std::signal(SIGINT, db_sigint_handler);

    FastaReader reader(infasta);
    int total_written = 0;
    const int BATCH = std::max(threads * 4, 32);

    struct Job { std::string name, seq; };

    while (!g_interrupted) {
        std::vector<Job> batch;
        batch.reserve(BATCH);
        {
            std::string rn, rs;
            while (static_cast<int>(batch.size()) < BATCH && reader.next(rn, rs)) {
                batch.push_back({std::move(rn), std::move(rs)});
            }
        }
        if (batch.empty()) break;
        const int bn = static_cast<int>(batch.size());

#pragma omp parallel for schedule(dynamic, 1) if (bn > 1)
        for (int i = 0; i < bn; i++) {
            try {
                if (g_interrupted) continue;
                auto &job = batch[i];

                GenomeContigs gc;
                gc.name = job.name;
                gc.contig_names.push_back(job.name);
                gc.seqs.push_back(std::move(job.seq));
                gc.lens.push_back(gc.seqs[0].size());
                gc.offsets.push_back(0);
                gc.nuc_length = gc.seqs[0].size();

                const GeneFinderParams prm{min_orf, coding_filter, map_window, threads};
                auto calls = call_genome(gc, prm);
                { std::string().swap(gc.seqs[0]); }
                write_genome_db(gc, calls, outfolder,
                                    verbose, coding_filter);
            } catch (const std::exception &e) {
#pragma omp critical
                std::cerr << "[call ERROR] " << batch[i].name << ": "
                          << e.what() << std::endl;
            }
        }
        total_written += bn;
    }
    std::signal(SIGINT, prev_handler);
    if (g_interrupted) throw std::runtime_error("Interrupted");
    return total_written;
}

static void write_chaai(
    const std::string &path, const std::string &name, uint64_t nuc_length,
    const std::vector<std::tuple<std::string, uint64_t, uint64_t>> &contigs,
    const std::vector<std::tuple<int, uint32_t, int, std::string>> &segs)
{
    int inv[256] = {};
    for (int i = 0; AA_LETTERS[i]; i++) inv[static_cast<uint8_t>(AA_LETTERS[i])] = i;

    std::vector<SegHeader> headers;
    headers.reserve(segs.size());
    std::vector<uint8_t> packed;
    for (const auto &s : segs) {
        SegHeader h;
        h.frame = static_cast<uint8_t>(std::get<0>(s));
        h.nuc_start = std::get<1>(s);
        h.aa_len = static_cast<uint16_t>(std::get<2>(s));
        headers.push_back(h);
        const std::string &prot = std::get<3>(s);
        std::vector<uint8_t> idx(h.aa_len);
        for (uint16_t i = 0; i < h.aa_len; i++) idx[i] = static_cast<uint8_t>(inv[static_cast<uint8_t>(prot[i])]);
        const size_t off = packed.size();
        packed.resize(off + packed_aa_bytes(h.aa_len));
        pack_aa_5bit(idx.data(), h.aa_len, packed.data() + off);
    }
    GenomeDb sk;
    sk.name = name;
    sk.nuc_length = nuc_length;
    sk.seg_headers = headers;
    sk.compute_coding_nuc();
    std::vector<ContigEntry> ctgs;
    for (const auto &c : contigs) ctgs.push_back({std::get<0>(c), std::get<1>(c), std::get<2>(c)});
    write_db(sk.name, sk.nuc_length, sk.coding_nuc, ctgs, sk.seg_headers, packed, path);
}

static std::vector<std::tuple<int, uint32_t, int, std::string>>
find_genes(const std::string &name, const std::string &seq, const GeneFinderParams &prm)
{
    GenomeContigs gc;
    gc.name = name;
    gc.contig_names.push_back(name);
    gc.seqs.push_back(seq);
    gc.lens.push_back(seq.size());
    gc.offsets.push_back(0);
    gc.nuc_length = seq.size();

    const GenomeCalls calls = call_genome(gc, prm);

    std::vector<std::tuple<int, uint32_t, int, std::string>> out;
    out.reserve(calls.headers.size());
    size_t pbyte = 0;
    std::vector<uint8_t> aa;
    for (const auto &h : calls.headers) {
        aa.resize(h.aa_len);
        unpack_aa_5bit(calls.packed.data() + pbyte, h.aa_len, aa.data());
        pbyte += packed_aa_bytes(h.aa_len);
        std::string prot(h.aa_len, 'X');
        for (uint16_t i = 0; i < h.aa_len; i++) prot[i] = AA_LETTERS[aa[i]];
        out.emplace_back(static_cast<int>(h.frame), h.nuc_start,
                         static_cast<int>(h.aa_len), std::move(prot));
    }
    return out;
}
