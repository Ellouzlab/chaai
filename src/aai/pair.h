#pragma once

#include "aai/align.h"
#include "aai/merge.h"
#include "aai/chain.h"
#include "aai/pairwise.h"
#include "db/types.h"

#include <cstdio>
#include <vector>

struct PairBufs {
    std::vector<SeedMatch> matches;
    std::vector<Chain> chains;
    AlignPairBufs align;

    void release_matches() {
        if (matches.capacity() > (1 << 20) && matches.capacity() > matches.size() * 8) {
            matches.clear();
            matches.shrink_to_fit();
        }
    }
};

static inline bool compare_pair(const GenomeDb &q, const GenomeDb &r,
                                const std::vector<QueryHashGroup> &qgroups,
                                const DistParams &p, PairBufs &b, DistResult &out,
                                bool screen = true,
                                std::vector<OrfPairHit> *orf_hits = nullptr)
{
    sorted_merge_matches(q.seed_span(), r.seed_span(), qgroups,
                         screen ? static_cast<size_t>(p.min_chain_seeds) : 0,
                         p.max_pair_matches, b.matches);

    if (screen) {
        const size_t n = b.matches.size();
        if (n < static_cast<size_t>(p.min_chain_seeds)) return false;
    }

    if (orf_hits) orf_hits->clear();
    out = align_pair(q, r, b.matches, b.chains, p, 1, b.align, orf_hits);
    b.release_matches();

    return !screen || (out.aai >= p.min_aai && out.total_aligned >= p.min_aligned_aa &&
                       out.matched_orfs >= p.min_matched_orfs);
}

static inline int format_dist_row(char *buf, size_t n, const DistResult &res)
{
    return snprintf(buf, n,
                    "%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t"
                    "%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%lu\t%lu\n",
                    res.query_name.c_str(), res.ref_name.c_str(),
                    res.aai, res.gcov_q, res.gcov_r, res.gcov,
                    res.orf_depth_q, res.orf_depth_r, res.orf_depth,
                    res.composite_score, res.matched_orfs, res.orf_hit,
                    static_cast<unsigned long>(res.query_nuc_len),
                    static_cast<unsigned long>(res.ref_nuc_len));
}

static inline const char *DIST_HEADER =
    "query\tref\taai\t"
    "chaining_genome_coverage_q\tchaining_genome_coverage_r\t"
    "chaining_genome_coverage\t"
    "avg_orf_chain_coverage_q\tavg_orf_chain_coverage_r\t"
    "avg_orf_chain_coverage\t"
    "composite_score\tmatched_orfs\torf_hit\tquery_len\tref_len\n";
