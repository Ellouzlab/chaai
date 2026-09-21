#pragma once
#include "db/types.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include "genes/period3.h"
#include <vector>
#include <cmath>
#include <csignal>
#include <omp.h>

static constexpr int N_DICODONS = 4096;

static inline int model_dim() { return N_DICODONS; }

static inline void accumulate_dicodons(
    const std::vector<AASegment> &segments,
    const ContigView &gv,
    const std::vector<int> &seg_ci,
    const std::vector<uint32_t> &members,
    double *counts, double &total)
{
    const long n = static_cast<long>(members.size());
#pragma omp parallel
    {
        double lcounts[N_DICODONS] = {};
        double lt = 0;
#pragma omp for schedule(dynamic, 64)
        for (long m = 0; m < n; m++) {
            const uint32_t si = members[m];
            const int ci = seg_ci[si];
            int prev = -1;
            for_each_codon(segments[si], gv.seqs[ci], gv.offsets[ci], [&](int cur) {
                if (prev >= 0 && cur >= 0) { lcounts[prev * 64 + cur]++; lt++; }
                prev = cur;
            });
        }
#pragma omp critical
        {
            for (int i = 0; i < N_DICODONS; i++) counts[i] += lcounts[i];
            total += lt;
        }
    }
}

static inline void accumulate_dicodons(
    const std::vector<AASegment> &segments,
    const ContigView &gv,
    const std::vector<int> &seg_ci,
    const std::vector<bool> &use,
    double *counts, double &total)
{
    std::vector<uint32_t> members;
    for (size_t i = 0; i < use.size(); i++) if (use[i]) members.push_back(static_cast<uint32_t>(i));
    accumulate_dicodons(segments, gv, seg_ci, members, counts, total);
}

static inline void build_log_odds(
    const double *coding, double coding_total,
    const double *bg, double bg_total,
    double *log_odds)
{
    constexpr double PSEUDO = 1.0;
    const int n = model_dim();
    for (int i = 0; i < n; i++) {
        const double p_c = (coding[i] + PSEUDO) / (coding_total + n * PSEUDO);
        const double p_b = (bg[i] + PSEUDO) / (bg_total + n * PSEUDO);
        log_odds[i] = std::log2(p_c / p_b);
    }
}

static inline double score_segment(
    const AASegment &seg, const char *contig_seq,
    uint64_t contig_offset, const double *log_odds)
{
    double score = 0;
    int count = 0;
    int prev = -1;
    for_each_codon(seg, contig_seq, contig_offset, [&](int cur) {
        if (prev >= 0 && cur >= 0) { score += log_odds[prev * 64 + cur]; count++; }
        prev = cur;
    });
    return (count > 0) ? score / count : 0;
}

static constexpr double MIN_BG_COUNT = 100.0;
static constexpr size_t MIN_TRAIN_SEGMENTS = 10;
static constexpr size_t BG_CHUNK = 500000;

struct BgChunk { int ci; size_t start, end; int blk_base; };

static inline void accumulate_background(
    const std::vector<BgChunk> &bg_chunks, const ContigView &gv,
    const std::vector<uint8_t> &coding_blk, int bsz, bool skip_coding,
    volatile std::sig_atomic_t &interrupted, int threads,
    double *bg, double &bg_total)
{
    const int n_bg_chunks = static_cast<int>(bg_chunks.size());
#pragma omp parallel num_threads(threads)
    {
        double lbg[N_DICODONS] = {};
        double lt = 0;
#pragma omp for schedule(dynamic, 1)
        for (int ch = 0; ch < n_bg_chunks; ch++) {
            if (interrupted) continue;
            const auto &c = bg_chunks[ch];
            const char *cseq = gv.seqs[c.ci];
            const int n_blk = static_cast<int>(gv.lens[c.ci] / bsz);
            size_t pos = c.start;
            while (pos < c.end) {
                size_t blk_end = c.end;
                if (skip_coding) {
                    const int bidx = static_cast<int>(pos / bsz);
                    if (bidx < n_blk && coding_blk[c.blk_base + bidx]) {
                        pos = static_cast<size_t>(bidx + 1) * bsz;
                        continue;
                    }
                    blk_end = std::min(static_cast<size_t>(bidx + 1) * bsz, c.end);
                }
                unsigned hx = 0; int run = 0;
                for (size_t k = pos; k < pos + 5; k++) {
                    const uint8_t b = static_cast<uint8_t>(cseq[k]);
                    if (b > 3) { run = 0; hx = 0; } else { run++; hx = ((hx << 2) | b) & 4095u; }
                }
                for (; pos < blk_end; pos++) {
                    const uint8_t b = static_cast<uint8_t>(cseq[pos + 5]);
                    if (b > 3) { run = 0; hx = 0; continue; }
                    hx = ((hx << 2) | b) & 4095u;
                    if (++run < 6) continue;
                    lbg[hx]++;
                    lt++;
                }
            }
        }
#pragma omp critical
        {
            for (int i = 0; i < N_DICODONS; i++) bg[i] += lbg[i];
            bg_total += lt;
        }
    }
}
