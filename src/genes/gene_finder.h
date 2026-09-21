#pragma once
#include "db/types.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include "genes/period3.h"
#include "genes/codon_model.h"
#include "genes/orf.h"
#include <vector>
#include <memory>
#include <set>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <csignal>
#include <omp.h>

static constexpr double PSI_REF[4][3] = {
    {0.274534, 0.314761, 0.207881},
    {0.219326, 0.240591, 0.282598},
    {0.331816, 0.171595, 0.257585},
    {0.174324, 0.273053, 0.251936}
};

static const double PSI_LOGREF[4][3] = {
    {std::log(PSI_REF[0][0]), std::log(PSI_REF[0][1]), std::log(PSI_REF[0][2])},
    {std::log(PSI_REF[1][0]), std::log(PSI_REF[1][1]), std::log(PSI_REF[1][2])},
    {std::log(PSI_REF[2][0]), std::log(PSI_REF[2][1]), std::log(PSI_REF[2][2])},
    {std::log(PSI_REF[3][0]), std::log(PSI_REF[3][1]), std::log(PSI_REF[3][2])}
};

static constexpr double CHAAI_TRAIN_PSI = 0.015;

static constexpr double CHAAI_OVLP = 0.6;
static constexpr int CHAAI_FIT_LO = 30, CHAAI_FIT_HI = 90;
static constexpr double CHAAI_SHRINK_K = 200.0;

template <class P>
static inline std::vector<uint32_t> collect_indices(size_t n, int threads, P &&pred)
{
    constexpr size_t BLK = 1u << 16;
    const long nb = static_cast<long>((n + BLK - 1) / BLK);
    std::vector<size_t> cnt(nb + 1, 0);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long b = 0; b < nb; b++) {
        const size_t beg = static_cast<size_t>(b) * BLK, end = std::min(beg + BLK, n);
        size_t c = 0;
        for (size_t i = beg; i < end; i++) if (pred(i)) c++;
        cnt[b + 1] = c;
    }
    for (long b = 0; b < nb; b++) cnt[b + 1] += cnt[b];
    std::vector<uint32_t> out(cnt[nb]);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long b = 0; b < nb; b++) {
        const size_t beg = static_cast<size_t>(b) * BLK, end = std::min(beg + BLK, n);
        size_t o = cnt[b];
        for (size_t i = beg; i < end; i++) if (pred(i)) out[o++] = static_cast<uint32_t>(i);
    }
    return out;
}

struct CsRow {
    const uint32_t *b, *e;
    const uint32_t *begin() const { return b; }
    const uint32_t *end() const { return e; }
    size_t size() const { return static_cast<size_t>(e - b); }
    uint32_t operator[](size_t i) const { return b[i]; }
};

struct GfDiag {
    bool capture_vectors = true;
    std::vector<AASegment> kept;
    std::vector<AASegment> scored;
    std::vector<double>    sc, psi, asym;
    std::vector<uint8_t>   train_used;
    double cod_dens = -1;
    int    min_train_aa = 0;
    size_t n_train = 0;
};

static inline double coding_budget(const std::vector<double> &obs, int maxl, double N, double p_stop,
                                   int min_orf, const std::vector<double> &pw)
{
    double budget = 0;
    for (int L = min_orf; L <= maxl; L++) {
        if (obs[L] <= 0) continue;
        const double exp_at = 2.0 * N * p_stop * pw[L] * p_stop;
        const double excess = obs[L] - exp_at;
        if (excess > 0) budget += excess * L * 3.0;
    }
    return budget;
}

struct FilterState {
    int bsz = 0;
    int min_orf = MIN_SEG_LEN_DEFAULT;
    int nc  = 0;
    std::vector<int>     contig_blk_offset;
    std::vector<uint8_t> coding_blk;
    double cod_dens = -1;

    size_t S  = 0;
    int64_t iS = 0;
    std::vector<int>    seg_ci;
    bool contig_ordered = false;
    std::vector<double> asym;

    std::vector<uint8_t>  asym_train_ok, rbs;
    int                   min_train_aa = 100;
    std::vector<uint32_t> train_idx;

    double bg[N_DICODONS]   = {};
    double cod1[N_DICODONS] = {};
    double lo[N_DICODONS]   = {};
    double bg_total   = 0;
    double cod1_total = 0;
    const double *BG = nullptr;
    double BG_T      = 0;

    std::vector<double>  sc, psi;
    std::vector<uint8_t> keep_flag;
    std::vector<uint64_t> run_hist;
};

static inline bool coding_map(const ContigView &gv, FilterState &st, int threads,
                                volatile std::sig_atomic_t &interrupted)
{
    const int bsz = st.bsz;
    std::vector<int>     &contig_blk_offset = st.contig_blk_offset;
    std::vector<uint8_t> &coding_blk        = st.coding_blk;

    int n_blocks = 0;
    contig_blk_offset.resize(gv.n_contigs());
    for (int ci = 0; ci < gv.n_contigs(); ci++) {
        contig_blk_offset[ci] = n_blocks;
        n_blocks += static_cast<int>(gv.lens[ci] / bsz);
    }
    if (n_blocks < 1) return false;

    const double blk_thr = CHAAI_ASYM_CUT;
    coding_blk.assign(static_cast<size_t>(n_blocks), 0);
    int n_cod_blk = 0;
    st.nc = gv.n_contigs();

#pragma omp parallel for schedule(dynamic, 64) num_threads(threads) reduction(+:n_cod_blk)
    for (int bi = 0; bi < n_blocks; bi++) {
        if (interrupted) continue;
        const int ci = static_cast<int>(
            std::upper_bound(contig_blk_offset.begin(),
                             contig_blk_offset.end(), bi) -
            contig_blk_offset.begin()) - 1;
        const int local_b = bi - contig_blk_offset[ci];
        const double snr = compute_asymmetry(gv.seqs[ci],
                                          static_cast<uint32_t>(local_b) * bsz, bsz);
        if (snr > blk_thr) {
            coding_blk[bi] = 1;
            n_cod_blk++;
        }
    }
    st.cod_dens = static_cast<double>(n_cod_blk) / n_blocks;
    return true;
}

static inline void locate_segments(const std::vector<AASegment> &segments, const ContigView &gv,
                                   FilterState &st, int threads)
{
    std::vector<int> &seg_ci = st.seg_ci;
    seg_ci.assign(st.S, 0);
    const int nc = gv.n_contigs();
    int ci = 0;
    st.contig_ordered = true;
    for (size_t i = 0; i < st.S; i++) {
        const uint64_t p = segments[i].nuc_start;
        while (ci + 1 < nc && p >= gv.offsets[ci + 1]) ci++;
        if (p < gv.offsets[ci]) { ci = gv.find_contig(p); st.contig_ordered = false; }
        seg_ci[i] = ci;
    }
    (void)threads;
}

static inline void orf_measures(const std::vector<AASegment> &segments, const ContigView &gv,
                                FilterState &st, int threads,
                                volatile std::sig_atomic_t &interrupted)
{
    st.asym.assign(st.S, 0.0);
    st.psi.assign(st.S, 0.0);
    st.iS = static_cast<int64_t>(st.S);
    const int64_t iS = st.iS;
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
    for (int64_t i = 0; i < iS; i++) {
        if (interrupted) continue;
        const int ci = st.seg_ci[i];
        const uint64_t local = segments[i].nuc_start - gv.offsets[ci];
        const int nlen = segments[i].aa_len * 3;
        if (local + nlen > gv.lens[ci] || nlen < 9) continue;

        const char *s = gv.seqs[ci] + local;
        const bool rev = segments[i].frame >= 3;
        uint32_t base[4][3] = {}, cod[4][3] = {};
        uint32_t cq[3][5];
        const bool bad = count_positions(s, nlen, cq);
        const uint32_t *c0 = cq[0], *c1 = cq[1], *c2 = cq[2];
        int n = 0;
        if (!bad) {
            for (int b = 0; b < 4; b++) { base[b][0] = c0[b]; base[b][1] = c1[b]; base[b][2] = c2[b]; }
            if (rev) for (int b = 0; b < 4; b++) { cod[b][0] = c2[3 - b]; cod[b][1] = c1[3 - b]; cod[b][2] = c0[3 - b]; }
            else     for (int b = 0; b < 4; b++) { cod[b][0] = c0[b];     cod[b][1] = c1[b];     cod[b][2] = c2[b];     }
        } else {
            for (n = 0; n + 2 < nlen; n += 3) {
                const uint8_t a0 = static_cast<uint8_t>(s[n]);
                const uint8_t a1 = static_cast<uint8_t>(s[n + 1]);
                const uint8_t a2 = static_cast<uint8_t>(s[n + 2]);
                if (a0 <= 3) base[a0][0]++;
                if (a1 <= 3) base[a1][1]++;
                if (a2 <= 3) base[a2][2]++;
                if ((a0 | a1 | a2) > 3) continue;
                if (rev) { cod[3 - a2][0]++; cod[3 - a1][1]++; cod[3 - a0][2]++; }
                else     { cod[a0][0]++;     cod[a1][1]++;     cod[a2][2]++;     }
            }
        }

        st.asym[i] = asym_from_counts(base);

        double tot = 0, row[4] = {};
        for (int b = 0; b < 4; b++)
            for (int q = 0; q < 3; q++) { row[b] += cod[b][q]; tot += cod[b][q]; }
        if (tot > 0) {
            double ll = 0;
            for (int b = 0; b < 4; b++) {
                ll += cod[b][0] * PSI_LOGREF[b][0] + cod[b][1] * PSI_LOGREF[b][1]
                    + cod[b][2] * PSI_LOGREF[b][2];
                if (row[b] > 0) ll -= row[b] * std::log((row[b] + 0.5) / (tot + 2.0));
            }
            st.psi[i] = ll / tot;
        }
    }
}

static inline bool training_set(const std::vector<AASegment> &segments, const ContigView &gv,
                                FilterState &st, int threads, GfDiag *diag)
{
    const size_t S  = st.S;
    const int64_t iS = st.iS;
    const std::vector<double> &asym          = st.asym;
    std::vector<uint8_t>      &asym_train_ok = st.asym_train_ok;
    std::vector<uint8_t>      &rbs           = st.rbs;

    asym_train_ok.assign(S, 0);
    rbs.assign(S, 0);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int64_t i = 0; i < iS; i++)
        asym_train_ok[i] = asym[i] >= CHAAI_ASYM_CUT_TRAIN ? 1 : 0;

    int rbs_floor = 1 << 30;
    auto score_rbs_down_to = [&](int floor) {
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
        for (int64_t i = 0; i < iS; i++) {
            if (!asym_train_ok[i] || st.psi[i] >= CHAAI_TRAIN_PSI) continue;
            const AASegment &g = segments[i];
            if (g.aa_len < floor || g.aa_len >= rbs_floor) continue;
            const int ci = st.seg_ci[i];
            const int64_t local = static_cast<int64_t>(g.nuc_start - gv.offsets[ci]);
            const bool rev = g.frame >= 3;
            rbs[i] = orf_has_rbs(gv.seqs[ci], gv.lens[ci], rev ? local + static_cast<int64_t>(g.aa_len) * 3 : local, rev);
        }
        rbs_floor = floor;
    };

    int min_train_aa = 100;
    size_t n_train = 0;
    for (int attempt = 0; attempt < 4; attempt++) {
        score_rbs_down_to(min_train_aa);
        size_t cnt = 0;
#pragma omp parallel for schedule(static) num_threads(threads) reduction(+ : cnt)
        for (int64_t i = 0; i < iS; i++)
            if (asym_train_ok[i] && segments[i].aa_len >= min_train_aa &&
                (st.psi[i] >= CHAAI_TRAIN_PSI || rbs[i])) cnt++;
        n_train = cnt;
        if (n_train >= 30) break;
        min_train_aa = std::max(30, min_train_aa / 2);
    }
    st.min_train_aa = min_train_aa;
    if (diag) { diag->min_train_aa = min_train_aa; diag->n_train = n_train; }
    if (n_train < MIN_TRAIN_SEGMENTS) return false;

    st.train_idx = collect_indices(S, threads, [&](size_t i) {
        return asym_train_ok[i] && segments[i].aa_len >= min_train_aa &&
               (st.psi[i] >= CHAAI_TRAIN_PSI || rbs[i]);
    });
    return true;
}

static inline bool noncoding_background(const std::vector<AASegment> &segments, const ContigView &gv,
                                     FilterState &st, int threads,
                                     volatile std::sig_atomic_t &interrupted)
{
    const int bsz = st.bsz;
    const int nc  = st.nc;
    double (&bg)[N_DICODONS]   = st.bg;
    double (&cod1)[N_DICODONS] = st.cod1;
    double &bg_total   = st.bg_total;
    double &cod1_total = st.cod1_total;

    std::vector<BgChunk> bg_chunks;
    for (int ci = 0; ci < nc; ci++) {
        const size_t valid = gv.lens[ci] > 5 ? gv.lens[ci] - 5 : 0;
        const int bb = st.contig_blk_offset[ci];
        for (size_t s = 0; s < valid; s += BG_CHUNK) {
            bg_chunks.push_back({ci, s, std::min(s + BG_CHUNK, valid), bb});
        }
    }

    accumulate_background(bg_chunks, gv, st.coding_blk, bsz, true, interrupted, threads, bg, bg_total);

    if (bg_total < 1000) {
        std::memset(bg, 0, sizeof(bg));
        bg_total = 0;
        accumulate_background(bg_chunks, gv, st.coding_blk, bsz, false, interrupted, threads, bg, bg_total);
    }
    if (bg_total < MIN_BG_COUNT) return false;

    accumulate_dicodons(segments, gv, st.seg_ci, st.train_idx, cod1, cod1_total);
    if (cod1_total < MIN_BG_COUNT) return false;
    return true;
}

static inline void dicodon_log_odds(const std::vector<AASegment> &segments, const ContigView &gv,
                                 FilterState &st, int threads)
{
    const size_t S  = st.S;
    const int64_t iS = st.iS;
    const int min_train_aa = st.min_train_aa;
    const std::vector<int>      &seg_ci      = st.seg_ci;
    const std::vector<uint8_t>  &asym_train_ok = st.asym_train_ok;
    const std::vector<uint32_t> &train_idx   = st.train_idx;
    double (&lo)[N_DICODONS] = st.lo;
    std::vector<double> &sc  = st.sc;

    st.BG   = st.bg;
    st.BG_T = st.bg_total;
    const double *BG = st.BG;
    const double BG_T = st.BG_T;
    build_log_odds(st.cod1, st.cod1_total, BG, BG_T, lo);

    sc.assign(S, 0.0);
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
    for (int64_t i = 0; i < iS; i++) {
        if (!asym_train_ok[i] || segments[i].aa_len < min_train_aa) continue;
        const int ci = seg_ci[i];
        sc[i] = score_segment(segments[i], gv.seqs[ci], gv.offsets[ci], lo);
    }

    std::vector<double> tsc(train_idx.size());
    const long nti = static_cast<long>(train_idx.size());
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long m = 0; m < nti; m++) tsc[m] = sc[train_idx[m]];
    if (tsc.size() >= 20) {
        std::sort(tsc.begin(), tsc.end());
        const double med = tsc[tsc.size() / 2];
        double cod2[N_DICODONS] = {};
        double cod2_total = 0;
        const std::vector<uint32_t> use2 = collect_indices(S, threads, [&](size_t i) {
            return asym_train_ok[i] && segments[i].aa_len >= min_train_aa && sc[i] >= med;
        });
        accumulate_dicodons(segments, gv, seg_ci, use2, cod2, cod2_total);
        if (cod2_total >= MIN_BG_COUNT) {
            build_log_odds(cod2, cod2_total, BG, BG_T, lo);
        }
    }
}

static inline void dicodon_score(const std::vector<AASegment> &segments, const ContigView &gv,
                                 FilterState &st, int threads)
{
    const int64_t iS = st.iS;
    const std::vector<int>  &seg_ci = st.seg_ci;
    const double (&lo)[N_DICODONS] = st.lo;
    std::vector<double> &sc = st.sc;

#pragma omp parallel for schedule(dynamic, 64) num_threads(threads)
    for (int64_t i = 0; i < iS; i++)
        sc[i] = score_segment(segments[i], gv.seqs[seg_ci[i]], gv.offsets[seg_ci[i]], lo);
}

static inline void capture_diag_vectors(const std::vector<AASegment> &segments,
                                        const FilterState &st, GfDiag *diag)
{
    const size_t S = st.S;

    if (diag && diag->capture_vectors) {
        diag->asym = st.asym;
        diag->sc = st.sc;
        diag->psi = st.psi;
        diag->scored = segments;
        diag->train_used.assign(S, 0);
        for (uint32_t i : st.train_idx) diag->train_used[i] = 1;
    }
}

static inline double geometric_q_band(const std::vector<uint64_t> &h, int lo, int hi)
{
    double N = 0, sL = 0;
    for (int L = lo; L <= hi; L++) { N += static_cast<double>(h[L]); sL += static_cast<double>(h[L]) * L; }
    if (N < 30.0) return -1.0;
    const double mb = std::min(std::max(sL / N, lo + 0.5), 0.5 * (lo + hi) - 0.5);
    double a = 1e-9, b = 1.0 - 1e-12;
    for (int it = 0; it < 60; it++) {
        const double r = 0.5 * (a + b);
        double num = 0, den = 0, p = 1;
        for (int L = lo; L <= hi; L++) { num += L * p; den += p; p *= r; }
        if (num / den < mb) a = r; else b = r;
    }
    return 1.0 - 0.5 * (a + b);
}

static inline void budget_walk(const double *keys, const double *keys2,
                               const std::vector<AASegment> &segments, const CsRow &idx,
                               double budget, uint64_t base, uint64_t len, std::vector<uint8_t> &keep)
{
    std::vector<uint32_t> ord(idx.begin(), idx.end());
    std::sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b){ return keys[a] > keys[b]; });
    std::set<std::pair<uint64_t, uint32_t>> acc;
    std::vector<uint32_t> beaten;
    uint64_t maxspan = 0;
    double used = 0;
    for (uint32_t i : ord) {
        if (used >= budget) break;
        const uint64_t gs = segments[i].nuc_start - base;
        const uint64_t ge = std::min(gs + static_cast<uint64_t>(segments[i].aa_len) * 3, len);
        if (ge <= gs) continue;
        beaten.clear();
        bool lose = false;
        for (auto it = acc.lower_bound({gs > maxspan ? gs - maxspan : 0, 0});
             it != acc.end() && it->first < ge; ++it) {
            const uint32_t j = it->second;
            const uint64_t js = it->first;
            const uint64_t je = std::min(js + static_cast<uint64_t>(segments[j].aa_len) * 3, len);
            if (je <= gs) continue;
            const uint64_t ov = std::min(ge, je) - std::max(gs, js);
            if (static_cast<double>(ov) > CHAAI_OVLP * static_cast<double>(std::min(ge - gs, je - js))) {
                if (keys2[i] <= keys2[j]) { lose = true; break; }
                beaten.push_back(j);
            }
        }
        if (lose) { used += static_cast<double>(ge - gs); continue; }
        for (uint32_t j : beaten) {
            keep[j] = 0;
            acc.erase({segments[j].nuc_start - base, j});
        }
        keep[i] = 1;
        acc.insert({gs, i});
        used += static_cast<double>(ge - gs);
        maxspan = std::max(maxspan, ge - gs);
    }
}

static inline void keep_by_score(const std::vector<AASegment> &segments, const ContigView &gv,
                                 FilterState &st, int threads)
{
    const size_t S = st.S;
    const int64_t iS = st.iS;
    const std::vector<double> &sc = st.sc, &asym = st.asym, &psi = st.psi;
    std::vector<uint8_t> &keep = st.keep_flag;
    keep.assign(S, 0);
    if (S == 0) return;

    constexpr int NF = 4, MAXL = RUN_HIST_MAX;
    std::unique_ptr<double[]> keys_buf(new double[S]), keys2_buf(new double[S]), loglen_buf(new double[S]);
    double *keys = keys_buf.get(), *keys2 = keys2_buf.get(), *loglen = loglen_buf.get();
    static const std::vector<double> LOGLEN = [] {
        std::vector<double> t(65536);
        for (int L = 0; L < 65536; L++) t[L] = std::log(std::max<double>(L, 1.0));
        return t;
    }();
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int64_t i = 0; i < iS; i++) loglen[i] = LOGLEN[segments[i].aa_len];

    double q = geometric_q_band(st.run_hist, CHAAI_FIT_LO, std::min(CHAAI_FIT_HI, MAXL));
    if (q <= 0) {
        uint64_t runs = 0;
        for (uint64_t v : st.run_hist) runs += v;
        q = runs ? std::min(0.5, static_cast<double>(runs) / (2.0 * static_cast<double>(gv.total_len))) : 0.05;
    }

    const int ncc = gv.n_contigs();
    std::vector<uint32_t> cs_off(ncc + 1, 0), cs_idx(S);
    for (size_t i = 0; i < S; i++) cs_off[st.seg_ci[i] + 1]++;
    for (int ci = 0; ci < ncc; ci++) cs_off[ci + 1] += cs_off[ci];
    if (st.contig_ordered) {
#pragma omp parallel for schedule(static) num_threads(threads)
        for (int64_t i = 0; i < iS; i++) cs_idx[i] = static_cast<uint32_t>(i);
    } else {
        std::vector<uint32_t> fill(cs_off.begin(), cs_off.end() - 1);
        for (size_t i = 0; i < S; i++) cs_idx[fill[st.seg_ci[i]]++] = static_cast<uint32_t>(i);
    }

    double gm[NF] = {}, gsd[NF] = {};
    for (size_t i = 0; i < S; i++) { gm[0] += sc[i]; gm[1] += asym[i]; gm[2] += psi[i]; gm[3] += loglen[i]; }
    for (int k = 0; k < NF; k++) gm[k] /= static_cast<double>(S);
    for (size_t i = 0; i < S; i++) {
        const double d[NF] = {sc[i] - gm[0], asym[i] - gm[1], psi[i] - gm[2], loglen[i] - gm[3]};
        for (int k = 0; k < NF; k++) gsd[k] += d[k] * d[k];
    }
    for (int k = 0; k < NF; k++) gsd[k] = std::sqrt(gsd[k] / static_cast<double>(S)) + 1e-9;

    std::vector<double> pw(MAXL + 1);
    for (int L = 0; L <= MAXL; L++) pw[L] = std::pow(1.0 - q, static_cast<double>(L));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int ci = 0; ci < ncc; ci++) {
        const CsRow idx{cs_idx.data() + cs_off[ci], cs_idx.data() + cs_off[ci + 1]};
        const size_t n = idx.size();
        if (n == 0) continue;
        const double lam = static_cast<double>(n) / (static_cast<double>(n) + CHAAI_SHRINK_K);
        double m[NF] = {}, sd[NF] = {};
        for (uint32_t i : idx) { m[0] += sc[i]; m[1] += asym[i]; m[2] += psi[i]; m[3] += loglen[i]; }
        for (int k = 0; k < NF; k++) m[k] /= static_cast<double>(n);
        for (uint32_t i : idx) {
            const double d[NF] = {sc[i] - m[0], asym[i] - m[1], psi[i] - m[2], loglen[i] - m[3]};
            for (int k = 0; k < NF; k++) sd[k] += d[k] * d[k];
        }
        double mm[NF], ss[NF];
        for (int k = 0; k < NF; k++) {
            mm[k] = lam * m[k] + (1.0 - lam) * gm[k];
            ss[k] = lam * std::sqrt(sd[k] / static_cast<double>(n)) + (1.0 - lam) * gsd[k] + 1e-9;
        }
        for (uint32_t i : idx) {
            const double v[NF] = {sc[i], asym[i], psi[i], loglen[i]};
            double z[NF];
            for (int k = 0; k < NF; k++) z[k] = (v[k] - mm[k]) / ss[k];
            keys[i]  = z[0] + z[1] + z[2] + z[3];
            keys2[i] = z[0] + z[2] + z[3];
        }
        std::vector<double> obs(MAXL + 1, 0.0);
        for (uint32_t i : idx) obs[std::min<int>(MAXL, segments[i].aa_len)] += 1.0;
        const double budget = coding_budget(obs, MAXL, static_cast<double>(gv.lens[ci]), q, st.min_orf, pw);
        budget_walk(keys, keys2, segments, idx, budget, gv.offsets[ci], gv.lens[ci], keep);
    }
}

static inline void compact_kept(std::vector<AASegment> &segments, const FilterState &st,
                                int threads)
{
    const size_t S = st.S;
    const std::vector<uint8_t> &keep_flag = st.keep_flag;

    constexpr size_t CBLK = 1u << 16;
    const long nblk = static_cast<long>((S + CBLK - 1) / CBLK);
    std::vector<size_t> bcnt(nblk + 1, 0);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long b = 0; b < nblk; b++) {
        const size_t beg = static_cast<size_t>(b) * CBLK, end = std::min(beg + CBLK, S);
        size_t c = 0;
        for (size_t i = beg; i < end; i++) c += keep_flag[i];
        bcnt[b + 1] = c;
    }
    for (long b = 0; b < nblk; b++) bcnt[b + 1] += bcnt[b];
    std::vector<AASegment> kept(bcnt[nblk]);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (long b = 0; b < nblk; b++) {
        const size_t beg = static_cast<size_t>(b) * CBLK, end = std::min(beg + CBLK, S);
        size_t o = bcnt[b];
        for (size_t i = beg; i < end; i++) {
            if (!keep_flag[i]) continue;
            kept[o++] = segments[i];
        }
    }
    segments = std::move(kept);
}

static inline void apply_coding_filter(
    std::vector<AASegment> &segments,
    const ContigView &gv,
    volatile std::sig_atomic_t &interrupted,
    int map_window = 300,
    int min_orf = MIN_SEG_LEN_DEFAULT,
    double *out_coding_density = nullptr,
    int threads = 1,
    GfDiag *diag = nullptr,
    const std::vector<uint64_t> *run_hist = nullptr)
{
    if (out_coding_density) *out_coding_density = -1;
    if (segments.empty() || gv.total_len < 9) return;

    FilterState st;
    st.bsz = std::max(map_window, 30);
    st.min_orf = min_orf;
    if (run_hist) st.run_hist = *run_hist;
    else st.run_hist.assign(RUN_HIST_MAX + 1, 0);

    if (!coding_map(gv, st, threads, interrupted)) {
        if (out_coding_density) *out_coding_density = 0;
        return;
    }
    if (out_coding_density) *out_coding_density = st.cod_dens;
    if (diag) diag->cod_dens = st.cod_dens;

    st.S = segments.size();

    locate_segments(segments, gv, st, threads);
    orf_measures(segments, gv, st, threads, interrupted);

    if (!training_set(segments, gv, st, threads, diag) ||
        !noncoding_background(segments, gv, st, threads, interrupted)) {
        st.sc.assign(st.S, 0.0);
        keep_by_score(segments, gv, st, threads);
        compact_kept(segments, st, threads);
        if (diag && diag->capture_vectors) diag->kept = segments;
        return;
    }

    dicodon_log_odds(segments, gv, st, threads);
    dicodon_score(segments, gv, st, threads);

    capture_diag_vectors(segments, st, diag);
    keep_by_score(segments, gv, st, threads);
    compact_kept(segments, st, threads);

    if (diag && diag->capture_vectors) diag->kept = segments;
}
