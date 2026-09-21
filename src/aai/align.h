#pragma once

#include "db/types.h"
#include "aai/chain.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <omp.h>

static constexpr int N_AMINO_ACIDS = 20;

static constexpr int8_t BLOSUM62[21][21] = {
    { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0, 0},
    {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, 0},
    {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 0},
    {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 0},
    { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1, 0},
    {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0},
    {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 0},
    { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, 0},
    {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0},
    {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3, 0},
    {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1, 0},
    {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0},
    {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1, 0},
    {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1, 0},
    {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2, 0},
    { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0},
    { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0},
    {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3, 0},
    {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1, 0},
    { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
};

static inline int xdrop_extend_forward(
    const uint8_t *q, const uint8_t *r,
    int q_start, int r_start, int q_len, int r_len,
    int xdrop)
{
    int score = 0, max_score = 0, best_len = 0;
    const int max_ext = std::min(q_len - q_start, r_len - r_start);
    for (int i = 0; i < max_ext; i++) {
        const uint8_t qa = q[q_start + i];
        const uint8_t ra = r[r_start + i];
        if (qa >= 20 || ra >= 20) break;
        score += BLOSUM62[qa][ra];
        if (score > max_score) {
            max_score = score;
            best_len = i + 1;
        }
        if (max_score - score > xdrop) break;
    }
    return best_len;
}

static inline int xdrop_extend_backward(
    const uint8_t *q, const uint8_t *r,
    int q_start, int r_start,
    int xdrop)
{
    int score = 0, max_score = 0, best_len = 0;
    const int max_ext = std::min(q_start, r_start);
    for (int i = 0; i < max_ext; i++) {
        const uint8_t qa = q[q_start - 1 - i];
        const uint8_t ra = r[r_start - 1 - i];
        if (qa >= 20 || ra >= 20) break;
        score += BLOSUM62[qa][ra];
        if (score > max_score) {
            max_score = score;
            best_len = i + 1;
        }
        if (max_score - score > xdrop) break;
    }
    return best_len;
}

static constexpr int SW_GAP_OPEN   = 11;
static constexpr int SW_GAP_EXTEND = 1;
static constexpr int SW_BAND       = 15;

static constexpr double KA_LAMBDA = 0.267;
static constexpr double KA_K      = 0.041;
static constexpr double KA_H      = 0.140;

static inline void ka_effective(int qlen, int rlen, double &meff, double &neff)
{
    const double m = static_cast<double>(qlen);
    const double n = static_cast<double>(rlen);
    double l = 0.0;
    for (int it = 0; it < 8; ++it) {
        const double mm = std::max(1.0, m - l);
        const double nn = std::max(1.0, n - l);
        const double next = std::log(KA_K * mm * nn) / KA_H;
        l = (next > 0.0) ? next : 0.0;
    }
    meff = std::max(1.0, m - l);
    neff = std::max(1.0, n - l);
}

static inline double compute_evalue_sum(const int *scores, int r, int qlen, int rlen)
{
    if (r <= 0) return 1e30;
    if (r > 10) r = 10;
    double meff, neff;
    ka_effective(qlen, rlen, meff, neff);
    const double lnKmn = std::log(KA_K * meff * neff);
    double T = 0.0;
    for (int i = 0; i < r; i++) T += KA_LAMBDA * static_cast<double>(scores[i]) - lnKmn;
    if (T <= 0.0) return 1e30;
    double lf_r = 0.0, lf_rm1 = 0.0;
    for (int i = 2; i <= r; i++) lf_r += std::log(static_cast<double>(i));
    for (int i = 2; i <= r - 1; i++) lf_rm1 += std::log(static_cast<double>(i));
    const double lnE = -T + (r - 1) * std::log(T) - lf_r - lf_rm1;
    return std::exp(lnE);
}

static inline double compute_evalue(int score, int qlen, int rlen)
{
    if (score <= 0) return 1e30;
    const double m = static_cast<double>(qlen);
    const double n = static_cast<double>(rlen);
    double l = 0.0;
    for (int it = 0; it < 8; ++it) {
        const double mm = std::max(1.0, m - l);
        const double nn = std::max(1.0, n - l);
        const double next = std::log(KA_K * mm * nn) / KA_H;
        l = (next > 0.0) ? next : 0.0;
    }
    const double meff = std::max(1.0, m - l);
    const double neff = std::max(1.0, n - l);
    return KA_K * meff * neff * std::exp(-KA_LAMBDA * static_cast<double>(score));
}

struct SWSeg {
    int matches, aligned, score;
    int q_beg, q_end, r_beg, r_end;
    int buf_start, buf_end;
};

struct SWBufs {
    std::vector<int16_t> h_prev, h_curr, e_curr, f_prev;
    std::vector<uint8_t> tb;
    std::vector<SWSeg> segs;
    std::vector<uint32_t> pos;

    void ensure(int max_qlen, int band) {
        const int W = 2 * band + 1;
        const size_t row = static_cast<size_t>(W);
        if (h_prev.size() < row) {
            h_prev.resize(row); h_curr.resize(row);
            e_curr.resize(row); f_prev.resize(row);
        }
        const size_t need = static_cast<size_t>(max_qlen + 1) * row;
        if (tb.size() < need) tb.resize(need);
    }
};

enum : uint8_t { TB_END = 0, TB_DIAG = 1, TB_UP = 2, TB_LEFT = 3 };

static constexpr int SW_MIN_SEG_SCORE = 5;

static inline void sw_traceback(
    const uint8_t *qa, int qs, const uint8_t *ra, int rs,
    int band, int W, const uint8_t *TB,
    int max_i, int max_k, int max_sc, bool track_query,
    SWBufs &sw)
{
    SWSeg g{};
    g.score = max_sc;
    g.buf_start = static_cast<int>(sw.pos.size());
    int matches = 0, aligned = 0;
    int q_min = INT32_MAX, q_max = -1, r_min = INT32_MAX, r_max = -1;

    int ti = max_i, tk = max_k;
    while (ti > 0 && tk >= 0 && tk < W) {
        const uint8_t op = TB[ti * W + tk];
        if (op == TB_END) break;
        if (op == TB_DIAG) {
            const int j  = ti + tk - band;
            const int qp = qs + ti - 1;
            const int rp = rs + j - 1;
            const uint8_t q_r = qa[qp];
            const uint8_t r_r = ra[rp];
            const bool eq = (q_r == r_r);
            aligned++;
            if (eq) matches++;
            if (qp < q_min) q_min = qp;
            if (qp > q_max) q_max = qp;
            if (rp < r_min) r_min = rp;
            if (rp > r_max) r_max = rp;
            const uint32_t pos = static_cast<uint32_t>(track_query ? qp : rp);
            sw.pos.push_back((pos << 1) | (eq ? 1u : 0u));
            ti--;
        } else if (op == TB_UP) {
            ti--; tk++;
        } else {
            tk--;
        }
    }
    if (aligned == 0) { sw.pos.resize(g.buf_start); return; }

    g.matches = matches; g.aligned = aligned;
    g.q_beg = q_min; g.q_end = q_max + 1;
    g.r_beg = r_min; g.r_end = r_max + 1;
    g.buf_end = static_cast<int>(sw.pos.size());
    sw.segs.push_back(g);
}

static inline void banded_sw_segments(
    const uint8_t *qa, int qs, int qe,
    const uint8_t *ra, int rs, int re,
    int band, bool track_query, SWBufs &sw)
{
    sw.segs.clear();
    sw.pos.clear();

    const int qlen = qe - qs;
    const int rlen = re - rs;
    if (qlen <= 0 || rlen <= 0) return;

    const int W = 2 * band + 1;
    sw.ensure(qlen, band);

    int16_t *H  = sw.h_prev.data();
    int16_t *Hc = sw.h_curr.data();
    int16_t *Ec = sw.e_curr.data();
    int16_t *Fp = sw.f_prev.data();
    uint8_t *TB = sw.tb.data();

    constexpr int16_t NEG = -30000;
    for (int k = 0; k < W; k++) { H[k] = 0; Fp[k] = NEG; }

    int blk_sc = 0, blk_i = 0, blk_k = 0;
    sw.pos.reserve(256);
    sw.segs.reserve(8);

    for (int i = 1; i <= qlen; i++) {
        const uint8_t q_r = qa[qs + i - 1];
        for (int k = 0; k < W; k++) { Hc[k] = 0; Ec[k] = NEG; }

        bool row_any = false;
        for (int k = 0; k < W; k++) {
            const int j = i + k - band;
            if (j < 1 || j > rlen) { TB[i * W + k] = TB_END; continue; }
            const uint8_t r_r = ra[rs + j - 1];

            const int16_t s_diag = H[k] + BLOSUM62[q_r][r_r];

            int16_t e_val = NEG;
            if (k > 0) {
                e_val = std::max(
                    static_cast<int16_t>(Hc[k - 1] - (SW_GAP_OPEN + SW_GAP_EXTEND)),
                    static_cast<int16_t>(Ec[k - 1] - SW_GAP_EXTEND));
            }
            Ec[k] = e_val;

            int16_t f_val = NEG;
            if (k + 1 < W) {
                f_val = std::max(
                    static_cast<int16_t>(H[k + 1] - (SW_GAP_OPEN + SW_GAP_EXTEND)),
                    static_cast<int16_t>(Fp[k + 1] - SW_GAP_EXTEND));
            }

            int16_t best = 0;
            uint8_t dir = TB_END;
            if (s_diag > best) { best = s_diag; dir = TB_DIAG; }
            if (e_val  > best) { best = e_val;  dir = TB_LEFT; }
            if (f_val  > best) { best = f_val;  dir = TB_UP;   }

            Hc[k] = best;
            sw.f_prev[k] = f_val;
            TB[i * W + k] = dir;

            if (best > 0) row_any = true;
            if (best > blk_sc) { blk_sc = best; blk_i = i; blk_k = k; }
        }

        std::swap(sw.h_prev, sw.h_curr);
        H  = sw.h_prev.data();
        Hc = sw.h_curr.data();

        if (!row_any) {
            if (blk_sc >= SW_MIN_SEG_SCORE)
                sw_traceback(qa, qs, ra, rs, band, W, TB, blk_i, blk_k, blk_sc, track_query, sw);
            blk_sc = 0;
        }
    }
    if (blk_sc >= SW_MIN_SEG_SCORE)
        sw_traceback(qa, qs, ra, rs, band, W, TB, blk_i, blk_k, blk_sc, track_query, sw);
}

struct DistParams
{
    int max_seed_occ;
    int chain_band, chain_max_qgap, min_chain_seeds;
    float min_aai;
    int min_aligned_aa;
    int min_chain_aligned_aa;
    int seed_weight;
    int shape_span;
    int min_matched_orfs;
    float min_orf_chain_cov;
    bool aggressive_filter;
    int xdrop;
    double max_evalue;
    size_t max_pair_matches = 0;
    bool skip_rbh = false;
    uint32_t gene_scaled = 1;

    ChainParams chain_params() const {
        return {chain_band, chain_max_qgap, min_chain_seeds, seed_weight};
    }
};

struct DistResult
{
    std::string query_name, ref_name;
    uint64_t query_nuc_len, ref_nuc_len;
    double aai;
    double gcov_q, gcov_r, gcov;
    double orf_depth_q, orf_depth_r, orf_depth;
    double composite_score;
    int total_aligned, total_aa_matches;
    int matched_orfs = 0;
    double orf_hit = 0;
};

struct OrfPairHit
{
    uint32_t q_seg, r_seg;
    uint8_t q_frame, r_frame;
    int64_t q_nuc_start, q_nuc_end;
    int64_t r_nuc_start, r_nuc_end;
    int n_seeds;
    double identity;
};

struct NucInterval { int64_t start, end; };

static inline uint64_t merge_intervals(std::vector<NucInterval> &ivs)
{
    if (ivs.empty()) return 0;
    std::sort(ivs.begin(), ivs.end(),
              [](const NucInterval &a, const NucInterval &b) {
                  return a.start < b.start;
              });
    uint64_t total = 0;
    int64_t ce = ivs[0].start;
    for (const auto &iv : ivs) {
        const int64_t s = std::max(iv.start, ce);
        if (iv.end > s) {
            total += static_cast<uint64_t>(iv.end - s);
            ce = iv.end;
        }
    }
    return total;
}

struct AnchorQR { int q, r; };
struct ClRange  { int start, end; };

struct AlignPairBufs {
    std::vector<int> order;
    std::vector<uint8_t> qa_buf, ra_buf;
    std::vector<NucInterval> qi_ivs, ri_ivs;
    SWBufs sw;
    std::vector<int> best_q, best_r;
    ChainDPBufs chain_bufs;
    std::vector<AnchorQR> anch_buf;
    std::vector<ClRange> cl_buf;
    std::vector<int> diag_buf;
    std::vector<uint64_t> match_bits, align_bits;
    std::vector<SWSeg> pending;
    std::vector<int> pend_scores;
    struct SegWin { int d, aband, q_lo, q_hi; };
    std::vector<SegWin> win_buf;
};

static inline DistResult align_pair(
    const GenomeDb &qsk, const GenomeDb &rsk,
    std::vector<SeedMatch> &matches,
    std::vector<Chain> &cbuf,
    const DistParams &p, int threads,
    AlignPairBufs &bufs,
    std::vector<OrfPairHit> *out_orf_hits = nullptr)
{
    DistResult res{};
    res.query_name = qsk.name;
    res.ref_name = rsk.name;
    res.query_nuc_len = qsk.nuc_length;
    res.ref_nuc_len = rsk.nuc_length;

    auto cp = p.chain_params();
    build_chains(matches, cbuf, cp, threads, bufs.chain_bufs);

    if (p.min_orf_chain_cov > 0) {
        size_t out = 0;
        for (size_t i = 0; i < cbuf.size(); i++) {
            const int q_span = static_cast<int>(cbuf[i].q_max) -
                               static_cast<int>(cbuf[i].q_min) + p.shape_span;
            const int r_span = static_cast<int>(cbuf[i].r_max) -
                               static_cast<int>(cbuf[i].r_min) + p.shape_span;
            const int orf_len = std::min(
                static_cast<int>(qsk.seg_headers[cbuf[i].q_seg].aa_len),
                static_cast<int>(rsk.seg_headers[cbuf[i].r_seg].aa_len));
            if (orf_len <= 0) continue;
            const float cov = static_cast<float>(std::max(q_span, r_span)) / orf_len;
            if (cov >= p.min_orf_chain_cov) {
                if (out != i) cbuf[out] = std::move(cbuf[i]);
                out++;
            }
        }
        cbuf.resize(out);
    }
    if (cbuf.empty()) return res;

    std::sort(cbuf.begin(), cbuf.end(), [](const Chain &a, const Chain &b) {
        if (a.q_seg != b.q_seg) return a.q_seg < b.q_seg;
        if (a.r_seg != b.r_seg) return a.r_seg < b.r_seg;
        return a.score > b.score;
    });

    struct PairInfo { int start, end; uint32_t q_seg, r_seg; float best_score; };
    std::vector<PairInfo> pair_info;

    for (int i = 0; i < static_cast<int>(cbuf.size()); ) {
        const uint32_t qs = cbuf[i].q_seg;
        const uint32_t rs = cbuf[i].r_seg;
        const int start = i;
        while (i < static_cast<int>(cbuf.size()) &&
               cbuf[i].q_seg == qs && cbuf[i].r_seg == rs) {
            i++;
        }
        pair_info.push_back({start, i, qs, rs, cbuf[start].score});
    }

    struct RBHPair { uint32_t q_seg, r_seg; int chain_start, chain_end; };
    std::vector<RBHPair> rbh_pairs;

    if (p.skip_rbh) {
        for (const auto &pi : pair_info) {
            rbh_pairs.push_back({pi.q_seg, pi.r_seg, pi.start, pi.end});
        }
    } else {
        uint32_t max_qs = 0, max_rs = 0;
        for (const auto &pi : pair_info) {
            if (pi.q_seg > max_qs) max_qs = pi.q_seg;
            if (pi.r_seg > max_rs) max_rs = pi.r_seg;
        }
        bufs.best_q.assign(max_qs + 1, -1);
        bufs.best_r.assign(max_rs + 1, -1);

        for (int i = 0; i < static_cast<int>(pair_info.size()); i++) {
            const int bq = bufs.best_q[pair_info[i].q_seg];
            if (bq < 0 || pair_info[i].best_score > pair_info[bq].best_score)
                bufs.best_q[pair_info[i].q_seg] = i;
            const int br = bufs.best_r[pair_info[i].r_seg];
            if (br < 0 || pair_info[i].best_score > pair_info[br].best_score)
                bufs.best_r[pair_info[i].r_seg] = i;
        }

        for (int i = 0; i < static_cast<int>(pair_info.size()); i++) {
            if (bufs.best_q[pair_info[i].q_seg] == i &&
                bufs.best_r[pair_info[i].r_seg] == i) {
                rbh_pairs.push_back({pair_info[i].q_seg, pair_info[i].r_seg,
                                     pair_info[i].start, pair_info[i].end});
            }
        }
    }
    if (rbh_pairs.empty()) return res;

    uint16_t max_aa = 0;
    for (const auto &rp : rbh_pairs) {
        if (qsk.seg_headers[rp.q_seg].aa_len > max_aa)
            max_aa = qsk.seg_headers[rp.q_seg].aa_len;
        if (rsk.seg_headers[rp.r_seg].aa_len > max_aa)
            max_aa = rsk.seg_headers[rp.r_seg].aa_len;
    }
    bufs.qa_buf.resize(max_aa);
    bufs.ra_buf.resize(max_aa);

    bufs.qi_ivs.clear();
    bufs.ri_ivs.clear();
    int tm = 0, ta = 0;
    int n_matched = 0;
    struct OrfEntry { double id; int weight; };
    std::vector<OrfEntry> orf_entries;

    struct SegSpan { uint32_t seg; int aa_start, aa_end; };
    std::vector<SegSpan> q_spans, r_spans;

    uint32_t last_q = UINT32_MAX, last_r = UINT32_MAX;

    auto &mset = bufs.match_bits;
    auto &aset = bufs.align_bits;
    auto &pend = bufs.pending;
    auto &wins = bufs.win_buf;

    for (const auto &rp : rbh_pairs) {
        const auto &qsh = qsk.seg_headers[rp.q_seg];
        const auto &rsh = rsk.seg_headers[rp.r_seg];
        const int ql = static_cast<int>(qsh.aa_len);
        const int rl = static_cast<int>(rsh.aa_len);
        if (ql <= 0 || rl <= 0) continue;

        const bool track_query = (ql <= rl);
        const int  short_len   = track_query ? ql : rl;
        const size_t nw = static_cast<size_t>((short_len + 63) / 64);
        if (mset.size() < nw) { mset.resize(nw); aset.resize(nw); }
        std::memset(mset.data(), 0, nw * sizeof(uint64_t));
        std::memset(aset.data(), 0, nw * sizeof(uint64_t));

        if (rp.q_seg != last_q) { qsk.unpack_segment(rp.q_seg, bufs.qa_buf.data()); last_q = rp.q_seg; }
        if (rp.r_seg != last_r) { rsk.unpack_segment(rp.r_seg, bufs.ra_buf.data()); last_r = rp.r_seg; }
        const uint8_t *qa = bufs.qa_buf.data();
        const uint8_t *ra = bufs.ra_buf.data();

        int total_seeds = 0;
        bool pair_valid = false;
        int64_t pair_q_nuc_s = 0, pair_q_nuc_e = 0;
        int64_t pair_r_nuc_s = 0, pair_r_nuc_e = 0;
        pend.clear();

        auto commit = [&](const SWSeg &g) {
            for (int b = g.buf_start; b < g.buf_end; b++) {
                const uint32_t v = bufs.sw.pos[b];
                const uint32_t pos = v >> 1;
                aset[pos >> 6] |= (1ull << (pos & 63));
                if (v & 1u) mset[pos >> 6] |= (1ull << (pos & 63));
            }
            const int64_t q_nuc_s = static_cast<int64_t>(qsh.nuc_start) + static_cast<int64_t>(g.q_beg) * 3;
            const int64_t q_nuc_e = static_cast<int64_t>(qsh.nuc_start) + static_cast<int64_t>(g.q_end) * 3;
            const int64_t r_nuc_s = static_cast<int64_t>(rsh.nuc_start) + static_cast<int64_t>(g.r_beg) * 3;
            const int64_t r_nuc_e = static_cast<int64_t>(rsh.nuc_start) + static_cast<int64_t>(g.r_end) * 3;
            bufs.qi_ivs.push_back({q_nuc_s, q_nuc_e});
            bufs.ri_ivs.push_back({r_nuc_s, r_nuc_e});
            q_spans.push_back({rp.q_seg, g.q_beg, g.q_end});
            r_spans.push_back({rp.r_seg, g.r_beg, g.r_end});
            if (!pair_valid) {
                pair_q_nuc_s = q_nuc_s; pair_q_nuc_e = q_nuc_e;
                pair_r_nuc_s = r_nuc_s; pair_r_nuc_e = r_nuc_e;
            } else {
                if (q_nuc_s < pair_q_nuc_s) pair_q_nuc_s = q_nuc_s;
                if (q_nuc_e > pair_q_nuc_e) pair_q_nuc_e = q_nuc_e;
                if (r_nuc_s < pair_r_nuc_s) pair_r_nuc_s = r_nuc_s;
                if (r_nuc_e > pair_r_nuc_e) pair_r_nuc_e = r_nuc_e;
            }
            pair_valid = true;
        };

        for (int ci = rp.chain_start; ci < rp.chain_end; ci++) {
            const auto &anch = cbuf[ci].anchors;
            const int na = static_cast<int>(anch.size());
            if (na == 0) continue;
            total_seeds += cbuf[ci].n_matches;

            wins.clear();
            int seg_start = 0;
            while (seg_start < na) {
                int dmin = static_cast<int>(anch[seg_start].q_pos) - static_cast<int>(anch[seg_start].r_pos);
                int dmax = dmin;
                int seg_end = seg_start + 1;
                while (seg_end < na) {
                    const int dg = static_cast<int>(anch[seg_end].q_pos) - static_cast<int>(anch[seg_end].r_pos);
                    const int lo = std::min(dmin, dg), hi = std::max(dmax, dg);
                    if (hi - lo > SW_BAND - 3) break;
                    dmin = lo; dmax = hi;
                    seg_end++;
                }

                auto &diags = bufs.diag_buf;
                diags.clear();
                int q_lo = ql, q_hi = 0;
                for (int i = seg_start; i < seg_end; i++) {
                    diags.push_back(static_cast<int>(anch[i].q_pos) - static_cast<int>(anch[i].r_pos));
                    const int qp = static_cast<int>(anch[i].q_pos);
                    if (qp < q_lo) q_lo = qp;
                    if (qp > q_hi) q_hi = qp;
                }
                const int cn = static_cast<int>(diags.size());
                std::nth_element(diags.begin(), diags.begin() + cn / 2, diags.end());
                const int d = diags[cn / 2];
                int max_dev = 0;
                for (int i = 0; i < cn; i++) {
                    const int dev = std::abs(diags[i] - d);
                    if (dev > max_dev) max_dev = dev;
                }
                assert(max_dev + 3 <= SW_BAND);
                const int aband = std::min(std::max(3, max_dev + 3), SW_BAND);

                q_hi = std::min(q_hi + p.shape_span - 1, ql);
                q_lo = std::max(q_lo, d);
                q_hi = std::min(q_hi, d + rl);
                const int da = static_cast<int>(anch[seg_start].q_pos) - static_cast<int>(anch[seg_start].r_pos);
                const int dz = static_cast<int>(anch[seg_end - 1].q_pos) - static_cast<int>(anch[seg_end - 1].r_pos);
                if (q_lo > 0 && (q_lo - da) > 0) {
                    int e = xdrop_extend_backward(qa, ra, q_lo, q_lo - da, p.xdrop);
                    e += std::max(0, da - d);
                    e = std::min(e, std::min(q_lo - d, q_lo));
                    q_lo -= e;
                }
                if (q_hi < ql && (q_hi - dz) < rl) {
                    int e = xdrop_extend_forward(qa, ra, q_hi, q_hi - dz, ql, rl, p.xdrop);
                    e += std::max(0, d - dz);
                    e = std::min(e, std::min(d + rl - q_hi, ql - q_hi));
                    q_hi += e;
                }

                wins.push_back({d, aband, q_lo, q_hi});
                seg_start = seg_end;
            }

            for (size_t t = 1; t < wins.size(); t++) {
                if (wins[t - 1].q_hi > wins[t].q_lo) {
                    const int mid = (wins[t - 1].q_hi + wins[t].q_lo) / 2;
                    wins[t - 1].q_hi = mid;
                    wins[t].q_lo = mid;
                }
            }

            for (const auto &wdw : wins) {
                const int d = wdw.d;
                const int aband = wdw.aband;
                const int q_lo = wdw.q_lo, q_hi = wdw.q_hi;
                const int r_lo = std::max(0, q_lo - d);
                const int r_hi = std::min(rl, q_hi - d);
                if (q_hi <= q_lo || r_hi <= r_lo) continue;
                if ((q_hi - q_lo) < p.min_chain_aligned_aa) continue;

                if (p.aggressive_filter && (q_hi - q_lo) >= 10) {
                    const uint8_t *sides[2] = {qa, ra};
                    const int starts[2] = {q_lo, r_lo};
                    const int ends[2] = {q_hi, r_hi};
                    bool reject = false;
                    for (int side = 0; side < 2; side++) {
                        const int len = ends[side] - starts[side];
                        if (len < 10) continue;
                        uint16_t freq[N_AMINO_ACIDS] = {};
                        for (int i = starts[side]; i < ends[side]; i++) {
                            const uint8_t aa = sides[side][i];
                            if (aa < N_AMINO_ACIDS) freq[aa]++;
                        }
                        uint16_t mx1 = 0, mx2 = 0, mx3 = 0;
                        double entropy = 0;
                        for (int a = 0; a < N_AMINO_ACIDS; a++) {
                            if (freq[a] > mx1)      { mx3 = mx2; mx2 = mx1; mx1 = freq[a]; }
                            else if (freq[a] > mx2) { mx3 = mx2; mx2 = freq[a]; }
                            else if (freq[a] > mx3) { mx3 = freq[a]; }
                            if (freq[a] > 0) {
                                const double pr = static_cast<double>(freq[a]) / len;
                                entropy -= pr * std::log2(pr);
                            }
                        }
                        if (static_cast<float>(mx1) / len > 0.30f ||
                            static_cast<float>(mx1 + mx2) / len > 0.50f ||
                            static_cast<float>(mx1 + mx2 + mx3) / len > 0.65f ||
                            entropy < 2.5) { reject = true; break; }
                    }
                    if (reject) continue;
                }

                banded_sw_segments(qa, q_lo, q_hi, ra, r_lo, r_hi, aband, track_query, bufs.sw);

                for (const auto &g : bufs.sw.segs) {
                    if (compute_evalue(g.score, ql, rl) <= p.max_evalue) commit(g);
                    else pend.push_back(g);
                }
                if (!pend.empty()) {
                    bufs.pend_scores.clear();
                    for (const auto &g : pend) bufs.pend_scores.push_back(g.score);
                    if (compute_evalue_sum(bufs.pend_scores.data(),
                                           static_cast<int>(bufs.pend_scores.size()), ql, rl) <= p.max_evalue) {
                        for (const auto &g : pend) commit(g);
                    }
                    pend.clear();
                }
            }
        }

        if (!pair_valid) continue;

        int pair_matches = 0, pair_aligned = 0;
        for (size_t w = 0; w < nw; w++) {
            pair_matches += __builtin_popcountll(mset[w]);
            pair_aligned += __builtin_popcountll(aset[w]);
        }
        if (pair_matches == 0) continue;

        tm += pair_matches;
        ta += pair_aligned;

        const int denom = std::min(ql, rl);
        if (denom > 0) {
            assert(pair_matches <= denom);
            const double oid = static_cast<double>(pair_matches) / denom;
            orf_entries.push_back({oid, denom});
            n_matched++;
        }
        if (out_orf_hits) {
            OrfPairHit hit;
            hit.q_seg = rp.q_seg;
            hit.r_seg = rp.r_seg;
            hit.q_frame = qsh.frame;
            hit.r_frame = rsh.frame;
            hit.q_nuc_start = pair_q_nuc_s;
            hit.q_nuc_end = pair_q_nuc_e;
            hit.r_nuc_start = pair_r_nuc_s;
            hit.r_nuc_end = pair_r_nuc_e;
            hit.n_seeds = total_seeds;
            hit.identity = static_cast<double>(pair_matches) / denom;
            out_orf_hits->push_back(hit);
        }
    }

    res.total_aa_matches = tm;
    res.total_aligned = ta;
    if (n_matched > 0) {
        double wsum = 0, wtotal = 0;
        for (const auto &e : orf_entries) {
            wsum += e.id * e.weight;
            wtotal += e.weight;
        }
        if (wtotal > 0) res.aai = wsum / wtotal;
    }

    const double gs2 = static_cast<double>(p.gene_scaled) * p.gene_scaled;
    const double n_matched_est = n_matched * gs2;
    res.matched_orfs = static_cast<int>(std::lround(n_matched_est));

    const uint64_t qc = merge_intervals(bufs.qi_ivs);
    const uint64_t rc = merge_intervals(bufs.ri_ivs);
    res.gcov_q = qsk.nuc_length > 0 ?
        std::min(static_cast<double>(qc) * gs2 / qsk.nuc_length, 1.0) : 0;
    res.gcov_r = rsk.nuc_length > 0 ?
        std::min(static_cast<double>(rc) * gs2 / rsk.nuc_length, 1.0) : 0;
    const double ccov_q = qsk.coding_nuc > 0 ?
        std::min(static_cast<double>(qc) * gs2 / qsk.coding_nuc, 1.0) : 0;
    const double ccov_r = rsk.coding_nuc > 0 ?
        std::min(static_cast<double>(rc) * gs2 / rsk.coding_nuc, 1.0) : 0;

    auto weighted_orf_cov = [](std::vector<SegSpan> &spans,
                               const Arr<SegHeader> &headers) -> double {
        if (spans.empty()) return 0;
        std::sort(spans.begin(), spans.end(),
                  [](const SegSpan &a, const SegSpan &b) { return a.seg < b.seg; });
        int total_aligned = 0, total_orf_len = 0;
        size_t i = 0;
        while (i < spans.size()) {
            const uint32_t seg = spans[i].seg;
            const int seg_len = static_cast<int>(headers[seg].aa_len);
            if (seg_len <= 0) { i++; continue; }
            std::vector<std::pair<int, int>> ivs;
            while (i < spans.size() && spans[i].seg == seg) {
                ivs.push_back({spans[i].aa_start, spans[i].aa_end});
                i++;
            }
            std::sort(ivs.begin(), ivs.end());
            int merged = 0;
            int ce = ivs[0].first;
            for (auto &iv : ivs) {
                const int s = std::max(iv.first, ce);
                if (iv.second > s) {
                    merged += iv.second - s;
                    ce = iv.second;
                }
            }
            total_aligned += std::min(merged, seg_len);
            total_orf_len += seg_len;
        }
        return total_orf_len > 0 ? static_cast<double>(total_aligned) / total_orf_len : 0;
    };
    res.orf_depth_q = weighted_orf_cov(q_spans, qsk.seg_headers);
    res.orf_depth_r = weighted_orf_cov(r_spans, rsk.seg_headers);

    res.gcov = std::max(res.gcov_q, res.gcov_r);
    res.orf_depth = std::max(res.orf_depth_q, res.orf_depth_r);

    const double q_n_orfs = static_cast<double>(
        qsk.n_orfs_full ? qsk.n_orfs_full : qsk.seg_headers.size());
    const double r_n_orfs = static_cast<double>(
        rsk.n_orfs_full ? rsk.n_orfs_full : rsk.seg_headers.size());
    const double orf_hit_q = q_n_orfs > 0 ?
        std::min(n_matched_est / q_n_orfs, 1.0) : 0;
    const double orf_hit_r = r_n_orfs > 0 ?
        std::min(n_matched_est / r_n_orfs, 1.0) : 0;
    res.orf_hit = std::max(orf_hit_q, orf_hit_r);

    res.composite_score = std::sqrt(std::max(ccov_q, ccov_r) * res.aai);

    return res;
}
