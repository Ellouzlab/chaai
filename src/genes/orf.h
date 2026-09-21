#pragma once
#include "db/types.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include <vector>
#include <algorithm>
#include <cstdint>

static constexpr int    CHAAI_RBS_MIN     = 5;
static constexpr int    CHAAI_RBS_SCAN    = 20;
static constexpr int    CHAAI_TRIM_START  = 30;
static constexpr uint16_t MAX_SEG_LEN = 65535;

static inline void trim_to_start(AASegment &seg, const char *seq, uint64_t coff, size_t clen, int scan,
                                 int min_len)
{
    const bool reverse = seg.frame >= 3;
    const int64_t five_prime = reverse
        ? (static_cast<int64_t>(seg.nuc_start - coff) + static_cast<int64_t>(seg.aa_len) * 3)
        : static_cast<int64_t>(seg.nuc_start - coff);
    auto rbase = [&](int64_t off) -> int {
        if (!reverse) { const int64_t p = five_prime + off;
            return (p >= 0 && p < static_cast<int64_t>(clen)) ? static_cast<uint8_t>(seq[p]) : 4; }
        const int64_t p = five_prime - 1 - off;
        if (p < 0 || p >= static_cast<int64_t>(clen)) return 4;
        return ENC_COMP[static_cast<uint8_t>(seq[p])];
    };
    const int lim = std::min(scan, static_cast<int>(seg.aa_len) - min_len);
    for (int c = 0; c < lim; c++) {
        const int b0 = rbase(c*3), b1 = rbase(c*3+1), b2 = rbase(c*3+2);
        if ((b1 == 3 && b2 == 2) && (b0 == 0 || b0 == 2 || b0 == 3)) {
            if (c > 0) { if (!reverse) seg.nuc_start += static_cast<uint64_t>(c) * 3;
                         seg.aa_len -= static_cast<uint16_t>(c); }
            return;
        }
    }
}

static inline int orf_rbs_score(const char *seq, size_t clen, int64_t five_prime, bool reverse)
{
    auto rbase = [&](int64_t off) -> int {
        if (!reverse) {
            const int64_t p = five_prime + off;
            return (p >= 0 && p < static_cast<int64_t>(clen)) ? static_cast<uint8_t>(seq[p]) : 4;
        }
        const int64_t p = five_prime - 1 - off;
        if (p < 0 || p >= static_cast<int64_t>(clen)) return 4;
        return ENC_COMP[static_cast<uint8_t>(seq[p])];
    };
    static const int SD[6] = {0, 2, 2, 0, 2, 2};
    int overall = 0;
    for (int c = 0; c < CHAAI_RBS_SCAN; c++) {
        const int b0 = rbase(c * 3), b1 = rbase(c * 3 + 1), b2 = rbase(c * 3 + 2);
        const bool is_start = (b1 == 3 && b2 == 2) && (b0 == 0 || b0 == 2 || b0 == 3);
        if (!is_start) continue;
        const int start_off = c * 3;
        int best = 0;
        {
            int up[16];
            for (int d = 1; d <= 15; d++) up[d] = rbase(start_off - d);
            for (int D = 6; D <= 14; D++) {
                int run = 0, bestrun = 0;
                for (int k = 0; k < 6; k++) {
                    if (up[D - k] == SD[k]) { if (++run > bestrun) bestrun = run; } else run = 0;
                }
                const int spacing = D - 5;
                const int score = bestrun + ((spacing >= 4 && spacing <= 10) ? 1 : 0);
                if (score > best) best = score;
            }
        }
        if (best > overall) overall = best;
    }
    return overall;
}

static inline bool orf_has_rbs(const char *seq, size_t clen, int64_t five_prime, bool reverse)
{
    return orf_rbs_score(seq, clen, five_prime, reverse) >= CHAAI_RBS_MIN;
}

static constexpr uint64_t FWD_STOPS = (1ULL << 48) | (1ULL << 50) | (1ULL << 56);
static constexpr uint64_t REV_STOPS = (1ULL << 60) | (1ULL << 28) | (1ULL << 52);

static void scan_window_set(
    const char *seq, size_t seq_len, uint64_t coff, uint16_t min_seg_len, int g,
    std::vector<AASegment> &fwd, std::vector<AASegment> &rev, uint64_t *hist = nullptr)
{
    if (seq_len < static_cast<size_t>(g) + 3) return;
    const uint64_t clen = seq_len;
    const int f = static_cast<int>((seq_len - g) % 3);
    const size_t nk = (seq_len - g) / 3;

    uint16_t fl = 0, rl = 0;
    uint64_t fs0 = 0;
    bool fopen = false, ropen = false;

    auto emit_fwd = [&](size_t k, bool by_stop) {
        if (hist && by_stop && fopen) hist[std::min<int>(fl, RUN_HIST_MAX)]++;
        if (fl >= min_seg_len) {
            AASegment seg;
            seg.frame = static_cast<uint8_t>(g);
            seg.aa_len = fl;
            seg.nuc_start = coff + g + fs0 * 3;
            fwd.push_back(seg);
        }
        fl = 0;
        fs0 = k + 1;
        fopen = by_stop;
    };
    auto emit_rev = [&](size_t i, bool by_stop) {
        if (hist && by_stop && ropen) hist[std::min<int>(rl, RUN_HIST_MAX)]++;
        if (rl >= min_seg_len) {
            const uint64_t s0 = i + 1;
            AASegment seg;
            seg.frame = static_cast<uint8_t>(3 + f);
            seg.aa_len = rl;
            seg.nuc_start = coff + clen - f - (s0 + rl) * 3;
            rev.push_back(seg);
        }
        rl = 0;
        ropen = by_stop;
    };

    const char *w = seq + g;
    for (size_t k = 0; k < nk; k++, w += 3) {
        const uint8_t n0 = static_cast<uint8_t>(w[0]);
        const uint8_t n1 = static_cast<uint8_t>(w[1]);
        const uint8_t n2 = static_cast<uint8_t>(w[2]);
        const size_t i = nk - 1 - k;
        if ((n0 | n1 | n2) > 3) {
            emit_fwd(k, false);
            emit_rev(i, false);
            continue;
        }
        const unsigned idx = n0 * 16u + n1 * 4u + n2;
        if ((FWD_STOPS >> idx) & 1) emit_fwd(k, true);
        else if (++fl == MAX_SEG_LEN) emit_fwd(k, false);
        if ((REV_STOPS >> idx) & 1) emit_rev(i, true);
        else if (++rl == MAX_SEG_LEN) emit_rev(i - 1, false);
    }
    emit_fwd(nk, false);
    emit_rev(static_cast<size_t>(-1), false);
    std::reverse(rev.begin(), rev.end());
}

static inline size_t fwd_close_idx(const AASegment &s, uint64_t coff, int f)
{
    return (s.nuc_start - coff - f) / 3 + s.aa_len;
}
static inline size_t rev_close_idx(const AASegment &s, uint64_t coff, size_t clen, int f)
{
    return (clen - f - (s.nuc_start - coff)) / 3;
}

static void merge_frame(const std::vector<AASegment> &fwd, const std::vector<AASegment> &rev,
                        uint64_t coff, size_t clen, int f, AASegment *out)
{
    size_t a = 0, b = 0, o = 0;
    while (a < fwd.size() && b < rev.size()) {
        if (fwd_close_idx(fwd[a], coff, f) <= rev_close_idx(rev[b], coff, clen, f))
            out[o++] = fwd[a++];
        else
            out[o++] = rev[b++];
    }
    while (a < fwd.size()) out[o++] = fwd[a++];
    while (b < rev.size()) out[o++] = rev[b++];
}

static void pack_segment_from_seq(
    const AASegment &seg, const char *seq,
    uint64_t contig_offset, uint8_t *out)
{
    uint64_t buf = 0;
    int bits = 0;
    uint8_t *dst = out;
    const uint64_t local = seg.nuc_start - contig_offset;

    for (uint16_t i = 0; i < seg.aa_len; i++) {
        uint8_t aa;
        int e0, e1, e2;
        if (seg.frame < 3) {
            const size_t pos = local + static_cast<size_t>(i) * 3;
            e0 = static_cast<uint8_t>(seq[pos]); e1 = static_cast<uint8_t>(seq[pos + 1]); e2 = static_cast<uint8_t>(seq[pos + 2]);
        } else {
            const size_t pos = local + static_cast<size_t>(seg.aa_len - 1 - i) * 3;
            e0 = ENC_COMP[static_cast<uint8_t>(seq[pos + 2])]; e1 = ENC_COMP[static_cast<uint8_t>(seq[pos + 1])]; e2 = ENC_COMP[static_cast<uint8_t>(seq[pos])];
        }
        aa = CODON_TABLE[e0 * 16 + e1 * 4 + e2];
        buf |= (static_cast<uint64_t>(aa & 0x1F)) << bits;
        bits += 5;
        while (bits >= 8) {
            *dst++ = static_cast<uint8_t>(buf & 0xFF);
            buf >>= 8;
            bits -= 8;
        }
    }
    if (bits > 0) *dst = static_cast<uint8_t>(buf & 0xFF);
}

static std::vector<uint8_t> translate_and_pack(
    const std::vector<AASegment> &segments, const ContigView &gv)
{
    const int n = static_cast<int>(segments.size());
    std::vector<size_t> offsets(n + 1);
    offsets[0] = 0;
    for (int i = 0; i < n; i++) {
        offsets[i + 1] = offsets[i] + packed_aa_bytes(segments[i].aa_len);
    }
    std::vector<uint8_t> packed(offsets[n], 0);

#pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < n; i++) {
        const int ci = gv.find_contig(segments[i].nuc_start);
        pack_segment_from_seq(segments[i], gv.seqs[ci], gv.offsets[ci],
                              packed.data() + offsets[i]);
    }
    return packed;
}

static constexpr char AA_LETTERS[] = "ARNDCQEGHILKMFPSTWYV*";
