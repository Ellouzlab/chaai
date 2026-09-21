#pragma once
#include "db/types.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include <vector>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
#include <omp.h>

template <class F>
static inline void for_each_codon(const AASegment &seg, const char *cseq, uint64_t coff, F &&f)
{
    const uint16_t len = seg.aa_len;
    if (len == 0) return;
    const size_t base = seg.nuc_start - coff;
    if (seg.frame < 3) {
        const char *p = cseq + base;
        for (uint16_t i = 0; i < len; i++, p += 3) {
            const uint8_t n0 = static_cast<uint8_t>(p[0]), n1 = static_cast<uint8_t>(p[1]), n2 = static_cast<uint8_t>(p[2]);
            f(((n0 | n1 | n2) > 3) ? -1 : (n0 * 16 + n1 * 4 + n2));
        }
    } else {
        const char *p = cseq + base + static_cast<size_t>(len - 1) * 3;
        for (uint16_t i = 0; i < len; i++, p -= 3) {
            const uint8_t r0 = static_cast<uint8_t>(p[0]), r1 = static_cast<uint8_t>(p[1]), r2 = static_cast<uint8_t>(p[2]);
            f(((r0 | r1 | r2) > 3) ? -1 : ((3 - r2) * 16 + (3 - r1) * 4 + (3 - r0)));
        }
    }
}

static inline bool count_positions(const char *s, int nlen, uint32_t cq[3][5])
{
    for (int q = 0; q < 3; q++) for (int b = 0; b < 5; b++) cq[q][b] = 0;
    int n = 0;
    bool bad = false;
#if defined(__AVX2__)
    static const uint32_t PM[3][3] = {
        {0x49249249u, 0x92492492u, 0x24924924u},
        {0x24924924u, 0x49249249u, 0x92492492u},
        {0x92492492u, 0x24924924u, 0x49249249u}};
    const __m256i v0 = _mm256_setzero_si256(), v1 = _mm256_set1_epi8(1), v2 = _mm256_set1_epi8(2), v3 = _mm256_set1_epi8(3);
    uint32_t badm = 0;
    int ph = 0;
    for (; n + 32 <= nlen; n += 32) {
        const __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(s + n));
        const uint32_t m0 = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v, v0)));
        const uint32_t m1 = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v, v1)));
        const uint32_t m2 = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v, v2)));
        const uint32_t m3 = static_cast<uint32_t>(_mm256_movemask_epi8(_mm256_cmpeq_epi8(v, v3)));
        badm |= ~(m0 | m1 | m2 | m3);
        const uint32_t *pm = PM[ph];
        cq[0][0] += __builtin_popcount(m0 & pm[0]); cq[1][0] += __builtin_popcount(m0 & pm[1]); cq[2][0] += __builtin_popcount(m0 & pm[2]);
        cq[0][1] += __builtin_popcount(m1 & pm[0]); cq[1][1] += __builtin_popcount(m1 & pm[1]); cq[2][1] += __builtin_popcount(m1 & pm[2]);
        cq[0][2] += __builtin_popcount(m2 & pm[0]); cq[1][2] += __builtin_popcount(m2 & pm[1]); cq[2][2] += __builtin_popcount(m2 & pm[2]);
        cq[0][3] += __builtin_popcount(m3 & pm[0]); cq[1][3] += __builtin_popcount(m3 & pm[1]); cq[2][3] += __builtin_popcount(m3 & pm[2]);
        ph = (ph + 2) % 3;
    }
    bad = badm != 0;
#endif
    int q = n % 3;
    for (; n < nlen; n++) {
        const uint8_t b = static_cast<uint8_t>(s[n]);
        cq[q][b > 3 ? 4 : b]++;
        q = (q == 2) ? 0 : q + 1;
    }
    return bad || cq[0][4] || cq[1][4] || cq[2][4];
}

static inline double asym_from_counts(const uint32_t ct[4][3])
{
    uint32_t nvi = 0;
    for (int b = 0; b < 4; b++) nvi += ct[b][0] + ct[b][1] + ct[b][2];
    const double nv = nvi;
    if (nv < 9) return 0;

    double power = 0, het = 0;
    for (int b = 0; b < 4; b++) {
        const double s0 = ct[b][0], s1 = ct[b][1], s2 = ct[b][2];
        const double tot  = s0 + s1 + s2;
        const double mean = tot / 3.0;
        const double d0 = s0 - mean, d1 = s1 - mean, d2 = s2 - mean;
        power += 1.5 * (d0 * d0 + d1 * d1 + d2 * d2);
        const double pb = tot / nv;
        het += pb * (1.0 - pb);
    }
    const double expected = nv * het;
    return (expected > 1.0) ? power / expected : 0;
}

static inline double compute_asymmetry(
    const char *seq, uint32_t start, int nuc_len)
{
    if (nuc_len < 9) return 0;

    uint32_t cq[3][5];
    count_positions(seq + start, nuc_len, cq);
    uint32_t ct[4][3];
    for (int b = 0; b < 4; b++) { ct[b][0] = cq[0][b]; ct[b][1] = cq[1][b]; ct[b][2] = cq[2][b]; }
    return asym_from_counts(ct);
}

static constexpr double CHAAI_ASYM_CUT = 1.6;
static constexpr double CHAAI_ASYM_CUT_TRAIN = 2.0;
