#pragma once

#include "db/types.h"
#include "aai/hash.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>

static constexpr int AA_ALPHA = 20;
static constexpr int MAX_SPACED_W = 12;

struct SpacedPattern {
    int weight;
    int span;
    int offsets[MAX_SPACED_W];
};
static_assert(MAX_SPACED_W <= 13, "20^13 < 2^56, must fit in uint64_t hash input");

static constexpr int MAX_PATSET = 4;
struct PatSet { int n; SpacedPattern p[MAX_PATSET]; };

static inline int patset_weight(const PatSet &ps) { return ps.n > 0 ? ps.p[0].weight : 0; }

static inline int patset_span(const PatSet &ps)
{
    int t = 0;
    for (int i = 0; i < ps.n; i++) if (ps.p[i].span > t) t = ps.p[i].span;
    return t;
}

static inline PatSet parse_seedspec(const char *spec)
{
    PatSet ps{}; ps.n = 0;
    if (!spec) return ps;
    const char *colon = std::strchr(spec, ':');
    if (!colon) return ps;
    const int w = std::atoi(spec);
    const char *cur = colon + 1;
    while (cur && *cur && ps.n < MAX_PATSET) {
        SpacedPattern sp{}; sp.weight = w;
        int k = 0, mx = 0;
        const char *e = cur;
        while (*e && *e != ';') {
            const int v = std::atoi(e);
            if (k < MAX_SPACED_W) { sp.offsets[k++] = v; if (v > mx) mx = v; }
            while (*e && *e != ',' && *e != ';') e++;
            if (*e == ',') e++;
        }
        sp.span = mx + 1;
        if (k == w) ps.p[ps.n++] = sp;
        cur = (*e == ';') ? e + 1 : nullptr;
    }
    return ps;
}

static inline uint16_t quant_nlogp(double log_e)
{
    const double v = -log_e * 256.0;
    if (v <= 0) return 0;
    return v >= 65535.0 ? 65535 : static_cast<uint16_t>(v + 0.5);
}

static inline void extract_pat(const uint8_t *aa, int len, uint32_t si,
                               const SpacedPattern &pat, uint64_t thr,
                               std::vector<Seed> &out,
                               const double *logp = nullptr, double log_npos = 0.0,
                               uint64_t thr32 = uint64_t{1} << 32)
{
    for (int q = 0; q + pat.span <= len; q++) {
        uint64_t kv = 0;
        double lp = 0.0;
        for (int i = 0; i < pat.weight; i++) {
            const uint8_t a = aa[q + pat.offsets[i]];
            kv = kv * AA_ALPHA + a;
            if (logp) lp += logp[a < AA_ALPHA ? a : 0];
        }
        const uint64_t h = hash64(kv);
        if (h <= thr && (h & 0xffffffffu) < thr32) {
            Seed sd{static_cast<uint32_t>(h), si, static_cast<uint16_t>(q), 0};
            if (logp) sd._pad = quant_nlogp(log_npos + lp);
            out.push_back(sd);
        }
    }
}

static inline void sort_and_mask_seeds(std::vector<Seed> &seeds, int max_occ)
{
    std::sort(seeds.begin(), seeds.end(),
              [](const Seed &a, const Seed &b) {
                  if (a.hash32 != b.hash32) return a.hash32 < b.hash32;
                  if (a.seg_idx != b.seg_idx) return a.seg_idx < b.seg_idx;
                  return a.aa_pos < b.aa_pos;
              });

    std::vector<Seed> filtered;
    filtered.reserve(seeds.size());
    size_t i = 0;
    while (i < seeds.size()) {
        size_t j = i;
        while (j < seeds.size() && seeds[j].hash32 == seeds[i].hash32) j++;
        int cut = max_occ;
        if (cut <= 0) {
            const double e = std::exp(-static_cast<double>(seeds[i]._pad) / 256.0);
            cut = std::max(2, static_cast<int>(e + 4.0 * std::sqrt(e) + 0.5));
        }
        if (static_cast<int>(j - i) <= cut) {
            for (size_t k = i; k < j; k++) {
                filtered.push_back(seeds[k]);
            }
        }
        i = j;
    }
    seeds = std::move(filtered);
}

static inline std::vector<Seed> extract_seeds_spec_from_db(
    const GenomeDb &sk, const PatSet &ps, uint32_t seed_scaled, int max_occ = 0,
    uint64_t thr32 = uint64_t{1} << 32)
{
    std::vector<Seed> seeds;
    const uint64_t thr = fmh_threshold(seed_scaled);
    std::vector<uint8_t> buf;

    double counts[AA_ALPHA] = {}, total = 0, npos = 0;
    for (uint32_t si = 0; si < static_cast<uint32_t>(sk.seg_headers.size()); si++) {
        const uint16_t len = sk.seg_headers[si].aa_len;
        buf.resize(len);
        sk.unpack_segment(si, buf.data());
        for (uint16_t k = 0; k < len; k++)
            if (buf[k] < AA_ALPHA) { counts[buf[k]] += 1.0; total += 1.0; }
        for (int i = 0; i < ps.n; i++)
            if (len >= ps.p[i].span) npos += len - ps.p[i].span + 1;
    }
    double logp[AA_ALPHA];
    for (int a = 0; a < AA_ALPHA; a++)
        logp[a] = std::log((counts[a] + 1.0) / (total + AA_ALPHA));
    const double log_npos = std::log(npos > 0 ? npos : 1.0);

    for (uint32_t si = 0; si < static_cast<uint32_t>(sk.seg_headers.size()); si++) {
        const uint16_t len = sk.seg_headers[si].aa_len;
        buf.resize(len);
        sk.unpack_segment(si, buf.data());
        for (int i = 0; i < ps.n; i++) {
            if (len < ps.p[i].span) continue;
            extract_pat(buf.data(), len, si, ps.p[i], thr, seeds, logp, log_npos, thr32);
        }
    }
    sort_and_mask_seeds(seeds, max_occ);
    return seeds;
}
