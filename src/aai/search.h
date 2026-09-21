#pragma once

#include "db/types.h"
#include "db/io.h"
#include "db/seedcache.h"
#include "aai/seeds.h"
#include "aai/chain.h"
#include "aai/align.h"
#include "aai/pair.h"
#include "aai/pairwise.h"
#include "run/progress.h"
#include "run/interrupt.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

static inline uint64_t marker_threshold(uint32_t factor)
{
    return factor <= 1 ? (uint64_t{1} << 32) : (uint64_t{1} << 32) / factor;
}

static inline void marker_hashes(SeedSpan seeds, uint64_t thr, std::vector<uint32_t> &out)
{
    out.clear();
    bool first = true;
    uint32_t last = 0;
    for (const auto &s : seeds) {
        if (static_cast<uint64_t>(s.hash32) >= thr) continue;
        if (!first && s.hash32 == last) continue;
        out.push_back(s.hash32);
        last = s.hash32;
        first = false;
    }
}

struct KeptSeed { Seed s; int pat; uint8_t res[MAX_SPACED_W]; };

#if defined(__AVX2__)
static inline __m256i mul64_by(__m256i a, uint64_t b)
{
    const __m256i bl = _mm256_set1_epi64x(static_cast<int64_t>(b & 0xffffffffu));
    const __m256i bh = _mm256_set1_epi64x(static_cast<int64_t>(b >> 32));
    const __m256i lo = _mm256_mul_epu32(a, bl);
    const __m256i cross = _mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(a, 32), bl),
                                           _mm256_mul_epu32(a, bh));
    return _mm256_add_epi64(lo, _mm256_slli_epi64(cross, 32));
}

static inline __m256i hash64_x4(__m256i k)
{
    k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 33));
    k = mul64_by(k, 0xff51afd7ed558ccdULL);
    k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 33));
    k = mul64_by(k, 0xc4ceb9fe1a85ec53ULL);
    k = _mm256_xor_si256(k, _mm256_srli_epi64(k, 33));
    return k;
}
#endif

static inline void keep_seed(const uint8_t *aa, int q, uint32_t si, int pi, const SpacedPattern &pat,
                             uint64_t thr64, uint64_t thr32, std::vector<KeptSeed> &kept)
{
    uint64_t kv = 0;
    for (int i = 0; i < pat.weight; i++) kv = kv * AA_ALPHA + aa[q + pat.offsets[i]];
    const uint64_t h = hash64(kv);
    if (h > thr64 || (h & 0xffffffffu) >= thr32) return;
    KeptSeed k{};
    k.s = Seed{static_cast<uint32_t>(h), si, static_cast<uint16_t>(q), 0};
    k.pat = pi;
    for (int i = 0; i < pat.weight; i++) k.res[i] = aa[q + pat.offsets[i]];
    kept.push_back(k);
}

static inline void extract_markers(const GenomeDb &g, const PatSet &ps, uint32_t seed_scaled,
                                   int max_occ, uint64_t thr32, std::vector<uint32_t> &out)
{
    const uint64_t thr64 = fmh_threshold(seed_scaled);
    std::vector<uint8_t> buf;
    std::vector<KeptSeed> kept;
    double counts[AA_ALPHA] = {}, total = 0, npos = 0;
#if defined(__AVX2__)
    const __m256i sign = _mm256_set1_epi64x(static_cast<int64_t>(0x8000000000000000ULL));
    const __m256i thr64v = _mm256_xor_si256(_mm256_set1_epi64x(static_cast<int64_t>(thr64)), sign);
    const __m256i thr32v = _mm256_set1_epi64x(static_cast<int64_t>(thr32));
    const __m256i low32 = _mm256_set1_epi64x(0xffffffffLL);
#endif
    for (uint32_t si = 0; si < static_cast<uint32_t>(g.seg_headers.size()); si++) {
        const int len = g.seg_headers[si].aa_len;
        buf.resize(len);
        g.unpack_segment(si, buf.data());
        const uint8_t *aa = buf.data();
        for (int k = 0; k < len; k++)
            if (aa[k] < AA_ALPHA) { counts[aa[k]] += 1.0; total += 1.0; }
        for (int pi = 0; pi < ps.n; pi++) {
            const SpacedPattern &pat = ps.p[pi];
            if (len < pat.span) continue;
            npos += len - pat.span + 1;
            const int last = len - pat.span;
            int q = 0;
#if defined(__AVX2__)
            for (; q + 3 <= last; q += 4) {
                __m256i kv = _mm256_setzero_si256();
                for (int i = 0; i < pat.weight; i++) {
                    int32_t four;
                    std::memcpy(&four, aa + q + pat.offsets[i], 4);
                    const __m256i a = _mm256_cvtepu8_epi64(_mm_cvtsi32_si128(four));
                    kv = _mm256_add_epi64(_mm256_add_epi64(_mm256_slli_epi64(kv, 4),
                                                           _mm256_slli_epi64(kv, 2)), a);
                }
                const __m256i h = hash64_x4(kv);
                const __m256i over = _mm256_cmpgt_epi64(_mm256_xor_si256(h, sign), thr64v);
                const __m256i under = _mm256_cmpgt_epi64(thr32v, _mm256_and_si256(h, low32));
                int m = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_andnot_si256(over, under)));
                while (m) {
                    const int lane = __builtin_ctz(m);
                    m &= m - 1;
                    keep_seed(aa, q + lane, si, pi, pat, thr64, thr32, kept);
                }
            }
#endif
            for (; q <= last; q++) keep_seed(aa, q, si, pi, pat, thr64, thr32, kept);
        }
    }
    double logp[AA_ALPHA];
    for (int a = 0; a < AA_ALPHA; a++)
        logp[a] = std::log((counts[a] + 1.0) / (total + AA_ALPHA));
    const double log_npos = std::log(npos > 0 ? npos : 1.0);
    std::vector<Seed> seeds;
    seeds.reserve(kept.size());
    for (const auto &k : kept) {
        double lp = 0.0;
        for (int i = 0; i < ps.p[k.pat].weight; i++) lp += logp[k.res[i] < AA_ALPHA ? k.res[i] : 0];
        Seed sd = k.s;
        sd._pad = quant_nlogp(log_npos + lp);
        seeds.push_back(sd);
    }
    sort_and_mask_seeds(seeds, max_occ);
    marker_hashes(SeedSpan(seeds), thr32, out);
}

struct MarkerIndex {
    std::vector<uint32_t> hashes;
    std::vector<uint32_t> offsets;
    std::vector<uint32_t> postings;
    std::vector<uint32_t> bucket;

    void build(const std::vector<std::vector<uint32_t>> &markers, int threads)
    {
        static constexpr size_t NB = 65536;
        const size_t nr = markers.size();
        std::vector<uint64_t> cnt(NB + 1, 0);
        for (const auto &m : markers)
            for (uint32_t h : m) cnt[(h >> 16) + 1]++;
        for (size_t b = 0; b < NB; b++) cnt[b + 1] += cnt[b];
        const uint64_t total = cnt[NB];
        std::vector<uint64_t> ent(total);
        {
            std::vector<uint64_t> pos(cnt.begin(), cnt.end() - 1);
            for (size_t r = 0; r < nr; r++)
                for (uint32_t h : markers[r])
                    ent[pos[h >> 16]++] = (static_cast<uint64_t>(h) << 32) | r;
        }
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
        for (int64_t b = 0; b < static_cast<int64_t>(NB); b++)
            std::sort(ent.begin() + cnt[b], ent.begin() + cnt[b + 1]);

        hashes.clear(); offsets.clear(); postings.clear();
        postings.resize(total);
        bucket.assign(NB + 1, 0);
        for (uint64_t i = 0; i < total; i++) {
            const uint32_t h = static_cast<uint32_t>(ent[i] >> 32);
            if (hashes.empty() || hashes.back() != h) {
                hashes.push_back(h);
                offsets.push_back(static_cast<uint32_t>(i));
            }
            postings[i] = static_cast<uint32_t>(ent[i] & 0xffffffffu);
        }
        offsets.push_back(static_cast<uint32_t>(total));
        size_t k = 0;
        for (size_t b = 0; b <= NB; b++) {
            while (k < hashes.size() && (hashes[k] >> 16) < b) k++;
            bucket[b] = static_cast<uint32_t>(k);
        }
    }

    std::pair<uint32_t, uint32_t> find(uint32_t h) const
    {
        const uint32_t b = h >> 16;
        const auto first = hashes.begin() + bucket[b];
        const auto last = hashes.begin() + bucket[b + 1];
        const auto it = std::lower_bound(first, last, h);
        if (it == last || *it != h) return {0, 0};
        const size_t k = static_cast<size_t>(it - hashes.begin());
        return {offsets[k], offsets[k + 1]};
    }
};

struct SearchCand {
    uint32_t rid;
    uint32_t shared;
    float excess;
    float cont;
    float score;
};

inline int run_search(
    const std::vector<std::string> &query_paths,
    const std::vector<std::string> &ref_paths,
    const std::string &output_path,
    int threads,
    int top_n,
    int candidates,
    const std::string &seed_shape, int seed_scaled,
    const std::string &marker_shape,
    int marker_factor,
    int max_posting,
    int batch_refs,
    bool containment_only,
    bool exclude_self,
    int max_seed_occ,
    int chain_band, int chain_max_qgap, int min_chain_seeds,
    float min_aai, int min_matched_orfs, float min_orf_cov,
    bool aggressive_filter, int xdrop, double max_evalue,
    size_t max_pair_matches,
    int gene_scaled = 1,
    const std::string &temp_dir = "")
{
    InterruptScope guard;
    set_temp_dir(temp_dir);

    const PatSet SHAPE = parse_seedspec(seed_shape.c_str());
    if (SHAPE.n == 0)
        throw std::invalid_argument("seed_shape must be \"w:off,off,...\"; got \"" + seed_shape + "\"");
    const PatSet MSHAPE = marker_shape.empty() ? SHAPE : parse_seedspec(marker_shape.c_str());
    if (MSHAPE.n == 0)
        throw std::invalid_argument("marker_shape must be \"w:off,off,...\"; got \"" + marker_shape + "\"");
    if (top_n < 1) top_n = 1;
    if (candidates < top_n) candidates = top_n;
    if (marker_factor < 1) marker_factor = 1;
    const uint32_t cap = max_posting > 0 ? static_cast<uint32_t>(max_posting) : UINT32_MAX;
    if (gene_scaled < 1) throw std::invalid_argument("gene_scaled must be >= 1");

    DistParams p{max_seed_occ,
                 chain_band, chain_max_qgap, min_chain_seeds,
                 min_aai, 0, 0,
                 patset_weight(SHAPE), patset_span(SHAPE), min_matched_orfs,
                 min_orf_cov, aggressive_filter, xdrop,
                 max_evalue, max_pair_matches};
    p.gene_scaled = static_cast<uint32_t>(gene_scaled);

    const int nq = static_cast<int>(query_paths.size());
    const int nr = static_cast<int>(ref_paths.size());
    const int bsz = (batch_refs <= 0 || batch_refs > nr) ? std::max(nr, 1) : batch_refs;
    const int nbatch = nr == 0 ? 0 : (nr + bsz - 1) / bsz;
    const uint64_t thr = marker_threshold(static_cast<uint32_t>(marker_factor));

    std::cerr << "[search] " << nq << " queries x " << nr << " refs, seeds " << seed_shape
              << " c=" << seed_scaled << ", markers " << (marker_shape.empty() ? seed_shape : marker_shape)
              << " c=" << seed_scaled << " thinned 1/" << marker_factor << ", posting cap " << max_posting << ", "
              << nbatch << " batch(es) of up to " << bsz << " refs, " << candidates
              << " candidate(s) aligned per query\n";

    const auto t_start = std::chrono::steady_clock::now();
    auto t_last = t_start;
    auto phase = [&](const char *what) {
        const auto now = std::chrono::steady_clock::now();
        std::cerr << "[search] " << what << ": "
                  << std::chrono::duration<double>(now - t_last).count() << " s\n";
        t_last = now;
    };
    auto load_one = [&](const std::vector<std::string> &paths, int i, uint32_t role) {
        std::vector<GenomeDb> one;
        load_and_seed(paths, i, i + 1, one, SHAPE, static_cast<uint32_t>(seed_scaled),
                      max_seed_occ, 1, p.gene_scaled, role);
        return std::move(one[0]);
    };
    auto markers_of = [&](const std::string &path, uint32_t role, std::string &name,
                          std::vector<uint32_t> &out, uint64_t *residues = nullptr) {
        GenomeDb g = read_db(path);
        subsample_orfs(g, p.gene_scaled, role);
        name = g.name;
        if (residues) {
            *residues = 0;
            for (const auto &h : g.seg_headers) *residues += h.aa_len;
        }
        extract_markers(g, MSHAPE, static_cast<uint32_t>(seed_scaled), max_seed_occ, thr, out);
    };

    std::vector<std::vector<uint32_t>> qm(nq);
    std::vector<std::string> qname(nq);
    std::vector<uint64_t> qres(nq, 0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int qi = 0; qi < nq; qi++) {
        if (interrupted()) continue;
        markers_of(query_paths[qi], 1, qname[qi], qm[qi], &qres[qi]);
    }
    if (interrupted()) throw Interrupted{};
    phase("query markers");

    std::unordered_map<std::string, int> query_of;
    if (p.gene_scaled == 1)
        for (int qi = 0; qi < nq; qi++) query_of.emplace(query_paths[qi], qi);
    auto ref_markers = [&](int ri, std::string &name, std::vector<uint32_t> &out) {
        const auto it = query_of.find(ref_paths[ri]);
        if (it != query_of.end()) { name = qname[it->second]; out = qm[it->second]; return; }
        markers_of(ref_paths[ri], 0, name, out);
    };

    const int ns = std::min(nr, 256);
    std::vector<std::vector<uint32_t>> sm(ns);
    std::vector<std::string> sname(ns);
    std::vector<int> sample_of(nr, -1);
    for (int k = 0; k < ns; k++)
        sample_of[std::min(static_cast<int>((static_cast<int64_t>(k) * nr + nr / 2) / ns), nr - 1)] = k;
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int ri = 0; ri < nr; ri++)
        if (sample_of[ri] >= 0) ref_markers(ri, sname[sample_of[ri]], sm[sample_of[ri]]);
    MarkerIndex sidx;
    sidx.build(sm, threads);
    auto chance_rate = [&](const std::vector<uint32_t> &m, const std::string &name,
                           std::vector<uint32_t> &counts, std::vector<double> &r) {
        std::fill(counts.begin(), counts.end(), 0);
        for (uint32_t h : m) {
            const auto pr = sidx.find(h);
            for (uint32_t e = pr.first; e < pr.second; e++) counts[sidx.postings[e]]++;
        }
        r.clear();
        for (int k = 0; k < ns; k++) {
            if (sm[k].empty() || sname[k] == name) continue;
            r.push_back(static_cast<double>(counts[k]) / sm[k].size());
        }
        if (r.empty()) return 0.0;
        std::nth_element(r.begin(), r.begin() + r.size() / 2, r.end());
        return r[r.size() / 2];
    };
    std::vector<double> lam(nq, 0.0), slam(ns, 0.0);
#pragma omp parallel num_threads(threads)
    {
        std::vector<uint32_t> counts(ns, 0);
        std::vector<double> r;
#pragma omp for schedule(dynamic, 4)
        for (int qi = 0; qi < nq; qi++) lam[qi] = chance_rate(qm[qi], qname[qi], counts, r);
#pragma omp for schedule(dynamic, 4)
        for (int k = 0; k < ns; k++) slam[k] = chance_rate(sm[k], sname[k], counts, r);
    }
    phase("chance calibration");

    std::vector<std::vector<SearchCand>> cand_a(nq), cand_b(nq);
    std::vector<uint32_t> rsize(nr, 0);
    std::vector<std::string> rname(nr);
    std::vector<float> rbg(nr, 0.0f);
    uint64_t total_markers = 0;
    auto keep_best = [&](std::vector<SearchCand> &scratch, std::vector<SearchCand> &keep) {
        if (scratch.size() > static_cast<size_t>(candidates)) {
            std::nth_element(scratch.begin(), scratch.begin() + candidates, scratch.end(),
                             [](const SearchCand &a, const SearchCand &x) { return a.score > x.score; });
            scratch.resize(candidates);
        }
        keep.insert(keep.end(), scratch.begin(), scratch.end());
        if (keep.size() > static_cast<size_t>(candidates)) {
            std::nth_element(keep.begin(), keep.begin() + candidates, keep.end(),
                             [](const SearchCand &a, const SearchCand &x) { return a.score > x.score; });
            keep.resize(candidates);
        }
        std::vector<SearchCand>(keep).swap(keep);
        scratch.clear();
    };
    for (int b = 0; b < nbatch; b++) {
        if (interrupted()) break;
        const int r0 = b * bsz, r1 = std::min(nr, r0 + bsz);
        const int nb = r1 - r0;
        std::vector<std::vector<uint32_t>> rm(nb);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (int k = 0; k < nb; k++) {
            const int ks = sample_of[r0 + k];
            if (ks >= 0) { rname[r0 + k] = sname[ks]; rm[k] = sm[ks]; }
            else ref_markers(r0 + k, rname[r0 + k], rm[k]);
            rsize[r0 + k] = static_cast<uint32_t>(rm[k].size());
        }
        MarkerIndex idx;
        idx.build(rm, threads);
        total_markers += idx.postings.size();
        rm.clear();
        rm.shrink_to_fit();

        {
            std::vector<float> cx(static_cast<size_t>(ns) * nb);
#pragma omp parallel num_threads(threads)
            {
                std::vector<uint32_t> counts(nb, 0);
#pragma omp for schedule(dynamic, 1)
                for (int k = 0; k < ns; k++) {
                    std::fill(counts.begin(), counts.end(), 0);
                    for (uint32_t h : sm[k]) {
                        const auto pr = idx.find(h);
                        if (pr.second - pr.first > cap) continue;
                        for (uint32_t e = pr.first; e < pr.second; e++) counts[idx.postings[e]]++;
                    }
                    const double nk = static_cast<double>(sm[k].size());
                    for (int j = 0; j < nb; j++) {
                        const double nj = rsize[r0 + j];
                        const double mn = std::min(nk, nj);
                        cx[static_cast<size_t>(j) * ns + k] =
                            mn > 0 ? static_cast<float>((counts[j] - slam[k] * nj) / mn) : 0.0f;
                    }
                }
#pragma omp for schedule(dynamic, 64)
                for (int j = 0; j < nb; j++) {
                    std::vector<float> v;
                    v.reserve(ns);
                    for (int k = 0; k < ns; k++)
                        if (sname[k] != rname[r0 + j]) v.push_back(cx[static_cast<size_t>(j) * ns + k]);
                    if (v.empty()) continue;
                    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
                    rbg[r0 + j] = v[v.size() / 2];
                }
            }
        }

#pragma omp parallel num_threads(threads)
        {
            std::vector<uint32_t> counts(nb, 0);
            std::vector<uint32_t> touched;
            std::vector<SearchCand> sa, sb;
#pragma omp for schedule(dynamic, 4)
            for (int qi = 0; qi < nq; qi++) {
                const auto &m = qm[qi];
                if (m.empty()) continue;
                touched.clear();
                for (uint32_t h : m) {
                    const auto pr = idx.find(h);
                    if (pr.second - pr.first > cap) continue;
                    for (uint32_t e = pr.first; e < pr.second; e++) {
                        const uint32_t k = idx.postings[e];
                        if (counts[k]++ == 0) touched.push_back(k);
                    }
                }
                const double nqm = static_cast<double>(m.size());
                for (uint32_t k : touched) {
                    const uint32_t s = counts[k];
                    counts[k] = 0;
                    const uint32_t rid = static_cast<uint32_t>(r0) + k;
                    if (exclude_self && rname[rid] == qname[qi]) continue;
                    const double nrm = static_cast<double>(rsize[rid]);
                    const double excess = s - lam[qi] * nrm;
                    if (excess <= 0) continue;
                    const float ex = static_cast<float>(excess);
                    const float cb = static_cast<float>(excess / std::min(nqm, nrm) - rbg[rid]);
                    sa.push_back({rid, s, ex, cb, ex});
                    sb.push_back({rid, s, ex, cb, cb});
                }
                keep_best(sa, cand_a[qi]);
                keep_best(sb, cand_b[qi]);
            }
        }
        std::cerr << "[search] batch " << (b + 1) << "/" << nbatch << ": " << nb << " refs, "
                  << idx.hashes.size() << " distinct markers\n";
    }
    phase("index and screen");
    std::vector<std::vector<SearchCand>> cand(nq);
    for (int qi = 0; qi < nq; qi++) {
        auto by_score = [](const SearchCand &a, const SearchCand &x) {
            return a.score != x.score ? a.score > x.score : a.rid < x.rid;
        };
        std::sort(cand_a[qi].begin(), cand_a[qi].end(), by_score);
        std::sort(cand_b[qi].begin(), cand_b[qi].end(), by_score);
        std::vector<uint32_t> seen;
        auto &cv = cand[qi];
        const size_t na = cand_a[qi].size(), nbl = cand_b[qi].size();
        for (size_t i = 0; i < std::max(na, nbl) && cv.size() < static_cast<size_t>(candidates); i++) {
            for (int l = 0; l < 2 && cv.size() < static_cast<size_t>(candidates); l++) {
                const auto &src = l == 0 ? cand_a[qi] : cand_b[qi];
                if (i >= src.size()) continue;
                if (std::find(seen.begin(), seen.end(), src[i].rid) != seen.end()) continue;
                seen.push_back(src[i].rid);
                cv.push_back(src[i]);
            }
        }
        cand_a[qi].clear(); cand_a[qi].shrink_to_fit();
        cand_b[qi].clear(); cand_b[qi].shrink_to_fit();
    }
    if (interrupted()) throw Interrupted{};

    std::ofstream out(output_path);
    if (!out) throw std::runtime_error("Cannot open: " + output_path);
    int total_hits = 0;

    if (containment_only) {
        out << "query\tref\tshared_markers\tmarker_excess\tmarker_containment\tquery_markers\tref_markers\n";
        for (int qi = 0; qi < nq; qi++) {
            const int k = std::min(static_cast<int>(cand[qi].size()), top_n);
            for (int j = 0; j < k; j++) {
                const auto &c = cand[qi][j];
                out << qname[qi] << '\t' << rname[c.rid] << '\t' << c.shared << '\t'
                    << c.excess << '\t' << c.cont << '\t' << qm[qi].size() << '\t' << rsize[c.rid] << '\n';
                total_hits++;
            }
        }
        std::cerr << "[search] " << total_hits << " hits, " << total_markers
                  << " reference markers indexed\n";
        return total_hits;
    }

    struct Hit { DistResult res; uint32_t shared; float excess, cont; };
    std::vector<std::vector<Hit>> hits(nq);
    std::vector<uint32_t> crefs;
    for (int qi = 0; qi < nq; qi++)
        for (const auto &c : cand[qi]) crefs.push_back(c.rid);
    std::sort(crefs.begin(), crefs.end());
    crefs.erase(std::unique(crefs.begin(), crefs.end()), crefs.end());
    const size_t csz = batch_refs > 0 ? static_cast<size_t>(batch_refs) : std::max<size_t>(crefs.size(), 1);
    static constexpr double QUERY_HEAP = 8e9;
    std::vector<int> slot(nr, -1);
    uint64_t n_aligned = 0;
    for (size_t c0 = 0; c0 < crefs.size(); c0 += csz) {
        if (interrupted()) break;
        const size_t c1 = std::min(crefs.size(), c0 + csz);
        std::vector<std::string> cpaths;
        for (size_t k = c0; k < c1; k++) {
            slot[crefs[k]] = static_cast<int>(k - c0);
            cpaths.push_back(ref_paths[crefs[k]]);
        }
        std::vector<GenomeDb> refs;
        load_and_seed(cpaths, 0, static_cast<int>(cpaths.size()), refs, SHAPE,
                      static_cast<uint32_t>(seed_scaled), max_seed_occ, threads, p.gene_scaled, 0);

        std::vector<std::pair<int, int>> work;
        for (int qi = 0; qi < nq; qi++)
            for (int j = 0; j < static_cast<int>(cand[qi].size()); j++)
                if (slot[cand[qi][j].rid] >= 0) work.push_back({qi, j});
        n_aligned += work.size();

        for (size_t w0 = 0; w0 < work.size() && !interrupted(); ) {
            std::vector<int> qids;
            size_t w1 = w0;
            double heap = 0;
            while (w1 < work.size()) {
                if (qids.empty() || work[w1].first != qids.back()) {
                    const double b = 0.625 * static_cast<double>(qres[work[w1].first]);
                    if (!qids.empty() && heap + b > QUERY_HEAP) break;
                    heap += b;
                    qids.push_back(work[w1].first);
                }
                w1++;
            }
            const int64_t npairs = static_cast<int64_t>(w1 - w0);
            const int chunk = static_cast<int>(std::max<int64_t>(1,
                std::min<int64_t>(candidates, npairs / (static_cast<int64_t>(threads) * 8))));
            std::vector<GenomeDb> qdb(qids.size());
            std::vector<int> qpos(w1 - w0);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
            for (int t = 0; t < static_cast<int>(qids.size()); t++)
                qdb[t] = load_one(query_paths, qids[t], 1);
            for (size_t w = w0, t = 0; w < w1; w++) {
                if (work[w].first != qids[t]) t++;
                qpos[w - w0] = static_cast<int>(t);
            }
            std::vector<uint8_t> ok(w1 - w0, 0);
            std::vector<DistResult> res(w1 - w0);
#pragma omp parallel num_threads(threads)
            {
                PairBufs bufs;
                std::vector<QueryHashGroup> qgroups;
                int last = -1;
#pragma omp for schedule(dynamic, chunk)
                for (int64_t x = 0; x < npairs; x++) {
                    if (interrupted()) continue;
                    const int t = qpos[x];
                    const GenomeDb &q = qdb[t];
                    if (q.seed_count() == 0) continue;
                    if (t != last) { group_query_seeds(q.seed_span(), qgroups); last = t; }
                    const auto &c = cand[work[w0 + x].first][work[w0 + x].second];
                    ok[x] = compare_pair(q, refs[slot[c.rid]], qgroups, p, bufs, res[x]) ? 1 : 0;
                }
            }
            for (size_t x = 0; x < w1 - w0; x++) {
                if (!ok[x]) continue;
                const int qi = work[w0 + x].first;
                const auto &c = cand[qi][work[w0 + x].second];
                hits[qi].push_back({std::move(res[x]), c.shared, c.excess, c.cont});
            }
            w0 = w1;
        }
        for (size_t k = c0; k < c1; k++) slot[crefs[k]] = -1;
        if (crefs.size() > csz)
            std::cerr << "[search] aligned " << c1 << "/" << crefs.size() << " candidate references\n";
    }
    if (interrupted()) throw Interrupted{};

    out << "query\tref\tshared_markers\tmarker_excess\tmarker_containment\t"
        << (DIST_HEADER + std::string(DIST_HEADER).find("aai"));
    char linebuf[1024];
    for (int qi = 0; qi < nq; qi++) {
        auto &hv = hits[qi];
        std::sort(hv.begin(), hv.end(), [](const Hit &a, const Hit &x) {
            return a.res.composite_score > x.res.composite_score;
        });
        const int k = std::min(static_cast<int>(hv.size()), top_n);
        for (int j = 0; j < k; j++) {
            const auto &h = hv[j];
            format_dist_row(linebuf, sizeof(linebuf), h.res);
            const char *rest = linebuf;
            for (int tabs = 0; tabs < 2 && *rest; rest++)
                if (*rest == '\t') tabs++;
            out << h.res.query_name << '\t' << h.res.ref_name << '\t' << h.shared << '\t'
                << h.excess << '\t' << h.cont << '\t' << rest;
            total_hits++;
        }
    }
    phase("alignment");
    std::cerr << "[search] " << total_hits << " hits from " << n_aligned << " aligned pairs, "
              << total_markers << " reference markers indexed\n";
    return total_hits;
}
