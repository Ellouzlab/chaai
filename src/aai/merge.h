#pragma once

#include "db/types.h"
#include "aai/chain.h"
#include <vector>
#include <algorithm>
#include <cmath>

struct QueryHashGroup {
    uint32_t hash;
    uint32_t start, end;
    float weight;
};

static inline void group_query_seeds(SeedSpan qseeds, std::vector<QueryHashGroup> &groups)
{
    groups.clear();
    for (size_t i = 0; i < qseeds.size(); ) {
        const uint32_t h = qseeds[i].hash32;
        size_t j = i;
        while (j < qseeds.size() && qseeds[j].hash32 == h) j++;
        groups.push_back({h, static_cast<uint32_t>(i), static_cast<uint32_t>(j), 1.0f});
        i = j;
    }
}

static constexpr size_t MAX_SEEDS_PER_SIDE = 8;

template <class F>
static inline void for_each_shared_hash(
    SeedSpan rseeds,
    const std::vector<QueryHashGroup> &groups, F &&fn)
{
    size_t ri = 0;
    const size_t rn = rseeds.size();

    for (const auto &g : groups) {
        if (g.weight == 0.0f) continue;

        if (ri < rn && rseeds[ri].hash32 < g.hash) {
            size_t step = 1;
            size_t pos = ri;
            while (pos + step < rn && rseeds[pos + step].hash32 < g.hash) {
                step *= 2;
            }
            size_t lo = pos + step / 2;
            size_t hi = std::min(pos + step, rn);
            while (lo < hi) {
                const size_t mid = lo + (hi - lo) / 2;
                if (rseeds[mid].hash32 < g.hash) lo = mid + 1; else hi = mid;
            }
            ri = lo;
        }

        if (ri >= rn || rseeds[ri].hash32 != g.hash) continue;

        size_t re = ri;
        while (re < rn && rseeds[re].hash32 == g.hash) re++;
        fn(g, ri, re);
        ri = re;
    }
}

static inline void sorted_merge_matches(
    SeedSpan qseeds,
    SeedSpan rseeds,
    const std::vector<QueryHashGroup> &groups,
    size_t min_matches, size_t max_matches,
    std::vector<SeedMatch> &out)
{
    out.clear();

    size_t cost[MAX_SEEDS_PER_SIDE + 1] = {};
    for_each_shared_hash(rseeds, groups,
        [&](const QueryHashGroup &g, size_t rb, size_t re) {
            const size_t nq = g.end - g.start, nr = re - rb;
            for (size_t c = 1; c <= MAX_SEEDS_PER_SIDE; c++)
                cost[c] += std::min(nq, c) * std::min(nr, c);
        });
    if (cost[MAX_SEEDS_PER_SIDE] < min_matches) return;

    size_t cap = MAX_SEEDS_PER_SIDE;
    if (max_matches > 0) {
        while (cap > 1 && cost[cap] > max_matches) cap--;
    }

    out.reserve(cost[cap]);
    for_each_shared_hash(rseeds, groups,
        [&](const QueryHashGroup &g, size_t rb, size_t re) {
            const size_t nq = g.end - g.start, nr = re - rb;
            const size_t nq_eff = std::min(nq, cap);
            const size_t nr_eff = std::min(nr, cap);
            const float w = 1.0f / std::sqrt(static_cast<float>(nq) * static_cast<float>(nr));
            for (size_t qi2 = 0; qi2 < nq_eff; qi2++) {
                const auto qidx = g.start + static_cast<uint32_t>(qi2 * nq / nq_eff);
                for (size_t ri2 = 0; ri2 < nr_eff; ri2++) {
                    const size_t ridx = rb + ri2 * nr / nr_eff;
                    out.push_back({qseeds[qidx].seg_idx, rseeds[ridx].seg_idx,
                                   qseeds[qidx].aa_pos, rseeds[ridx].aa_pos,
                                   w});
                }
            }
        });
}
