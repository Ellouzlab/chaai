#pragma once

#include "db/types.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <omp.h>

struct SeedMatch {
    uint32_t q_seg, r_seg;
    uint16_t q_pos, r_pos;
    float weight;
};

struct AnchorPos { uint16_t q_pos, r_pos; };

struct Chain
{
    uint32_t q_seg, r_seg;
    uint16_t q_min, q_max, r_min, r_max;
    int n_matches;
    int diag_center;
    float score;
    std::vector<AnchorPos> anchors;
};

struct ChainParams
{
    int chain_band;
    int chain_max_qgap;
    int min_chain_seeds;
    int seed_weight;
};

struct ChainDPBufs {
    std::vector<float> f;
    std::vector<int> pred;
    std::vector<int> root;
    std::vector<int> tip;
    std::vector<int> order;
    std::vector<int> path;
    std::vector<uint32_t> radix;

    void ensure(int n) {
        if (static_cast<int>(f.size()) < n) {
            f.resize(n); pred.resize(n);
            root.resize(n); tip.resize(n);
        }
    }
};

static inline void dp_chain_orf_pair(
    const SeedMatch *anchors, int N,
    std::vector<Chain> &out, const ChainParams &cp,
    ChainDPBufs &db)
{
    if (N == 0) return;

    const float C = static_cast<float>(cp.chain_band);
    const int B = cp.chain_max_qgap;
    const int A = std::max(B / std::max(1, cp.seed_weight), 50);

    db.ensure(N);
    float *f   = db.f.data();
    int *pred  = db.pred.data();
    int *root  = db.root.data();
    int *tip   = db.tip.data();

    for (int i = 0; i < N; i++) {
        f[i] = anchors[i].weight;
        pred[i] = -1;

        const int x_i = static_cast<int>(anchors[i].q_pos);
        const int y_i = static_cast<int>(anchors[i].r_pos);

        const int j_start = std::max(0, i - A);
        for (int j = i - 1; j >= j_start; j--) {
            const int x_j = static_cast<int>(anchors[j].q_pos);
            if (x_i - x_j > B) break;

            const int y_j = static_cast<int>(anchors[j].r_pos);
            if (y_j >= y_i) continue;

            const int diag_shift = std::abs((x_i - x_j) - (y_i - y_j));
            const float s = C - static_cast<float>(diag_shift);
            if (s <= 0) continue;

            const float candidate = f[j] + s + anchors[i].weight;
            if (candidate > f[i]) {
                f[i] = candidate;
                pred[i] = j;
            }
        }
    }
    for (int i = 0; i < N; i++) {
        root[i] = pred[i] < 0 ? i : root[pred[i]];
        tip[i] = -1;
    }
    for (int i = 0; i < N; i++) {
        const int r = root[i];
        if (tip[r] < 0 || f[i] > f[tip[r]]) tip[r] = i;
    }
    auto &order = db.order;
    order.clear();
    for (int i = 0; i < N; i++) {
        if (tip[i] >= 0) order.push_back(tip[i]);
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) { return f[a] > f[b]; });

    for (const int best : order) {
        auto &path = db.path;
        path.clear();
        for (int c = best; c >= 0; c = pred[c]) path.push_back(c);

        const int n_anchors = static_cast<int>(path.size());
        if (n_anchors < cp.min_chain_seeds) continue;

        Chain ch;
        ch.q_seg = anchors[path.back()].q_seg;
        ch.r_seg = anchors[path.back()].r_seg;
        ch.q_min = anchors[path.back()].q_pos;
        ch.q_max = anchors[path.front()].q_pos;
        ch.r_min = anchors[path.back()].r_pos;
        ch.r_max = anchors[path.front()].r_pos;
        ch.n_matches = n_anchors;
        ch.score = f[best];

        const int mid = path[n_anchors / 2];
        ch.diag_center = static_cast<int>(anchors[mid].q_pos) -
                         static_cast<int>(anchors[mid].r_pos);

        ch.anchors.resize(n_anchors);
        for (int k = 0; k < n_anchors; k++) {
            ch.anchors[k] = {anchors[path[k]].q_pos, anchors[path[k]].r_pos};
        }
        std::sort(ch.anchors.begin(), ch.anchors.end(),
                  [](const AnchorPos &a, const AnchorPos &b) {
                      return a.q_pos < b.q_pos;
                  });

        out.push_back(std::move(ch));
    }
}

static inline void build_chains(
    std::vector<SeedMatch> &matches, std::vector<Chain> &chains,
    const ChainParams &cp, int threads,
    ChainDPBufs &dpbufs)
{
    chains.clear();
    const size_t M = matches.size();
    if (M == 0) return;

    uint32_t max_qs = 0;
    for (const auto &m : matches) {
        if (m.q_seg > max_qs) max_qs = m.q_seg;
    }
    const uint32_t nb = max_qs + 1;

    auto &rc = dpbufs.radix;
    rc.resize(2 * nb + 2);
    uint32_t *starts = rc.data();
    uint32_t *wp     = rc.data() + nb + 1;

    std::memset(starts, 0, (nb + 1) * sizeof(uint32_t));
    for (const auto &m : matches) starts[m.q_seg + 1]++;
    for (uint32_t i = 1; i <= nb; i++) starts[i] += starts[i - 1];
    for (uint32_t i = 0; i < nb; i++) wp[i] = starts[i];

    for (uint32_t b = 0; b < nb; b++) {
        while (wp[b] < starts[b + 1]) {
            SeedMatch m = matches[wp[b]];
            while (m.q_seg != b) {
                SeedMatch &dst = matches[wp[m.q_seg]++];
                std::swap(dst, m);
            }
            matches[wp[b]++] = m;
        }
    }

    for (uint32_t qs = 0; qs < nb; qs++) {
        if (starts[qs + 1] - starts[qs] >= 2) {
            std::sort(matches.data() + starts[qs],
                      matches.data() + starts[qs + 1],
                      [](const SeedMatch &a, const SeedMatch &b) {
                          if (a.r_seg != b.r_seg) return a.r_seg < b.r_seg;
                          if (a.q_pos != b.q_pos) return a.q_pos < b.q_pos;
                          return a.r_pos < b.r_pos;
                      });
        }
    }
    const int min_sz = cp.min_chain_seeds;

    if (threads <= 1) {
        for (uint32_t qs = 0; qs < nb; qs++) {
            uint32_t ps = starts[qs];
            while (ps < starts[qs + 1]) {
                uint32_t pe = ps + 1;
                while (pe < starts[qs + 1] &&
                       matches[pe].r_seg == matches[ps].r_seg) pe++;
                if (static_cast<int>(pe - ps) >= min_sz)
                    dp_chain_orf_pair(&matches[ps], static_cast<int>(pe - ps),
                                      chains, cp, dpbufs);
                ps = pe;
            }
        }
    } else {
        struct Partition { uint32_t start, end; };
        std::vector<Partition> parts;
        for (uint32_t qs = 0; qs < nb; qs++) {
            uint32_t ps = starts[qs];
            while (ps < starts[qs + 1]) {
                uint32_t pe = ps + 1;
                while (pe < starts[qs + 1] &&
                       matches[pe].r_seg == matches[ps].r_seg) pe++;
                if (static_cast<int>(pe - ps) >= min_sz) parts.push_back({ps, pe});
                ps = pe;
            }
        }
        const int np = static_cast<int>(parts.size());
        std::vector<std::vector<Chain>> pc(np);
#pragma omp parallel num_threads(std::max(1, threads))
        {
            ChainDPBufs local_bufs;
#pragma omp for schedule(dynamic, 1)
            for (int pi = 0; pi < np; pi++) {
                dp_chain_orf_pair(&matches[parts[pi].start],
                                  static_cast<int>(parts[pi].end - parts[pi].start),
                                  pc[pi], cp, local_bufs);
            }
        }
        for (auto &v : pc) {
            chains.insert(chains.end(),
                          std::make_move_iterator(v.begin()),
                          std::make_move_iterator(v.end()));
        }
    }
}
