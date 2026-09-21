#pragma once
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <cstdlib>

#include "aai/chain.h"
#include "aai/seeds.h"
#include "db/io.h"
#include "db/seedcache.h"
#include "db/types.h"
#include "aai/hash.h"

#include <cstring>
#include <string>
#include <vector>
#include <omp.h>

inline constexpr int PROGRESS_INTERVAL = 256;

inline void subsample_orfs(GenomeDb &sk, uint32_t gene_scaled, uint32_t role = 0)
{
    sk.n_orfs_full = sk.seg_headers.size();
    if (gene_scaled <= 1 || sk.seg_headers.empty()) return;
    const uint64_t thr = UINT64_MAX / gene_scaled;

    uint64_t salt = 0xcbf29ce484222325ULL;
    for (unsigned char ch : sk.name) salt = (salt ^ ch) * 0x100000001b3ULL;
    salt = hash64(salt ^ role);

    std::vector<SegHeader> kept_h;
    std::vector<uint64_t> kept_off;
    kept_h.reserve(sk.seg_headers.size() / gene_scaled + 16);
    kept_off.reserve(kept_h.capacity());
    for (size_t si = 0; si < sk.seg_headers.size(); si++) {
        const uint8_t *pk = sk.packed_aa_ptr + sk.seg_aa_offsets[si];
        const size_t nb = packed_aa_bytes(sk.seg_headers[si].aa_len);
        uint64_t h = salt ^ sk.seg_headers[si].aa_len;
        size_t b = 0;
        for (; b + 8 <= nb; b += 8) {
            uint64_t chunk;
            std::memcpy(&chunk, pk + b, 8);
            h = hash64(h ^ chunk);
        }
        if (b < nb) {
            uint64_t chunk = 0;
            std::memcpy(&chunk, pk + b, nb - b);
            h = hash64(h ^ chunk);
        }
        if (h < thr) {
            kept_h.push_back(sk.seg_headers[si]);
            kept_off.push_back(sk.seg_aa_offsets[si]);
        }
    }
    sk.seg_headers.swap(kept_h);
    sk.seg_aa_offsets.swap(kept_off);
}

inline void load_and_seed(
    const std::vector<std::string> &paths, int start, int end,
    std::vector<GenomeDb> &out, const PatSet &shape,
    uint32_t seed_scaled, int max_seed_occ, int threads,
    uint32_t gene_scaled = 1, uint32_t role = 0)
{
    const int n = end - start;
    out.resize(n);

    std::string shape_key;
    for (int pi = 0; pi < shape.n; pi++) {
        shape_key += std::to_string(shape.p[pi].weight) + ":";
        for (int k = 0; k < shape.p[pi].weight; k++)
            shape_key += std::to_string(shape.p[pi].offsets[k]) + ",";
        shape_key += ";";
    }

#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (int i = 0; i < n; i++) {
        const std::string &path = paths[start + i];
        uint64_t db_size = 0, db_mtime = 0;
        const bool cache = db_identity(path, db_size, db_mtime);
        const std::string cpath =
            cache ? seed_cache_path(path, shape_key, seed_scaled, gene_scaled, max_seed_occ, role)
                  : std::string();

        if (cache) {
            GenomeDb g = read_db(path, true);
            if (map_seed_cache(g, cpath, db_size, db_mtime)) { out[i] = std::move(g); continue; }
        }

        out[i] = read_db(path);
        subsample_orfs(out[i], gene_scaled, role);
        out[i].seeds = extract_seeds_spec_from_db(out[i], shape, seed_scaled, max_seed_occ);
        if (cache) store_seed_cache(out[i], cpath, db_size, db_mtime);
    }
}
