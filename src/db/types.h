#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/mman.h>

static inline size_t packed_aa_bytes(uint16_t len)
{
    return (static_cast<size_t>(len) * 5 + 7) / 8;
}

static inline void pack_aa_5bit(const uint8_t *aa, uint16_t len, uint8_t *out)
{
    uint64_t buf = 0;
    int bits = 0;
    uint8_t *dst = out;
    for (uint16_t i = 0; i < len; i++) {
        buf |= (static_cast<uint64_t>(aa[i] & 0x1F)) << bits;
        bits += 5;
        while (bits >= 8) {
            *dst++ = static_cast<uint8_t>(buf & 0xFF);
            buf >>= 8;
            bits -= 8;
        }
    }
    if (bits > 0) {
        *dst = static_cast<uint8_t>(buf & 0xFF);
    }
}

static inline void unpack_aa_5bit(const uint8_t *packed, uint16_t len, uint8_t *aa)
{
    uint64_t buf = 0;
    int bits = 0;
    const uint8_t *src = packed;
    for (uint16_t i = 0; i < len; i++) {
        while (bits < 5) {
            buf |= (static_cast<uint64_t>(*src++)) << bits;
            bits += 8;
        }
        aa[i] = static_cast<uint8_t>(buf & 0x1F);
        buf >>= 5;
        bits -= 5;
    }
}

struct Seed
{
    uint32_t hash32;
    uint32_t seg_idx;
    uint16_t aa_pos;
    uint16_t _pad = 0;
};
static_assert(sizeof(Seed) == 12, "Seed layout");

struct SegHeader
{
    uint8_t frame;
    uint64_t nuc_start;
    uint16_t aa_len;
};

struct ContigEntry {
    std::string name;
    uint64_t offset;
    uint64_t length;
};

template <class T>
struct Arr {
    std::vector<T> own;
    const T *p = nullptr;
    size_t n = 0;

    bool mapped() const { return p != nullptr; }
    size_t size() const { return p ? n : own.size(); }
    bool empty() const { return size() == 0; }
    const T *data() const { return p ? p : own.data(); }
    const T &operator[](size_t i) const { return p ? p[i] : own[i]; }
    T &operator[](size_t i) { return own[i]; }
    const T *begin() const { return data(); }
    const T *end() const { return data() + size(); }
    void resize(size_t m) { own.resize(m); }
    void reserve(size_t m) { own.reserve(m); }
    void push_back(const T &v) { own.push_back(v); }
    void swap(std::vector<T> &v) { own.swap(v); p = nullptr; n = 0; }
    void view(const T *q, size_t m) { std::vector<T>().swap(own); p = q; n = m; }
    Arr &operator=(std::vector<T> &&v) { own = std::move(v); p = nullptr; n = 0; return *this; }
    Arr &operator=(const std::vector<T> &v) { own = v; p = nullptr; n = 0; return *this; }
    Arr() = default;
    Arr(Arr &&) = default;
    Arr &operator=(Arr &&) = default;
};

struct SeedSpan {
    const Seed *p = nullptr;
    size_t n = 0;
    SeedSpan() = default;
    SeedSpan(const std::vector<Seed> &v) : p(v.data()), n(v.size()) {}
    SeedSpan(const Seed *q, size_t m) : p(q), n(m) {}
    const Seed &operator[](size_t i) const { return p[i]; }
    size_t size() const { return n; }
    bool empty() const { return n == 0; }
    const Seed *begin() const { return p; }
    const Seed *end() const { return p + n; }
};

struct GenomeDb
{
    std::string name;
    uint64_t nuc_length = 0;
    uint64_t coding_nuc = 0;
    std::vector<ContigEntry> contigs;

    std::vector<Seed> seeds;
    Arr<SegHeader> seg_headers;
    uint64_t n_orfs_full = 0;
    uint64_t n_segs_file = 0;

    std::vector<uint8_t> packed_aa_owned;
    const uint8_t *packed_aa_ptr = nullptr;
    Arr<uint64_t> seg_aa_offsets;

    void *mmap_base = nullptr;
    size_t mmap_len = 0;

    void *seed_map_base = nullptr;
    size_t seed_map_len = 0;
    size_t seeds_n = 0;

    const Seed *seeds_view = nullptr;
    SeedSpan seed_span() const {
        return seeds_view ? SeedSpan(seeds_view, seeds_n) : SeedSpan(seeds);
    }
    size_t seed_count() const { return seeds_view ? seeds_n : seeds.size(); }

    void unpack_segment(uint32_t idx, uint8_t *out) const {
        unpack_aa_5bit(packed_aa_ptr + seg_aa_offsets[idx],
                       seg_headers[idx].aa_len, out);
    }

    void compute_coding_nuc()
    {
        coding_nuc = 0;
        if (seg_headers.empty()) return;

        struct Iv { uint64_t s, e; };
        std::vector<Iv> ivs;
        ivs.reserve(seg_headers.size());
        for (const auto &sh : seg_headers) {
            ivs.push_back({sh.nuc_start,
                           sh.nuc_start + static_cast<uint64_t>(sh.aa_len) * 3});
        }
        std::sort(ivs.begin(), ivs.end(),
                  [](const Iv &a, const Iv &b) { return a.s < b.s; });

        uint64_t cs = ivs[0].s;
        uint64_t ce = ivs[0].e;
        for (size_t i = 1; i < ivs.size(); i++) {
            if (ivs[i].s <= ce) {
                ce = std::max(ce, ivs[i].e);
            } else {
                coding_nuc += ce - cs;
                cs = ivs[i].s;
                ce = ivs[i].e;
            }
        }
        coding_nuc += ce - cs;
    }

    GenomeDb() = default;
    GenomeDb(GenomeDb &&o) noexcept { move_from(std::move(o)); }
    GenomeDb &operator=(GenomeDb &&o) noexcept {
        if (this != &o) { cleanup(); move_from(std::move(o)); }
        return *this;
    }
    GenomeDb(const GenomeDb &) = delete;
    GenomeDb &operator=(const GenomeDb &) = delete;
    ~GenomeDb() { cleanup(); }

private:
    void cleanup() {
        if (mmap_base) {
            munmap(mmap_base, mmap_len);
            mmap_base = nullptr;
        }
        if (seed_map_base) {
            munmap(seed_map_base, seed_map_len);
            seed_map_base = nullptr;
        }
    }

    void move_from(GenomeDb &&o) {
        name = std::move(o.name);
        nuc_length = o.nuc_length;
        coding_nuc = o.coding_nuc;
        contigs = std::move(o.contigs);
        seeds = std::move(o.seeds);
        seg_headers = std::move(o.seg_headers);
        n_orfs_full = o.n_orfs_full;
        n_segs_file = o.n_segs_file;
        packed_aa_owned = std::move(o.packed_aa_owned);
        seg_aa_offsets = std::move(o.seg_aa_offsets);
        mmap_base = o.mmap_base;
        mmap_len = o.mmap_len;
        o.mmap_base = nullptr;
        o.mmap_len = 0;
        seeds_view = o.seeds_view; o.seeds_view = nullptr;
        seed_map_base = o.seed_map_base;
        seed_map_len = o.seed_map_len;
        seeds_n = o.seeds_n;
        o.seed_map_base = nullptr;
        o.seed_map_len = 0;
        o.seeds_n = 0;
        packed_aa_ptr = mmap_base ? o.packed_aa_ptr : packed_aa_owned.data();
        o.packed_aa_ptr = nullptr;
    }
};
