#pragma once

#include "db/types.h"
#include <fstream>
#include <stdexcept>
#include <cstring>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

static constexpr char CHAAI_MAGIC[8] = {'C', 'H', 'A', 'A', 'I', 'D', 'B', '2'};

static inline void write_db(
    const std::string &name, uint64_t nuc_length, uint64_t coding_nuc,
    const std::vector<ContigEntry> &contigs,
    const Arr<SegHeader> &seg_headers,
    const std::vector<uint8_t> &packed_aa,
    const std::string &path)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot open for writing: " + path);

    out.write(CHAAI_MAGIC, 8);

    auto name_len = static_cast<uint16_t>(name.size());
    out.write(reinterpret_cast<const char *>(&name_len), 2);
    out.write(name.data(), name_len);
    out.write(reinterpret_cast<const char *>(&nuc_length), 8);
    out.write(reinterpret_cast<const char *>(&coding_nuc), 8);

    auto n_contigs = static_cast<uint32_t>(contigs.size());
    out.write(reinterpret_cast<const char *>(&n_contigs), 4);
    for (const auto &c : contigs) {
        auto cname_len = static_cast<uint16_t>(c.name.size());
        out.write(reinterpret_cast<const char *>(&cname_len), 2);
        out.write(c.name.data(), cname_len);
        out.write(reinterpret_cast<const char *>(&c.offset), 8);
        out.write(reinterpret_cast<const char *>(&c.length), 8);
    }

    auto n_segs = static_cast<uint64_t>(seg_headers.size());
    out.write(reinterpret_cast<const char *>(&n_segs), 8);
    for (const auto &sh : seg_headers) {
        out.write(reinterpret_cast<const char *>(&sh.frame), 1);
        out.write(reinterpret_cast<const char *>(&sh.nuc_start), 8);
        out.write(reinterpret_cast<const char *>(&sh.aa_len), 2);
    }

    if (!packed_aa.empty()) {
        out.write(reinterpret_cast<const char *>(packed_aa.data()), packed_aa.size());
    }
}

static inline GenomeDb read_db(const std::string &path, bool lazy_headers = false)
{
    GenomeDb sk;

    char magic[8];
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open: " + path);
    in.read(magic, 8);
    if (std::memcmp(magic, CHAAI_MAGIC, 8) != 0) {
        static constexpr char OLD_MAGIC[8] = {'C', 'H', 'A', 'A', 'I', 'D', 'B', '\0'};
        if (std::memcmp(magic, OLD_MAGIC, 8) == 0)
            throw std::runtime_error(
                "Database was built by an older chaai with 32-bit coordinates and must be "
                "rebuilt: " + path);
        throw std::runtime_error("Not a chaai database: " + path);
    }

    uint16_t name_len;
    in.read(reinterpret_cast<char *>(&name_len), 2);
    sk.name.resize(name_len);
    in.read(&sk.name[0], name_len);
    in.read(reinterpret_cast<char *>(&sk.nuc_length), 8);
    in.read(reinterpret_cast<char *>(&sk.coding_nuc), 8);

    uint32_t n_contigs;
    in.read(reinterpret_cast<char *>(&n_contigs), 4);
    sk.contigs.resize(n_contigs);
    for (uint32_t i = 0; i < n_contigs; i++) {
        uint16_t cname_len;
        in.read(reinterpret_cast<char *>(&cname_len), 2);
        sk.contigs[i].name.resize(cname_len);
        in.read(&sk.contigs[i].name[0], cname_len);
        in.read(reinterpret_cast<char *>(&sk.contigs[i].offset), 8);
        in.read(reinterpret_cast<char *>(&sk.contigs[i].length), 8);
    }

    uint64_t n_segs;
    in.read(reinterpret_cast<char *>(&n_segs), 8);
    sk.n_segs_file = n_segs;

    if (lazy_headers) {
        in.seekg(static_cast<std::streamoff>(n_segs) * 11, std::ios::cur);
    } else {
        sk.seg_headers.resize(n_segs);
        for (uint64_t i = 0; i < n_segs; i++) {
            in.read(reinterpret_cast<char *>(&sk.seg_headers[i].frame), 1);
            in.read(reinterpret_cast<char *>(&sk.seg_headers[i].nuc_start), 8);
            in.read(reinterpret_cast<char *>(&sk.seg_headers[i].aa_len), 2);
        }
    }

    auto aa_file_offset = static_cast<size_t>(in.tellg());

    if (!lazy_headers) {
        sk.seg_aa_offsets.resize(n_segs);
        uint64_t off = 0;
        for (uint64_t i = 0; i < n_segs; i++) {
            sk.seg_aa_offsets[i] = off;
            off += packed_aa_bytes(sk.seg_headers[i].aa_len);
        }
    }

    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("Cannot mmap: " + path);
    struct stat st;
    fstat(fd, &st);
    sk.mmap_len = static_cast<size_t>(st.st_size);
    sk.mmap_base = mmap(nullptr, sk.mmap_len, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (sk.mmap_base == MAP_FAILED) {
        sk.mmap_base = nullptr;
        throw std::runtime_error("mmap failed: " + path);
    }
    madvise(static_cast<char *>(sk.mmap_base) + aa_file_offset,
            sk.mmap_len - aa_file_offset, MADV_SEQUENTIAL);
    sk.packed_aa_ptr = static_cast<const uint8_t *>(sk.mmap_base) + aa_file_offset;

    return sk;
}
