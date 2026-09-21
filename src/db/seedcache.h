#pragma once

#include "db/types.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <climits>
#include <csignal>
#include <dirent.h>
#include <cerrno>
#include <csignal>
#include <dirent.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

struct SeedCacheHead {
    uint64_t magic, db_size, db_mtime;
    uint64_t n_seeds, n_segs, n_segs_file, n_orfs_full;
};
static constexpr uint64_t SEEDCACHE_MAGIC = 0x43414149'53454544ULL;

static inline void sweep_stale_tmp(const std::string &dir)
{
    DIR *d = opendir(dir.c_str());
    if (!d) return;
    while (const dirent *e = readdir(d)) {
        const std::string name = e->d_name;
        const size_t at = name.find(".seedcache.tmp");
        if (at == std::string::npos) continue;
        const long pid = strtol(name.c_str() + at + 14, nullptr, 10);
        if (pid > 0 && kill(static_cast<pid_t>(pid), 0) != 0 && errno == ESRCH) {
            unlink((dir + "/" + name).c_str());
        }
    }
    closedir(d);
}

static inline std::string &temp_dir_ref()
{
    static std::string dir;
    return dir;
}

static inline void set_temp_dir(const std::string &dir)
{
    std::string d = dir;
    if (d.empty()) {
        char cwd[PATH_MAX];
        d = std::string(getcwd(cwd, sizeof cwd) ? cwd : ".") + "/chaai_temp";
    }
    mkdir(d.c_str(), 0777);
    sweep_stale_tmp(d);
    temp_dir_ref() = d;
}

static inline const std::string &temp_dir()
{
    if (temp_dir_ref().empty()) set_temp_dir("");
    return temp_dir_ref();
}

static inline std::string seed_cache_path(const std::string &db_path, const std::string &shape,
                                          uint32_t scaled, uint32_t gene_scaled, int max_occ,
                                          uint32_t role = 0)
{
    char real[PATH_MAX];
    const char *abs = realpath(db_path.c_str(), real) ? real : db_path.c_str();

    uint64_t h = 1469598103934665603ULL;
    auto mix = [&h](const void *p, size_t n) {
        const auto *b = static_cast<const uint8_t *>(p);
        for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    };
    mix(abs, std::strlen(abs));
    mix(shape.data(), shape.size());
    mix(&scaled, sizeof scaled);
    mix(&gene_scaled, sizeof gene_scaled);
    mix(&max_occ, sizeof max_occ);
    if (gene_scaled > 1) mix(&role, sizeof role);

    char name[64];
    snprintf(name, sizeof name, "/chaai-%016llx.seedcache", static_cast<unsigned long long>(h));
    return temp_dir() + name;
}

static inline bool db_identity(const std::string &db_path, uint64_t &size, uint64_t &mtime)
{
    struct stat st{};
    if (stat(db_path.c_str(), &st) != 0) return false;
    size = static_cast<uint64_t>(st.st_size);
    mtime = static_cast<uint64_t>(st.st_mtime);
    return true;
}

static inline bool map_seed_cache(GenomeDb &sk, const std::string &path,
                                  uint64_t db_size, uint64_t db_mtime)
{
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;
    struct stat st{};
    if (fstat(fd, &st) != 0 || st.st_size < static_cast<off_t>(sizeof(SeedCacheHead))) {
        close(fd);
        return false;
    }
    const size_t len = static_cast<size_t>(st.st_size);
    void *m = mmap(nullptr, len, PROT_READ, MAP_PRIVATE, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return false;

    const auto *h = static_cast<const SeedCacheHead *>(m);
    const size_t want = sizeof(SeedCacheHead) + h->n_seeds * sizeof(Seed) +
                        h->n_segs * (sizeof(SegHeader) + sizeof(uint64_t));
    if (h->magic != SEEDCACHE_MAGIC || want != len || h->db_size != db_size ||
        h->db_mtime != db_mtime || h->n_segs_file != sk.n_segs_file) {
        munmap(m, len);
        return false;
    }

    const auto *base = static_cast<const uint8_t *>(m) + sizeof(SeedCacheHead);
    const auto *hdr = reinterpret_cast<const SegHeader *>(base + h->n_seeds * sizeof(Seed));
    const auto *off = reinterpret_cast<const uint64_t *>(hdr + h->n_segs);

    sk.seed_map_base = m;
    sk.seed_map_len = len;
    sk.seeds_view = reinterpret_cast<const Seed *>(base);
    sk.seeds_n = h->n_seeds;
    std::vector<Seed>().swap(sk.seeds);
    sk.seg_headers.view(hdr, h->n_segs);
    sk.seg_aa_offsets.view(off, h->n_segs);
    sk.n_orfs_full = h->n_orfs_full;
    return true;
}

static inline void store_seed_cache(GenomeDb &sk, const std::string &path,
                                    uint64_t db_size, uint64_t db_mtime)
{
    if (sk.seeds.empty() || sk.seg_headers.mapped()) return;

    const std::string tmp = path + ".tmp" + std::to_string(getpid());
    int fd = open(tmp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) return;

    const SeedCacheHead h{SEEDCACHE_MAGIC, db_size, db_mtime, sk.seeds.size(),
                          sk.seg_headers.size(), sk.n_segs_file, sk.n_orfs_full};
    auto put = [fd](const void *p, size_t n) {
        return n == 0 || write(fd, p, n) == static_cast<ssize_t>(n);
    };
    const bool ok = put(&h, sizeof h) &&
                    put(sk.seeds.data(), sk.seeds.size() * sizeof(Seed)) &&
                    put(sk.seg_headers.data(), sk.seg_headers.size() * sizeof(SegHeader)) &&
                    put(sk.seg_aa_offsets.data(), sk.seg_aa_offsets.size() * sizeof(uint64_t));
    close(fd);

    if (!ok || rename(tmp.c_str(), path.c_str()) != 0) {
        unlink(tmp.c_str());
        return;
    }
    map_seed_cache(sk, path, db_size, db_mtime);
}
