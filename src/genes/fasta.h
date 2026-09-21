#pragma once
#include <cstdint>
#include <string>
#include <omp.h>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <cstring>
#include <cctype>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <omp.h>

static bool is_fasta_ext(const std::string &path)
{
    const auto dot = path.rfind('.');
    if (dot == std::string::npos) return false;
    std::string ext = path.substr(dot);
    for (auto &c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext == ".fa" || ext == ".fna" || ext == ".fasta" ||
           ext == ".faa" || ext == ".fas";
}

class FastaReader
{
    std::ifstream fin;
    std::string pending_name;
    bool has_pending = false;

public:
    explicit FastaReader(const std::string &path) : fin(path)
    {
        if (!fin) throw std::runtime_error("Cannot open fasta: " + path);
    }

    bool next(std::string &name, std::string &seq)
    {
        name.clear();
        seq.clear();

        if (has_pending) {
            name = std::move(pending_name);
            has_pending = false;
        }
        std::string line;
        while (name.empty() && std::getline(fin, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty() && line[0] == '>') {
                const size_t end = line.find_first_of(" \t", 1);
                name = line.substr(1, end == std::string::npos
                                      ? std::string::npos : end - 1);
            }
        }
        if (name.empty()) return false;

        while (std::getline(fin, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line[0] == '>') {
                const size_t end = line.find_first_of(" \t", 1);
                pending_name = line.substr(1, end == std::string::npos
                                              ? std::string::npos : end - 1);
                has_pending = true;
                break;
            }
            seq.append(line);
        }
        return seq.size() >= 3;
    }
};

static std::string safe_filename(const std::string &name)
{
    std::string s = name;
    for (char &c : s) {
        if (c == '/' || c == '\\' || c == ':' || c == '*' ||
            c == '?' || c == '"' || c == '<' || c == '>' || c == '|') {
            c = '_';
        }
    }
    return s;
}

struct GenomeContigs
{
    std::string name;
    std::vector<std::string> seqs;
    std::vector<std::string> contig_names;
    std::vector<size_t> lens;
    std::vector<uint64_t> offsets;
    uint64_t nuc_length = 0;
    bool encoded = false;
};

static inline void encode_in_place(GenomeContigs &gc)
{
    if (gc.encoded) return;
    const int nr = static_cast<int>(gc.seqs.size());
#pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel())
    for (int i = 0; i < nr; i++) {
        char *s = gc.seqs[i].empty() ? nullptr : &gc.seqs[i][0];
        const size_t n = gc.seqs[i].size();
        encode_bytes(s, s, n);
    }
    gc.encoded = true;
}

static GenomeContigs read_genome_contigs(
    const std::string &path, const std::string &genome_name, bool encode = false)
{
    GenomeContigs gc;
    gc.name = genome_name;
    gc.encoded = encode;

    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("Cannot open: " + path);
    struct stat st;
    fstat(fd, &st);
    const auto fsize = static_cast<size_t>(st.st_size);
    if (fsize == 0) { close(fd); return gc; }
    const auto *data = static_cast<const char *>(
        mmap(nullptr, fsize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0));
    close(fd);
    if (data == MAP_FAILED) throw std::runtime_error("mmap failed: " + path);
    madvise(const_cast<char *>(data), fsize, MADV_SEQUENTIAL);

    struct Range { size_t seq_start, seq_end; };
    std::vector<Range> ranges;
    std::vector<std::string> header_names;
    const char *ptr = data;
    const char *fend = data + fsize;
    while (ptr < fend) {
        ptr = static_cast<const char *>(memchr(ptr, '>', fend - ptr));
        if (!ptr) break;
        const char *hdr_start = ptr + 1;
        const char *nl = static_cast<const char *>(memchr(ptr, '\n', fend - ptr));
        if (!nl) break;
        const char *name_end = hdr_start;
        while (name_end < nl && *name_end != ' ' && *name_end != '\t'
               && *name_end != '\r') name_end++;
        header_names.emplace_back(hdr_start, name_end);
        ptr = nl + 1;
        const auto seq_start = static_cast<size_t>(ptr - data);
        const auto *next = static_cast<const char *>(memchr(ptr, '>', fend - ptr));
        const size_t seq_end = next ? static_cast<size_t>(next - data) : fsize;
        ranges.push_back({seq_start, seq_end});
        ptr = next ? next : fend;
    }

    const int nr = static_cast<int>(ranges.size());
    gc.seqs.resize(nr);
    gc.lens.resize(nr);
    gc.contig_names = std::move(header_names);
    struct Chunk { int rec; size_t raw_start, raw_end, n_nl, n_cr, out_off; };
    std::vector<Chunk> chunks;
    constexpr size_t CHUNK = 1u << 22;
    for (int i = 0; i < nr; i++) {
        size_t p = ranges[i].seq_start;
        while (p < ranges[i].seq_end) {
            size_t e = std::min(p + CHUNK, ranges[i].seq_end);
            if (e < ranges[i].seq_end) {
                const char *nl = static_cast<const char *>(memchr(data + e, '\n', ranges[i].seq_end - e));
                e = nl ? static_cast<size_t>(nl - data) + 1 : ranges[i].seq_end;
            }
            chunks.push_back({i, p, e, 0, 0, 0});
            p = e;
        }
    }
    const long nch = static_cast<long>(chunks.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (long c = 0; c < nch; c++) {
        Chunk &ch = chunks[c];
        const char *q = data + ch.raw_start, *e = data + ch.raw_end;
        size_t nl = 0, cr = 0;
        while (q < e) {
            const char *n = static_cast<const char *>(memchr(q, '\n', e - q));
            if (!n) break;
            nl++;
            if (n > data + ch.raw_start && n[-1] == '\r') cr++;
            q = n + 1;
        }
        ch.n_nl = nl; ch.n_cr = cr;
    }
    {
        int cur = -1; size_t off = 0;
        for (auto &ch : chunks) {
            if (ch.rec != cur) { if (cur >= 0) gc.lens[cur] = off; cur = ch.rec; off = 0; }
            ch.out_off = off;
            off += (ch.raw_end - ch.raw_start) - ch.n_nl - ch.n_cr;
        }
        if (cur >= 0) gc.lens[cur] = off;
        for (int i = 0; i < nr; i++) gc.seqs[i].reserve(gc.lens[i]);
#if defined(MADV_POPULATE_WRITE)
        {
            const size_t pg = static_cast<size_t>(sysconf(_SC_PAGESIZE));
            std::vector<std::pair<char *, size_t>> spans;
            for (int i = 0; i < nr; i++) {
                if (gc.lens[i] < (1u << 20)) continue;
                char *b = &gc.seqs[i][0];
                const uintptr_t lo = reinterpret_cast<uintptr_t>(b) & ~(pg - 1);
                const uintptr_t hi = (reinterpret_cast<uintptr_t>(b) + gc.lens[i] + pg - 1) & ~(pg - 1);
                for (uintptr_t a = lo; a < hi; a += (64u << 20))
                    spans.push_back({reinterpret_cast<char *>(a), std::min<size_t>(64u << 20, hi - a)});
            }
            const long nsp = static_cast<long>(spans.size());
#pragma omp parallel for schedule(dynamic, 1)
            for (long k = 0; k < nsp; k++) madvise(spans[k].first, spans[k].second, MADV_POPULATE_WRITE);
        }
#endif
        for (int i = 0; i < nr; i++) gc.seqs[i].resize(gc.lens[i]);
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (long c = 0; c < nch; c++) {
        const Chunk &ch = chunks[c];
        char *dst = &gc.seqs[ch.rec][0] + ch.out_off;
        size_t p = ch.raw_start;
        while (p < ch.raw_end) {
            const auto *nl = static_cast<const char *>(memchr(data + p, '\n', ch.raw_end - p));
            size_t part = nl ? static_cast<size_t>(nl - (data + p)) : (ch.raw_end - p);
            if (part > 0 && data[p + part - 1] == '\r') part--;
            if (part > 0) {
                if (encode) encode_bytes(data + p, dst, part);
                else memcpy(dst, data + p, part);
                dst += part;
            }
            p += (nl ? static_cast<size_t>(nl - (data + p)) + 1 : (ch.raw_end - p));
        }
    }

    munmap(const_cast<char *>(data), fsize);

    gc.offsets.resize(nr);
    uint64_t total = 0;
    for (int i = 0; i < nr; i++) {
        gc.offsets[i] = total;
        total += gc.lens[i];
    }
    gc.nuc_length = total;
    return gc;
}
