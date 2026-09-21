#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstddef>

struct ContigView {
    std::vector<const char *> seqs;
    std::vector<size_t> lens;
    std::vector<uint64_t> offsets;
    size_t total_len = 0;

    ContigView() = default;

    ContigView(const char *s, size_t l) : total_len(l) {
        seqs.push_back(s);
        lens.push_back(l);
        offsets.push_back(0);
    }

    ContigView(const std::vector<std::string> &contigs,
               const std::vector<uint64_t> &offs, size_t tot) : total_len(tot) {
        seqs.resize(contigs.size());
        lens.resize(contigs.size());
        offsets = offs;
        for (size_t i = 0; i < contigs.size(); i++) {
            seqs[i] = contigs[i].c_str();
            lens[i] = contigs[i].size();
        }
    }

    int n_contigs() const {
        return static_cast<int>(seqs.size());
    }

    int find_contig(uint64_t pos) const {
        auto it = std::upper_bound(offsets.begin(), offsets.end(), pos);
        return std::max(0, static_cast<int>(it - offsets.begin()) - 1);
    }
};
