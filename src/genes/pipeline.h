#pragma once

#include <vector>
#include <csignal>
#include <omp.h>

#include "db/types.h"
#include "genes/segment.h"
#include "genes/nuc.h"
#include "genes/contigs.h"
#include "genes/orf.h"
#include "genes/fasta.h"
#include "genes/gene_finder.h"

static volatile std::sig_atomic_t g_interrupted = 0;

struct GenomeCalls
{
    std::vector<SegHeader> headers;
    std::vector<uint8_t> packed;
};

struct GeneFinderParams
{
    int  min_orf       = MIN_SEG_LEN_DEFAULT;
    bool coding_filter = true;
    int  map_window     = 300;
    int  threads       = 1;
};

struct PipelineStages
{
    bool capture_all = true;
    bool capture_diag = true;
    std::vector<AASegment> orfs, coding_filter, select, final_;
    GfDiag diag;
    double cod_dens = -1;
};

static GenomeCalls call_genome_impl(const GenomeContigs &gc_in, const GeneFinderParams &prm,
                                    PipelineStages *st)
{
    if (!gc_in.encoded) encode_in_place(const_cast<GenomeContigs &>(gc_in));
    const GenomeContigs &gc = gc_in;
    const int min_orf = prm.min_orf;
    const bool coding_filter = prm.coding_filter;
    const int map_window = prm.map_window;
    const int threads = prm.threads;

    const int n_c = static_cast<int>(gc.seqs.size());

    std::vector<AASegment> all_segs;
    std::vector<size_t> contig_beg(n_c + 1, 0);
    std::vector<uint64_t> run_hist(RUN_HIST_MAX + 1, 0);
    {
        std::vector<int> tci;
        for (int ci = 0; ci < n_c; ci++) if (gc.lens[ci] >= 3) tci.push_back(ci);
        const size_t nt = tci.size();
        std::vector<std::vector<AASegment>> fwd(nt * 3), rev(nt * 3);
        const long ntask = static_cast<long>(nt * 3);
        std::vector<std::vector<uint64_t>> hists(
            static_cast<size_t>(std::max(1, omp_get_max_threads())),
            std::vector<uint64_t>(RUN_HIST_MAX + 1, 0));
#pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel())
        for (long t = 0; t < ntask; t++) {
            if (g_interrupted) continue;
            const size_t ti = static_cast<size_t>(t) / 3;
            const int g = static_cast<int>(t % 3), ci = tci[ti];
            std::vector<AASegment> r;
            scan_window_set(gc.seqs[ci].c_str(), gc.lens[ci], gc.offsets[ci],
                            static_cast<uint16_t>(min_orf), g, fwd[ti * 3 + g], r,
                            hists[static_cast<size_t>(omp_get_thread_num()) % hists.size()].data());
            rev[ti * 3 + (gc.lens[ci] - g) % 3] = std::move(r);
        }
        std::vector<size_t> off(nt * 3 + 1, 0);
        for (size_t x = 0; x < nt * 3; x++) off[x + 1] = off[x] + fwd[x].size() + rev[x].size();
        all_segs.resize(off[nt * 3]);
#pragma omp parallel for schedule(dynamic, 1) if (!omp_in_parallel())
        for (long x = 0; x < ntask; x++) {
            const size_t ti = static_cast<size_t>(x) / 3;
            const int f = static_cast<int>(x % 3), ci = tci[ti];
            merge_frame(fwd[x], rev[x], gc.offsets[ci], gc.lens[ci], f, all_segs.data() + off[x]);
        }
        size_t ti = 0, cur = 0;
        for (int ci = 0; ci < n_c; ci++) {
            contig_beg[ci] = cur;
            if (ti < nt && tci[ti] == ci) { cur = off[(ti + 1) * 3]; ti++; }
        }
        contig_beg[n_c] = cur;
        for (const auto &h : hists)
            for (int L = 0; L <= RUN_HIST_MAX; L++) run_hist[L] += h[L];
    }
    ContigView gv(gc.seqs, gc.offsets, gc.nuc_length);
    if (st && st->capture_all) st->orfs = all_segs;

    if (coding_filter) {
        double cdens = -1.0;
        if (st) st->diag.capture_vectors = st->capture_all || st->capture_diag;
        apply_coding_filter(all_segs, gv, g_interrupted, map_window, min_orf, &cdens, threads,
                         st ? &st->diag : nullptr, &run_hist);
        if (st) {
            st->cod_dens = cdens;
            if (st->capture_all) st->coding_filter = std::move(st->diag.kept);
            std::vector<AASegment>().swap(st->diag.kept);
            if (!st->capture_diag) {
                std::vector<AASegment>().swap(st->diag.scored);
                std::vector<double>().swap(st->diag.asym);
                std::vector<uint8_t>().swap(st->diag.train_used);
            }
        }

        if (st && st->capture_all) st->select = all_segs;
    }

    {
        const long ns = static_cast<long>(all_segs.size());
#pragma omp parallel for schedule(static) num_threads(threads) if (!omp_in_parallel())
        for (long i = 0; i < ns; i++) {
            AASegment &seg = all_segs[i];
            const int ci = gv.find_contig(seg.nuc_start);
            trim_to_start(seg, gv.seqs[ci], gv.offsets[ci], gv.lens[ci], CHAAI_TRIM_START, min_orf);
        }
    }
    if (st) st->final_ = all_segs;

    GenomeCalls out;
    out.headers.resize(all_segs.size());
    for (size_t i = 0; i < all_segs.size(); i++) {
        out.headers[i].frame = all_segs[i].frame;
        out.headers[i].nuc_start = all_segs[i].nuc_start;
        out.headers[i].aa_len = all_segs[i].aa_len;
    }
    out.packed = translate_and_pack(all_segs, gv);
    return out;
}

static inline GenomeCalls call_genome(const GenomeContigs &gc, const GeneFinderParams &prm)
{
    return call_genome_impl(gc, prm, nullptr);
}

static inline std::vector<std::string> segs_to_proteins(
    const std::vector<AASegment> &segs, const ContigView &gv)
{
    std::vector<std::string> out(segs.size());
    const int n = static_cast<int>(segs.size());
#pragma omp parallel for schedule(dynamic, 256) if (!omp_in_parallel())
    for (int i = 0; i < n; i++) {
        const auto &s = segs[i];
        const int ci = gv.find_contig(s.nuc_start);
        std::vector<uint8_t> packed(packed_aa_bytes(s.aa_len), 0);
        pack_segment_from_seq(s, gv.seqs[ci], gv.offsets[ci], packed.data());
        std::vector<uint8_t> aa(s.aa_len);
        unpack_aa_5bit(packed.data(), s.aa_len, aa.data());
        std::string prot(s.aa_len, 'X');
        for (uint16_t k = 0; k < s.aa_len; k++) prot[k] = AA_LETTERS[aa[k]];
        out[i] = std::move(prot);
    }
    return out;
}
