
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "genes/build.h"

namespace py = pybind11;

static py::dict find_genes_stages(const std::string &name, const std::string &seq,
                                  const GeneFinderParams &prm)
{
    GenomeContigs gc;
    gc.name = name;
    gc.contig_names.push_back(name);
    gc.seqs.push_back(seq);
    gc.lens.push_back(seq.size());
    gc.offsets.push_back(0);
    gc.nuc_length = seq.size();

    PipelineStages st;
    call_genome_impl(gc, prm, &st);
    ContigView gv(gc.seqs, gc.offsets, gc.nuc_length);

    py::dict out;
    out["orfs"]          = py::cast(segs_to_proteins(st.orfs, gv));
    out["coding_filter"] = py::cast(segs_to_proteins(st.coding_filter, gv));
    out["select"]        = py::cast(segs_to_proteins(st.select, gv));
    out["final"]         = py::cast(segs_to_proteins(st.final_, gv));
    out["coding_density"] = st.cod_dens;
    out["min_train_aa"]  = st.diag.min_train_aa;
    out["n_train"]       = st.diag.n_train;

    py::dict f;
    const size_t n = st.diag.scored.size();
    std::vector<double> score(n), psi(n), gc3(n), asymv(n);
    std::vector<int> alen(n), frame(n), tr(n);
    std::vector<uint32_t> ns(n);
    for (size_t i = 0; i < n; i++) {
        const auto &g = st.diag.scored[i];
        score[i] = i < st.diag.sc.size() ? st.diag.sc[i] : 0.0; psi[i] = i < st.diag.psi.size() ? st.diag.psi[i] : PSI_UNSCORED; gc3[i] = 0.0;
        asymv[i] = i < st.diag.asym.size() ? st.diag.asym[i] : 0.0;
        alen[i] = g.aa_len; frame[i] = g.frame; ns[i] = g.nuc_start;
        tr[i] = i < st.diag.train_used.size() ? st.diag.train_used[i] : 0;
    }
    f["protein"]    = py::cast(segs_to_proteins(st.diag.scored, gv));
    f["score"]      = py::cast(score);
    f["psi"]  = py::cast(psi);
    f["gc3m1"]      = py::cast(gc3);
    f["asym"]         = py::cast(asymv);
    f["aa_len"]     = py::cast(alen);
    f["frame"]      = py::cast(frame);
    f["nuc_start"]  = py::cast(ns);
    f["train_used"] = py::cast(tr);
    out["features"] = f;
    return out;
}

PYBIND11_MODULE(genecall, m)
{
    m.doc() = "chaai gene calling: 6-frame translate, stop-split, segment selection";

    py::class_<GeneFinderParams>(m, "GeneFinderParams")
        .def(py::init<>())
        .def_readwrite("min_orf", &GeneFinderParams::min_orf)
        .def_readwrite("coding_filter", &GeneFinderParams::coding_filter)
        .def_readwrite("map_window", &GeneFinderParams::map_window)
        .def_readwrite("threads", &GeneFinderParams::threads);

    m.def("find_genes", &find_genes,
          py::arg("name"), py::arg("seq"), py::arg("params") = GeneFinderParams(),
          py::call_guard<py::gil_scoped_release>(),
          "Find genes in one sequence. Returns [(frame, nuc_start, aa_len, protein)].");

    m.def("find_genes_stages", &find_genes_stages,
          py::arg("name"), py::arg("seq"), py::arg("params") = GeneFinderParams(),
          "Diagnostics: same pipeline as find_genes but returns the ORF set after every stage "
          "{orfs, coding_filter, cleave, select, final, coding_density}.");

    m.def("build_from_folder", &build_from_folder,
          py::arg("infolder"), py::arg("outfolder"),
          py::arg("threads") = 1,
          py::arg("coding_filter") = false,
          py::arg("min_orf") = MIN_SEG_LEN_DEFAULT,
          py::arg("map_window") = 300,
          py::arg("verbose") = false,
          py::call_guard<py::gil_scoped_release>(),
          "Call genes in every FASTA in a folder, writing one .chaai each.");

    m.def("build_from_fasta", &build_from_fasta,
          py::arg("infasta"), py::arg("outfolder"),
          py::arg("threads") = 1,
          py::arg("coding_filter") = false,
          py::arg("min_orf") = MIN_SEG_LEN_DEFAULT,
          py::arg("map_window") = 300,
          py::arg("verbose") = false,
          py::call_guard<py::gil_scoped_release>(),
          "Call genes in a multi-FASTA, one .chaai per record.");

    m.def("write_chaai", &write_chaai,
          py::arg("path"), py::arg("name"), py::arg("nuc_length"),
          py::arg("contigs"), py::arg("segs"),
          "Write a .chaai from external ORFs. contigs=[(name,offset,len)], segs=[(frame,nuc_start,aa_len,protein)].");

    py::class_<SegHeader>(m, "SegHeader")
        .def_readonly("frame", &SegHeader::frame)
        .def_readonly("nuc_start", &SegHeader::nuc_start)
        .def_readonly("aa_len", &SegHeader::aa_len);

    py::class_<GenomeDb>(m, "GenomeDb")
        .def_readonly("name", &GenomeDb::name)
        .def_readonly("nuc_length", &GenomeDb::nuc_length)
        .def_readonly("coding_nuc", &GenomeDb::coding_nuc)
        .def_readonly("seg_headers", &GenomeDb::seg_headers);
}
