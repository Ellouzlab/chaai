#include "aai/dist.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(dist, m)
{
    m.doc() = "chaai pairwise distance calculation";
    m.def("run_dist", &run_dist,
          py::arg("query_paths"), py::arg("ref_paths"),
          py::arg("output_path"), py::arg("threads") = 1,
          py::arg("is_self") = false,
          py::arg("seed_shape") = "",
          py::arg("seed_scaled") = 10,
          py::arg("max_seed_occ") = 0,
          py::arg("chain_band") = 25, py::arg("chain_max_qgap") = 500,
          py::arg("min_chain_seeds") = 2,
          py::arg("min_aai") = 0.01f, py::arg("min_aligned_aa") = 30,
          py::arg("min_chain_aligned_aa") = 20,
          py::arg("min_matched_orfs") = 0,
          py::arg("min_orf_cov") = 0.20f,
          py::arg("aggressive_filter") = false,
          py::arg("xdrop") = 15,
          py::arg("max_evalue") = 0.1,
          py::arg("max_pair_matches") = 0,
          py::arg("gene_scaled") = 1,
          py::arg("temp_dir") = "",
          py::arg("no_triangle") = false,
          py::call_guard<py::gil_scoped_release>());

    py::register_exception_translator([](std::exception_ptr p) {
        try { if (p) std::rethrow_exception(p); }
        catch (const Interrupted &) { PyErr_SetNone(PyExc_KeyboardInterrupt); }
    });
}
