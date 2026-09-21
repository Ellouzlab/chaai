#include "aai/search.h"
#include "run/interrupt.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <exception>

namespace py = pybind11;

PYBIND11_MODULE(search, m)
{
    m.doc() = "chaai search: marker index over the references, dist on the closest candidates";
    m.def("run_search", &run_search,
          py::arg("query_paths"), py::arg("ref_paths"),
          py::arg("output_path"), py::arg("threads") = 1,
          py::arg("top_n") = 1,
          py::arg("candidates") = 10,
          py::arg("seed_shape") = "",
          py::arg("seed_scaled") = 10,
          py::arg("marker_shape") = "",
          py::arg("marker_factor") = 1,
          py::arg("max_posting") = 400,
          py::arg("batch_refs") = 0,
          py::arg("containment_only") = false,
          py::arg("exclude_self") = false,
          py::arg("max_seed_occ") = 0,
          py::arg("chain_band") = 25,
          py::arg("chain_max_qgap") = 500,
          py::arg("min_chain_seeds") = 2,
          py::arg("min_aai") = 0.01f,
          py::arg("min_matched_orfs") = 0,
          py::arg("min_orf_cov") = 0.20f,
          py::arg("aggressive_filter") = false,
          py::arg("xdrop") = 15,
          py::arg("max_evalue") = 0.1,
          py::arg("max_pair_matches") = 0,
          py::arg("gene_scaled") = 1,
          py::arg("temp_dir") = "",
          py::call_guard<py::gil_scoped_release>());

    py::register_exception_translator([](std::exception_ptr p) {
        try { if (p) std::rethrow_exception(p); }
        catch (const Interrupted &) { PyErr_SetNone(PyExc_KeyboardInterrupt); }
    });
}
