#include "../../Torus/src/utils/PJMebf.cc"

namespace py = pybind11;

void init_pjmebf(py::module_ &torus) {
	torus.def("get_tag_names", [](const std::string& filename) {
        std::vector<std::string> tagnames;
        if(Ebf_GetTagNames(filename, tagnames) == 0) {
            return tagnames;
        } else {
            throw py::value_error("Failed to get tag names from file");
        }
    }, "Get a list of tag names from an EBF file");

    torus.def("get_tori_names", [](const std::string& filename) {
        std::vector<std::string> tori;
        if(Ebf_GetToriNames(filename, tori) == 0) {
            return tori;
        } else {
            throw py::value_error("Failed to get tori names from file");
        }
    }, "Get a list of tori names from an EBF file");
}
