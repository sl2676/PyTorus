#include "../../Torus/src/utils/Err.h"
namespace py = pybind11;

void init_err(py::module_ &torus) {
	py::class_<TorusExcept>(torus, "TorusExcept")
        .def(py::init<string, int>(), py::arg("m"), py::arg("i"))
        .def_readwrite("msgs", &TorusExcept::msgs)
        .def_readwrite("index", &TorusExcept::index);
	torus.attr("toruserrno") = py::cast(&toruserrno);
	torus.def("TorusError", &TorusError);
}
