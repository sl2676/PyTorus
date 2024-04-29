#include "../../Torus/src/utils/Err.cc"
/*
void TorusError(const string& m, int i) {
    toruserrno = i;
    if (i < 0) {
        throw py::value_error("ERROR: " + m); // Use PyBind11's value_error
    } else {
        py::print("WARNING:", m); // Use PyBind11 to print the warning to stderr
    }
}
*/
namespace py = pybind11;

void init_err(py::module_ &torus) {
	py::class_<TorusExcept>(torus, "TorusExcept")
        .def(py::init<string, int>(), py::arg("m"), py::arg("i"))
        .def_readwrite("msgs", &TorusExcept::msgs)
        .def_readwrite("index", &TorusExcept::index);
	torus.attr("toruserrno") = py::cast(&toruserrno);
	torus.def("TorusError", &TorusError);
}
