#include "../../Torus/src/utils/jjb_utils.cc"


namespace py = pybind11;
void init_jjb_utils(py::module_ &torus) {
	torus.def("dmatrix", py::overload_cast<int>(&dmatrix));
	torus.def("dmatrix", py::overload_cast<int, int>(&dmatrix));
	torus.def("dmatrix", py::overload_cast<int, int, int>(&dmatrix));
	torus.def("dmatrix", py::overload_cast<int, int, int, int>(&dmatrix));
	torus.def("d3array", py::overload_cast<int, int, int>(&d3array));
	torus.def("free_d3array", (&free_d3array));
	torus.def("delmatrix", py::overload_cast<double, int>(&delmatrix));
	torus.def("delmatrix", py::overload_cast<double, int, int>(&delmatrix));
	torus.def("delmatrix", py::overload_cast<double, int, int, int>(&delmatrix));
	torus.def("quadratic", (&quadratic));
	torus.def("lntp", (&lntp));
	torus.def("fourn", (&fourn));
	torus.def("rlft3", (&rlft3));
}
