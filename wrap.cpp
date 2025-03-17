#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string>

#include "./Torus/src/utils/FreeMemory.h"


#include "./src/utils/wrap_pi.cpp"
#include "./src/utils/wrap_units.cpp"
#include "./src/utils/wrap_constants.cpp"
#include "./src/utils/wrap_compress.cpp"
#include "./src/utils/wrap_err.cpp"
#include "./src/utils/wrap_chb.cpp"
#include "./src/utils/wrap_vector.cpp"
#include "./src/utils/wrap_matrix.cpp"
#include "./src/utils/wrap_numerics.cpp"
#include "./src/utils/wrap_wdmath.cpp"
#include "./src/utils/wrap_PJMNum.cpp"
#include "./src/utils/wrap_numerics_templates.cpp"
//#include "./src/utils/wrap_PJMCoords.cpp"
//#include "./src/utils/wrap_pspline.cpp"

void init_pi(py::module_ &);
void init_units(py::module_ &);
void init_constants(py::module_ &);
void init_compress(py::module_ &);
void init_err(py::module_ &);
void init_chb(py::module_ &);
void init_vector(py::module_ &);
void init_matrix(py::module_ &);
void init_pjmebf(py::module_ &);
void init_numerics(py::module_ &);
void init_wdmath(py::module_ &);
void init_pjmnum(py::module_ &);
void init_pjmcoords(py::module_ &);
//void init_pspline(py::module_ &);

namespace py = pybind11;
using namespace std;
using namespace pybind11::literals;

PYBIND11_MODULE(_torus, torus) {
	torus.doc() = "pytorus torus library";
	init_pi(torus);
	init_units(torus);
	init_constants(torus);
	init_compress(torus);
	init_err(torus);
	init_vector(torus);
	init_matrix(torus);
  //init_chb(torus);
	init_numerics(torus);
	init_wdmath(torus);
	init_pjmnum(torus);
	init_numerics_templates(torus);
//	init_pjmcoords(torus);
	//init_pspline(torus);
}


