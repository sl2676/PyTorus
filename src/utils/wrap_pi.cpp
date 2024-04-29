#include "../../Torus/src/utils/Pi.h"

namespace py = pybind11;

void init_pi(py::module_ &torus) {
	torus.attr("Pi") = Pi;
	torus.attr("Pih") = Pih;
	torus.attr("Piq") = Piq;
	torus.attr("Pi3h") = Pi3h;
	torus.attr("TPi") = TPi;
	torus.attr("FPi") = FPi;
	torus.attr("iTPi") = iTPi;
	torus.attr("SPi") = SPi;
	torus.attr("STPi") = STPi;
	torus.attr("iSTPi") = iSTPi;
}
