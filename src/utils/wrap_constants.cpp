#include "../../Torus/src/utils/Constants.h"
using namespace GalactoConstants;
namespace py = pybind11;

void init_constants(py::module_ &torus) {
	torus.attr("EtoPB1950_0") = py::make_tuple(EtoPB1950_0[0], EtoPB1950_0[1], EtoPB1950_0[2]);
	torus.attr("EtoPB1950_1") = py::make_tuple(EtoPB1950_1[0], EtoPB1950_1[1], EtoPB1950_1[2]);
	torus.attr("EtoPB1950_2") = py::make_tuple(EtoPB1950_2[0], EtoPB1950_2[1], EtoPB1950_2[2]);
	torus.attr("EtoPB1950") = py::make_tuple(
        py::array_t<double>(3, EtoPB1950_0),
        py::array_t<double>(3, EtoPB1950_1),
        py::array_t<double>(3, EtoPB1950_2)
    );
	torus.attr("EtoPJ1991_0") = py::make_tuple(EtoPJ1991_0[0], EtoPJ1991_0[1], EtoPJ1991_0[2]);
	torus.attr("EtoPJ1991_1") = py::make_tuple(EtoPJ1991_1[0], EtoPJ1991_1[1], EtoPJ1991_1[2]);
	torus.attr("EtoPJ1991_2") = py::make_tuple(EtoPJ1991_2[0], EtoPJ1991_2[1], EtoPJ1991_2[2]);
	torus.attr("EtoPJ1991") = py::make_tuple(
		py::array_t<double>(3, EtoPJ1991_0),
		py::array_t<double>(3, EtoPJ1991_1),
		py::array_t<double>(3, EtoPJ1991_2)
	);
	torus.attr("EtoPJ2000_0") = py::make_tuple(EtoPJ2000_0[0], EtoPJ2000_0[1], EtoPJ2000_0[2]);
	torus.attr("EtoPJ2000_1") = py::make_tuple(EtoPJ2000_1[0], EtoPJ2000_1[1], EtoPJ2000_1[2]);
	torus.attr("EtoPJ2000_2") = py::make_tuple(EtoPJ2000_2[0], EtoPJ2000_2[1], EtoPJ2000_2[2]);
	torus.attr("EtoPJ2000") = py::make_tuple(
		py::array_t<double>(3, EtoPJ2000_0),
		py::array_t<double>(3, EtoPJ2000_1),
		py::array_t<double>(3, EtoPJ2000_2)
	);

	torus.attr("Rsun_in_kpc") = Rsun_in_kpc;
	torus.attr("zsun_in_kpc") = zsun_in_kpc;
	torus.attr("vcsun_in_kms") = vcsun_in_kms;
	torus.attr("usun_in_kms") = usun_in_kms;
	torus.attr("vsun_in_kms") = vsun_in_kms;
	torus.attr("wsun_in_kms") = wsun_in_kms;
	torus.attr("Zsun") = Zsun;
	torus.attr("Rsun") = Rsun;
	torus.attr("zsun") = zsun;
	torus.attr("vcsun") = vcsun;
	torus.attr("usun") = usun;
	torus.attr("vsun") = vsun;
	torus.attr("wsun") = wsun; 

}
