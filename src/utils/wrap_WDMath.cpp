#include "../../Torus/src/utils/WDMath.h"
namespace py = pybind11;
void init_constants(py::module_ &torus) {
	torus.attr("EulerGamma") = EulerGamma;
	torus.attr("LogofTwo") = LogofTwo;
	torus.attr("LogofTwoInv") = LogofTwoInv;
	torus.attr("LogofTen") = LogofTen;
	torus.attr("LogofTenInv") = LogofTenInv;
		// misc func
	torus.def("MathError", &MathError);
	torus.def("MathWarning", &MathWarning);
	torus.def("SphVol", &SphVol);
		// logs and exps func
	torus.def("ln", &ln);
	torus.def("ld", &ld);
	torus.def("lg", &lg);
	torus.def("Tento", &Tento);
	torus.def("Twoto", &Twoto);
		// logs of complex trig and hyperbolic func
	#ifdef __COMPLEX__
	torus.def("lnsin", &lnsin);
	torus.def("lncos", &lncos);
	torus.def("lnsinh", &lnsinh);
	torus.def("lncosh", &lncosh);
	#endif
		// gamma functions
	torus.def("LogGamma", py::overload_cast<const double> (&LogGamma));
	torus.def("GammaP", &GammaP);
	torus.def("LogGamma", py::overload_cast<const double, const double>(&LogGamma));
	torus.def("Loggamma", &Loggamma);
	#ifdef __COMPLEX__
	torus.def("LogGamma", py::overload_cast<const complex<double>>(&LogGamma));
	#endif
		// exponential integrals
	torus.def("En", &En);
	torus.def("Ei", &Ei);
		// Bessel func
	torus.def("J0", &J0);
	torus.def("J1", &J1);
	torus.def("Jn", &Jn);
	torus.def("Y0", &Y0);
	torus.def("Y1", &Y1);
	torus.def("Yn", &Yn);
	torus.def("I0", &I0);
	torus.def("I1", &I1);
	torus.def("In", &In);
	torus.def("K0", &K0);
	torus.def("K1", &K1);
	torus.def("Kn", &Kn);
		// ortho polynomials
		// hermite polynomials
	torus.def("HermiteH", py::overload_cast<const int, const double>(&HermiteH));
	torus.def("HermiteH", py::overload_cast<const int, const double, double*>(&HermiteH));
	torus.def("NormSqHermite", &NormSqHermite);
	torus.def("HermiteH_normalized", py::overload_cast<const int, const double>(&HermiteH_normalized));
	torus.def("HermiteH_normalized", py::overload_cast<const int, const double, double*>(&HermiteH_normalized));
}
