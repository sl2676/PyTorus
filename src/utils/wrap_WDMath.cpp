#include "../../Torus/src/utils/WDMath.h"

namespace py = pybind11;
const int    maxit  = 100;
const double fpmin  = 1.e-40,
             eps    = 1.e-9,
	     logeps =-20.72326583694641115616192309216;


static void gser(double& gamser, const double a, const double x, double& lng) {
    lng=LogGamma(a);
    if(x<=0.) {
	if(x<0.) MathError("x<0 in gser()");
	gamser=0.;
	return;
    }
    int    n;
    double sum,del,ap;
    ap  = a;
    del = sum = 1.0/a;
    for(n=1; n<=maxit; n++) {
	++ap;
	del *= x/ap;
	sum += del;
	if(WDabs(del) < WDabs(sum)*eps) {
	    gamser = sum*exp(-x+a*log(x)-lng);
	    return;
	}
    }
    MathError("a too large, maxit too small in gser()");
}

static void gcf(double& gammcf, double a, double x, double& lng) {
    int i;
    double an,b,c,d,del,h;

    lng = LogGamma(a);
    b   = x+1.-a;
    c   = 1./fpmin;
    d   = 1./b;
    h   = d;
    for(i=1; i<=maxit; i++) {
	an =-i*(i-a);
	b += 2.;
	d  = an*d+b; if(WDabs(d)<fpmin) d=fpmin;
	c  = b+an/c; if(WDabs(c)<fpmin) c=fpmin;
	d  = 1./d;
	del= d*c;
	h *= del;    if(WDabs(del-1.)<eps) break;
    }
    if(i>maxit) MathError("a too large, maxit too small in gcf()");
    gammcf = exp(-x+a*log(x)-lng) * h;
}

void init_wdmath(py::module_ &torus) {
	torus.attr("EulerGamma") = EulerGamma;
	torus.attr("LogofTwo") = LogofTwo;
	torus.attr("LogofTwoInv") = LogofTwoInv;
	torus.attr("LogofTen") = LogofTen;
	torus.attr("LogofTenInv") = LogofTenInv;
		// misc func
	torus.def("gser", &gser);
	torus.def("gcf", &gcf);
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
