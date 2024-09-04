#include <iostream>
#include <algorithm>
#include <cstdio>
#include <complex>

namespace py = pybind11;


void Numerics_error(const char* msgs) {
	std::cerr << "Error in Numerics: " << msgs << std::endl;
	std::exit(1);
}

py::object WDabs(const py::object& x) {
    if (py::isinstance<py::int_>(x)) {
        int val = x.cast<int>();
        return py::cast(std::abs(val)); 
    } else if (py::isinstance<py::float_>(x)) {
        double val = x.cast<double>();
        return py::cast(std::abs(val)); 
    } else if (py::isinstance<py::object>(x) && py::hasattr(x, "real")) {
        std::complex<double> c = x.cast<std::complex<double>>();
        return py::cast(std::abs(c)); 
    } else {
        throw std::runtime_error("Unsupported type for WDabs");
    }
}

py::object sign(const py::object& x) {
	if (py::isinstance<py::int_>(x)) {
		int val = x.cast<int>();
		return py::cast((val<0)?-1:((val>0)? 1:0));
	} else if (py::isinstance<py::float_>(x)) {
		double val = x.cast<double>();
		return py::cast((val<0)?-1:((val>0)?1:0));
	} else if (py::isinstance<py::object>(x) && py::hasattr(x, "real")) {
		std::complex<double> val = x.cast<std::complex<double>>();
		return py::cast((val.real()<0)?-1:((val.real()>0)?1:0));
	} else {
		throw std::runtime_error("Unsupported type for sign");
	}
}

py::object sign(const py::object& x, const py::object& s) {
    if (py::isinstance<py::int_>(x) && py::isinstance<py::int_>(s)) {
        int val = x.cast<int>();
        int sign_val = s.cast<int>();
        return py::cast((sign_val > 0) ? WDabs(val) : -WDabs(val));
    } else if (py::isinstance<py::float_>(x) && py::isinstance<py::float_>(s)) {
        double val = x.cast<double>();
        double sign_val = s.cast<double>();
        return py::cast((sign_val > 0) ? WDabs(val) : -WDabs(val));
    } else if (py::isinstance<py::object>(x) && py::isinstance<py::object>(s) && py::hasattr(x, "real")) {
        std::complex<double> val = x.cast<std::complex<double>>();
        double sign_val = s.cast<std::complex<double>>().real();
        return py::cast((sign_val > 0) ? WDabs(val.real()) : -WDabs(val.real()));
    } else {
        throw std::runtime_error("Unsupported types for sign");
    }
}

py::object rtsafe(py::object funcd, const py::object x1, const py::object x2, const py::object xacc) {
    const int maxit = 100;
    py::object xl, xh, dx, dxo, f, df, fh, fl, rts, swap, temp;
    double dx_val, dxo_val, f_val, df_val, rts_val;

    double x1_val = x1.cast<double>();
    double x2_val = x2.cast<double>();
    double xacc_val = xacc.cast<double>();

    py::tuple result1 = py::cast<py::tuple>(funcd(x1));
    fl = result1[0];
    df = result1[1];
    
    py::tuple result2 = py::cast<py::tuple>(funcd(x2));
    fh = result2[0];
    df = result2[1];

    if (py::cast<double>(fl) * py::cast<double>(fh) >= 0.) 
        Numerics_error("root must be bracketed in rtsafe");
    
    if (py::cast<double>(fl) < 0.) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
        swap = fl;
        fl = fh;
        fh = swap;
    }

    rts = py::cast(0.5 * (x1_val + x2_val));
    dxo = py::cast(std::abs(x2_val - x1_val));
    dx = dxo;
    py::tuple result3 = py::cast<py::tuple>(funcd(rts));
    f = result3[0];
    df = result3[1];
    
    for (int j = 0; j < maxit; ++j) {
        dx_val = py::cast<double>(dx);
        dxo_val = py::cast<double>(dxo);
        f_val = py::cast<double>(f);
        df_val = py::cast<double>(df);
        rts_val = py::cast<double>(rts);
        
        if ((((rts_val - py::cast<double>(xh)) * df_val - f_val) * ((rts_val - py::cast<double>(xl)) * df_val - f_val) >= 0.) ||
            (std::abs(2. * f_val) > std::abs(dxo_val * df_val))) {
            dxo = py::cast(dx_val);
            dx = py::cast(0.5 * (py::cast<double>(xh) - py::cast<double>(xl)));
            rts = py::cast(py::cast<double>(xl) + py::cast<double>(dx));
            if (py::cast<double>(xl) == py::cast<double>(rts)) return rts;
        } else {
            dxo = py::cast(dx_val);
            dx = py::cast(f_val / df_val);
            temp = rts;
            rts = py::cast(py::cast<double>(rts) - py::cast<double>(dx));
            if (py::cast<double>(temp) == py::cast<double>(rts)) return rts;
        }

        if (std::abs(py::cast<double>(dx)) < xacc_val) return rts;

        py::tuple result4 = py::cast<py::tuple>(funcd(rts));
        f = result4[0];
        df = result4[1];
        
        if (py::cast<double>(f) < 0.) {
            xl = rts;
            fl = f;
        } else {
            xh = rts;
            fh = f;
        }
    }

    char msg[200];
    sprintf(msg, "maximum number of iterations exceeded in rtsafe:\n%g %g %g %g\n", py::cast<double>(xl), py::cast<double>(xh), py::cast<double>(fl), py::cast<double>(fh));
    Numerics_error(msg);
    return rts;
}

py::object rtsafe(const py::object& o, py::object funcd, const py::object x1, const py::object x2, const py::object xacc) {
    const int maxit = 100;
    py::object xl, xh, dx, dxo, f, df, fh, fl, rts, swap, temp;

    py::object result1 = funcd(x1);
    fl = result1.attr("f");
    df = result1.attr("df");
    
    py::object result2 = funcd(x2);
    fh = result2.attr("f");
    df = result2.attr("df");

    if (py::cast<double>(fl) * py::cast<double>(fh) >= 0.) 
        Numerics_error("root must be bracketed in rtsafe");
    
    if (py::cast<double>(fl) < 0.) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
        swap = fl;
        fl = fh;
        fh = swap;
    }

    rts = py::cast(0.5 * (py::cast<double>(x1) + py::cast<double>(x2)));
    dxo = py::cast(WDabs(py::cast<double>(x2) - py::cast<double>(x1)));
    dx = dxo;
    
    result1 = funcd(rts);
    f = result1.attr("f");
    df = result1.attr("df");
    
    for (int j = 0; j < maxit; ++j) {
        if ((((py::cast<double>(rts) - py::cast<double>(xh)) * py::cast<double>(df) - py::cast<double>(f)) * 
            ((py::cast<double>(rts) - py::cast<double>(xl)) * py::cast<double>(df) - py::cast<double>(f)) >= 0.) || 
            (WDabs(2. * py::cast<double>(f)) > WDabs(py::cast<double>(dxo) * py::cast<double>(df)))) {
            dxo = dx;
            dx = py::cast(0.5 * (py::cast<double>(xh) - py::cast<double>(xl)));
            rts = py::cast(py::cast<double>(xl) + py::cast<double>(dx));
            if (py::cast<double>(xl) == py::cast<double>(rts)) return rts;
        } else {
            dxo = dx;
            dx = py::cast(py::cast<double>(f) / py::cast<double>(df));
            temp = rts;
            rts = py::cast(py::cast<double>(rts) - py::cast<double>(dx));
            if (py::cast<double>(temp) == py::cast<double>(rts)) return rts;
        }

        if (WDabs(py::cast<double>(dx)) < py::cast<double>(xacc)) return rts;

        result1 = funcd(rts);
        f = result1.attr("f");
        df = result1.attr("df");
        
        if (py::cast<double>(f) < 0.) {
            xl = rts;
            fl = f;
        } else {
            xh = rts;
            fh = f;
        }
    }

    char msg[200];
    sprintf(msg, "maximum number of iterations exceeded in rtsafe:\n%g %g %g %g\n", py::cast<double>(xl), py::cast<double>(xh), py::cast<double>(fl), py::cast<double>(fh));
    Numerics_error(msg);
    return rts;
}

int hunt(const PyVector& xarr, const int n, const py::object x, const int j) {
	int jm, jlo=j, jhi, l=n-1;
	int ascnd=(xarr.__getitem__(l)>xarr.__getitem__(0));
	
	if (!ascnd && xarr.__getitem__(l) == xarr.__getitem__(0)) return -1;
	if ( (ascnd && x < xarr.__getitem__(0) || (!ascnd && x > xarr.__getitem__(0)))) return -1;
	if ( (ascnd && x >xarr.__getitem__(l) || (!ascnd && x < xarr.__getitem__(l))))	 return n;
	
	if (jlo < 0 || jlo > l) {
		jlo = -1;
		jhi = n;
	} else {
		int inc = 1;
		if (x >= xarr.__getitem__(jlo) == ascnd) {
			if (jlo == l) return (x == xarr.__getitem__(l)) ? l : n;
			jhi = jlo+1;
			while (x >= xarr.__getitem__(jhi) == ascnd) {
				jlo = jhi;
				inc += inc;
				jhi = jlo + inc;
				if (jhi > l) {
					jhi = n;
					break;
				}
			}
		} else {
			if (jlo == 0) return -1;
			jhi = jlo;
			jlo-=1;
			while (x < xarr.__getitem__(jlo) == ascnd) {
				jhi = jlo;
				inc += inc;
				jlo = jhi - inc;
				if (jlo < 0) {
					jlo = 0;
					break;
				}
			}
		}
	}
	while (jhi - jlo != 1) {
		jm = (jhi + jlo) >> 1;
		if (x >= xarr.__getitem__(jm) == ascnd) jlo = jm;
		else jhi = jm;
	}
	return jlo;
}

void find(int& klo, const int n, const PyVector& x, const py::object xi) {
	if (klo < 0 || klo > n-1 || x.__getitem__(klo) > xi || x.__getitem__(klo + 1) < xi) {
		klo = int( py::cast<int>(xi - x.__getitem__(0)) / py::cast<int>((x.__getitem__(n-1)-x.__getitem__(0))) * (n-1) );
		if (klo < 0 || klo >=n) {
			std::cerr << ' ' << xi << ' ' << x.__getitem__(0) << ' ' << x.__getitem__(n-1);
			Numerics_error("Find: x out of range");		
		}	
	}
}

int find_for_polev(int& j, const int n, const int m, const PyVector& x, const py::object xi) {
	int M = m;
	j = int( py::cast<int>(xi - x.__getitem__(0)) / py::cast<int>(x.__getitem__(n-1) - x.__getitem__(0)) * (n-1) );
	j = hunt(x, n, xi, j) - (m+1)/2 + 1;
	if (j >= 0 && j < n && x.__getitem__(j) == xi) M = 1;
	else if (j < 0) j = 0;
	else if (j >n-M) j= n-M;
	return M;
}

// auto convertedVector = PyVector::convertToFloatIfNeeded(p);
// p = PyVector(std::move(convertedVector));

py::object polint(const PyVector& xa, const PyVector& ya, const int n, const py::object x) {
	int i, m;
	py::object y;
	std::vector<py::object> p_data(n);
	PyVector p(py::cast(p_data));

	for (i = 0; i < n; i++) p.__setitem__(py::cast<int>(p.__getitem__(i)), ya.__getitem__(i));
	for (m = 1; m < n; m++)
	for (i = 0; i < n - m; i++) {
		if (xa.__getitem__(i) == xa.__getitem__(i+m)) Numerics_error("x's not distinct in polint");
		p.__setitem__(i, ( (x-xa.__getitem__(i+m) * p.__getitem__(i) + 
							(xa.__getitem__(i) - x) * (p.__getitem__(i+1)) /
							(xa.__getitem__(i) - xa.__getitem__(i+m))
						)));
	}
	y = p.__getitem__(0);
	return y;	
}
py::object polev(py::object x, const PyVector& xarr, const int n, const int m=4) {
	int j, M=find_for_polev(j, n, m, xarr, x);
	return polint(xarr+j, yarr+j, M, x);
}


void init_numerics_templates(py::module_ &torus) {
	torus.def("find_for_polev", [](int& j, const int n, const int m, const PyVector& x, const py::object xi) {
		return find_for_polev(j, n, m, x, xi);
	});
	torus.def("polev", [](py::object x, const PyVector& xarr, const int n, const int m=4) {
		return polev(x, xarr, n, m);
	});
	torus.def("polint", [](const PyVector& xa, const PyVector& ya, const int n, const py::object x) {
		return polint(xa, ya, n, x);
	});
	torus.def("WDabs", [](const py::object& x) {
        return WDabs(x);
    });
	torus.def("sign", [](const py::object& x) {
		return sign(x);
	});
	torus.def("rtsafe", [](py::object funcd, const py::object x1, const py::object x2, const py::object xacc) {
        return rtsafe(funcd, x1, x2, xacc);
    });
	torus.def("rtsafe", [](const py::object& o, py::object funcd, const py::object x1, const py::object x2, const py::object xacc) {
        return rtsafe(o, funcd, x1, x2, xacc);
    });
	torus.def("sign", [](const py::object& x, const py::object& s) {
        return sign(x, s);
    });
	torus.def("hunt", [](const PyVector& xarr, const int n, const py::object x, const int j) {
		return hunt(xarr, n, x, j);
	});
	torus.def("find", [](int& klo, const int n, const PyVector& x, py::object xi) {
		return find(klo, n, x, xi);
	});
}
