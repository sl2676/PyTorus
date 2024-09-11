#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
namespace py = pybind11;

double roundIfCloseToInteger(double value, double epsilon = 1e-9) {
    double nearestInt = std::round(value);
    if (std::abs(value - nearestInt) < epsilon) {
        return nearestInt;
    }
    return value;
}

double sanitizePrecision(double value, double threshold = 1e-9) {
	if (std::abs(value) < threshold) {
		return 0.0;
	}
	return value;
}

int Numerics_message(const char* msgs) {
	std::cerr << "ERROR in Numerics: " << msgs << '\n';
	return -1;
}


int GaussJordan(PyMatrix& a, int n, PyMatrix& b, int m) {
	a.promoteMatrixVariantIfNeeded<double>();
	b.promoteMatrixVariantIfNeeded<double>();

    std::vector<int> ipiv(n, 0);
    std::vector<int> indxr(n), indxc(n);

    for (int i = 0; i < n; ++i) {
        double big = 0.0;
        int irow = -1, icol = -1;

        for (int j = 0; j < n; ++j) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < n; ++k) {
                    if (ipiv[k] == 0) {
                        double aij = std::abs(py::cast<double>(a.getValueAt(j, k)));
                        if (aij >= big) {
                            big = aij;
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) {
                        return -1;  
                    }
                }
            }
        }

        if (irow < 0 || icol < 0) return -1;  

        ++(ipiv[icol]);

        if (irow != icol) {
            for (int l = 0; l < n; ++l) {
                py::object temp = a.getValueAt(irow, l);
                a.set_value(irow, l, a.getValueAt(icol, l));
                a.set_value(icol, l, temp);
            }
            for (int l = 0; l < m; ++l) {
                py::object temp = b.getValueAt(irow, l);
                b.set_value(irow, l, b.getValueAt(icol, l));
                b.set_value(icol, l, temp);
            }
        }

        indxr[i] = irow;
        indxc[i] = icol;

        double piv = py::cast<double>(a.getValueAt(icol, icol));
        if (piv == 0.0) return -1;  

        double pivinv = 1.0 / piv;
		if (std::abs(pivinv) < 1e-10) {
			pivinv = 0.0;
		}
        for (int l = 0; l < n; ++l) {
            auto value = py::float_(py::cast<double>(a.getValueAt(icol, l)) * pivinv);
			
			a.set_value(icol, l, py::float_(sanitizePrecision(py::cast<double>(a.getValueAt(icol, l)) * pivinv)));
        }
        for (int l = 0; l < m; ++l) {
            b.set_value(icol, l, py::float_(sanitizePrecision(py::cast<double>(b.getValueAt(icol, l)) * pivinv)));
        }

        for (int ll = 0; ll < n; ++ll) {
            if (ll != icol) {
                double dum = py::cast<double>(a.getValueAt(ll, icol));
                a.set_value(ll, icol, py::float_(0.0));
                for (int l = 0; l < n; ++l) {
                    double value = py::cast<double>(a.getValueAt(ll, l)) - py::cast<double>(a.getValueAt(icol, l)) * dum;
                    a.set_value(ll, l, py::float_(value));
                }
                for (int l = 0; l < m; ++l) {
                    double value = py::cast<double>(b.getValueAt(ll, l)) - py::cast<double>(b.getValueAt(icol, l)) * dum;
                    b.set_value(ll, l, py::float_(value));
                }
            }
        }
    }

    for (int l = n - 1; l >= 0; --l) {
        if (indxr[l] != indxc[l]) {
            for (int k = 0; k < n; ++k) {
                py::object temp = a.getValueAt(k, indxr[l]);
                a.set_value(k, indxr[l], a.getValueAt(k, indxc[l]));
                a.set_value(k, indxc[l], temp);
            }
        }
    }

    return 0;  
}


int GaussJordanVec(PyMatrix& a, const int n, PyVector& b) {

    a.promoteMatrixVariantIfNeeded<double>();
    auto convertedVector = PyVector::convertToFloatIfNeeded(b.getBaseVector().get());
    b = PyVector(std::move(convertedVector));    

    std::vector<int> ipiv(n, 0);
    std::vector<int> indxr(n), indxc(n);

    for (int i = 0; i < n; ++i) {
        double big = 0.0;
        int irow = -1, icol = -1;

        for (int j = 0; j < n; ++j) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < n; ++k) {
                    if (ipiv[k] == 0) {
                        double aij = std::abs(py::cast<double>(a.getValueAt(j, k)));
                        if (aij > big) {
                            big = aij;
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1) {
                        return Numerics_message("Singular Matrix 1");
                    }
                }
            }
        }

        if (irow < 0 || icol < 0) return Numerics_message("No valid pivot found");

        ++(ipiv[icol]);

        if (irow != icol) {
            for (int l = 0; l < n; ++l) {
                py::object temp = a.getValueAt(irow, l);
                a.set_value(irow, l, a.getValueAt(icol, l));
                a.set_value(icol, l, temp);
            }
            py::object temp = b.__getitem__(irow);
            b.__setitem__(irow, b.__getitem__(icol));
            b.__setitem__(icol, temp);
        }

        double piv = py::cast<double>(a.getValueAt(icol, icol));
        if (std::abs(piv) < 1e-10) {
            return Numerics_message("Matrix is singular or nearly singular");
        }

        double pivinv = 1.0 / piv;
        if (std::abs(pivinv) < 1e-10) {
            pivinv = 0.0;  
        }

        a.set_value(icol, icol, py::float_(1.0));
        b.__setitem__(icol, py::float_(sanitizePrecision(py::cast<double>(b.__getitem__(icol)) * pivinv)));
        for (int l = 0; l < n; ++l) {
            a.set_value(icol, l, py::float_(sanitizePrecision(py::cast<double>(a.getValueAt(icol, l)) * pivinv)));
        }

        for (int ll = 0; ll < n; ++ll) {
            if (ll != icol) {
                double dum = py::cast<double>(a.getValueAt(ll, icol));
                a.set_value(ll, icol, py::float_(0.0));
                double b_val = sanitizePrecision(py::cast<double>(b.__getitem__(ll)) - py::cast<double>(b.__getitem__(icol)) * dum);
                b.__setitem__(ll, py::float_(b_val));
                for (int l = 0; l < n; ++l) {
                    double value1 = sanitizePrecision(py::cast<double>(a.getValueAt(ll, l)) - py::cast<double>(a.getValueAt(icol, l)) * dum);
                    a.set_value(ll, l, py::float_(value1));
                }
            }
        }
    }

    for (int l = n - 1; l >= 0; --l) {
        if (indxr[l] != indxc[l]) {
            for (int k = 0; k < n; ++k) {
                py::object temp = a.getValueAt(k, indxr[l]);
                a.set_value(k, indxr[l], a.getValueAt(k, indxc[l]));
                a.set_value(k, indxc[l], temp);
            }
        }
    }
    return 0;
}

int GaussBack(PyMatrix& a, const int n, PyVector& b) {

	
	a.promoteMatrixVariantIfNeeded<double>();
	auto convertedVector = PyVector::convertToFloatIfNeeded(b.getBaseVector().get());
	b = PyVector(std::move(convertedVector));
    int i, irow, j, l;
    double pivinv, temp;

    for (i = 0; i < n; i++) {
        irow = i;
        for (j = i + 1; j < n; j++) {
            if (std::abs(py::cast<double>(a.getValueAt(j, i))) > std::abs(py::cast<double>(a.getValueAt(irow, i)))) {
                irow = j;
            }
        }
        if (py::cast<double>(a.getValueAt(irow, i)) == 0.0) {
            return Numerics_message("GaussBack: Singular Matrix");
        }

        if (irow != i) {
            for (l = i; l < n; l++) {
                py::object temp = a.getValueAt(i, l);
                a.set_value(i, l, a.getValueAt(irow, l));
                a.set_value(irow, l, temp);
            }
            py::object temp = b.__getitem__(i);
            b.__setitem__(i, b.__getitem__(irow));
            b.__setitem__(irow, temp);
        }

        pivinv = 1.0 / py::cast<double>(a.getValueAt(i, i));
        a.set_value(i, i, py::float_(1.0));
        b.__setitem__(i, py::float_(py::cast<double>(b.__getitem__(i)) * pivinv));
        for (l = i + 1; l < n; l++) {
            a.set_value(i, l, py::float_(py::cast<double>(a.getValueAt(i, l)) * pivinv));
        }

        for (j = i + 1; j < n; j++) {
            if (py::cast<double>(a.getValueAt(j, i)) != 0.0) {
                temp = py::cast<double>(a.getValueAt(j, i));
                a.set_value(j, i, py::float_(0.0));
                b.__setitem__(j, py::float_(py::cast<double>(b.__getitem__(j)) - temp * py::cast<double>(b.__getitem__(i))));
                for (l = i + 1; l < n; l++) {
                    a.set_value(j, l, py::float_(py::cast<double>(a.getValueAt(j, l)) - temp * py::cast<double>(a.getValueAt(i, l))));
                }
            }
        }
    }

    for (i = n - 2; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            b.__setitem__(i, py::float_(py::cast<double>(b.__getitem__(i)) - py::cast<double>(a.getValueAt(i, j)) * py::cast<double>(b.__getitem__(j))));
        }
    }

    return 0;
}

double qbulir(std::function<double(double)> func, double a, double b, double eps_, double& err) {
    double ba = b - a;
    if (ba == 0.0) return 0.0;

    int n = 2, nn = 3, mx = 25, m, mr, bo, bu = 0, odd = 1;
    double c, d1, ddt, den, e, eps, eta = std::numeric_limits<double>::epsilon() * 10, gr, hm, nt, sm, t, t1, t2, t2a, ta, tab = 0.0, tb, v = 0.0, w;
    double d[7], dt[7];

    eps = std::max(eps_, eta);
    sm = 0.0;
    gr = 0.0;
    t1 = 0.0;
    t2 = 0.5 * (func(a) + func(b));
    t2a = t2;
    tb = std::fabs(t2a);
    c = t2 * ba;
    dt[0] = c;

    for (m = 1; m <= mx; ++m) {  
        bo = (m >= 7);
        hm = ba / n;
        if (odd) {
            for (int i = 1; i <= n; i += 2) {
                w = func(a + i * hm);
                t2 += w;
                tb += std::fabs(w);
            }
            nt = t2;
            tab = tb * std::fabs(hm);
            d[1] = 16.0 / 9.0;
            d[3] = 64.0 / 9.0;
            d[5] = 256.0 / 9.0;
        } else {
            for (int i = 1; i <= n; i += 6) {
                w = i * hm;
                t1 += func(a + w) + func(b - w);
            }
            nt = t1 + t2a;
            t2a = t2;
            d[1] = 9.0 / 4.0;
            d[3] = 9.0;
            d[5] = 36.0;
        }
        ddt = dt[0];
        t = nt * hm;
        dt[0] = t;
        nt = dt[0];
        if (bo) {
            mr = 6;
            d[6] = 64.0;
            w = 144.0;
        } else {
            mr = m;
            d[m] = n * n;
            w = d[m];
        }
        for (int i = 1; i <= mr; ++i) {
            d1 = d[i] * ddt;
            den = d1 - nt;
            e = nt - ddt;
            if (den != 0.0) {
                e /= den;
                v = nt * e;
                nt = d1 * e;
                t += v;
            } else {
                nt = 0.0;
                v = 0.0;
            }
            ddt = dt[i];
            dt[i] = v;
        }
        ta = c;
        c = t;
        if (!bo) t -= v;
        v = t - ta;
        t += v;
        err = std::fabs(v);
        if (ta < t) {
            d1 = ta;
            ta = t;
            t = d1;
        }
        bo = bo || (ta < gr && t > sm);
        if (bu && bo && err < eps * tab * w) break;
        gr = ta;
        sm = t;
        odd = !odd;
        int i = n;
        n = nn;
        nn = i + i;
        bu = bo;
        d[2] = 4.0;
        d[4] = 16.0;
    }
    v = tab * eta;
    if (err < v) err = v;
    if (m == mx) Numerics_message("qbulir exceeding maximum of iterations");

    c = sanitizePrecision(c);
    return roundIfCloseToInteger(c);
}

void GaussLegendre(PyVector& x, PyVector& w, const int n) {
    double eps = 1e-10;

    auto x_converted = PyVector::convertToFloatIfNeeded(x.getBaseVector().get());
    auto w_converted = PyVector::convertToFloatIfNeeded(w.getBaseVector().get());
    x = PyVector(std::move(x_converted));
    w = PyVector(std::move(w_converted));

    int m = (n + 1) / 2;
    double z1, z, pp, p3, p2, p1;

    for (int i = 0; i < m; ++i) {
        z = cos(M_PI * (i + 0.75) / (n + 0.5));
        do {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 0; j < n; ++j) {
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j + 1) * z * p2 - j * p3) / double(j + 1);
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > eps);

        x.__setitem__(i, py::float_(-z));
        x.__setitem__(n - 1 - i, py::float_(z));
        w.__setitem__(i, py::float_(2.0 / ((1.0 - z * z) * pp * pp)));
        w.__setitem__(n - 1 - i, w.__getitem__(i));
    }
}

void LegendrePeven(PyVector& p, const double x, const int np) {
    double x2 = x * x;

    auto p_converted = PyVector::convertToFloatIfNeeded(p.getBaseVector().get());
    p = PyVector(std::move(p_converted));

    p.__setitem__(0, py::float_(1.0));
    if (np > 1) {
        p.__setitem__(1, py::float_(1.5 * x2 - 0.5));

        for (int n = 2; n < np; ++n) {
            int l = 2 * (n - 1);
            int l2 = 2 * l;

            double term1 = -py::cast<double>(p.__getitem__(n - 2)) * l * (l - 1) / double((l2 + 1) * (l2 - 1));
            double term2 = py::cast<double>(p.__getitem__(n - 1)) * (x2 - (l2 * l + l2 - 1) / double((l2 - 1) * (l2 + 3)));

            double result = (term1 + term2) * (l2 + 1) * (l2 + 3) / double((l + 1) * (l + 2));
            p.__setitem__(n, py::float_(result));
        }
    }
}

void dLegendrePeven(PyVector& p, PyVector& d, const double x, const int np) {
    double x2 = x * x;

    auto p_converted = PyVector::convertToFloatIfNeeded(p.getBaseVector().get());
    auto d_converted = PyVector::convertToFloatIfNeeded(d.getBaseVector().get());
    p = PyVector(std::move(p_converted));
    d = PyVector(std::move(d_converted));

    p.__setitem__(0, py::float_(1.0));
    d.__setitem__(0, py::float_(0.0));
    if (np > 1) {
        p.__setitem__(1, py::float_(1.5 * x2 - 0.5));
        d.__setitem__(1, py::float_(1.5));

        for (int n = 2; n < np; ++n) {
            int l = 2 * (n - 1);
            int l2 = 2 * l;

            double p_term1 = -py::cast<double>(p.__getitem__(n - 2)) * l * (l - 1) / double((l2 + 1) * (l2 - 1));
            double p_term2 = py::cast<double>(p.__getitem__(n - 1)) * (x2 - (l2 * l + l2 - 1) / double((l2 - 1) * (l2 + 3)));

            double p_result = (p_term1 + p_term2) * (l2 + 1) * (l2 + 3) / double((l + 1) * (l + 2));
            p.__setitem__(n, py::float_(p_result));

            double d_term1 = -py::cast<double>(d.__getitem__(n - 2)) * l * (l - 1) / double((l2 + 1) * (l2 - 1));
            double d_term2 = py::cast<double>(d.__getitem__(n - 1)) * (x2 - (l2 * l + l2 - 1) / double((l2 - 1) * (l2 + 3)));
            double d_term3 = py::cast<double>(p.__getitem__(n - 1));

            double d_result = (d_term1 + d_term2 + d_term3) * (l2 + 1) * (l2 + 3) / double((l + 1) * (l + 2));
            d.__setitem__(n, py::float_(d_result));
        }
    }

    x2 = 2 * x;
    for (int n = 0; n < np; ++n) {
        d.__setitem__(n, py::float_(py::cast<double>(d.__getitem__(n)) * x2));
    }
}

double zbrent(std::function<double(double)> func, const double x1, const double x2, const double tol) {
    int iter = 0, itmax = 100;
    double a, b, c = 0., d = 0., e = 0., fa, fb, fc, tol1, tolh = 0.5 * tol, eps = 6.e-8, s, p, q, r, xm;

    fa = func(a = x1);
    fb = func(b = x2);
    if ((fa > 0. && fb > 0.) || (fa < 0. && fb < 0.))
        throw std::runtime_error("zbrent: root must be bracketed");

    fc = fb;
    while (iter++ < itmax) {
        if ((fc > 0. && fb > 0.) || (fc < 0. && fb < 0.)) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        tol1 = eps * fabs(b) + tolh;
        xm = 0.5 * (c - b);
        if (fabs(xm) < tol1 || fb == 0.) return b;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s = fb / fa;
            if (a == c) {
                p = 2. * xm * s;
                q = 1. - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.));
                q = (q - 1.) * (r - 1.) * (s - 1.);
            }
            if (p > 0.) q = -q;
            else p = -p;
            if (2. * p < std::min(3. * xm * q - fabs(tol1 * q), fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (fabs(d) > tol1) b += d;
        else b += (xm < 0.) ? -tol1 : tol1;
        fb = func(b);
    }
    throw std::runtime_error("zbrent exceeding iterations");
    return b;
}

void int_heap_index(PyVector& A, const int n, PyVector& indx) {
    int l, j, ir, indxt, i;
    double q;

    for (j = 0; j < n; j++) {
        indx.__setitem__(j, py::int_(j));
    }

    l = n >> 1;
    ir = n - 1;

    for (;;) {
        if (l > 0) {
            q = py::cast<double>(A.__getitem__(indxt = py::cast<int>(indx.__getitem__(--l))));
        } else {
            q = py::cast<double>(A.__getitem__(indxt = py::cast<int>(indx.__getitem__(ir))));
            indx.__setitem__(ir, indx.__getitem__(0));
            if (--ir == 0) {
                indx.__setitem__(0, py::int_(indxt));
                return;
            }
        }
        i = l;
        j = (l << 1) + 1;
        while (j <= ir) {
            if (j < ir && py::cast<double>(A.__getitem__(py::cast<int>(indx.__getitem__(j)))) < py::cast<double>(A.__getitem__(py::cast<int>(indx.__getitem__(j + 1))))) {
                j++;
            }
            if (q < py::cast<double>(A.__getitem__(py::cast<int>(indx.__getitem__(j))))) {
                indx.__setitem__(i, indx.__getitem__(j));
                j += 1 + (i = j);
            } else {
                j = ir + 1;
            }
        }
        indx.__setitem__(i, py::int_(indxt));
    }
}

void int_func_heap_index(std::function<int(int)> func, const int n, PyVector& indx) {
    int l, j, ir, indxt, i;
    int q;

    for (j = 0; j < n; j++) {
        indx.__setitem__(j, py::int_(j));
    }

    l = n >> 1;
    ir = n - 1;

    for (;;) {
        if (l > 0) {
            q = func(indxt = py::cast<int>(indx.__getitem__(--l)));
        } else {
            q = func(indxt = py::cast<int>(indx.__getitem__(ir)));
            indx.__setitem__(ir, indx.__getitem__(0));
            if (--ir == 0) {
                indx.__setitem__(0, py::int_(indxt));
                return;
            }
        }
        i = l;
        j = (l << 1) + 1;
        while (j <= ir) {
            if (j < ir && func(py::cast<int>(indx.__getitem__(j))) < func(py::cast<int>(indx.__getitem__(j + 1)))) {
                j++;
            }
            if (q < func(py::cast<int>(indx.__getitem__(j)))) {
                indx.__setitem__(i, indx.__getitem__(j));
                j += 1 + (i = j);
            } else {
                j = ir + 1;
            }
        }
        indx.__setitem__(i, py::int_(indxt));
    }
}

void double_func_heap_index(std::function<double(const int)> func, const int n, PyVector& indx) {
    int l, j, ir, indxt, i;
    double q;

    for (j = 0; j < n; j++) {
        indx.__setitem__(j, py::int_(j));
    }

    l = n >> 1;
    ir = n - 1;

    for (;;) {
        if (l > 0) {
            q = func(indxt = py::cast<int>(indx.__getitem__(--l)));
        } else {
            q = func(indxt = py::cast<int>(indx.__getitem__(ir)));
            indx.__setitem__(ir, indx.__getitem__(0));
            if (--ir == 0) {
                indx.__setitem__(0, py::int_(indxt));
                return;
            }
        }
        i = l;
        j = (l << 1) + 1;
        while (j <= ir) {
            if (j < ir && func(py::cast<int>(indx.__getitem__(j))) < func(py::cast<int>(indx.__getitem__(j + 1)))) {
                j++;
            }
            if (q < func(py::cast<int>(indx.__getitem__(j)))) {
                indx.__setitem__(i, indx.__getitem__(j));
                j += 1 + (i = j);
            } else {
                j = ir + 1;
            }
        }
        indx.__setitem__(i, py::int_(indxt));
    }
}

int hunt(const PyVector& x, const int n, const double xx, int jlo) {
    int jm, jhi, inc = 1;
    bool ascend = (x.__getitem__(n - 1).cast<double>() >= x.__getitem__(0).cast<double>());
    if (jlo < 0 || jlo > n - 1) {
        jlo = 0;
        jhi = n - 1;
    } else {
        if ((ascend && xx >= x.__getitem__(jlo).cast<double>()) || (!ascend && xx <= x.__getitem__(jlo).cast<double>())) {
            for (;;) {
                jhi = jlo + inc;
                if (jhi >= n) {
                    jhi = n;
                    break;
                } else if ((ascend && xx < x.__getitem__(jhi).cast<double>()) || (!ascend && xx > x.__getitem__(jhi).cast<double>())) {
                    break;
                } else {
                    jlo = jhi;
                    inc += inc;
                }
            }
        } else {
            jhi = jlo;
            for (;;) {
                jlo = jhi - inc;
                if (jlo <= 0) {
                    jlo = 0;
                    break;
                } else if ((ascend && xx >= x.__getitem__(jlo).cast<double>()) || (!ascend && xx <= x.__getitem__(jlo).cast<double>())) {
                    break;
                } else {
                    jhi = jlo;
                    inc += inc;
                }
            }
        }
    }
    while (jhi - jlo != 1) {
        jm = (jhi + jlo) >> 1;
        if ((ascend && xx >= x.__getitem__(jm).cast<double>()) || (!ascend && xx <= x.__getitem__(jm).cast<double>())) {
            jlo = jm;
        } else {
            jhi = jm;
        }
    }
    return std::max(0, std::min(n - 2, jlo));
}

double qsplin(PyVector& x, PyVector& y, PyVector& y2, const int n, const double al, const double x1, const double x2) {
    if (x1 == x2) return 0.0;
    if ((x.__getitem__(n-1).cast<double>() - x1) * (x1 - x.__getitem__(0).cast<double>()) < 0.0) Numerics_message("qsplin: x1 not in range");
    if ((x.__getitem__(n-1).cast<double>() - x2) * (x2 - x.__getitem__(0).cast<double>()) < 0.0) Numerics_message("qsplin: x2 not in range");
    if (x1 > x2) Numerics_message("qsplin: x1 > x2");
    if (x1 < 0.0 && al < 0.0) Numerics_message("qsplin: integral complex");
    if (x1 == 0.0 && al <= -1.0) Numerics_message("qsplin: integral diverging");

    int i, k;
    double q = 0.0, ali, ali1, h, h2, t, xl, xh;
    double a[4];

    k = hunt(x, n, x1, int((x1 - x.__getitem__(0).cast<double>()) / (x.__getitem__(n-1).cast<double>() - x.__getitem__(0).cast<double>()) * (n - 1)));

    while (x.__getitem__(k).cast<double>() < x2) {
        xl = std::max(x.__getitem__(k).cast<double>(), x1);
        xh = std::min(x.__getitem__(k+1).cast<double>(), x2);
        h = x.__getitem__(k+1).cast<double>() - x.__getitem__(k).cast<double>();
        if (h == 0.0) Numerics_message("qsplin: bad x input");
        h2 = h * h;
        a[0] = x.__getitem__(k+1).cast<double>() * y.__getitem__(k).cast<double>() - x.__getitem__(k).cast<double>() * y.__getitem__(k+1).cast<double>() + x.__getitem__(k).cast<double>() * x.__getitem__(k+1).cast<double>() / 6.0 *
               ((x.__getitem__(k+1).cast<double>() + h) * y2.__getitem__(k).cast<double>() - (x.__getitem__(k).cast<double>() - h) * y2.__getitem__(k+1).cast<double>());
        a[1] = y.__getitem__(k+1).cast<double>() - y.__getitem__(k).cast<double>() + ((h2 - 3.0 * x.__getitem__(k+1).cast<double>() * x.__getitem__(k+1).cast<double>()) * y2.__getitem__(k).cast<double>()
               - (h2 - 3.0 * x.__getitem__(k).cast<double>() * x.__getitem__(k).cast<double>()) * y2.__getitem__(k+1).cast<double>()) / 6.0;
        a[2] = 0.5 * (x.__getitem__(k+1).cast<double>() * y2.__getitem__(k).cast<double>() - x.__getitem__(k).cast<double>() * y2.__getitem__(k+1).cast<double>());
        a[3] = (y2.__getitem__(k+1).cast<double>() - y2.__getitem__(k).cast<double>()) / 6.0;

        for (i = 0, ali = al, t = 0.0; i < 4; i++, ali += 1.0) {
            if (ali == -1.0)
                t += a[i] * log(xh / xl);
            else {
                ali1 = ali + 1.0;
                t += a[i] * (pow(xh, ali1) - pow(xl, ali1)) / ali1;
            }
        }
        q += t / h;
        k++;
    }
    return q;
}

int CholeskyDecomposition(PyMatrix& a, const int n) {
    int i, j, k;
    double sum;

    a.promoteMatrixVariantIfNeeded<double>();

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            sum = py::cast<double>(a.getValueAt(i, j));
            for (k = i - 1; k >= 0; k--) {
                sum -= py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(a.getValueAt(j, k));
            }
            if (i == j) {
                if (sum <= 0.0) {
                    return Numerics_message("CholeskyDecomposition: Matrix not positive definite");
                }
                a.set_value(i, i, py::float_(sqrt(sum)));
            } else {
                a.set_value(j, i, py::float_(sum / py::cast<double>(a.getValueAt(i, i))));
                a.set_value(i, j, py::float_(0.0));  
            }
        }
    }
    return 0;
}

void CholeskySolution(PyMatrix& a, const int n, PyVector& b) {
    int i, k;
    double sum;

    a.promoteMatrixVariantIfNeeded<double>();
    auto convertedVector = PyVector::convertToFloatIfNeeded(b.getBaseVector().get());
    b = PyVector(std::move(convertedVector));

    for (i = 0; i < n; i++) {
        sum = py::cast<double>(b.__getitem__(i));
        for (k = 0; k < i; k++) {
            sum -= py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(b.__getitem__(k));
        }
        b.__setitem__(i, py::float_(sum / py::cast<double>(a.getValueAt(i, i))));
    }

    for (i = n - 1; i >= 0; i--) {
        sum = py::cast<double>(b.__getitem__(i));
        for (k = i + 1; k < n; k++) {
            sum -= py::cast<double>(a.getValueAt(k, i)) * py::cast<double>(b.__getitem__(k));
        }
        b.__setitem__(i, py::float_(sum / py::cast<double>(a.getValueAt(i, i))));
    }
}
void CholeskyInvertL(PyMatrix& a, const int n) {
    int i, j, k;
    double sum;

    a.promoteMatrixVariantIfNeeded<double>();

    for (i = 0; i < n; i++) {
        a.set_value(i, i, py::float_(1.0 / py::cast<double>(a.getValueAt(i, i))));
        for (j = i + 1; j < n; j++) {
            sum = 0.0;
            for (k = i; k < j; k++) {
                sum -= py::cast<double>(a.getValueAt(j, k)) * py::cast<double>(a.getValueAt(k, i));
            }
            a.set_value(j, i, py::float_(sum / py::cast<double>(a.getValueAt(j, j))));
        }
    }
}


void CholeskyInvertF(PyMatrix& a, const int n) {
    int i, j, k;
    double sum;

    a.promoteMatrixVariantIfNeeded<double>();

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) { 
            sum = 0.0;
            for (k = j; k < n; k++) { 
                sum += py::cast<double>(a.getValueAt(k, i)) * py::cast<double>(a.getValueAt(k, j));
            }
            a.set_value(i, j, py::float_(sum));
        }
    }

    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            a.set_value(j, i, py::float_(py::cast<double>(a.getValueAt(i, j))));
        }
    }
}

int LUDecomposition(PyMatrix& a, const int n, PyVector& indx, int& d) {
    const double tiny = 1.e-5;
    int i, imax, j, k;
    double big, dum, sum;

    a.promoteMatrixVariantIfNeeded<double>();

    std::vector<double> vv(n);

    d = 1;

    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            big = std::max(big, std::abs(py::cast<double>(a.getValueAt(i, j))));
        }
        if (big == 0.0) {
            return Numerics_message("LUDecomposition: singular matrix");
        }
        vv[i] = 1.0 / big;
    }

    for (j = 0; j < n; j++) {
        imax = j;

        for (i = 0; i < j; i++) {
            sum = py::cast<double>(a.getValueAt(i, j));
            for (k = 0; k < i; k++) {
                sum -= py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(a.getValueAt(k, j));
            }
            a.set_value(i, j, py::float_(sum));
        }

        big = 0.0;
        for (i = j; i < n; i++) {
            sum = py::cast<double>(a.getValueAt(i, j));
            for (k = 0; k < j; k++) {
                sum -= py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(a.getValueAt(k, j));
            }
            a.set_value(i, j, py::float_(sum));
            if ((dum = vv[i] * std::fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) {
            for (k = 0; k < n; k++) {
                dum = py::cast<double>(a.getValueAt(imax, k));
                a.set_value(imax, k, a.getValueAt(j, k));
                a.set_value(j, k, py::float_(dum));
            }
            d = -d;
            vv[imax] = vv[j];
        }

        indx.__setitem__(j, py::int_(imax));
        if (py::cast<double>(a.getValueAt(j, j)) == 0.0) {
            a.set_value(j, j, py::float_(tiny));
        }

        if (j < n - 1) {
            dum = 1.0 / py::cast<double>(a.getValueAt(j, j));
            for (i = j + 1; i < n; i++) {
                a.set_value(i, j, py::float_(py::cast<double>(a.getValueAt(i, j)) * dum));
            }
        }
    }
    return 0;
}

double LUDet3(PyMatrix& A) {
    const int n = 3;
    int d;
    PyVector indx(py::make_tuple(0, 0, 0)); 
    
    A.promoteMatrixVariantIfNeeded<double>();

    PyMatrix AA(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            AA.set_value(i, j, A.getValueAt(i, j));
        }
    }

    LUDecomposition(AA, n, indx, d);

    double D = static_cast<double>(d);
		for (int i = 0; i < n; ++i) {
        D *= py::cast<double>(AA.getValueAt(i, i));
    }

    return roundIfCloseToInteger(D);
}

void LUSolution(PyMatrix& a, const int n, PyVector& indx, PyVector& b) {
    int i, ii = -1, ip, j;
    double sum;
	auto convertedVector = PyVector::convertToFloatIfNeeded(b.getBaseVector().get());
	b = PyVector(std::move(convertedVector));
    a.promoteMatrixVariantIfNeeded<double>();
    for (i = 0; i < n; i++) {
        ip = py::cast<int>(indx.__getitem__(i));
        sum = py::cast<double>(b.__getitem__(ip));
        b.__setitem__(ip, b.__getitem__(i));
        if (ii >= 0) {
            for (j = ii; j < i; j++) {
                sum -= py::cast<double>(a.getValueAt(i, j)) * py::cast<double>(b.__getitem__(j));
            }
        } else if (sum != 0.0) {
            ii = i;
        }
        b.__setitem__(i, py::float_(sum));
    }

    for (i = n - 1; i >= 0; i--) {
        sum = py::cast<double>(b.__getitem__(i));
        for (j = i + 1; j < n; j++) {
            sum -= py::cast<double>(a.getValueAt(i, j)) * py::cast<double>(b.__getitem__(j));
        }
        b.__setitem__(i, py::float_(sum / py::cast<double>(a.getValueAt(i, i))));
    }
}

void LUInvert(PyMatrix& a, PyMatrix& y, const int n, PyVector& indx) {
    int i, j;
    std::vector<double> col_data(n, 0.0);  
    PyVector col(py::cast(col_data));
    y.promoteMatrixVariantIfNeeded<double>();

    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            col.__setitem__(i, py::float_(0.0));
        }
        col.__setitem__(j, py::float_(1.0));
        LUSolution(a, n, indx, col);
        for (i = 0; i < n; i++) {
            y.set_value(i, j, col.__getitem__(i));
        }
    }
}

void tred2(PyMatrix& a, const int n, PyVector& d, PyVector& e, const char EV) {
    int l, k, j, i;
    double scale, hh, h, g, f;

    a.promoteMatrixVariantIfNeeded<double>();
	auto convertedDVector = PyVector::convertToFloatIfNeeded(d.getBaseVector().get());
	auto convertedEVector = PyVector::convertToFloatIfNeeded(e.getBaseVector().get());
	d = PyVector(std::move(convertedDVector));
	e = PyVector(std::move(convertedEVector));
	

    for (i = n - 1; i > 0; i--) {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0) {
            for (k = 0; k <= l; k++) {
                scale += fabs(py::cast<double>(a.getValueAt(i, k)));
            }
            if (scale == 0.0) {
                e.__setitem__(i, a.getValueAt(i, l));
            } else {
                for (k = 0; k <= l; k++) {
                    a.set_value(i, k, py::float_(py::cast<double>(a.getValueAt(i, k)) / scale));
                    h += py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(a.getValueAt(i, k));
                }
                f = py::cast<double>(a.getValueAt(i, l));
                g = (f >= 0.0) ? -sqrt(h) : sqrt(h);
                e.__setitem__(i, py::float_(scale * g));
                h -= f * g;
                a.set_value(i, l, py::float_(f - g));
                f = 0.0;
                for (j = 0; j <= l; j++) {
                    if (EV) a.set_value(j, i, py::float_(py::cast<double>(a.getValueAt(i, j)) / h));
                    g = 0.0;
                    for (k = 0; k <= j; k++) {
                        g += py::cast<double>(a.getValueAt(j, k)) * py::cast<double>(a.getValueAt(i, k));
                    }
                    for (k = j + 1; k <= l; k++) {
                        g += py::cast<double>(a.getValueAt(k, j)) * py::cast<double>(a.getValueAt(i, k));
                    }
                    e.__setitem__(j, py::float_(g / h));
                    f += py::cast<double>(e.__getitem__(j)) * py::cast<double>(a.getValueAt(i, j));
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++) {
                    f = py::cast<double>(a.getValueAt(i, j));
                    e.__setitem__(j, py::float_(py::cast<double>(e.__getitem__(j)) - hh * f));
                    g = py::cast<double>(e.__getitem__(j));
                    for (k = 0; k <= j; k++) {
                        a.set_value(j, k, py::float_(py::cast<double>(a.getValueAt(j, k)) - (f * py::cast<double>(e.__getitem__(k)) + g * py::cast<double>(a.getValueAt(i, k)))));
                    }
                }
            }
        } else {
            e.__setitem__(i, a.getValueAt(i, l));
        }
        d.__setitem__(i, py::float_(h));
    }
    d.__setitem__(0, py::float_(0.0));
    e.__setitem__(0, py::float_(0.0));

    if (EV) {
        for (i = 0; i < n; i++) {
            l = i - 1;
            if (py::cast<double>(d.__getitem__(i)) != 0.0) {
                for (j = 0; j <= l; j++) {
                    g = 0.0;
                    for (k = 0; k <= l; k++) {
                        g += py::cast<double>(a.getValueAt(i, k)) * py::cast<double>(a.getValueAt(k, j));
                    }
                    for (k = 0; k <= l; k++) {
                        a.set_value(k, j, py::float_(py::cast<double>(a.getValueAt(k, j)) - g * py::cast<double>(a.getValueAt(k, i))));
                    }
                }
            }
            d.__setitem__(i, a.getValueAt(i, i));
            a.set_value(i, i, py::float_(1.0));
            for (j = 0; j <= l; j++) {
                a.set_value(j, i, py::float_(0.0));
                a.set_value(i, j, py::float_(0.0));
            }
        }
    } else {
        for (i = 0; i < n; i++) {
            d.__setitem__(i, a.getValueAt(i, i));
        }
    }
}

void tqli(PyVector& d, PyVector& e, const int n, PyMatrix& z, const char EV) {
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

    auto convertedDVector = PyVector::convertToFloatIfNeeded(d.getBaseVector().get());
	auto convertedEVector = PyVector::convertToFloatIfNeeded(e.getBaseVector().get());
	d = PyVector(std::move(convertedDVector));
	e = PyVector(std::move(convertedEVector));

	z.promoteMatrixVariantIfNeeded<double>();

    for (i = 1; i < n; i++) {
        e.__setitem__(i - 1, e.__getitem__(i));
    }
    e.__setitem__(n - 1, py::float_(0.0));

    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n - 1; m++) {
                dd = fabs(py::cast<double>(d.__getitem__(m))) + fabs(py::cast<double>(d.__getitem__(m + 1)));
                if (fabs(py::cast<double>(e.__getitem__(m))) + dd == dd) break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    throw std::runtime_error("tqli: too many iterations");
                }
                g = (py::cast<double>(d.__getitem__(l + 1)) - py::cast<double>(d.__getitem__(l))) / (2.0 * py::cast<double>(e.__getitem__(l)));
                r = (g == 0.0) ? 1.0 : hypot(g, 1.0);
                g = py::cast<double>(d.__getitem__(m)) - py::cast<double>(d.__getitem__(l)) + py::cast<double>(e.__getitem__(l)) / (g + (g >= 0.0 ? fabs(r) : -fabs(r)));
                s = c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--) {
                    f = s * py::cast<double>(e.__getitem__(i));
                    b = c * py::cast<double>(e.__getitem__(i));
                    e.__setitem__(i + 1, py::float_(r = hypot(f, g)));
                    if (r == 0.0) {
                        d.__setitem__(i + 1, py::float_(py::cast<double>(d.__getitem__(i + 1)) - p));
                        e.__setitem__(m, py::float_(0.0));
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = py::cast<double>(d.__getitem__(i + 1)) - p;
                    r = (py::cast<double>(d.__getitem__(i)) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d.__setitem__(i + 1, py::float_(g + p));
                    g = c * r - b;
                    if (EV) {
                        for (k = 0; k < n; k++) {
                            f = py::cast<double>(z.getValueAt(k, i + 1));
                            z.set_value(k, i + 1, py::float_(s * py::cast<double>(z.getValueAt(k, i)) + c * f));
                            z.set_value(k, i, py::float_(c * py::cast<double>(z.getValueAt(k, i)) - s * f));
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d.__setitem__(l, py::float_(py::cast<double>(d.__getitem__(l)) - p));
                e.__setitem__(l, py::float_(g));
                e.__setitem__(m, py::float_(0.0));
            }
        } while (m != l);
    }
}

void balanc(PyMatrix& a, const int n) {
    const double radix = 2.0, sqrdx = radix * radix;
    int last = 0, j, i;
    double s, r, g, f, c;

    a.promoteMatrixVariantIfNeeded<double>();

    while (last == 0) {
        last = 1;
        for (i = 0; i < n; i++) {
            r = c = 0.0;
            for (j = 0; j < n; j++) {
                if (j != i) {
                    c += fabs(py::cast<double>(a.getValueAt(j, i)));
                    r += fabs(py::cast<double>(a.getValueAt(i, j)));
                }
            }
            if (c != 0.0 && r != 0.0) {
                g = r / radix;
                f = 1.0;
                s = c + r;
                while (c < g) {
                    f *= radix;
                    c *= sqrdx;
                }
                g = r * radix;
                while (c > g) {
                    f /= radix;
                    c /= sqrdx;
                }
                if ((c + r) / f < 0.95 * s) {
                    last = 0;
                    g = 1.0 / f;
                    for (j = 0; j < n; j++) {
                        a.set_value(i, j, py::float_(py::cast<double>(a.getValueAt(i, j)) * g));
                        a.set_value(j, i, py::float_(py::cast<double>(a.getValueAt(j, i)) * f));
                    }
                }
            }
        }
    }
}

void elmhes(PyMatrix& a, const int n) {
    int m, j, i;
    double y, x;

    a.promoteMatrixVariantIfNeeded<double>();


    for (m = 1; m < n - 1; m++) {
        x = 0.0;
        i = m;
        for (j = m; j < n; j++) {
            if (fabs(py::cast<double>(a.getValueAt(j, m - 1))) > fabs(x)) {
                x = py::cast<double>(a.getValueAt(j, m - 1));
                i = j;
            }
        }
        if (i != m) {
            for (j = m - 1; j < n; j++) {
                double temp = py::cast<double>(a.getValueAt(i, j));
                a.set_value(i, j, a.getValueAt(m, j));
                a.set_value(m, j, py::float_(temp));
            }
            for (j = 0; j < n; j++) {
                double temp = py::cast<double>(a.getValueAt(j, i));
                a.set_value(j, i, a.getValueAt(j, m));
                a.set_value(j, m, py::float_(temp));
            }
        }
        if (x != 0.0) {
            for (i = m + 1; i < n; i++) {
                y = py::cast<double>(a.getValueAt(i, m - 1));
                if (y != 0.0) {
                    y /= x;
                    a.set_value(i, m - 1, py::float_(y));
                    for (j = m; j < n; j++) {
                        a.set_value(i, j, py::float_(py::cast<double>(a.getValueAt(i, j)) - y * py::cast<double>(a.getValueAt(m, j))));
                    }
                    for (j = 0; j < n; j++) {
                        a.set_value(j, m, py::float_(py::cast<double>(a.getValueAt(j, m)) + y * py::cast<double>(a.getValueAt(j, i))));
                    }
                }
            }
        }
    }
}

void hqr(PyMatrix& a, const int n, PyVector& wr, PyVector& wi) {
    int nn, m, l, k, j, its, i, mmin;
    double z = 0.0, y = 0.0, x = 0.0, w = 0.0, v = 0.0, u = 0.0, t = 0.0, s = 0.0, r = 0.0, q = 0.0, p = 0.0, anrm;

    a.promoteMatrixVariantIfNeeded<double>();
	auto convertedWrVector = PyVector::convertToFloatIfNeeded(wr.getBaseVector().get());
	auto convertedWiVector = PyVector::convertToFloatIfNeeded(wi.getBaseVector().get());
	wr = PyVector(std::move(convertedWrVector));
	wi = PyVector(std::move(convertedWiVector));	


    anrm = fabs(py::cast<double>(a.getValueAt(0, 0)));
    for (i = 1; i < n; i++) {
        for (j = i - 1; j < n; j++) {
            anrm += fabs(py::cast<double>(a.getValueAt(i, j)));
        }
    }
    nn = n - 1;
    t = 0.0;
    while (nn >= 0) {
        its = 0;
        do {
            for (l = nn; l >= 1; l--) {
                s = fabs(py::cast<double>(a.getValueAt(l - 1, l - 1))) + fabs(py::cast<double>(a.getValueAt(l, l)));
                if (s == 0.0) s = anrm;
                if (fabs(py::cast<double>(a.getValueAt(l, l - 1))) + s == s) break;
            }
            x = py::cast<double>(a.getValueAt(nn, nn));
            if (l == nn) {
                wr.__setitem__(nn, py::float_(x + t));
                wi.__setitem__(nn--, py::float_(0.0));
            } else {
                y = py::cast<double>(a.getValueAt(nn - 1, nn - 1));
                w = py::cast<double>(a.getValueAt(nn, nn - 1)) * py::cast<double>(a.getValueAt(nn - 1, nn));
                if (l == (nn - 1)) {
                    p = 0.5 * (y - x);
                    q = p * p + w;
                    z = sqrt(fabs(q));
                    x += t;
                    if (q >= 0.0) {
                        z = p + ((p >= 0.0) ? fabs(z) : -fabs(z));
                        wr.__setitem__(nn - 1, py::float_(x + z));
                        wr.__setitem__(nn, py::float_(x - w / z));
                        wi.__setitem__(nn - 1, py::float_(0.0));
                        wi.__setitem__(nn, py::float_(0.0));
                    } else {
                        wr.__setitem__(nn - 1, py::float_(x + p));
                        wr.__setitem__(nn, py::float_(x + p));
                        wi.__setitem__(nn - 1, py::float_(-(fabs(z))));
                        wi.__setitem__(nn, py::float_(fabs(z)));
                    }
                    nn -= 2;
                } else {
                    if (its == 30) throw std::runtime_error("hqr: exceeding iterations");
                    if (its == 10 || its == 20) {
                        t += x;
                        for (i = 0; i <= nn; i++) a.set_value(i, i, py::float_(py::cast<double>(a.getValueAt(i, i)) - x));
                        s = fabs(py::cast<double>(a.getValueAt(nn, nn - 1))) + fabs(py::cast<double>(a.getValueAt(nn - 1, nn - 2)));
                        y = x = 0.75 * s;
                        w = -0.4375 * s * s;
                    }
                    ++its;
                    for (m = nn - 2; m >= l; m--) {
                        z = py::cast<double>(a.getValueAt(m, m));
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / py::cast<double>(a.getValueAt(m + 1, m)) + py::cast<double>(a.getValueAt(m, m + 1));
                        q = py::cast<double>(a.getValueAt(m + 1, m + 1)) - z - r - s;
                        r = py::cast<double>(a.getValueAt(m + 2, m + 1));
                        s = fabs(p) + fabs(q) + fabs(r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if (m == l) break;
                        u = fabs(py::cast<double>(a.getValueAt(m, m - 1))) * (fabs(q) + fabs(r));
                        v = fabs(p) * (fabs(py::cast<double>(a.getValueAt(m - 1, m - 1))) + fabs(z) + fabs(py::cast<double>(a.getValueAt(m + 1, m + 1))));
                        if (u + v == v) break;
                    }
                    for (i = m + 2; i <= nn; i++) {
                        a.set_value(i, i - 2, py::float_(0.0));
                        if (i != (m + 2)) a.set_value(i, i - 3, py::float_(0.0));
                    }
                    for (k = m; k <= nn - 1; k++) {
                        if (k != m) {
                            p = py::cast<double>(a.getValueAt(k, k - 1));
                            q = py::cast<double>(a.getValueAt(k + 1, k - 1));
                            r = 0.0;
                            if (k != (nn - 1)) r = py::cast<double>(a.getValueAt(k + 2, k - 1));
                            if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        if ((s = ((p >= 0.0) ? fabs(p) : -fabs(p)) * sqrt(p * p + q * q + r * r)) != 0.0) {
                            if (k == m) {
                                if (l != m) a.set_value(k, k - 1, py::float_(-py::cast<double>(a.getValueAt(k, k - 1))));
                            } else {
                                a.set_value(k, k - 1, py::float_(-s * x));
                            }
                            p += s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q /= p;
                            r /= p;
                            for (j = k; j < n; j++) {
                                p = py::cast<double>(a.getValueAt(k, j)) + q * py::cast<double>(a.getValueAt(k + 1, j));
                                if (k != (nn - 1)) {
                                    p += r * py::cast<double>(a.getValueAt(k + 2, j));
                                    a.set_value(k + 2, j, py::float_(py::cast<double>(a.getValueAt(k + 2, j)) - p * z));
                                }
                                a.set_value(k + 1, j, py::float_(py::cast<double>(a.getValueAt(k + 1, j)) - p * y));
                                a.set_value(k, j, py::float_(py::cast<double>(a.getValueAt(k, j)) - p * x));
                            }
                            mmin = std::min(nn, k + 3);
                            for (i = l; i <= mmin; i++) {
                                p = x * py::cast<double>(a.getValueAt(i, k)) + y * py::cast<double>(a.getValueAt(i, k + 1));
                                if (k != (nn - 1)) {
                                    p += z * py::cast<double>(a.getValueAt(i, k + 2));
                                    a.set_value(i, k + 2, py::float_(py::cast<double>(a.getValueAt(i, k + 2)) - p * r));
                                }
                                a.set_value(i, k + 1, py::float_(py::cast<double>(a.getValueAt(i, k + 1)) - p * q));
                                a.set_value(i, k, py::float_(py::cast<double>(a.getValueAt(i, k)) - p));
                            }
                        }
                    }
                }
            }
        } while (l < nn - 1);
    }
}


static double LevCof(PyVector& x, PyVector& y, PyVector& sig, const int N, PyVector& a, PyVector& fit, const int M, PyMatrix& A, PyVector& B, std::function<double(const double, PyVector&, PyVector&, const int)> func) {
    int i, j, k, l, mf, n;
    double si, dy, wt, cq = 0.0;
    std::vector<double> dyda_data(M, 0.0);
    PyVector dyda(py::cast(dyda_data));

    auto convertedXVector = PyVector::convertToFloatIfNeeded(x.getBaseVector().get());
    auto convertedYVector = PyVector::convertToFloatIfNeeded(y.getBaseVector().get());
    auto convertedSigVector = PyVector::convertToFloatIfNeeded(sig.getBaseVector().get());
    auto convertedaVector = PyVector::convertToFloatIfNeeded(a.getBaseVector().get());
    auto convertedFitVector = PyVector::convertToFloatIfNeeded(fit.getBaseVector().get());
    auto convertedBVector = PyVector::convertToFloatIfNeeded(B.getBaseVector().get());

    x = PyVector(std::move(convertedXVector));
    y = PyVector(std::move(convertedYVector));
    sig = PyVector(std::move(convertedSigVector));
    a = PyVector(std::move(convertedaVector));
    fit = PyVector(std::move(convertedFitVector));
    B = PyVector(std::move(convertedBVector));

    A.promoteMatrixVariantIfNeeded<double>();

    for (i = mf = 0; i < M; i++) if (py::cast<int>(fit.__getitem__(i))) mf++;
    for (j = 0; j < mf; j++) {
        B.__setitem__(j, py::float_(0.0));
        for (l = 0; l < mf; l++) {
            A.set_value(j, l, py::float_(0.0));
        }
    }

    for (n = 0; n < N; n++) {
        dy = py::cast<double>(y.__getitem__(n)) - func(py::cast<double>(x.__getitem__(n)), a, dyda, M);
        si = 1.0 / (py::cast<double>(sig.__getitem__(n)) * py::cast<double>(sig.__getitem__(n)));
        cq += pow(dy / py::cast<double>(sig.__getitem__(n)), 2);

        for (i = j = 0; i < M; i++) {
            if (py::cast<int>(fit.__getitem__(i))) {
                wt = py::cast<double>(dyda.__getitem__(i)) * si;
                for (k = l = 0; k < M; k++) {
                    if (py::cast<int>(fit.__getitem__(k))) {
                        A.set_value(j, l, py::float_(py::cast<double>(A.getValueAt(j, l)) + wt * py::cast<double>(dyda.__getitem__(k))));
                        l++;
                    }
                }
                B.__setitem__(j, py::float_(py::cast<double>(B.__getitem__(j)) + dy * wt));
                j++;
            }
        }
    }

    return cq;
}

double LevMar(PyVector& x, PyVector& y, PyVector& sig, const int N, PyVector& a, PyVector& fit, const int M, std::function<double(const double, PyVector&, PyVector&, const int)> func, const double dcmax, const int itmax) {
    int it = 0, i, j, mf;
    int mm[2];
    double dc, lam = 0.125, tm, cq, cqo;
    	
	std::vector<double> zero_mf(mf, 0.0);
	PyMatrix A(mf, mf);
    PyMatrix Ay(mf, mf);

	PyVector B(py::cast(zero_mf));
	PyVector By(py::cast(zero_mf));
	PyVector ay(py::cast(a.extractDataAs<double>()));

    for (i = mf = 0; i < M; i++) {
        if (py::cast<int>(fit.__getitem__(i))) mf++;
    }
    mm[0] = mm[1] = mf;

    cqo = LevCof(x, y, sig, N, a, fit, M, A, B, func);
    for (dc = 0., i = 0; i < mf; i++) dc += py::cast<double>(B.__getitem__(i)) * py::cast<double>(B.__getitem__(i));
    dc = sqrt(dc) / double(N);

    while (dc > dcmax && it++ < itmax) {
        tm = 1. + lam;
        for (i = 0; i < mf; i++) {
            for (j = 0; j < mf; j++) Ay.set_value(i, j, py::float_(py::cast<double>(A.getValueAt(i, j))));
            Ay.set_value(i, i, py::float_(py::cast<double>(A.getValueAt(i, i)) * tm));
            By.__setitem__(i, B.__getitem__(i));
        }
        GaussBack(Ay, mf, By);
        for (i = j = 0; i < M; i++) if (py::cast<int>(fit.__getitem__(i))) ay.__setitem__(i, py::float_(py::cast<double>(a.__getitem__(i)) + py::cast<double>(By.__getitem__(j++))));
        if (cqo > (cq = LevCof(x, y, sig, N, ay, fit, M, Ay, By, func))) {
            lam *= 0.125;
            cqo = cq;
            for (dc = 0., i = 0; i < mf; i++) {
                B.__setitem__(i, By.__getitem__(i));
                dc += py::cast<double>(B.__getitem__(i)) * py::cast<double>(B.__getitem__(i));
                for (j = 0; j < mf; j++) A.set_value(i, j, Ay.getValueAt(i, j));
            }
            dc = sqrt(dc) / double(N);
        } else
            lam *= 8;
    }

    return cqo;
}

inline double gauss_fit(const double x, PyVector& p, PyVector& df, const int D) {
    if (D != 3) Numerics_message("FitGauss: D != 3");

    double fg, tm;
    tm = (x - py::cast<double>(p.__getitem__(1))) / py::cast<double>(p.__getitem__(2));
    df.__setitem__(0, py::float_(exp(-0.5 * pow(tm, 2))));
    fg = py::cast<double>(p.__getitem__(0)) * py::cast<double>(df.__getitem__(0));
    df.__setitem__(1, py::float_(tm / py::cast<double>(p.__getitem__(2)) * fg));
    df.__setitem__(2, py::float_(tm * py::cast<double>(df.__getitem__(1))));
    
    return fg;
}


double FitGauss(PyVector& x, PyVector& y, PyVector& dy, const int N, PyVector& p, PyVector& f) {

	return LevMar(x, y, dy, N, p, f, 3, &gauss_fit, 1.e-6, 100);
}

void init_numerics(py::module_ &torus) {
	torus.def("GaussJordan", &GaussJordan);
	torus.def("GaussJordan", &GaussJordanVec);
	torus.def("GaussBack", &GaussBack);
	torus.def("qbulir", &qbulir, py::arg("func"), py::arg("a"), py::arg("b"), py::arg("eps_"), py::arg("err"));
	torus.def("GaussLegendre", &GaussLegendre, py::arg("x"), py::arg("w"), py::arg("n"));
	torus.def("LegendrePeven", &LegendrePeven, py::arg("p"), py::arg("x"), py::arg("np"));
	torus.def("dLegendrePeven", &dLegendrePeven, py::arg("p"), py::arg("d"), py::arg("x"), py::arg("np"));
	torus.def("zbrent", &zbrent, py::arg("func"), py::arg("x1"), py::arg("x2"), py::arg("tol"));
	torus.def("heap_index", &int_heap_index, py::arg("A"), py::arg("n"), py::arg("indx"));
	torus.def("heap_index", &int_func_heap_index, py::arg("func"), py::arg("n"), py::arg("indx"));
	torus.def("heap_index", [](std::function<double(const int)> func, int n, PyVector& pyindx) {
    	double_func_heap_index(func, n, pyindx);
	}, py::arg("func"), py::arg("n"), py::arg("indx"));
	torus.def("qsplin", &qsplin, py::arg("x"), py::arg("y"), py::arg("y2"), py::arg("n"), py::arg("al"), py::arg("x1"), py::arg("x2"));
	torus.def("CholeskyDecomposition", &CholeskyDecomposition, py::arg("a"), py::arg("n"));
	torus.def("CholeskySolution", &CholeskySolution, py::arg("a"), py::arg("n"), py::arg("b"));
	torus.def("CholeskyInvertL", &CholeskyInvertL, py::arg("a"), py::arg("b"));
	torus.def("CholeskyInvertF", &CholeskyInvertF, py::arg("a"), py::arg("n"));
	torus.def("LUDecomposition", &LUDecomposition, py::arg("a"), py::arg("n"), py::arg("indx"), py::arg("d"));
	torus.def("LUDet3", &LUDet3, py::arg("A"));
	torus.def("LUSolution", &LUSolution, py::arg("a"), py::arg("n"), py::arg("indx"), py::arg("b"));
	torus.def("LUInvert", &LUInvert, py::arg("a"), py::arg("y"), py::arg("n"), py::arg("indx"));
	torus.def("tred2", &tred2, py::arg("a"), py::arg("n"), py::arg("d"), py::arg("e"), py::arg("EV"));
	torus.def("tqli", &tqli, py::arg("d"), py::arg("e"), py::arg("n"), py::arg("z"), py::arg("EV"));
	torus.def("balanc", &balanc, py::arg("a"), py::arg("n"));
	torus.def("elmhes", &elmhes, py::arg("a"), py::arg("n"));
	torus.def("hqr", &hqr, py::arg("a"), py::arg("n"), py::arg("wr"), py::arg("wi"));
	torus.def("LevCof", &LevCof, py::arg("x"), py::arg("y"), py::arg("sig"), py::arg("N"), py::arg("a"), py::arg("fit"), py::arg("M"), py::arg("A"), py::arg("B"), py::arg("func"));
	torus.def("LevMar", &LevMar, py::arg("x"), py::arg("y"), py::arg("sig"), py::arg("N"), py::arg("a"), py::arg("fit"), py::arg("M"), py::arg("func"), py::arg("dcmax"), py::arg("itmax"));
	torus.def("gauss_fit", &gauss_fit, py::arg("x"), py::arg("p"), py::arg("df"), py::arg("D"));
	torus.def("FitGauss", &FitGauss, py::arg("x"), py::arg("y"), py::arg("dy"), py::arg("N"), py::arg("p"), py::arg("f"));
}


