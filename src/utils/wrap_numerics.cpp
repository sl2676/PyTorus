//#include "../../Torus/src/utils/Numerics.cc"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>
namespace py = pybind11;

int Numerics_message(const char* msgs) {
	std::cerr << "ERROR in Numerics: " << msgs << '\n';
	return -1;
}

int GaussJordan(PyMatrix& a, size_t n, PyMatrix& b, size_t m) {
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
        for (int l = 0; l < n; ++l) {
            a.set_value(icol, l, py::float_(py::cast<double>(a.getValueAt(icol, l)) * pivinv));
        }
        for (int l = 0; l < m; ++l) {
            b.set_value(icol, l, py::float_(py::cast<double>(b.getValueAt(icol, l)) * pivinv));
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



void init_numerics(py::module_ &torus) {
	torus.def("GaussJordan", &GaussJordan);
//	torus.def("GaussJordan", &GaussJordanVec);
}
