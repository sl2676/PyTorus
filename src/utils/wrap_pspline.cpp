#include <cmath>


namespace py = pybind11;
/*
template<class S>
inline void find_for_Pspline(int& klo, const int n, PyVector& x, const S xi) {
    auto x_vec = x.extractDataAs<S>();
    find(klo, n, x_vec.data(), xi);
    if (klo == n - 1) {
        klo--;
    }
}
void find_for_PsplineWrapper(int& klo, const int n, PyVector& x, const double xi) {
    find_for_Pspline<double>(klo, n, x, xi);
}
*/

template<class S, class T>
void Pspline(
    PyVector& x,   
    PyVector& Y,   
    PyVector& Y1,  
    int n,         
    PyVector& Y3)  
{
    auto x_vec = x.extractDataAs<S>();
    auto Y_vec = Y.extractDataAs<T>();
    auto Y1_vec = Y1.extractDataAs<T>();
    std::vector<T> Y3_vec(n);

    const S zero = 0., one = 1., three = 3., seven = 7., ten = 10., twelve = 12.;
    int i;
    S p, sig, dx, dx1, dx2;
    T dy = Y_vec[1] - Y_vec[0], dy1 = dy;
    std::vector<S> v(n - 1);
    dx = x_vec[1] - x_vec[0];
    Y3_vec[0] = v[0] = zero;
    for (i = 1; i < n - 1; i++) {
        dx1 = x_vec[i + 1] - x_vec[i];
        dx2 = x_vec[i + 1] - x_vec[i - 1];
        dy1 = Y_vec[i + 1] - Y_vec[i];
        sig = dx / dx2;
        p = sig * v[i - 1] - three;
        v[i] = (sig - one) / p;
        Y3_vec[i] = twelve * (seven * Y1_vec[i] * dx2 / (dx * dx1)
                             + three * (Y1_vec[i - 1] / dx + Y1_vec[i + 1] / dx1)
                             - ten * (dy / (dx * dx) + dy1 / (dx1 * dx1))) / dx2;
        Y3_vec[i] = (Y3_vec[i] - sig * Y3_vec[i - 1]) / p;
        dx = dx1;
        dy = dy1;
    }
    Y3_vec[n - 1] = zero;
    for (i = n - 2; i >= 0; i--)
        Y3_vec[i] += v[i] * Y3_vec[i + 1];

    Y3 = PyVector(std::make_unique<TypedVector<T>>(Y3_vec));
}

void PsplineWrapper(PyVector& x, PyVector& Y, PyVector& Y1, const int n, PyVector& Y3) {
    if (x.getType() == "double" && Y.getType() == "double" && Y1.getType() == "double") {
        Pspline<double, double>(x, Y, Y1, n, Y3);
    } else {
        throw std::runtime_error("Unsupported vector types for Pspline");
    }
}



void init_pspline(py::module_ &torus) {
	torus.def("Pspline", &PsplineWrapper, "Compute the third derivative of y using Pspline",
              py::arg("x"), py::arg("Y"), py::arg("Y1"), py::arg("n"), py::arg("Y3"));	
//	torus.def("find_for_Pspline", &find_for_PsplineWrapper, "Find for Pspline with PyVector",
//          py::arg("klo"), py::arg("n"), py::arg("x"), py::arg("xi"));
}

