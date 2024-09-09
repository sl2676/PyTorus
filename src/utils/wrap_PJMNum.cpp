#include <vector>
#include <cmath>
namespace py = pybind11;

double trapzd(double(*func)(double), const double a, const double b, const int n) {
    static double s;

    if (n == 1)
        return (s = 0.5 * (b - a) * func(a) + func(b));
    else {
        double x, tmp, sum, del;
        int it, j;
        for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
        tmp = 1. / double(it);
        del = (b - a) * tmp;
        x = a + 0.5 * del;
        for (sum = 0., j = 0; j < it; j++, x += del) sum += func(x);
        s = 0.5 * (s + (b - a) * sum * tmp);
        return s;
    }
}

void polint(PyVector& xa, PyVector& ya, const int n, const double x, double &y, double &dy) {
    int ns = 0;
    std::vector<double> c(n), d(n);
    double dift, dif = fabs(x - xa.__getitem__(0).cast<double>());
    for (int i = 0; i != n; i++) {
        if ((dift = fabs(x - xa.__getitem__(i).cast<double>())) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya.__getitem__(i).cast<double>();
        d[i] = ya.__getitem__(i).cast<double>();
    }
    y = ya.__getitem__(ns).cast<double>();
    ns--;
    for (int m = 1; m != n; m++) {
        for (int i = 0; i < n - m; i++) {
            double ho = xa.__getitem__(i).cast<double>() - x;
            double hp = xa.__getitem__(i + m).cast<double>() - x;
            double w = c[i + 1] - d[i];
            double den = ho - hp;
            if (den == 0.0) {
                std::cerr << "Error in polint\n";
                return; 
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        dy = (2 * (ns + 1) < (n - m)) ? c[ns + 1] : d[ns--];
        y += dy;
    }
}

double qromb(double(*func)(double), const double a, const double b, const double EPS) {
    const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
    double ss, dss, s[JMAX], h[JMAXP], s_t[K], h_t[K];

    h[0] = 1.;
    for (int j = 1; j <= JMAX; j++) {
        s[j - 1] = trapzd(func, a, b, j);
        if (j >= K) {
            for (int i = 0; i < K; i++) {
                h_t[i] = h[j - K + i];
                s_t[i] = s[j - K + i];
            }
            PyVector h_t_py = PyVector(py::cast(std::vector<double>(h_t, h_t + K)));
            PyVector s_t_py = PyVector(py::cast(std::vector<double>(s_t, s_t + K)));
            polint(h_t_py, s_t_py, K, 0.0, ss, dss);
            if (fabs(dss) <= EPS * fabs(ss)) return ss;
        }
        h[j] = 0.25 * h[j - 1];
    }
    return ss;
}

double probks(const double alam) {
    const double EPS1 = 1.e-6, EPS2 = 1.e-16;
    int j;
    double a2, fac = 2., sum = 0., term, termbf = 0.;

    a2 = -2. * alam * alam;
    for (j = 1; j != 100; j++) {
        term = fac * exp(a2 * j * j);
        sum += term;
        if (fabs(term) <= EPS1 * termbf || fabs(term) <= EPS2 * sum) return sum;
        fac = -fac;
        termbf = fabs(term);
    }
    return 1.; 
}

void kstwo(PyVector& data1, int n1, PyVector& data2, int n2, double &d, double &prob) {
    auto vec1 = data1.extractDataAs<double>();
    auto vec2 = data2.extractDataAs<double>();

    int j1 = 0, j2 = 0;
    double d1, d2, dt, en1, en2, en, fn1 = 0., fn2 = 0.;

    std::sort(vec1.begin(), vec1.end());
    std::sort(vec2.begin(), vec2.end());

    en1 = n1;
    en2 = n2;
    d = 0.;

    while (j1 < n1 && j2 < n2) {
        if ((d1 = vec1[j1]) <= (d2 = vec2[j2])) fn1 = j1++ / en1;
        if (d2 <= d1) fn2 = j2++ / en2;
        if ((dt = fabs(fn2 - fn1)) > d) d = dt;
    }

    en = sqrt(en1 * en2 / (en1 + en2));
    prob = probks((en + 0.12 + 0.11 / en) * d);
}

template <class C>
class stored_qromb {
public:
    C* o;
    double (C::*func)(double) const;
    double a, b, *table;
    int ntab, ntabMAX, tabMAX, nmax;
    double sout;

    stored_qromb(C* obj, double (C::*f)(double) const) : o(obj), func(f) {
        ntab = 0;
        ntabMAX = 20;
        table = new double[ntabMAX];
        nmax = 20;
    }

    void sort_table() {
        std::sort(table, table + ntab);
    }

    double trapzd_store(const int n) {
        static double s;

        if (n == 1) {
            return (s = 0.5 * (b - a) * ((o->*func)(a) + (o->*func)(b)));
        } else {
            double x, tmp, sum, del;
            int it, j;
            for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
            tmp = 1. / double(it);
            del = (b - a) * tmp;
            x = a + 0.5 * del;
            for (sum = 0., j = 0; j < it; j++, x += del) sum += (o->*func)(x);
            s = 0.5 * (s + (b - a) * sum * tmp);
            return s;
        }
    }

    double qromb_store(const double aa, const double bb) {
        a = aa;
        b = bb;
        const double EPS = 1.e-6;
        const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
        double ss, dss, s[JMAX], h[JMAXP], s_t[K], h_t[K];

        h[0] = 1.;
        for (int j = 1; j <= JMAX; j++) {
            s[j - 1] = trapzd_store(j);
            if (j >= K) {
                for (int i = 0; i < K; i++) {
                    h_t[i] = h[j - K + i];
                    s_t[i] = s[j - K + i];
                }
                PyVector h_t_vec = PyVector(py::cast(std::vector<double>(h_t, h_t + K)));
				PyVector s_t_vec = PyVector(py::cast(std::vector<double>(s_t, s_t + K)));
				polint(h_t_vec, s_t_vec, K, 0.0, ss, dss);
                if (fabs(dss) <= EPS * fabs(ss)) return ss;
            }
            h[j] = 0.25 * h[j - 1];
        }
        return ss;
    }

    ~stored_qromb() {
        delete[] table;
    }
};

class PyStoredQromb {
public:
    py::object py_obj;
    py::function py_func;
    double a, b, *table;
    int ntab, ntabMAX, tabMAX, nmax;
    double sout;

    PyStoredQromb(py::object obj, py::function func) : py_obj(obj), py_func(func) {
        ntab = 0;
        ntabMAX = 20;
        table = new double[ntabMAX];
        nmax = 20;
    }

	void sort_table() {
        //cerr << ntab <<  ' ' << ntabMAX << '\n';

        double *integrands, *cumulative;
        integrands = new double[ntab];
        cumulative = new double[ntabMAX];
        integrands[0] = table[0];
        integrands[ntab-1] = table[ntab-1];
        int oldend = ntab;
        //for(int i=0;i!=ntab;i++) cerr << table[i] << '\n';
        for(int i=0;i!=nmax-1;i++) { // ??
            int n = nmax-i, step = pow(2,i+1), point=pow(2,i); 
            //cerr <<"point step "<< point << ' ' << step << '\n';
            int it=1;
            for(int j=1;j<nmax-i-1;j++) it <<= 1;
            //cerr << oldend << ' ' << it << '\n';
            for(int i=0;i!=it;i++) {
                //cerr << table[oldend-it+i] << '\n';
                integrands[point] = table[oldend-it+i];
                point += step;
            }
            oldend -= it;
        }
        //for(int i=0;i!=ntab;i++) cerr << integrands[i] << '\n';
        cumulative[0] = 0.;
        for(int i=1;i!=ntab;i++) {
            cumulative[i] = cumulative[i-1]+(0.5*(integrands[i-1]+integrands[i]))/double(ntab-1);
            //cerr << cumulative[i] << '\n';
        }
        for(int i=ntab;i!=ntabMAX;i++)  cumulative[i] = 0.;
        delete[] integrands;
        delete[] table;
        table = cumulative;
    }    

	double trapzd_store(const int n) {
        static double s;

        if (n == 1) {
            return (s = 0.5 * (b - a) * (py_func(py_obj, a).cast<double>() + py_func(py_obj, b).cast<double>()));
        } else {
            double x, tmp, sum, del;
            int it, j;
            for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
            tmp = 1. / double(it);
            del = (b - a) * tmp;
            x = a + 0.5 * del;
            for (sum = 0., j = 0; j < it; j++, x += del) sum += py_func(py_obj, x).cast<double>();
            s = 0.5 * (s + (b - a) * sum * tmp);
            return s;
        }
    }

    double qromb_store(const double aa, const double bb) {
        a = aa;
        b = bb;
        const double EPS = 1.e-6;
        const int JMAX = 20, JMAXP = JMAX + 1, K = 5;
        double ss, dss, s[JMAX], h[JMAXP], s_t[K], h_t[K];

        h[0] = 1.;
        for (int j = 1; j <= JMAX; j++) {
            s[j - 1] = trapzd_store(j);
            if (j >= K) {
                for (int i = 0; i < K; i++) {
                    h_t[i] = h[j - K + i];
                    s_t[i] = s[j - K + i];
                }
				PyVector h_t_vec = PyVector(py::cast(std::vector<double>(h_t, h_t + K)));
				PyVector s_t_vec = PyVector(py::cast(std::vector<double>(s_t, s_t + K)));
                polint(h_t_vec, s_t_vec, K, 0.0, ss, dss);
                if (fabs(dss) <= EPS * fabs(ss)) return ss;
            }
            h[j] = 0.25 * h[j - 1];
        }
        return ss;
    }

    ~PyStoredQromb() {
        delete[] table;
    }
};


void init_pjmnum(py::module_ &torus) {

	py::class_<PyStoredQromb>(torus, "PyStoredQromb")
        .def(py::init<py::object, py::function>())
        .def("trapzd_store", &PyStoredQromb::trapzd_store)
        .def("qromb_store", &PyStoredQromb::qromb_store);

	torus.def("polint", &polint, "A function that performs polynomial interpolation",
          py::arg("xa"), py::arg("ya"), py::arg("n"), py::arg("x"), py::arg("y"), py::arg("dy"));
	torus.def("trapzd", &trapzd, "A function that performs numerical integration using the trapezoidal rule",
          py::arg("func"), py::arg("a"), py::arg("b"), py::arg("n"));	
	torus.def("qromb", &qromb, "A function that performs integration using the Romberg method",
          py::arg("func"), py::arg("a"), py::arg("b"), py::arg("EPS"));	
	torus.def("probks", &probks, "A function that computes the Kolmogorov-Smirnov probability",
          py::arg("alam"));
	torus.def("kstwo", &kstwo, "A function that computes the Kolmogorov-Smirnov two-sample test",
			py::arg("data1"), py::arg("n1"), py::arg("data2"), py::arg("n2"), py::arg("d"), py::arg("prob"));
}
