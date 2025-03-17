#include <pybind11/stl.h>

//#include "../../Torus/src/utils/CHB.cc"
#include "../../Torus/src/utils/CHB.h"

namespace py = pybind11;

void init_chb(py::module_ &torus) {
	py::class_<Cheby>(torus, "Cheby")
        .def(py::init<>())
        .def(py::init<const Cheby&>())
        .def(py::init<const int>())
        .def(py::init<double*, const int>(), py::arg("inp"), py::arg("N"))
        .def(py::init<double*, double*, const int, const int>(), py::arg("x"), py::arg("y"), py::arg("np"), py::arg("N"))
        .def("writecoeffs", [](const Cheby& self) {
            std::ostringstream oss;
            self.writecoeffs(oss);
            return oss.str();
        })
        .def("setcoeffs", &Cheby::setcoeffs, py::arg("inp"), py::arg("N"))
        .def("chebyfit", &Cheby::chebyfit, py::arg("y"), py::arg("x"), py::arg("np"), py::arg("NC") = 0)
        .def("__add__", [](const Cheby &self, const Cheby &other) -> Cheby {
            return const_cast<Cheby&>(self) + other;
        }, py::is_operator())
        .def("__mul__", [](const Cheby &self, double scalar) -> Cheby {
            return const_cast<Cheby&>(self) * scalar;
        }, py::is_operator())
		.def("unfit1", &Cheby::unfit1)
        .def("unfitn", [](const Cheby &self, const std::vector<double> &x) {
  			std::vector<double> y(x.size());
    		std::vector<double> x_copy = x;
    		self.unfitn(x_copy.data(), y.data(), static_cast<int>(x.size()));
    		return y;
		})
		.def("unfitderiv", [](const Cheby &self, double x) {
            double y, dy, d2y, d3y; 
            self.unfitderiv(x, y, dy, d2y, d3y);
            return py::make_tuple(y, dy, d2y, d3y);
        })
        .def("unfitderiv", [](const Cheby &self, double x) {
            double y, dy, d2y; 
            self.unfitderiv(x, y, dy, d2y);
            return py::make_tuple(y, dy, d2y);
        })
        .def("unfitderiv", [](const Cheby &self, double x) {
            double y, dy; 
            self.unfitderiv(x, y, dy);
            return py::make_tuple(y, dy);
        });
}
