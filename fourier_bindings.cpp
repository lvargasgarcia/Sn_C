#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "FourierTransform.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fourier_transform, m) {
    py::class_<FourierTransform<double>>(m, "FourierTransform")
        .def(py::init<int, std::map<int, double>, std::string>())
        .def("inverseFT", &FourierTransform<double>::inverseFT)
        .def_readwrite("coefficients", &FourierTransform<double>::coefficients)
        .def_readwrite("irreps", &FourierTransform<double>::irreps)
        .def_readwrite("f", &FourierTransform<double>::f)
        .def_readwrite("invFT", &FourierTransform<double>::invFT)
        .def_readwrite("williams_seq", &FourierTransform<double>::williams_seq);
};
