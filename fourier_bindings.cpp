// #include <pybind11/pybind11.h>
// #include <pybind11/operators.h>
// #include <pybind11/stl.h>
// #include "FourierTransform.hpp"
// #include "rationalscalar.hpp"

// // namespace py = pybind11;

// // PYBIND11_MODULE(fourier_transform, m) {
// //     py::class_<FourierTransform<double>>(m, "FourierTransform")
// //         .def(py::init<int, std::map<int, double>, std::string>())
// //         .def("inverseFT", &FourierTransform<double>::inverseFT)
// //         .def_readwrite("coefficients", &FourierTransform<double>::coefficients)
// //         .def_readwrite("irreps", &FourierTransform<double>::irreps)
// //         .def_readwrite("f", &FourierTransform<double>::f)
// //         .def_readwrite("invFT", &FourierTransform<double>::invFT)
// //         .def_readwrite("williams_seq", &FourierTransform<double>::williams_seq);
// // };

// namespace py = pybind11;

// PYBIND11_MODULE(fourier_transform, m) {
//     py::class_<Rational>(m, "Rational")
//         .def(py::init<long, long>(), py::arg("numerator"), py::arg("denominator"))
//         .def(py::init<double>())  // Constructor desde double
//         .def("__repr__", [](const Rational &r) {
//             return "<Rational(" + std::to_string(r.numerator) + "/" + std::to_string(r.denominator) + ")>";
//         })
//         .def("__str__", &Rational::toString)
//         .def("to_double", &Rational::toDouble)
//         .def(py::self + py::self)
//         .def(py::self - py::self)
//         .def(py::self * py::self)
//         .def(py::self / py::self)
//         .def(py::self += py::self)
//         .def(py::self *= py::self)
//         .def(py::self == py::self);

//     // 📌 Exponer FourierTransform para tipo double
//     py::class_<FourierTransform<double>>(m, "FourierTransform_double")
//         .def(py::init<int, std::map<int, double>, std::string, int>())
//         .def("inverseFT", &FourierTransform<double>::inverseFT)
//         .def_readwrite("coefficients", &FourierTransform<double>::coefficients)
//         .def_readwrite("irreps", &FourierTransform<double>::irreps)
//         .def_readwrite("f", &FourierTransform<double>::f)
//         .def_readwrite("invFT", &FourierTransform<double>::invFT)
//         .def_readwrite("williams_seq", &FourierTransform<double>::williams_seq);

//     // 📌 Exponer FourierTransform para tipo Rational
//     py::class_<FourierTransform<Rational>>(m, "FourierTransform_Rational")
//         .def(py::init<int, std::map<int, Rational>, std::string, int>())
//         .def("inverseFT", &FourierTransform<Rational>::inverseFT)
//         .def_readwrite("coefficients", &FourierTransform<Rational>::coefficients)
//         .def_readwrite("irreps", &FourierTransform<Rational>::irreps)
//         .def_readwrite("f", &FourierTransform<Rational>::f)
//         .def_readwrite("invFT", &FourierTransform<Rational>::invFT)
//         .def_readwrite("williams_seq", &FourierTransform<Rational>::williams_seq);
// };


#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include "FourierTransform.hpp"
#include "rationalscalar.hpp"
#include "Irrep.hpp"

namespace py = pybind11;

PYBIND11_MODULE(fourier_transform, m) {
    // Binding para la clase Rational
    py::class_<Rational>(m, "Rational")
        .def(py::init<long, long>(), py::arg("numerator"), py::arg("denominator"))
        .def(py::init<double>())  // Constructor desde double
        .def("__repr__", [](const Rational &r) {
            return "<Rational(" + std::to_string(r.numerator) + "/" + std::to_string(r.denominator) + ")>";
        })
        .def("__str__", &Rational::toString)
        .def("to_double", &Rational::toDouble)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self / py::self)
        .def(py::self += py::self)
        .def(py::self *= py::self)
        .def(py::self == py::self)
        .def_readwrite("numerator", &Rational::numerator)
        .def_readwrite("denominator", &Rational::denominator);

    // Binding para la clase FourierTransform con tipo double
    py::class_<FourierTransform<double>>(m, "FourierTransform_double")
        .def(py::init<int, std::map<int, double>, std::string, int>())
        .def("inverseFT", &FourierTransform<double>::inverseFT)
        .def_readwrite("coefficients", &FourierTransform<double>::coefficients)
        .def_readwrite("irreps", &FourierTransform<double>::irreps)
        .def_readwrite("f", &FourierTransform<double>::f)
        .def_readwrite("invFT", &FourierTransform<double>::invFT)
        .def_readwrite("williams_seq", &FourierTransform<double>::williams_seq);

    // Binding para la clase FourierTransform con tipo Rational
    py::class_<FourierTransform<Rational>>(m, "FourierTransform_Rational")
        .def(py::init<int, std::map<int, Rational>, std::string, int>())
        .def("inverseFT", &FourierTransform<Rational>::inverseFT)
        .def_readwrite("coefficients", &FourierTransform<Rational>::coefficients)
        .def_readwrite("irreps", &FourierTransform<Rational>::irreps)
        .def_readwrite("f", &FourierTransform<Rational>::f)
        .def_readwrite("invFT", &FourierTransform<Rational>::invFT)
        .def_readwrite("williams_seq", &FourierTransform<Rational>::williams_seq);

    // Binding para la clase Irrep
    py::class_<Irrep<Rational>>(m, "Irrep_Rational")
        .def(py::init<std::vector<int>, std::string>())
        .def("evaluate", &Irrep<Rational>::evaluate)
        .def_readwrite("partition", &Irrep<Rational>::partition)
        .def_readwrite("n", &Irrep<Rational>::n)
        .def_readwrite("mode", &Irrep<Rational>::mode)
        .def_readwrite("d_lambda", &Irrep<Rational>::d_lambda)
        .def_readwrite("matrices", &Irrep<Rational>::matrices)
        .def("__repr__", [](const Irrep<Rational> &irrep) {
            return "<Irrep_Rational(partition=" + std::to_string(irrep.partition[0]) + ", mode=" + irrep.mode + ")>";
        })
        .def("__str__", [](const Irrep<Rational> &irrep) {
            return "Irrep_Rational with partition: " + std::to_string(irrep.partition[0]) + " and mode: " + irrep.mode;
        });

    // Documentación del módulo
    m.doc() = R"pbdoc(
        Fourier Transform Module
        -----------------------
        .. currentmodule:: fourier_transform
        .. autosummary::
           :toctree: _generate

           Rational
           FourierTransform_double
           FourierTransform_Rational
           Irrep_Rational
    )pbdoc";
};