#ifndef RATIONAL_CASTER_HPP
#define RATIONAL_CASTER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "rationalscalar.hpp"
#include <pybind11/cast.h>
#include <tuple>

namespace pybind11 { namespace detail {

    // Basic type caster for individual Rational values
    template <>
    struct type_caster<Rational> {
    public:
        PYBIND11_TYPE_CASTER(Rational, _("Rational"));

        // Convert from Python to C++
        bool load(handle src, bool) {
            if (!isinstance<tuple>(src) || len(src) != 2) {
                return false;
            }
            auto t = reinterpret_borrow<tuple>(src);
            value = Rational(t[0].cast<long>(), t[1].cast<long>());
            return true;
        }

        // Convert from C++ to Python
        static handle cast(const Rational& src, return_value_policy, handle) {
            return pybind11::make_tuple(src.numerator, src.denominator).release();
        }
    };

    // Explicitly disable direct numpy integration for Rational
    template <>
    struct is_pod_struct<Rational> : std::false_type {};

    // Specialization of npy_format_descriptor for Rational
    template <>
    struct npy_format_descriptor<Rational> {
        static pybind11::dtype dtype() {
            return pybind11::dtype::of<pybind11::object>();
        }
        static constexpr auto name = _("Rational");
        static std::string format() {
            return "O";  // Python object format
        }
        static constexpr bool is_vectorizable = false;
        static constexpr bool is_array = false;
        static constexpr bool has_get_item = false;
        static constexpr auto value = pybind11::detail::npy_api::NPY_OBJECT_;
    };

    // Override the Eigen type caster for Rational matrices
    template <int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    struct type_caster<Eigen::Matrix<Rational, _Rows, _Cols, _Options, _MaxRows, _MaxCols>> {
    private:
        using MatrixType = Eigen::Matrix<Rational, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
        using NumMatrix = Eigen::Matrix<long, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
        using DenMatrix = Eigen::Matrix<long, _Rows, _Cols, _Options, _MaxRows, _MaxCols>;
        
    public:
        PYBIND11_TYPE_CASTER(MatrixType, _("Tuple[numpy.ndarray, numpy.ndarray]"));
        
        bool load(handle src, bool convert) {
            if (!isinstance<tuple>(src) || len(src) != 2) {
                return false;
            }
            
            auto t = reinterpret_borrow<tuple>(src);
            handle num_handle = t[0];
            handle den_handle = t[1];
            
            // Use existing Eigen type casters for the integer matrices
            type_caster<NumMatrix> num_caster;
            type_caster<DenMatrix> den_caster;
            
            if (!num_caster.load(num_handle, convert) || !den_caster.load(den_handle, convert)) {
                return false;
            }
            
            const NumMatrix& num = num_caster.value;
            const DenMatrix& den = den_caster.value;
            
            // Check dimensions
            if (num.rows() != den.rows() || num.cols() != den.cols()) {
                return false;
            }
            
            // Create the Rational matrix
            value.resize(num.rows(), num.cols());
            for (int i = 0; i < num.rows(); ++i) {
                for (int j = 0; j < num.cols(); ++j) {
                    value(i, j) = Rational(num(i, j), den(i, j));
                }
            }
            
            return true;
        }
        
        static handle cast(const MatrixType& src, return_value_policy policy, handle parent) {
            // Create numerator and denominator matrices
            NumMatrix num(src.rows(), src.cols());
            DenMatrix den(src.rows(), src.cols());
            
            for (int i = 0; i < src.rows(); ++i) {
                for (int j = 0; j < src.cols(); ++j) {
                    num(i, j) = src(i, j).numerator;
                    den(i, j) = src(i, j).denominator;
                }
            }
            
            // Cast integer matrices to Python numpy arrays
            auto num_obj = type_caster<NumMatrix>::cast(num, policy, parent);
            auto den_obj = type_caster<DenMatrix>::cast(den, policy, parent);
            
            if (!num_obj || !den_obj) {
                return handle();
            }
            
            return pybind11::make_tuple(num_obj, den_obj).release();
        }
    };

}} // namespace pybind11::detail

#endif // RATIONAL_CASTER_HPP