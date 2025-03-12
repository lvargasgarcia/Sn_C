#ifndef RATIONALSCALAR_HPP
#define RATIONALSCALAR_HPP
#include <iostream>
#include "eigen-3.4.0/Eigen/Core"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <numeric>

// // Clase Racional
class Rational {
public:
    long numerator, denominator;

    // Constructor
    Rational(long num, long denom) : numerator(num), denominator(denom) {
        // Asegurarnos de que el denominador sea siempre positivo
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
        // Simplificar la fracción usando el máximo común divisor (gcd)
        simplify();
    }

    Rational(double d){
        long precision = 1000000;
        long num = (long)(d*precision);
        long denom = precision;
        numerator = num;
        denominator = denom;
        simplify();
    }

    // Constructor por defecto
    Rational(){
        numerator = 0;
        denominator = 1;
    }

    // Método para simplificar la fracción
    void simplify() {
        long gcd = std::gcd(numerator, denominator);
        numerator /= gcd;
        denominator /= gcd;
    }

    // // Operadores de suma
    // Rational operator+(const Rational& other) const {
    //     long num = numerator * other.denominator + denominator * other.numerator;
    //     long denom = denominator * other.denominator;
    //     Rational result(num, denom);
    //     result.simplify();
    //     return result;
    // }

    // // Operadores de resta
    // Rational operator-(const Rational& other) const {
    //     long num = numerator * other.denominator - denominator * other.numerator;
    //     long denom = denominator * other.denominator;
    //     Rational result(num, denom);
    //     result.simplify();
    //     return result;
    // }

    Rational operator+(const Rational& other) const {
        long gcd_d = std::gcd(denominator, other.denominator);
        long lcm_d = (denominator / gcd_d) * other.denominator; // MCM sin overflow
    
        long num = numerator * (lcm_d / denominator) + other.numerator * (lcm_d / other.denominator);
        
        return Rational(num, lcm_d);
    }
    
    Rational operator-(const Rational& other) const {
        long gcd_d = std::gcd(denominator, other.denominator);
        long lcm_d = (denominator / gcd_d) * other.denominator;
    
        long num = numerator * (lcm_d / denominator) - other.numerator * (lcm_d / other.denominator);
        
        return Rational(num, lcm_d);
    }
    

    Rational operator-() const {
        return Rational(-numerator, denominator);
    }

    // Rational operator/(const Rational& other) const {
    //     long num = numerator * other.denominator;
    //     long denom = denominator * other.numerator;
    //     Rational result(num, denom);
    //     result.simplify();
    //     return result;
    // }

    Rational operator/(const Rational& other) const {
        long gcd_a = std::gcd(numerator, other.numerator);
        long gcd_b = std::gcd(denominator, other.denominator);
    
        long num = (numerator / gcd_a) * (other.denominator / gcd_b);
        long denom = (denominator / gcd_b) * (other.numerator / gcd_a);
    
        return Rational(num, denom);
    }    

    // // Operadores de multiplicación
    // Rational operator*(const Rational& other) const {
    //     long num = numerator * other.numerator;
    //     long denom = denominator * other.denominator;
    //     Rational result(num, denom);
    //     result.simplify();
    //     return result;
    // }

    Rational operator*(const Rational& other) const {
        long gcd_a = std::gcd(numerator, other.denominator);
        long gcd_b = std::gcd(denominator, other.numerator);
    
        long num = (numerator / gcd_a) * (other.numerator / gcd_b);
        long denom = (denominator / gcd_b) * (other.denominator / gcd_a);
    
        return Rational(num, denom);
    }    

    // Operadores de igualdad
    bool operator==(const Rational& other) const {
        return numerator == other.numerator && denominator == other.denominator;
    }

    // Operadores de asignación con suma
    Rational& operator+=(const Rational& other) {
        *this = *this + other;
        return *this;
    }

    // Operadores de asignación con multiplicación
    Rational& operator*=(const Rational& other) {
        *this = *this * other;
        return *this;
    }

    // Operador de asignación
    Rational& operator=(const Rational& other) {
        if (this != &other) {
            numerator = other.numerator;
            denominator = other.denominator;
        }
        return *this;
    }

    // Método para imprimir el número racional
    string toString() const {
        return to_string(numerator) + "/" + to_string(denominator);
    }

    // Método para convertir a flotante (double)
    double toDouble() const {
        return static_cast<double>(numerator) / denominator;
    }

};

// Especialización de Eigen::NumTraits para la clase Rational
namespace Eigen {
    template<>
    struct NumTraits<Rational> : public NumTraits<long> {
        typedef Rational Scalar;
        typedef Rational Real;
        typedef Rational NonInteger;
        typedef Rational Nested;

        // Métodos necesarios para trabajar con Eigen
        static inline Scalar epsilon() { return Rational(1, 1000000); }
        static inline Scalar dummy_precision() { return Rational(0, 1); }
        static inline Scalar highest() { return Rational(std::numeric_limits<long>::max(), 1); }
        static inline Scalar lowest() { return Rational(std::numeric_limits<long>::min(), 1); }

        static inline int digits10() { return 0; }
        
        static inline bool isInteger(const Scalar& a) { return a.denominator == 1; }
        static inline bool isZero(const Scalar& a) { return a.numerator == 0; }
    };
}

Rational sqrt(const Rational& r) {
    return Rational(0,1);
}

Rational abs(const Rational& r) {
    return Rational(std::abs(r.numerator), r.denominator);
}

std::ostream& operator<<(std::ostream& os, const Rational& r) {
    os << r.numerator << "/" << r.denominator;
    return os;
}

// // // Clase Racional
// class Rational {
//     public:
//         __int128_t numerator, denominator;
    
//         // Constructor
//         Rational(__int128_t num, __int128_t denom) : numerator(num), denominator(denom) {
//             // Asegurarnos de que el denominador sea siempre positivo
//             if (denominator < 0) {
//                 numerator = -numerator;
//                 denominator = -denominator;
//             }
//             // Simplificar la fracción usando el máximo común divisor (gcd)
//             simplify();
//         }
    
//         Rational(double d){
//             __int128_t precision = 1000000;
//             __int128_t num = (long)(d*precision);
//             __int128_t denom = precision;
//             numerator = num;
//             denominator = denom;
//             simplify();
//         }
    
//         // Constructor por defecto
//         Rational(){
//             numerator = 0;
//             denominator = 1;
//         }
    
//         // Método para simplificar la fracción
//         void simplify() {
//             __int128_t gcd = std::gcd(numerator, denominator);
//             numerator /= gcd;
//             denominator /= gcd;
//         }
    
//         // Operadores de suma
//         Rational operator+(const Rational& other) const {
//             __int128_t num = numerator * other.denominator + denominator * other.numerator;
//             __int128_t denom = denominator * other.denominator;
//             Rational result(num, denom);
//             result.simplify();
//             return result;
//         }
    
//         // Operadores de resta
//         Rational operator-(const Rational& other) const {
//             __int128_t num = numerator * other.denominator - denominator * other.numerator;
//             __int128_t denom = denominator * other.denominator;
//             Rational result(num, denom);
//             result.simplify();
//             return result;
//         }
    
//         Rational operator-() const {
//             return Rational(-numerator, denominator);
//         }
    
//         Rational operator/(const Rational& other) const {
//             __int128_t num = numerator * other.denominator;
//             __int128_t denom = denominator * other.numerator;
//             Rational result(num, denom);
//             result.simplify();
//             return result;
//         }
    
//         // Operadores de multiplicación
//         Rational operator*(const Rational& other) const {
//             __int128_t num = numerator * other.numerator;
//             __int128_t denom = denominator * other.denominator;
//             Rational result(num, denom);
//             result.simplify();
//             return result;
//         }
    
//         // Operadores de igualdad
//         bool operator==(const Rational& other) const {
//             return numerator == other.numerator && denominator == other.denominator;
//         }
    
//         // Operadores de asignación con suma
//         Rational& operator+=(const Rational& other) {
//             *this = *this + other;
//             return *this;
//         }
    
//         // Operadores de asignación con multiplicación
//         Rational& operator*=(const Rational& other) {
//             *this = *this * other;
//             return *this;
//         }
    
//         // Operador de asignación
//         Rational& operator=(const Rational& other) {
//             if (this != &other) {
//                 numerator = other.numerator;
//                 denominator = other.denominator;
//             }
//             return *this;
//         }
    
//         // Método para imprimir el número racional
//         string toString() const {
//             return to_string((long)numerator) + "/" + to_string((long)denominator);
//         }
    
//         // Método para convertir a flotante (double)
//         double toDouble() const {
//             return static_cast<double>(numerator) / denominator;
//         }
    
//     };
    
//     // Especialización de Eigen::NumTraits para la clase Rational
//     namespace Eigen {
//         template<>
//         struct NumTraits<Rational> : public NumTraits<__int128_t> {
//             typedef Rational Scalar;
//             typedef Rational Real;
//             typedef Rational NonInteger;
//             typedef Rational Nested;
    
//             // Métodos necesarios para trabajar con Eigen
//             static inline Scalar epsilon() { return Rational(1, 1000000); }
//             static inline Scalar dummy_precision() { return Rational(0, 1); }
//             static inline Scalar highest() { return Rational(std::numeric_limits<__int128_t>::max(), 1); }
//             static inline Scalar lowest() { return Rational(std::numeric_limits<__int128_t>::min(), 1); }
    
//             static inline int digits10() { return 0; }
            
//             static inline bool isInteger(const Scalar& a) { return a.denominator == 1; }
//             static inline bool isZero(const Scalar& a) { return a.numerator == 0; }
//         };
//     }
    
//     Rational sqrt(const Rational& r) {
//         return Rational(0,1);
//     }
    
//     Rational abs(const Rational& r) {
//         return Rational(std::abs((long)r.numerator), r.denominator);
//     }
    
//     std::ostream& operator<<(std::ostream& os, const Rational& r) {
//         os << (long)r.numerator << "/" << (long)r.denominator;
//         return os;
//     }


#endif