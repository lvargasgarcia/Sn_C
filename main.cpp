#include "FourierTransform.hpp"
#include "Irrep.hpp"
#include <random>
#include "rationalscalar.hpp"
#include <chrono>

using namespace std;

double normal_random(double mean, double stddev) {
    // Create a random number generator
    static std::random_device rd;
    static std::mt19937 gen(rd());
    
    // Create a normal distribution with the specified mean and standard deviation
    std::normal_distribution<double> distribution(mean, stddev);
    
    // Generate and return a random number from the distribution
    return distribution(gen);
}


// mpq_class sqrt(mpq_class x) {
//     return mpq_class(0);
// }

// // Sobrecarga del operador de salida
// std::ostream& operator<<(std::ostream& os, const mpq_class& r) {
//     os << r.get_str();
//     return os;
// }

// mpf_class sqrt(mpf_class x) {
//     return mpf_class(0);
// }

// // Sobrecarga del operador de salida
// std::ostream& operator<<(std::ostream& os, const mpf_class& r) {
//     os << r.get_d();
//     return os;
// }

int main() {

    using MpqScalar = Rational;

    double t_0 = std::chrono::system_clock::now().time_since_epoch().count();
    
    map<int,MpqScalar> f;
    map<int,double> f1;
    int n = 8;
    int mean = 100;
    int stddev = 1;
    vector<vector<int>> perms = permutations(n);
    for(auto permutation : perms){
        f[to_int(permutation)] = MpqScalar((int)(normal_random(400,200)*1000),1);
        //f1[to_int(permutation)] = f[to_int(permutation)].toDouble();
        cout << "f(" << to_int(permutation) << ") = " << f[to_int(permutation)] << endl;
    }

    auto ft1 = FourierTransform<MpqScalar>(n, f, "YKR", 14);

    ft1.build_coefficients();

    ft1.inverse_fourier_transform();
    // auto ft2 = FourierTransform<double>(n, f1, "YKR", 4);


    vector<double> errores_maximos;
    
    for(int order = 0; order < n; order++){  
        double max_error = 0;
        for(auto permutation : perms){
            
            // cout << to_int(permutation) << endl;
            // cout << f[to_int(permutation)].toDouble() << endl;
            // cout << ft1.inverseFT(permutation, order).toDouble() << endl;
            // cout << "--------------------------------------" << endl;
            
            double error = (abs(f[to_int(permutation)] - ft1.inverseFT(permutation, order))).toDouble();
            if(error > max_error){
                max_error = error;
            }
        }
        errores_maximos.push_back(max_error);
    }

    cout << "Maximos errores: " << endl;
    for(int i = 0; i < errores_maximos.size(); i++){
        cout << "Orden " << i << " : " << errores_maximos[i] << endl;
    }

    // auto irrep = Irrep<double>({4,2}, "YKR");

    // for(auto pi : perms){
    //     cout << "Permutation: " << to_int(pi) << endl;
    //     cout << irrep.evaluate(pi) << endl;
    // }

    // for(auto [key, value] : ft1.coefficients){
    //     cout << "Partition: " << key << endl;
    //     cout << value << endl;
    // }

    // for(auto [key, value] : ft2.coefficients){
    //     cout << "Partition: " << key << endl;
    //     cout << value << endl;
    // }

    // auto irrep1 = Irrep<MpqScalar>({2,2}, "YKR");
    // auto irrep2 = Irrep<double>({2,2}, "YKR");

    // for (auto m : irrep1.matrices){
    //     cout << m << endl;
    // }

    // for (auto m : irrep2.matrices){
    //     cout << m << endl;
    // }

    // cout << irrep1.evaluate({1,2,3,4}) << endl;
    // cout << irrep2.evaluate({1,2,3,4}) << endl;

    // cout << irrep1.evaluate({2,3,1,4}) << endl;
    // cout << irrep2.evaluate({2,3,1,4}) << endl;

    double t_1 = std::chrono::system_clock::now().time_since_epoch().count();

    cout << "Execution time: " << (t_1 - t_0) / 1e9 << " seconds" << endl;

    return 0;
}