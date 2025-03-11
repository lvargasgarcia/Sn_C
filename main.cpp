#include "FourierTransform.hpp"
#include "Irrep.hpp"
#include <random>
#include "rationalscalar.hpp"

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

int main() {

    using MpqScalar = Rational;
    
    map<int,MpqScalar> f;
    map<int,double> f1;
    int n = 8;
    int mean = 1000;
    int stddev = 1;
    vector<vector<int>> perms = permutations(n);
    for(auto permutation : perms){
        f[to_int(permutation)] = MpqScalar((int)normal_random(10,10)*100000,1);
        f1[to_int(permutation)] = f[to_int(permutation)].toDouble();
        cout << "f(" << to_int(permutation) << ") = " << f[to_int(permutation)].toString() << endl;
    }

    auto ft1 = FourierTransform<MpqScalar>(n, f, "YKR");
    auto ft2 = FourierTransform<double>(n, f1, "YKR");


    vector<double> errores_maximos;
    
    for(int order = 0; order < n; order++){  
        double max_error = 0;
        for(auto permutation : perms){
            
            cout << to_int(permutation) << endl;
            cout << f[to_int(permutation)].toDouble() << endl;
            cout << ft1.inverseFT(permutation, order).toDouble() << endl;
            cout << "--------------------------------------" << endl;
            
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

    for(auto [key, value] : ft1.coefficients){
        cout << "Partition: " << key << endl;
        cout << value << endl;
    }

    for(auto [key, value] : ft2.coefficients){
        cout << "Partition: " << key << endl;
        cout << value << endl;
    }

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

    return 0;
}