#ifndef FOURIERTRANSFORM_HPP
#define FOURIERTRANSFORM_HPP


#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <tuple>
#include <chrono>
#include <future>  

#include <mutex>  
#include "eigen-3.4.0/Eigen/Core"
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/Dense"
#include "YoungTableau.hpp"

#include "Irrep.hpp"
#include "permutation_utils.hpp"

using namespace std;

int getFirstDigit(int num) {
    num = std::abs(num);  // Asegurar que sea positivo
    while (num >= 10) {
        num /= 10;  // Eliminar dígitos hasta que quede solo el primero
    }
    return num;
}


template <typename T> class FourierTransform {
    
    public:

    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<T>;
    using SnVector = Eigen::Vector<T, Eigen::Dynamic>;

    int n;
    int nfact;
    string mode;
    vector<vector<int>> partitions;
    map<int, Matrix> coefficients;
    map<int, Irrep<T>> irreps;
    map<int, T> f;
    map<int, T> invFT;
    map<int,SnVector> inverse_ft;
    vector<SnVector> inv_ft_orders;
    int num__threads;

    FourierTransform(int n, map<int,T> f, string mode, int nthreads) {
        
        this->partitions = integer_partitions(n);
        this->n = n;
        this->nfact = 1;
        for(size_t i = 1; i <= n; i++){
            this->nfact *= i;
        }
        this->mode = mode;
        this->f = f;
        this->num__threads = nthreads;

        for(auto p : partitions){
            auto irrep = Irrep<T>(p, this->mode);
            inverse_ft.emplace(to_int(p), SnVector::Zero(nfact));
            this->irreps.emplace(to_int(p), irrep);
            this->coefficients.emplace(to_int(p), Matrix::Zero(irrep.matrices[0].rows(), irrep.matrices[0].cols()));
        }

        for(size_t i = 0; i < n; i++){
            inv_ft_orders.push_back(SnVector::Zero(nfact));
        }

    }

    void build_coefficients(){
        
        #pragma omp parallel for num_threads(1) schedule(dynamic)
        for(int i = 0; i < partitions.size(); i++) {
            auto partition = partitions[i];
            
            double t_0 = std::chrono::system_clock::now().time_since_epoch().count();
            string task = "---- Task " + to_string(to_int(partition));
            
            // Thread-safe output
            #pragma omp critical(cout)
            {
                cout << task << " started" << endl;
            }
            
            auto irrep = this->irreps[to_int(partition)];
            
            #pragma omp critical(cout)
            {
                cout << task << " d_lambda: " << irrep.d_lambda << endl;
            }
            
            auto coef = calculate_coefficient(irrep, this->f, this->nfact, this->num__threads);
            
            #pragma omp critical(maps)
            {
                this->coefficients[to_int(partition)] = coef;
            }
            
            double t_1 = std::chrono::system_clock::now().time_since_epoch().count();
            
            #pragma omp critical(cout)
            {
                cout << task << " finished in " << (t_1 - t_0) / 1e9 << " seconds" << endl;
            }
        }

    }

    void inverse_fourier_transform(){
        
        #pragma omp parallel for num_threads(1) schedule(dynamic)
        for(int i = 0; i < partitions.size(); i++) {
            auto partition = partitions[i];
            
            double t_0 = std::chrono::system_clock::now().time_since_epoch().count();
            string task = "---- Task " + to_string(to_int(partition));
            
            // Thread-safe output
            #pragma omp critical(cout)
            {
                cout << task << " started" << endl;
            }
            
            auto irrep = this->irreps[to_int(partition)];
            
            #pragma omp critical(cout)
            {
                cout << task << " d_lambda: " << irrep.d_lambda << endl;
            }
            
            auto coef = this->coefficients[to_int(partition)];
            if(!(coef == Matrix::Zero(irrep.matrices[0].rows(), irrep.matrices[0].cols()))){
                
                auto vec = invFT_coef(irrep, this->f, this->nfact, coef, this->num__threads);
            
                // Thread-safe update of the coefficients and inverse_ft maps
                #pragma omp critical(maps)
                {
                    this->inverse_ft[to_int(partition)] = vec;
                }

            }
            
            double t_1 = std::chrono::system_clock::now().time_since_epoch().count();
            
            #pragma omp critical(cout)
            {
                cout << task << " finished in " << (t_1 - t_0) / 1e9 << " seconds" << endl;
            }
        }


        for(auto partition : partitions){
            int order = n - partition[0];
            for(int o = order; o < n; o++){
                inv_ft_orders[o] += inverse_ft[to_int(partition)];
            }
        }

        for(size_t i = 0; i < n; i++){
            build_inv_ft(inv_ft_orders[i], i);
        }

    }

    T inverseFT(vector<int> pi, int order){
        return this->invFT[to_int(pi)*10 + order];
    }

    void build_inv_ft(SnVector v, size_t order){
        vector<int> tau;
        tau.push_back(2);
        tau.push_back(1);
        for(size_t i = 3; i <= n; i++){
            tau.push_back(i);
        }
    
        vector<int> sigma;
        for(size_t i = 1; i < n; i++){
            sigma.push_back(i+1);
        }
        sigma.push_back(1);
    
        vector<int> q;
        for(size_t i = 1; i <= n; i++){
            q.push_back(n-i+1);
        }
    
        auto qtau = compose(q, tau);
        auto qsigma = compose(q, sigma);
        auto qsigmatau = compose(qsigma, tau);
        auto invSigma = inverse(sigma);
    
        auto p = compose(qsigma, tau);
        int k = 0;
        this->invFT.emplace(to_int(inverse(p))*10 + order, v[k]);
        k++;

        while(p != qtau){
            if(p != qsigmatau){
                if(williamsCondition(p,n) && p != qsigma){
                    p = compose(p, tau);
                }else{
                    p = compose(p, invSigma);
                }
            }else{
                p = compose(p, invSigma);
            }
            this->invFT.emplace(to_int(inverse(p))*10 + order, v[k]);
            k++;
        }
    }

    void set_coefficient(vector<vector<T>> &coef, vector<int> partition){

        Matrix new_matrix = Matrix(coef.size(), coef[0].size());
        for(size_t i = 0; i < coef.size(); i++){
            for(size_t j = 0; j < coef[i].size(); j++){
                new_matrix(i,j) = coef[i][j];
            }
        }

        this->coefficients[to_int(partition)] = new_matrix;
    }

    // static Matrix calculate_coefficient(Irrep<T> irrep, map<int, T> f, int nfact){

    //     int n = irrep.n;

    //     vector<int> tau;
    //     tau.push_back(2);
    //     tau.push_back(1);
    //     for(size_t i = 3; i <= n; i++){
    //         tau.push_back(i);
    //     }
    
    //     vector<int> sigma;
    //     for(size_t i = 1; i < n; i++){
    //         sigma.push_back(i+1);
    //     }
    //     sigma.push_back(1);
    
    //     vector<int> q;
    //     for(size_t i = 1; i <= n; i++){
    //         q.push_back(n-i+1);
    //     }
    
    //     auto qtau = compose(q, tau);
    //     auto qsigma = compose(q, sigma);
    //     auto qsigmatau = compose(qsigma, tau);
    //     auto invSigma = inverse(sigma);
    
    //     auto p = compose(qsigma, tau);

    //     auto p_matrix = irrep.evaluate(p);
    //     auto tau_matrix = irrep.evaluate(tau);
    //     auto invSigma_matrix = irrep.evaluate(invSigma);

    //     Matrix coefficient = f[to_int(p)]*p_matrix;

    //     while(p != qtau){
    //         if(p != qsigmatau){
    //             if(williamsCondition(p,n) && p != qsigma){
    //                 p = compose(p, tau);
    //                 p_matrix = p_matrix * tau_matrix;
    //             }else{
    //                 p = compose(p, invSigma);
    //                 p_matrix = p_matrix * invSigma_matrix;
    //             }
    //         }else{
    //             p = compose(p, invSigma);
    //             p_matrix = p_matrix * invSigma_matrix;
    //         }
    //         coefficient = coefficient + f[to_int(p)]*p_matrix;
    //     }

    //     cout << coefficient << endl;

    //     return coefficient;
    // }

    template <typename U>
    static inline vector<std::vector<U>> split_vector(const vector<U> vec, int k) {
        vector<vector<U>> result;
        int n = vec.size();
        int chunk = n / k;

        int start = 0;
        for (int i = 0; i < k; i++) {
            vector<U> currVector(chunk);
            for (int j = start; j < chunk*(i + 1); j++) {
                currVector[j - start] = (vec[j]);
            }
            result.push_back(currVector);
            start += chunk;
        }
        return result;
    }
    
    static Matrix calculate_coefficient(Irrep<T> irrep, map<int, T> f, int nfact, int num__threads){

        auto perms = permutations(irrep.n);

        auto parts = split_vector(perms, num__threads);

        auto coef = Matrix(irrep.matrices[0].rows(), irrep.matrices[0].cols()).setZero();

        
        #pragma omp parallel for num_threads(num__threads) 
        for(int i = 0; i < parts.size(); i++){
            
            auto local_coef = (f[to_int(parts[i][0])]*irrep.evaluate(parts[i][0])).toDense();
            for(int j = 1; j < parts[i].size(); j++){
                local_coef += f[to_int(parts[i][j])]*irrep.evaluate(parts[i][j]);
            }

            #pragma omp critical
            {
                coef += local_coef;
            }
        }

        return coef;
    }

    static SnVector invFT_coef(Irrep<T> irrep, map<int,T> f, int nfact, Matrix coefficient, int num__threads){

        auto perms = permutations(irrep.n);
        auto chunk = nfact/num__threads;

        coefficient.transposeInPlace();

        vector<T> invFT(nfact);

        auto parts = split_vector(perms, num__threads);
        T factor = T(coefficient.rows()) / T(nfact);
        
        #pragma omp parallel for num_threads(num__threads) 
        for(int i = 0; i < parts.size(); i++){
            
            vector<T> local_vec;
            for(int j = 0; j < parts[i].size(); j++){
                local_vec.push_back(frobenius_inner_product(coefficient, irrep.evaluate(parts[i][j])));
            }

            #pragma omp critical
            {
                for(int k = i*chunk; k<(i+1)*chunk; k++){
                    invFT[k] = local_vec[k - i*chunk];
                }
            }
        }

        coefficient.transposeInPlace();

        return factor*Eigen::Map<SnVector>(invFT.data(), invFT.size());
    }

    // static SnVector invFT_coef(Irrep<T> irrep, map<int,T> f, int nfact, Matrix coefficient){
        

        
    //     vector<T> invFT;

    //     int n = irrep.n;

    //     vector<int> tau;
    //     tau.push_back(2);
    //     tau.push_back(1);
    //     for(size_t i = 3; i <= n; i++){
    //         tau.push_back(i);
    //     }
    
    //     vector<int> sigma;
    //     for(size_t i = 1; i < n; i++){
    //         sigma.push_back(i+1);
    //     }
    //     sigma.push_back(1);
    
    //     vector<int> q;
    //     for(size_t i = 1; i <= n; i++){
    //         q.push_back(n-i+1);
    //     }
    
    //     auto qtau = compose(q, tau);
    //     auto qsigma = compose(q, sigma);
    //     auto qsigmatau = compose(qsigma, tau);
    //     auto invSigma = inverse(sigma);
        
    //     auto p = compose(qsigma, tau);

    //     auto p_matrix_ = irrep.evaluate(p);
    //     auto tau_matrix_ = irrep.evaluate(tau);
    //     auto invSigma_matrix_ = irrep.evaluate(invSigma);

    //     T factor = T(coefficient.rows()) / T(nfact);

    //     coefficient.transposeInPlace();

    //     T coef = frobenius_inner_product(coefficient, p_matrix_);
    //     invFT.push_back(coef);

    //     while(p != qtau){
    //         if(p != qsigmatau){
    //             if(williamsCondition(p,n) && p != qsigma){
    //                 p = compose(p, tau);
    //                 p_matrix_ = p_matrix_ * tau_matrix_;
    //             }else{
    //                 p = compose(p, invSigma);
    //                 p_matrix_ = p_matrix_ * invSigma_matrix_;
    //             }
    //         }else{
    //             p = compose(p, invSigma);
    //             p_matrix_ = p_matrix_ * invSigma_matrix_;
    //         }
    //         T coef = frobenius_inner_product(coefficient, (p_matrix_));
    //         invFT.push_back(coef);
    //     }

    //     coefficient.transposeInPlace();

    //     auto eigen_vec = Eigen::Map<SnVector>(invFT.data(), invFT.size());
    //     return factor*eigen_vec;
    // }

    static inline T frobenius_inner_product(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, 
                                            const Eigen::SparseMatrix<T, Eigen::ColMajor>& B) {
        T result = 0;
        for (int k = 0; k < B.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<T, Eigen::ColMajor>::InnerIterator it(B, k); it; ++it) {
                result += A.coeff(it.row(), it.col()) * it.value();
            }
        }
        return result;
    }

};

#endif