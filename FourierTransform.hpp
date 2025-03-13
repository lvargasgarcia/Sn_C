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
    string williams_seq;

    FourierTransform(int n, map<int,T> f, string mode, int nthreads) {
        
        this->partitions = integer_partitions(n);
        this->n = n;
        this->nfact = 1;
        for(size_t i = 1; i <= n; i++){
            this->nfact *= i;
        }
        this->mode = mode;
        this->f = f;
        this->williams_seq = williams_sequence(n);
        
        map<int,SnVector> inverse_ft;
        vector<SnVector> inv_ft_orders;

        for(size_t i = 0; i < n; i++){
            inv_ft_orders.push_back(SnVector::Zero(nfact));
        }
        
        // Forma secuencial

        // for(auto partition : partitions){
            
        //     cout << "Calculating coefficient for partition: " << to_int(partition) << endl;
        //     auto irrep = Irrep<T>(partition,mode);
        //     this->irreps.emplace(to_int(partition), irrep);
        //     cout << "D_lambda" << irrep.d_lambda << endl;
        //     auto info = calculate_coefficient(this->irreps[to_int(partition)], f, this->nfact);
        //     this->coefficients.emplace(to_int(partition), info.first);
        //     inverse_ft.emplace(to_int(partition), info.second);

        // }

        // Using mutex for thread-safe updates to shared data structures
        // std::mutex mtx;
        // std::mutex cout_mutex;

        // // Parallel execution using std::async
        // std::vector<std::future<void>> futures;
        
        // // Launch a task for each partition
        // for(auto partition : partitions){
        //     futures.push_back(std::async(std::launch::async, [this, &mtx, &cout_mutex, &inverse_ft, partition]() {
                
        //         double t_0 = std::chrono::system_clock::now().time_since_epoch().count();
        //         string task = "---- Task " + to_string(to_int(partition));
                
        //         {
        //             std::lock_guard<std::mutex> lock(cout_mutex);
        //             cout << task << " started" << endl;
                    
        //         }

                
        //         auto irrep = Irrep<T>(partition, this->mode);

        //         {
        //             std::lock_guard<std::mutex> lock(cout_mutex);
        //             cout << task << " d_lambda: " << irrep.d_lambda << endl;
        //         }
                
        //         // Thread-safe update of the irreps map
        //         {
        //             std::lock_guard<std::mutex> lock(mtx);
        //             this->irreps.emplace(to_int(partition), irrep);
        //         }
                
        //         auto info = calculate_coefficient(irrep, this->f, this->nfact);
                
        //         // Thread-safe update of the coefficients and inverse_ft maps
        //         {
        //             std::lock_guard<std::mutex> lock(mtx);
        //             this->coefficients.emplace(to_int(partition), info.first);
        //             inverse_ft.emplace(to_int(partition), info.second);
        //         }

        //         double t_1 = std::chrono::system_clock::now().time_since_epoch().count();

        //         {
        //             std::lock_guard<std::mutex> lock(cout_mutex);
        //             cout << task << " finished in " << (t_1 - t_0) / 1e9 << " seconds" << endl;
        //         }

        //     }));
        // }

        // Parallel loop using OpenMP with inline thread count specification
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
        for(int i = 0; i < partitions.size(); i++) {
            auto partition = partitions[i];
            
            double t_0 = std::chrono::system_clock::now().time_since_epoch().count();
            string task = "---- Task " + to_string(to_int(partition));
            
            // Thread-safe output
            #pragma omp critical(cout)
            {
                cout << task << " started" << endl;
            }
            
            auto irrep = Irrep<T>(partition, this->mode);
            
            #pragma omp critical(cout)
            {
                cout << task << " d_lambda: " << irrep.d_lambda << endl;
            }
            
            // Thread-safe update of the irreps map
            #pragma omp critical(maps)
            {
                this->irreps.emplace(to_int(partition), irrep);
            }
            
            auto info = calculate_coefficient(irrep, this->f, this->nfact);
            
            // Thread-safe update of the coefficients and inverse_ft maps
            #pragma omp critical(maps)
            {
                this->coefficients.emplace(to_int(partition), info.first);
                inverse_ft.emplace(to_int(partition), info.second);
            }
            
            double t_1 = std::chrono::system_clock::now().time_since_epoch().count();
            
            #pragma omp critical(cout)
            {
                cout << task << " finished in " << (t_1 - t_0) / 1e9 << " seconds" << endl;
            }
        }
        
        // // Wait for all tasks to complete
        // for(auto& future : futures) {
        //     future.wait();
        // }

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

    static pair<Matrix,SnVector> calculate_coefficient(Irrep<T> irrep, map<int, T> f, int nfact){

        vector<T> invFT;

        int n = irrep.n;

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

        auto p_matrix = irrep.evaluate(p);
        auto tau_matrix = irrep.evaluate(tau);
        auto invSigma_matrix = irrep.evaluate(invSigma);

        Matrix coefficient = f[to_int(p)]*p_matrix;

        while(p != qtau){
            if(p != qsigmatau){
                if(williamsCondition(p,n) && p != qsigma){
                    p = compose(p, tau);
                    p_matrix = p_matrix * tau_matrix;
                }else{
                    p = compose(p, invSigma);
                    p_matrix = p_matrix * invSigma_matrix;
                }
            }else{
                p = compose(p, invSigma);
                p_matrix = p_matrix * invSigma_matrix;
            }
            coefficient = coefficient + f[to_int(p)]*p_matrix;
        }

        p = compose(qsigma, tau);

        auto p_matrix_ = (irrep.evaluate(p));
        auto tau_matrix_ = (tau_matrix);
        auto invSigma_matrix_ = (invSigma_matrix);

        T factor = T(coefficient.rows()) / T(nfact);
        int k = 0;

        coefficient.transposeInPlace();

        T coef = frobenius_inner_product(coefficient, p_matrix_);
        invFT.push_back(coef);
        k++;

        while(p != qtau){
            if(p != qsigmatau){
                if(williamsCondition(p,n) && p != qsigma){
                    p = compose(p, tau);
                    p_matrix_ = p_matrix_ * tau_matrix_;
                }else{
                    p = compose(p, invSigma);
                    p_matrix_ = p_matrix_ * invSigma_matrix_;
                }
            }else{
                p = compose(p, invSigma);
                p_matrix_ = p_matrix_ * invSigma_matrix_;
            }
            T coef = frobenius_inner_product(coefficient, (p_matrix_));
            invFT.push_back(coef);
            k++;
        }

        coefficient.transposeInPlace();
        auto eigen_vec = Eigen::Map<SnVector>(invFT.data(), invFT.size());
        return make_pair(coefficient, factor*eigen_vec);
    }

    static inline T frobenius_inner_product(Matrix& A, Matrix& B){
        Eigen::Map<SnVector> vec_A(A.data(), A.cols()*A.rows());
        Eigen::Map<SnVector> vec_B(B.data(), B.cols()*B.rows());
        return vec_A.dot(vec_B);
    }

    // static inline T frobenius_inner_product(Matrix& A, SparseMatrix& B) {
    //     // Convert dense matrix A to a vector
    //     Eigen::Map<SnVector> vec_A(A.data(), A.cols()*A.rows());
        
    //     // For sparse matrix B, we need to use a different approach
    //     // Create a temporary sparse vector by copying the data
    //     SnVector vec_B(B.cols()*B.rows());
        
    //     // Fill vec_B using B's non-zero entries
    //     for (int k = 0; k < B.outerSize(); ++k) {
    //         for (typename SparseMatrix::InnerIterator it(B, k); it; ++it) {
    //             int linear_index = it.row() + it.col() * B.rows(); // Column-major layout
    //             vec_B.coeffRef(linear_index) = it.value();
    //         }
    //     }
        
    //     // Now compute the dot product
    //     return vec_A.dot(vec_B);
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