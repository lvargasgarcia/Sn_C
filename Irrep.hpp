#ifndef IRRREP_HPP
#define IRRREP_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <tuple>
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/Dense"
#include "YoungTableau.hpp"
#include "permutation_utils.hpp"

using namespace std;



template <typename T> class Irrep {

    public:

    using BlockMatrix = vector<vector<pair<pair<int, int>, T>>>;
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using SparseMatrix = Eigen::SparseMatrix<T>;

    vector<int> partition;
    int n;
    string mode;
    int d_lambda;
    vector<SparseMatrix> matrices;

    Irrep(vector<int> partition, string mode) {
        this->partition = partition;
        this->n = accumulate(partition.begin(), partition.end(), 0);
        this->mode = mode;
        this->d_lambda = frame_robinson_thrall(partition, n);
        for(int i = 2; i <= n; i++) {
            calc_matrix(i);
        }
        this->d_lambda = (this->matrices)[0].rows();
    }

    Irrep(){
        this->partition = {-1};
        this->n = -1;
        this->mode = "YKR";
        this->d_lambda = 1;
    }

    void calc_matrix(int i) {
        (this->matrices).push_back(decompress_matrix(generate_matrix(i, this->partition, this->mode)).sparseView());
    }

    SparseMatrix evaluate(vector<int> pi){
        auto transpositions = express_into_adjacent_transpositions(pi);

        if(transpositions.size() == 0){
            auto result = Matrix::Identity(matrices[0].rows(), matrices[0].cols());
            return result.sparseView();
        }
        auto result = matrices[transpositions[0] - 2];
        for(int i = 1; i < transpositions.size(); i++){
            result = result * matrices[transpositions[i] - 2];
        }
        //cout << result << endl;
        return result;
    }

    static BlockMatrix generate_matrix(int i, vector<int> partition, string mode) {
        int length = accumulate(partition.begin(), partition.end(), 0);
        if(i == length){
            return clausen_matrix(partition, mode);
        }else{
            vector<YoungTableau> tableaux = YoungTableau::young_tableaux_1stlevel(partition);
            vector<BlockMatrix> block_matrices;
            for(auto yt: tableaux) {
                block_matrices.push_back(generate_matrix(i, yt.subpartition, mode));
            }
            return direct_sum(block_matrices);
        }
    }

    static BlockMatrix direct_sum(vector<BlockMatrix> matrices) {
        // Calculate total rows and columns
        size_t total_rows = 0;
        size_t total_cols = 0;

        for (const auto& matrix : matrices) {
            total_rows += matrix.size();
            total_cols += matrix[0].size();
        }

        // Initialize the result BlockMatrix with default values
        BlockMatrix result(total_rows, vector<pair<pair<int, int>, T>>(total_cols, {{0, 0}, T(0)}));

        size_t row_start = 0;
        size_t col_start = 0;

        for (const auto& matrix : matrices) {
            size_t rows = matrix.size();
            size_t cols = matrix[0].size();

            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                    result[row_start + i][col_start + j] = matrix[i][j];
                }
            }

            row_start += rows;
            col_start += cols;
        }

        for (size_t i = 0; i < total_rows; ++i) {
            for (size_t j = 0; j < total_cols; ++j) {
                if (result[i][j] == make_pair(make_pair(0, 0), T(0))) {
                    result[i][j] = {{result[i][i].first.first, result[j][j].first.second}, T(0)};
                }
            }
        }

        return result;
    }


    static Matrix decompress_matrix(BlockMatrix matrix) {

        size_t dim = accumulate(matrix.begin(), matrix.end(), 0, [](size_t a, vector<pair<pair<int, int>, T>> b) {return a + b[0].first.first;});
        vector<Matrix> rows(matrix.size());

        for (size_t i = 0; i < matrix.size(); i++) {
            // Inicializa la primera submatriz escalada
            Matrix row = matrix[i][0].second * Matrix::Identity(matrix[i][0].first.first, matrix[i][0].first.second);
            for (int j = 1; j < matrix.size(); j++) {
                Matrix aux = matrix[i][j].second * Matrix::Identity(matrix[i][j].first.first, matrix[i][j].first.second);
                Matrix new_row(row.rows(), row.cols() + aux.cols());
                new_row << row, aux;
                row = new_row;
            }

            rows[i] = row;
        }

        Matrix result = rows[0];
        for (size_t i = 1; i < rows.size(); i++) {
            Matrix aux = rows[i];
            Matrix new_result(result.rows() + aux.rows(), result.cols());
            new_result << result, aux;
            result = new_result;
        }

        return result;
    }

    static BlockMatrix clausen_matrix(vector<int> partition, string mode) {
        vector<YoungTableau> tableaux = YoungTableau::young_tableaux_2ndlevel(partition);
        int n = tableaux.size();
        BlockMatrix matrix = BlockMatrix(n, vector<pair<pair<int, int>, T>>(n, {{-1, -1}, T(0)}));
        //diagonal
        for(int i = 0; i < n; i++) {
            YoungTableau yt = tableaux[i];
            int a = yt.y_symbol[0].first;
            int b = yt.y_symbol[1].second;
            int c = yt.y_symbol[1].first;
            int d = yt.y_symbol[1].second;
            T coef = 0;
            if (a == c) {
                coef = T(1);
            } else if(b == d){
                coef = T(-1);
            }
            matrix[i][i] = {{yt.d, yt.d}, coef};
        }
        //off-diagonal
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){

                if (i < j){

                    YoungTableau yt1 = tableaux[i];
                    YoungTableau yt2 = tableaux[j];

                    int d1 = yt1.d;
                    int d2 = yt2.d;

                    if(yt1.subpartition == yt2.subpartition){

                        int a = yt1.y_symbol[0].first;
                        int b = yt1.y_symbol[0].second;
                        int c = yt1.y_symbol[1].first;
                        int d = yt1.y_symbol[1].second;

                        T chi = T(abs(a - c) + abs(b - d));
                        T chi_sqr = chi * chi;
                        T inv_chi = T(1) / chi;
                        T q_sqr_chi = T(1) - T(1) / (chi_sqr);

                        if (mode == "YKR"){
                            matrix[i][i] = {{d1, d1}, inv_chi};
                            matrix[j][j] = {{d2, d2}, -inv_chi};
                            matrix[i][j] = {{d1, d2}, T(1)};
                            matrix[j][i] = {{d2, d1}, q_sqr_chi};
                        }else if (mode == "YSR"){
                            matrix[i][i] = {{d1, d1}, inv_chi};
                            matrix[j][j] = {{d2, d2}, -inv_chi};
                            matrix[i][j] = {{d1, d2}, q_sqr_chi};
                            matrix[j][i] = {{d2, d1}, T(1)};
                        }else{
                            matrix[i][i] = {{d1, d1}, inv_chi};
                            matrix[j][j] = {{d2, d2}, -inv_chi};
                            matrix[i][j] = {{d1, d2}, sqrt(q_sqr_chi)};
                            matrix[j][i] = {{d2, d1}, sqrt(q_sqr_chi)};
                        }

                    }else{
                        matrix[i][j] = {{d1, d2}, T(0)};
                        matrix[j][i] = {{d2, d1}, T(0)};
                    }


                }
            }
        }
        return matrix;
    }

};

#endif