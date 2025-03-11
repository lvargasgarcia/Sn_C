#ifndef YOUNG_TABLEAU_HPP
#define YOUNG_TABLEAU_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <tuple>

using namespace std;


int cells_below_cell(vector<int> partition, int i, int j) {
    int resp = 1;
    for (size_t k = i+1; k < partition.size(); ++k) {
        if (partition[k] > j) {
            resp += 1;
        }else{
            break;
        }
    }
    return resp;
}

int hook_length_prod(vector<int> partition) {
    int resp = 1;
    for (size_t i = 0; i < partition.size(); ++i) {
        for (size_t j = 0; j < partition[i]; ++j) {
            int cells_in_row = partition[i] - j;
            int cells_below = cells_below_cell(partition, i, j);
            resp *= (cells_in_row + cells_below - 1);
        }
    }
    return resp;
}

int frame_robinson_thrall(vector<int> partition, int nfact) {
    return (int) (nfact / hook_length_prod(partition));
}

class YoungTableau {
    
    public:
        
        vector<int> shape;
        int n;
        int nfact;
        vector<vector<int>> tableau;
        vector<pair<int, int>> y_symbol;
        vector<int> subpartition;
        int d;

        YoungTableau(const vector<int> partition) {
            
            shape = partition;
            n = accumulate(shape.begin(), shape.end(), 0);
            nfact = 1;
            for (int i = 2; i <= n; ++i) {
                nfact *= i;
            }
            tableau.resize(shape.size(), vector<int>(n, 0));

            for (size_t i = 0; i < shape.size(); ++i) {
                tableau[i].resize(shape[i], 0);
            }

        }

        void set_subpartition() {
            vector<int> partition(shape.size(), 0);
            for (size_t i = 0; i < shape.size(); ++i) {
                for (size_t j = 0; j < shape[i]; ++j) {
                    if (tableau[i][j] == 0) {
                        partition[i] += 1;
                    }
                }
            }

            partition.erase(remove(partition.begin(), partition.end(), 0), partition.end());
            this->subpartition = partition;
        }

        void set_y_symbol() {
            
            pair<int, int> pair_n (-1,-1);
            pair<int, int> pair_n1 (-1,-1); 
            vector<pair<int, int>> y_symbol (2);

            y_symbol[0] = pair_n;
            y_symbol[1] = pair_n1;

            for(size_t i = 0; i < tableau.size(); ++i) {
                for(size_t j = 0; j < tableau[i].size(); ++j) {
                    if (tableau[i][j] == n) {
                        y_symbol[0] = make_pair(i, j);
                    }else if (tableau[i][j] == n - 1) {
                        y_symbol[1] = make_pair(i, j);
                    }
                }
            }
            this->y_symbol = y_symbol;
        }
    
        // Function to generate first-level Young tableaux
    static vector<YoungTableau> young_tableaux_1stlevel(const vector<int> partition) {
        vector<YoungTableau> result;
        for (size_t i = 0; i < partition.size(); ++i) {
            int j = partition[i] - 1;
            bool no_cells_below = (cells_below_cell(partition, i, j) == 1);
            if (no_cells_below) {
                YoungTableau yt(partition);
                yt.tableau[i][j] = yt.n;
                yt.set_y_symbol();
                yt.set_subpartition();
                yt.d = frame_robinson_thrall(partition, yt.nfact);
                result.push_back(yt);
            }
        }
        return result;
    }

    // Function to generate second-level Young tableaux

    static vector<YoungTableau> young_tableaux_2ndlevel(const vector<int> partition) {
        vector<YoungTableau> first_level = young_tableaux_1stlevel(partition);
        vector<YoungTableau> second_level;
        for(size_t i = 0; i < first_level.size(); ++i) {
            YoungTableau yt = first_level[i];
            vector<int> subpartition = first_level[i].subpartition;
            vector<YoungTableau> tableaux_to_append = young_tableaux_1stlevel(subpartition);
            for(size_t j = 0; j < tableaux_to_append.size(); ++j) {
                if (j > 0) {
                    yt.tableau[yt.y_symbol[1].first][yt.y_symbol[1].second] = 0;
                }
                YoungTableau yt2 = tableaux_to_append[j];
                int a = yt2.y_symbol[0].first;
                int b = yt2.y_symbol[0].second;
                yt.tableau[a][b] = yt2.n;
                yt.y_symbol[1] = yt2.y_symbol[0];
                yt.subpartition = yt2.subpartition;
                yt.d = frame_robinson_thrall(yt.subpartition, (int)(yt.nfact / (yt.n*(yt.n-1))));
                second_level.push_back(yt);
            }
        }
        return second_level;
    }

    string to_string() {
        string result = "Tableau:\n";
        for (size_t i = 0; i < tableau.size(); ++i) {
            for (size_t j = 0; j < tableau[i].size(); ++j) {
                result += std::to_string(tableau[i][j]) + " ";
            }
            result += "\n";
        }

        result += "Y-symbol: ";
        for (size_t i = 0; i < y_symbol.size(); ++i) {
            result += "(" + std::to_string(y_symbol[i].first) + ", " + std::to_string(y_symbol[i].second) + ") ";
        }
        result += "\n";

        result += "Subpartition: ";
        for (size_t i = 0; i < subpartition.size(); ++i) {
            result += std::to_string(subpartition[i]) + " ";
        }
        result += "\n";

        result += "d_lambda: " + std::to_string(d) + "\n";

        return result;
    }


};

#endif  // YOUNG_TABLEAU_HPP
