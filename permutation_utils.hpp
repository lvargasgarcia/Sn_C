#ifndef PERMUTATION_UTILS_HPP
#define PERMUTATION_UTILS_HPP


#include <vector>
#include <iostream>

using namespace std;


vector<int> compose(vector<int> a, vector<int> b) {
    vector<int> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[b[i]-1];
    }
    return result;
}

vector<int> inverse(vector<int> a) {
    vector<int> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[a[i]-1] = i+1;
    }
    return result;
}

vector<vector<int>> cycles(vector<int> pi){
    vector<vector<int>> result;
    vector<int> visited(pi.size(), 0);
    for (int i = 0; i < pi.size(); i++) {
        if (visited[i] == 0) {
            vector<int> cycle;
            int j = i;
            while (visited[j] == 0) {
                cycle.push_back(j+1);
                visited[j] = 1;
                j = pi[j]-1;
            }
            result.push_back(cycle);
        }
    }
    return result;
}

vector<pair<int,int>> transpositions (vector<int> pi) {
    auto cycles_pi = cycles(pi);
    vector<pair<int,int>> result;
    for (auto cycle : cycles_pi) {
        for (int i = cycle.size() - 1; i >= 1; i--) {
            result.push_back(make_pair(cycle[0], cycle[i]));
        }
    }
    return result;
}

vector<int> adjacent_transpositions_left (pair<int,int> transposition) {
    if((transposition.first + 1) == transposition.second) {
        return {transposition.second};
    }else{
        auto t1 = adjacent_transpositions_left(make_pair(transposition.first+1, transposition.second));
        t1.push_back(transposition.first+1);
        return t1;
    }
}

vector<int> adjacent_transpositions (pair<int,int> transposition) {
    auto v = adjacent_transpositions_left(transposition);
    vector<int> aux1;
    aux1.insert(aux1.end(), v.rbegin(), v.rend());
    vector<int> aux2;
    aux2.insert(aux2.end(), v.begin() + 1, v.end());
    vector<int> result;
    result.insert(result.end(), aux1.begin(), aux1.end());
    result.insert(result.end(), aux2.begin(), aux2.end());
    return result;
}

vector<int> clean(vector<int> a) {
    vector<int> result;
    for(auto x : a){
        result.push_back(x);
    }
    auto iterator = result.begin();
    bool is_clean = true;
    while(iterator != result.end()) {
        if (iterator + 1 != result.end() && *iterator == *(iterator + 1)) {
            result.erase(iterator);
            result.erase(iterator);
            is_clean = false;
        } else {
            iterator++;
        }
    }
    return is_clean ? result : clean(result);
}

vector<int> express_into_adjacent_transpositions(vector<int> pi) {
    vector<int> result;
    auto transps = transpositions(pi);
    for(auto t : transps) {
        auto adj_t = adjacent_transpositions(t);
        result.insert(result.end(), adj_t.begin(), adj_t.end());
    }
    return clean(result);
}


string print_vec(vector<int> a) {
    string result = "[";
    for (int i = 0; i < a.size(); i++) {
        result += to_string(a[i]);
        if (i < a.size()-1) {
            result += ", ";
        }
    }
    result += "]";
    return result;
}

bool williamsCondition(vector<int> p, int n){
    int index = -1;
    int j = 0;
    while(index==-1){
        if(p[j] == n){
            index = j;
        }
        j++;
    }
    int r = p[(index % (n-1)) + 1];
    return r % (n-1) + 1 == p[0];
}

string williams_sequence(int n){
    
    string sequence = "";

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

    while(p != qtau){
        if(p != qsigmatau){
            if(williamsCondition(p,n) && p != qsigma){
                p = compose(p, tau);
                sequence += "t";
            }else{
                p = compose(p, invSigma);
                sequence += "i";
            }
        }else{
            p = compose(p, invSigma);
            sequence += "i";
        }
    }
    return sequence;
}

vector<vector<int>> permutations(int n){
    
    vector<vector<int>> perms;

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
    perms.push_back(p);

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
        perms.push_back(p);
    }
    return perms;
}

int to_int(vector<int> vec) {
    int number = 0;
    for (int i = 0; i < vec.size(); ++i) {
        number = number * 10 + vec[i];  // Shift current number left by one digit and add the next digit
    }
    return number;
}

vector<vector<int>> integer_partitions(int n) { 
    int p[n]; // An array to store a partition
    for(size_t i = 0; i < n; i++){
        p[i] = 0;
    }
    int k = 0; // Index of last element in a partition 
    p[k] = n; // Initialize first partition as number itself 
    vector<vector<int>> result;
    // This loop first prints current partition then generates next 
    // partition. The loop stops when the current partition has all 1s 
    while (true) 
    { 
        // print current partition 
        vector<int> partition;
        for (int i = 0; i <= k; i++){
            if(p[i] != 0){
                partition.push_back(p[i]);
            }
        }
        result.push_back(partition);
        // Generate next partition 
 
        // Find the rightmost non-one value in p[]. Also, update the 
        // rem_val so that we know how much value can be accommodated 
        int rem_val = 0; 
        while (k >= 0 && p[k] == 1) 
        { 
            rem_val += p[k]; 
            k--; 
        } 
 
        // if k < 0, all the values are 1 so there are no more partitions 
        if (k < 0) break; 
 
        // Decrease the p[k] found above and adjust the rem_val 
        p[k]--; 
        rem_val++; 
 
 
        // If rem_val is more, then the sorted order is violated. Divide 
        // rem_val in different values of size p[k] and copy these values at 
        // different positions after p[k] 
        while (rem_val > p[k]) 
        { 
            p[k+1] = p[k]; 
            rem_val = rem_val - p[k]; 
            k++; 
        } 
 
        // Copy rem_val to next position and increment position 
        p[k+1] = rem_val; 
        k++; 
    }
    return result;
} 

#endif