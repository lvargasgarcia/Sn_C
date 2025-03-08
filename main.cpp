#include "YoungTableau.hpp"


using namespace std;

int main() {
    vector<int> partition = {4,2,1,1};
    vector<YoungTableau> second_level_tableaux = YoungTableau::young_tableaux_2ndlevel(partition);
    for (size_t i = 0; i < second_level_tableaux.size(); ++i) {
        cout << second_level_tableaux[i].to_string() << endl;
    }
    return 0;
}