#include "tools.hpp"

void fill(std::vector<int>& vec, int n){
    vec.resize(n);
    for(int i = 0; i < n; i++){
        std::cin >> vec[i];
    }
}

void fill(std::vector<std::vector<float> >& vec, int n, int k){
    vec.resize(n);
    for(int i = 0; i < n; i++){
        vec[i].resize(k);
        for(int j = 0; j < k; j++){
            std::cin >> vec[i][j];
        }
    }
}

int findMax(std::vector<float>& vec) {
    float max = vec[0];
    int k = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] > max) {
            max = vec[i];
            k = i;
        }
    }

    return k;
}

int findMinPositive(std::vector<float>& vec) {
    float min = vec[findMax(vec)];
    int k = -1;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i] <= min && vec[i] > 0) {
            min = vec[i];
            k = i;
            
        }
    }
    return k;
}

void transpose(std::vector< std::vector<float> >& vec, int size) {
    int t;
    for(int i = 0; i < size; ++i)
    {
        for(int j = i; j < size; ++j)
        {
            t = vec[i][j];
            vec[i][j] = vec[j][i];
            vec[j][i] = t;
        }
    }
}

void printTask() {
    std::cout << R"(Primary task of LP:
    CX -> max(min)
    AX <= B
    X >= 0
                    )"<<std::endl;
}
