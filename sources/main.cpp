#include "Simplex.hpp"
#include "tools.hpp"

int main(int argc, const char * argv[]) {
    int variablesNumber, restrictionNumber;
    func funcTendention;
    std::string funcTendentionString;
    
    printTask();
    
    std::cout<<"Enter tendention of function(min or max): ";
    std::getline(std::cin, funcTendentionString);
    for(auto& c : funcTendentionString)
    {
       c = std::tolower(c);
    }
    if (funcTendentionString[1] == 'a') {
        funcTendention = MAX;
    } else {
        funcTendention = MIN;
    }
    
    std::cout<<"Enter number of variables: "; std::cin >> variablesNumber;
    std::cout<<"Enter number of restrictions: "; std::cin >> restrictionNumber;
    
    std::vector<int>                  B;
    std::vector<int>                  C;
    std::vector< std::vector<float> > A;
    
    std::cout<<"Enter matrix A: ";
    fill(A,restrictionNumber,variablesNumber);
    
    std::cout<<"Enter matrix B: ";
    fill(B,restrictionNumber);
    
    std::cout<<"Enter matrix C: ";
    fill(C,variablesNumber);
    
    Simplex lab1(A, B, C, variablesNumber, restrictionNumber, funcTendention);
    
    lab1.resolve();
    lab1.printOptimalSolution();
    
    return 0;
}
//3 3 4 1 1 1 2 0 0 0.5 1 4 3 2 7 5 3 F = 13 max
//3 3 2 1 1 1 1 0 0 0.5 2 3 6 3 5 3 8 F = 15.75 max
//2 3 1 -2 -2 1 1 1 2 -2 5 1 -1 F = -3 min
//4 2 3 1 -4 -1 -2 -4 -1 1 -3 -3 -4 -18 -30 -5 F = -36 max
// Проверка на отсутствие решений
//2 2 1 2 1 2 -1 1 1 1
// Проверка на отсутсвие оптимального решения
//3 3 1 -2 3 4 -5 6 7 -8 9 31 32 33 -4 2 1
