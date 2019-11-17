#ifndef Simplex_hpp
#define Simplex_hpp

#include "tools.hpp"

enum func {MIN, MAX};

class Simplex {
protected:
    int _iterationNumber = 1;
    int _r;
    int _k;
    int _rows;
    int _columns;
    int _variablesNumber;
    int _restrictionNumber;
    func _funcTendention;
    std::vector< std::vector<float> > _A;
    std::vector<int>                  _B;
    std::vector<int>                  _C;
    std::vector<std::vector<float> > _table;
    std::vector<int>                 _variablesRow;
    std::vector<int>                 _variablesColumn;
public:
    Simplex(std::vector<std::vector<float> >& A,
            std::vector<int>&                 B,
            std::vector<int>&                 C,
            int variablesNumber,
            int restrictionNumber,
            func funcTendention
    );
    
    void makeTable();
    
    void makeIteration();
    void makeIterationBasis();
    
    int findSolvingColumn();
    int findSolvingRow();
    
    void findSolvingElementBasis();
    
    float getSolvingElement();
    
    bool isFRowHasPositiveValues();
    bool isSColumnHasNegativeValues();
    
    void findOptimalSolution();
    void printOptimalSolution();
    
    void findBasisSolution();
    void printBasisSolition();
    
    void resolve();
    
    void printSolvingElement();
    void print();
};

#endif
