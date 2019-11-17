#include <iomanip>

#include "Simplex.hpp"

Simplex::Simplex(std::vector<std::vector<float> >& A,
                 std::vector<int>&                 B,
                 std::vector<int>&                 C,
                 int                 variablesNumber,
                 int               restrictionNumber,
                 func                 funcTendention,
                 bool                           flag) {
    
    _B = B;
    _C = C;
    _r = 0;
    _k = 0;
    dual = flag;
    _funcTendention = funcTendention;
    _variablesNumber = variablesNumber;
    _restrictionNumber = restrictionNumber;
    if (dual) {
        _A = transpose(A);
        _columns = _restrictionNumber + 1;
        _rows    = _variablesNumber + 1;
        _table.resize(_rows);
        for (int i = 0; i < _rows; i++)
            _table[i].resize(_columns);
        // Ряд и строка с номерами переменных
        _variablesRow.resize(_restrictionNumber);
        _variablesColumn.resize(_variablesNumber);
        
        for(int i = 0; i < variablesNumber; i++)
            _variablesRow[i] = i + 1;
        for(int j = 0; j < restrictionNumber; j++)
            _variablesColumn[j] = j + 1 + variablesNumber;
    } else {
        _A = A;
        _columns = _variablesNumber + 1;
        _rows    = _restrictionNumber + 1;
        _table.resize(_rows);
        for (int i = 0; i < _rows; i++)
            _table[i].resize(_columns);
        // Ряд и строка с номерами переменных
        _variablesRow.resize(variablesNumber);
        _variablesColumn.resize(restrictionNumber);
        
        for(int i = 0; i < variablesNumber; i++)
            _variablesRow[i] = i + 1;
        for(int j = 0; j < restrictionNumber; j++)
            _variablesColumn[j] = j + 1 + variablesNumber;
    }
    if (dual) {
        _variable = 'y';
        _functionChar = "F'";
    } else {
        _variable = 'x';
        _functionChar = "F ";
    }
}
// Составление таблицы
void Simplex::makeTable(){
    if (dual) {
        for(int i = 0; i < _rows - 1; i++)
            _table[i][0] = static_cast<float>(-_C[i]);
        for(int i = 0; i < _rows - 1; i++){
            for(int j = 1; j < _columns; j++){
                _table[i][j] = -_A[i][j-1];
            }
        }
        
        _table[_rows-1][0] = 0;
        for(int i = 1; i < _columns; i++)
            _table[_rows-1][i] = static_cast<float>(-_B[i-1]);
        if (_funcTendention == MIN)
            _funcTendention = MAX;
        else
            _funcTendention = MIN;
    } else {
        for(int i = 0; i < _rows - 1; i++)
            _table[i][0] = static_cast<float>(_B[i]);
        
        for(int i = 0; i < _rows - 1; i++){
            for(int j = 1; j < _columns; j++){
                _table[i][j] = _A[i][j-1];
            }
        }
        
        _table[_rows-1][0] = 0;
        for(int i = 1; i < _columns; i++)
            _table[_rows-1][i] = static_cast<float>(_C[i-1]);
    }
    
}
// Проверки для базисного и оптимального решений
bool Simplex::isFRowHasPositiveValues(){
    for(int i = 1; i < _columns; i++)
        if (_table[_rows-1][i] > 0) return true;
    return false;
}

bool Simplex::isSColumnHasNegativeValues() {
    for(int i = 0; i < _rows - 1; i++)
        if (_table[i][0] < 0) return true;
    return false;
}
// Получение разрешающего элемента
float Simplex::getSolvingElement() {
    return _table[_r][_k];
}
// Поиск разрешающего элемента для опорного решения
void Simplex::findSolvingElementBasis() {
    int k = -1;
    int r = 0;
    for (int i = 0; i < _rows - 1; i++) {
        if (_table[i][0] < 0) {
            r = i;
            for (int j = 1; j < _columns; j++) {
                if(_table[r][j] < 0) {
                    k = j;
                    break;
                }
            }
            break;
        }
    }
    try {
        if (k == -1)
            throw r;
    } catch (int k) {
        std::cerr<<"No negative values in row "<<r + 1<<", No solutions for LPP"<<std::endl;
        exit(0);
    }
    std::vector<float> buf(_rows);
    for(int i = 0; i < _rows - 1; i++){
        if(_table[i][k] != 0)
            buf[i] = _table[i][0] / _table[i][k];
        else buf[i] = -1;
        if(buf[r] > buf[i] && buf[i] > 0)
            r = i;
    }
    
     _k = k;
    _r = r;
}
// Итерация поиска опорного решения
void Simplex::makeIterationBasis() {
    std::vector<std::vector<float> > buffer_Table = _table;
    findSolvingElementBasis();
    int r = _r;
    int k = _k;
    printBasisSolition();
    printSolvingElement();
    float solvingElement = getSolvingElement();
    //new solving element in table
    buffer_Table[r][k] = 1 / solvingElement;
    //solving row
    for(int i = 0; i < _columns; i++) {
        if(i!=k)
            buffer_Table[r][i] = buffer_Table[r][i] / solvingElement;
        else
            continue;
    }
    //solving column
    for(int i = 0; i < _rows; i++) {
        if(i!=r)
            buffer_Table[i][k] = (buffer_Table[i][k] / solvingElement) * (-1);
        else
            continue;
    }
    //rest of table
    for(int i = 0; i < _rows; i++){
        for(int j = 0; j < _columns; j++){
            if(j!=k && i!=r){
                double buf = buffer_Table[i][j] - ( (_table[i][k]*_table[r][j])/solvingElement );
                buffer_Table[i][j] = buf;
            }
        }
    }
    _table = buffer_Table;
    
    std::swap(_variablesRow[k - 1] , _variablesColumn[r]);
    print();
}
// Поиск разрешающего элемента для оптимального решения
int Simplex::findSolvingColumn(){
    int k = 1;
    std::vector<float> buf;
    for(int i = 1; i < _columns; i++) {
        buf.push_back(_table[_rows-1][i]);
    }
    k = findMax(buf) + 1;
    _k = k;
    return k;
}
int Simplex::findSolvingRow() {
    int k = findSolvingColumn();
    int r = 0;
    std::vector<float> buf;
    for(int i = 0; i < _rows - 1; i++){
        if(_table[i][k] != 0)
            buf.push_back(_table[i][0] / _table[i][k]);
        else buf.push_back(-1);
    }
    
    try {
        r = findMinPositive(buf);
        if (r == -1)
            throw k;
    } catch (...) {
        std::cerr<<"No positive values in column "<< k <<", No optimal solution for LPP"<<std::endl;
        exit(0);
    }
    
    _r = r;
    
    return r;
}
// Итерация нахождения оптимального решения
void Simplex::makeIteration() {
    std::vector<std::vector<float> > buffer_Table = _table;
    int r = findSolvingRow();
    int k = findSolvingColumn();
    printBasisSolition();
    printSolvingElement();
    float solvingElement = getSolvingElement();
    //new solving element in table
    buffer_Table[r][k] = 1 / solvingElement;
    //solving row
    for(int i = 0; i < _columns; i++) {
        if(i!=k)
            buffer_Table[r][i] = buffer_Table[r][i] / solvingElement;
        else
            continue;
    }
    //solving column
    for(int i = 0; i < _rows; i++) {
        if(i!=r)
            buffer_Table[i][k] = (buffer_Table[i][k] / solvingElement) * (-1);
        else
            continue;
    }
    //rest of table
    for(int i = 0; i < _rows; i++){
        for(int j = 0; j < _columns; j++){
            if(j!=k && i!=r){
                double buf = buffer_Table[i][j] - ( (_table[i][k]*_table[r][j])/solvingElement );
                buffer_Table[i][j] = buf;
            }
        }
    }
    _table = buffer_Table;
    
    std::swap(_variablesRow[k - 1] , _variablesColumn[r]);
    print();
}
// Нахождение опорного решения
void Simplex::findBasisSolution() {
    while(isSColumnHasNegativeValues()){
        std::cout<<"######################################################################################################"<<std::endl;
        std::cout<<"Iteration #"<<_iterationNumber<<std::endl;
        makeIterationBasis();
        _iterationNumber++;
    }
}
// Нахождение оптимального решения
void Simplex::findOptimalSolution() {
    while(isFRowHasPositiveValues()){
        std::cout<<"######################################################################################################"<<std::endl;
        std::cout<<"Iteration #"<<_iterationNumber<<std::endl;
        makeIteration();
        _iterationNumber++;
    }
}
// Решение симплекс таблицы
void Simplex::resolve() {
    makeTable();
    print();
    if (isSColumnHasNegativeValues()) {
        findBasisSolution();
    }
    if (isFRowHasPositiveValues()) {
        findOptimalSolution();
    }
}
// Вывод информации
void Simplex::printSolvingElement(){
    std::cout<<"Solving Element = "<<getSolvingElement()<< "(row = " << _r + 1 << ", column = " << _k <<")"<< std::endl;
}
void Simplex::print(){
    std::cout<<std::endl;
    
    std::cout<<"+--+"<<"----------+";
    for (int i = 0; i < _variablesNumber; i++)
        std::cout<<"----------+";
    std::cout<<std::endl;
    
    std::cout<<"|  |"<<std::setw(11)<<"Si0|";
    for (int i = 0; i < _variablesNumber; i++)
        std::cout<<std::setw(9)<<_variable<<_variablesRow[i]<<"|";
    std::cout<<std::endl;
    
    for (int i = 0; i < _restrictionNumber; i++) {
        std::cout<<"+--+"<<"----------+";
        for (int i = 0; i < _variablesNumber; i++)
            std::cout<<"----------+";
        std::cout<<std::endl;
        std::cout<<"|"<<_variable<<_variablesColumn[i]<<"|";
        for (int j = 0; j < _columns; j++)
            std::cout<<std::setw(10)<<_table[i][j]<<"|";
        std::cout<<std::endl;
    }
    std::cout<<"+--+"<<"----------+";
    for (int i = 0; i < _variablesNumber; i++)
        std::cout<<"----------+";
    std::cout<<std::endl;
    std::cout<<"|"<<_functionChar<<"|";
    for (int i = 0; i < _columns; i++)
        std::cout<<std::setw(10)<<_table[_rows-1][i]<<"|";
    std::cout<<std::endl;
    std::cout<<"+--+"<<"----------+";
    for (int i = 0; i < _variablesNumber; i++)
        std::cout<<"----------+";
    std::cout<<std::endl;
}
void Simplex::printBasisSolition() {
    if (isSColumnHasNegativeValues() == false) {
        std::cout<<"(Basis solution) ";
        for(int i = 0; i < _columns-1; i++) {
            std::cout<<_variable<<_variablesRow[i]<<" = 0, ";
        }
        for(int i = 0; i < _rows-1; i++) {
            std::cout<<_variable<<_variablesColumn[i]<<" = "<<_table[i][0]<<", ";
        }
        if (_funcTendention == MAX)
            std::cout<< std::endl <<_functionChar<<"(max) = "<<-_table[_rows-1][0] << std::endl;
        else
            std::cout<< std::endl <<_functionChar<<"(min) = "<<_table[_rows-1][0] << std::endl;
    }
    else {
        std::cout<<"No basis solution"<<std::endl;
    }
}
void Simplex::printOptimalSolution(){
    if (isFRowHasPositiveValues() == false) {
        std::cout<<"(Optimal solution) ";
        for(int i = 0; i < _columns-1; i++) {
            std::cout<<_variable<<_variablesRow[i]<<" = 0, ";
        }
        for(int i = 0; i < _rows-1; i++) {
            std::cout<<_variable<<_variablesColumn[i]<<" = "<<_table[i][0]<<", ";
        }
        if (_funcTendention == MAX)
            std::cout<< std::endl <<_functionChar<<"(max) = "<<-_table[_rows-1][0] << std::endl;
        else
            std::cout<< std::endl <<_functionChar<<"(min) = "<<_table[_rows-1][0] << std::endl;
    }
    else {
        std::cout<<"No optimal solution"<<std::endl;
    }
}
