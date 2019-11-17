#ifndef Tools_hpp
#define Tools_hpp

#include <vector>
#include <iostream>

int findMax(std::vector<float>& vec);
int findMinPositive(std::vector<float>& vec);

void transpose(std::vector< std::vector<float> >& vec, int size);

void fill(std::vector<int>& vec, int n);
void fill(std::vector<std::vector<float> >& vec, int n, int k);

void printTask();
#endif /* tools_hpp */
