#ifndef PAWS_H
#define PAWS_H

#include <string>
#include <iostream>
#include <ctime>
#include <filesystem>

#include "wish.h"

// use ILOG's STL namespace
ILOSTLBEGIN

class SInstance : public WishInstance {
  public:
    int nbvar_x;
    int nbvar_y;
    IloIntVarArray *vars_y;
    std::vector <IloNumArray> solutions;
    
    SInstance();
    SInstance(WishInstance w_ins);
    void addYVariables(int l, int b);
    void addSConstraints(IloNum M, int l, int b);
    void addNaiveOptimizationGoal();
    int findPSolutions(int P);
};

int computeK(SInstance &S_base, int n, double delta, int P);
IloNumArray paws(char *file_path, int l, int b, double delta, int P, int alpha, unsigned long seed, filesystem::path cplex_log_path);
IloNumArray pawsWrapper(char *file_path, int n, double epsilon, int b, double delta0, int P, int alpha, unsigned long seed);

#endif