#ifndef WISH_H
#define WISH_H

#include <new>
#include <set>
#include <cmath>
#include <bitset>
#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>

#include <stdio.h>
#include <sys/time.h>
#include <ilcp/cpext.h>
#include <ilcplex/ilocplex.h>

// use ILOG's STL namespace
ILOSTLBEGIN


typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

class WishInstance {
  public:
    char pbname[1024];

    // From UAI file
    int nbvar,nbval,nbconstr;

    // CPLEX SOLVER
    IloEnv env;
    IloModel *model;
    IloTimer *timer;
    IloCplex *cplex;
    IloInt timelimit;

    // From arguments 
    IloInt parity_number;
    bool yannakis;    // yannakis encoding
    unsigned long seed;
    bool use_given_seed;

    IloIntVarArray *vars;
    IloNumExpr *objexpr;
    
    IloArray<IloNumArray> *cost;
    std::vector < std::vector< int> > scopes;

    // for removal
    IloConstraintArray *xor_constr_array;
    IloIntVarArray *xor_int_vars_array;
    IloBoolVarArray *xor_bool_vars_array;

    bool sparsify_coeff;
    std::vector < std::vector <bool> > coeffA;
    std::vector <bool> feasiblesol;

    char instance_name[1024];

    WishInstance();
    ~WishInstance();

    void loadInstance();
    void readInstance(char *file_path);
    void getObjectExpression();
    void addOptimizationGoal();
    void sampleHashCoeffA();
    bool getFeasibleSolution();
    void extractXorConstraints();
    bool solveInstance();
    void removeXorConstraints();
};

unsigned long get_seed(void);
void parseParityArgs(WishInstance & ins, int & argc, char **argv);
void parseArgs(WishInstance & ins, int argc, char **argv) ;
double median(vector<double> &v);
int sparsify(vector <vector <bool> > & A);
void print_matrix (vector <vector <bool> > A);
vector <vector <bool> > generate_matrix(int m, int n);
vector <vector <bool> > generate_Toeplitz_matrix(int m, int n);

#endif