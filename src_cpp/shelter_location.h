#ifndef SHELTER_LOCATION_H
#define SHELTER_LOCATION_H

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

#include "utils.h"
#include "graph.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


class ShelterLocation {
    public:
         // Parameters
        int _N;
        int _T;
        int _M;
        Graph *graph;
        vector<int> _src;
        vector<int> _q;

        // CPLEX SOLVER
        // IloEnv env;
        // IloModel *model;
        // IloTimer *timer;
        // IloCplex *cplex;
        // IloInt timelimit;
        // IloNumExpr *objexpr;
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloInt timelimit;


        // Variables       
        // The following variables appears in the optimization
        IloBoolVarArray shelter_assign;
        vector<vector<IloBoolVarArray>> bvars;
        vector<vector<vector<vector<int>>>> flows;
        vector<vector<vector<int>>> shelters;


        // Generate XOR  
        bool sparsify_coeff;
        bool yannakis;    // yannakis encoding by default
        std::vector < std::vector <bool> > coeffA;
        std::vector <bool> feasiblesol;
        
        // XOR Constraints
        int n_vars; 
        IloConstraintArray *xor_constr_array;
        IloIntVarArray *xor_int_vars_array;
        IloBoolVarArray *xor_bool_vars_array;


        ShelterLocation();
        ~ShelterLocation();

        void loadParameters(int N, int T, int M, vector<int> src, vector<int> q);
        void addFlowConstraints();
        void addShelterConstraints();
        void addXORConstraints();
        void parseAllConstraints();
        bool solveInstance();

};





#endif