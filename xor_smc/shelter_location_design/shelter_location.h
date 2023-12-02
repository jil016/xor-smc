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

#include <filesystem>

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
        bool allow_bf;

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

        char _output_dir[1024];


        // Variables       
        // The following variables appears in the optimization
        IloBoolVarArray shelter_assign;
        vector<vector<IloBoolVarArray>> bvars;
        IloBoolVarArray bvars_maj;


        // Index vector for variables
        vector<vector<vector<vector<int>>>> flows;  // idx = flow[t][s][i][j]; bvars[t][s][idx] represents flow from i to j for (t,s)
        vector<vector<vector<int>>> shelters;   // idx = shelters[t][s][i][j]; bvars[t][s][idx] represents the shelter selection of (t,s) 

        // Constraints
        vector<vector<IloConstraintArray>> const_flow;
        vector<vector<IloConstraint>> const_flow_union; // added constraints for majority 

        vector<vector<IloConstraintArray>> const_is_shelter;
        vector<vector<IloConstraint>> const_is_shelter_union;   // added constraints for majority 

        vector<vector<IloConstraintArray>> const_nbf;   // no back and forth

        IloConstraint const_max_shelter;
        IloConstraint const_majority;

        // Constraints for XOR
        bool sparsify_coeff;
        bool yannakis;    // yannakis encoding by default
        vector<vector<IloConstraintArray>> const_xor;
        vector<vector<IloIntVarArray>> ivars_xor;
        vector<vector<IloBoolVarArray>> bvars_xor;


        ShelterLocation();
        ~ShelterLocation();

        // running pipeline
        void loadParameters(char graph_file[], int N, int T, int M, 
                            vector<int> src, vector<int> q, char output_dir[]);
        void genAllConstraints();
        bool solveInstance();


        // separate functions
        void genFlowConstraints();
        void genShelterConstraints();
        void genMajorityConstraints();
        void genNBFConstraints();
        void genXORConstraints();     
        void prepareModel();   

        // utils
        bool makeHashFuncSolvable(vector<vector<bool>> &coeffA);
        void extractXorVarConst(vector<vector<bool>> coeffA, int t, int s);
};





#endif