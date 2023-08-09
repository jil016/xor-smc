#ifndef SUPPLY_CHAIN_H
#define SUPPLY_CHAIN_H

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


class SupplyChain {
    public:
         // Parameters
        int _N;
        int _T;
        int _M;
        Graph *graph;
        Graph *graph_theta;
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
        IloBoolVarArray factory_assign;

        vector<vector<IloBoolVarArray>> bvars;  // stacks theta, supplier, and discretization variable
                                                // T * material 


        vector<vector<IloBoolVarArray>> bvars_flow;
        IloBoolVarArray bvars_maj;      // ignore majority for now


        // supplier list (read from input)
        vector<vector<int>> supplier_locations;  // suppliers[m]: suppliers for material i
                                                // suppliers[m][j]: node j is one of the suppliers for material i


        // Index vector for variables
        vector<vector<vector<int>>> supplier_idx;   // idx = supplier[t][s][i][j]; bvars[t][s][idx] represents the shelter selection of (t,s) 
        // flows
        vector<vector<vector<vector<int>>>> flow_idx;  // idx = flow[t][s][i][j]; bvars[t][s][idx] represents flow from i to j for (t,s)
        vector<vector<vector<vector<int>>>> theta_idx; 
        

        // theta list -- for natural disaster
        // y list -- for the P(theta)
        // supplier identifier

        

         

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


        SupplyChain();
        ~SupplyChain();

        // running pipeline
        void loadParameters(char graph_file[], int N, int T, int M, 
                            vector<int> src, vector<int> q, char output_dir[]);
        void genAllConstraints();
        bool solveInstance();


        // separate functions
        void initializeVariables();
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