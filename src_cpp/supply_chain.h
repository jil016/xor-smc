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
#include "supply_net.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


class SupplyChain {
    public:
         // Parameters
        int _N;
        int _M;
        int _T;
        SupplyNet _network;
        char _output_dir[1024];

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
        
        IloBoolVarArray supplier_selection;

        vector<vector<vector<IloBoolVarArray>>> bvars;  // T * N * M * IloBoolVarArray
                                                        // some materials are empty
                                                        // only those required
                                                        // capacity variables
                                                        // theta variables
                                                        // ...............

        
        // contract variables
        vector<IloBoolVarArray> contract; 

        IloBoolVarArray bvars_maj;      // May ignore majority for now

        vector<vector<int>> index_selection;

        vector<vector<vector<int>>> index_supplies;     // Index of supplies in bvars
                                                        // _T * _N * _N the spply from j to i is 
                                                        // bvars[t][i][m][index_supplies[t][j][i]]. 
            
        vector<vector<vector<vector<vector<int>>>>> index_disasters;    // Index of disaster variables
                                                                            // _T * _N * _M * _Nd * dim(disaster)
                                                                            // look up disasters in bvars.
        
        vector<vector<vector<vector<int>>>> index_capacities;   // _T * _N * _M * dim(capacity)



        // Constraints

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
        void loadParameters(string net_folder, int n_disaster, int _T,  char output_dir[]);
        void genAllConstraints();
        bool solveInstance();


        // separate functions
        void initializeVariables();
        void genSupplyConstraints();
        
        void genCapacityConstraints();
        void genBudgetConstraints();

        void genMajorityConstraints();
        void genXORConstraints();     
        void prepareModel();   

        // utils
        bool makeHashFuncSolvable(vector<vector<bool>> &coeffA);
        void extractXorVarConst(vector<vector<bool>> coeffA, int t, int s);
};





#endif