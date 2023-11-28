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

        ////////NEW VERSION/////////
        // CPLEX SOLVER
        IloEnv env;
        IloModel model;
        IloCplex cplex;
        IloInt timelimit;

        // Parameters
        int _N;
        int _N_edges;
        int _N_dedges;
        int _N_end;
        int _T;
        int _target;
        int _N_inedge;
        string _net_folder;
        SupplyNet _network;
        char _output_dir[1024];
        vector<int> _end_nodes;
        vector<vector<int>> _edges;
        vector<vector<int>> _edge_map;

        // Variables
        IloBoolVarArray _var_select;  // edge selection
        vector<IloBoolVarArray> _var_disaster; // disaster edges

        vector<IloBoolVarArray> _var_prob_dis;
        vector<IloBoolVarArray> _var_cap_dis;
        vector<IloBoolVarArray> _var_inedge;    //

        // extra variables
        vector<IloBoolVarArray> _var_node_conn;    // indicator of node connection
        vector<IloBoolVarArray> _var_prob_add;  // single


        // Constraints
        IloConstraintArray _const_budget;


        // XOR Constraints

        // functions

        void genConnectionConstraints();
        void genProbConstraints();
        void genCapacityConstraints();


        ////////END NEW VERSION/////////

        int _M;

        // Variables
        vector<vector<IloBoolVarArray>> supply_selection;
        vector<vector<vector<IloBoolVarArray>>> bvars;
        vector<vector<vector<IloConstraintArray>>> constraints;
        IloConstraintArray const_budget;

        // Constraints for XOR
        bool sparsify_coeff;
        bool yannakis;    // yannakis encoding by default
        vector<vector<vector<IloConstraintArray>>> const_xor;
        vector<vector<vector<IloIntVarArray>>> ivars_xor;
        vector<vector<vector<IloBoolVarArray>>> bvars_xor;


        SupplyChain();
        ~SupplyChain();

        // running pipeline
        void loadParameters(const string& net_folder, int target, int _T, char output_dir[]);
        void genAllConstraints();
        bool solveInstance();

        // separate functions
        void initializeVariables();
        void genSupplyConstraints();

        void genXORConstraints();     
        void prepareModel();   

        // utils
        bool makeHashFuncSolvable(vector<vector<bool>> &coeffA);
        void extractXorVarConst(vector<vector<bool>> coeffA, int t, int n, int m);


        // Ignore majority in this application
        // IloBoolVarArray bvars_maj;      
        // IloConstraint const_majority;
        // void genMajorityConstraints();
};





#endif