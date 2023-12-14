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
        bool yannakis;    // yannakis encoding by default
        bool sparsify_coeff;

        // Parameters
        int _N;
        int _N_edges;
        int _N_dedges;
        int _N_end;
        int _T;
        int _threshold;

        string _net_folder;
        SupplyNet _network;
        char _output_dir[1024];
        vector<int> _end_nodes;
        vector<vector<int>> _edges;
        vector<vector<int>> _edge_map;
        vector<vector<int>> _dedges;
        vector<vector<int>> _dedge_map;

        int _prec_prob, _prec_cap;    // precision represented by the number of bits
        int _prec_cst, _prec_bgt;

        // Objective Variables
        IloBoolVarArray _var_select;  // edge selection
        IloConstraintArray _const_budget;

        // Disaster Probability
        vector<IloBoolVarArray> _var_disaster; // disaster edges
        vector<vector<IloBoolVarArray>> _var_prob_dis;   // [T, N_factors, prec_prob] discretization vars
        vector<vector<IloBoolVarArray>> _var_prob_add;   // [T, N_factors, table_size] added vars
        vector<IloConstraintArray> _const_prob_dis; // constraints on discretization variables
        vector<vector<IloConstraintArray>> _const_prob_add;  // constraints on added variables, T, N_factors, vars
        vector<vector<IloNumExpr>> _expr_prob;           // [T, N_factors] probability expression


        // Connection
        vector<IloBoolVarArray> _var_node_conn;    // indicator of node connection
        vector<IloBoolVarArray> _var_edge_conn;    // indicator of node connection
        vector<IloConstraintArray> _const_connect;


        // Capacity
        vector<IloBoolVarArray> _var_cap_dis;   //
        vector<IloConstraint> _const_cap_dis;
        vector<IloNumExpr> _expr_cap;   // should be related to connection variables


        // XOR Constraints
        vector<IloBoolVarArray> _var_all_inxor;     // all variables appear in XOR
        vector<IloIntVarArray> ivars_xor;
        vector<IloBoolVarArray> bvars_xor;
        vector<IloConstraintArray> const_xor;

        // Majority
        IloBoolVarArray _vars_maj;
        IloConstraint _const_majority;
        vector<IloConstraint> _const_all_union;  // all constraints including majority

        SupplyChain();
        ~SupplyChain();

        // Functions
        void initializeVariables();
        void genProbConstraints();
        void genConnectionConstraints();
        void genCapacityConstraints();
        void genBudgetConstraints();
        void genXORConstraints();
        void genMajorityConstraints();

        // utils
        bool makeHashFuncSolvable(vector<vector<bool>> &coeffA);
        void extractXorVarConst(vector<vector<bool>> coeffA, int t);

        // running pipeline:
        void loadParameters(const string& net_folder, int threshold, int _T, char output_dir[]);
        void genAllConstraints();
        void prepareModel();
        bool solveInstance();

};

#endif