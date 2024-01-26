#ifndef XOR_SMC_SMC_H
#define XOR_SMC_SMC_H

#include <new>
#include <set>
#include <cmath>
#include <bitset>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>

#include <filesystem>

#include <cstdio>
#include <sys/time.h>
#include <ilcp/cpext.h>
#include <ilcplex/ilocplex.h>

#include "utils.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

class SMC{
public:
    SMC();
    ~SMC();

    // running pipeline:
    void loadParameters(const string& cnf_file,
                        const string& uai_file,
                        const string& mar_file,
                        int threshold,
                        int T,
                        char output_dir[]);

    bool solveInstance();

protected:
    // CPLEX SOLVER
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloInt timelimit;
    bool yannakis;    // yannakis encoding by default
    bool sparsify_coeff;
    int _T;
    int _threshold;
    char _output_dir[1024];

    // CNF
    string _cnf_file;
    int _N_vars_cnf;
    vector<vector<int>> _clauses;

    // UAI
    string _uai_file;
    string _mar_file;
    int _N_vars_uai;
    vector<int> _scope_evi;
    vector<int> _scope_mar;
    int _prec_n_digits;    // precision represented by the number of bits
    int _prec_n_int;
    int _N_factors;
    vector<vector<int>> _factors;
    vector<vector<double>> _tables;

    // variables for cnf
    IloBoolVarArray _var_cnf;  // edge selection
    IloConstraintArray _const_cnf;

    // variables for uai
    vector<IloBoolVarArray> _var_uai; // disaster edges
    vector<vector<IloBoolVarArray>> _var_prob_dis;   // [T, N_factors, prec_prob] discretization vars
    vector<vector<IloBoolVarArray>> _var_prob_add;   // [T, N_factors, table_size] added vars
    vector<IloConstraintArray> _const_prob_dis; // constraints on discretization variables
    vector<vector<IloConstraintArray>> _const_prob_add;  // constraints on added variables, T, N_factors, vars
    vector<IloConstraintArray> _const_uai_evi;  // uai shares variables with cnf
    vector<vector<IloNumExpr>> _expr_prob;           // [T, N_factors] probability expression

    // XOR Constraints
    vector<IloBoolVarArray> _var_all_inxor;     // all variables appear in XOR
    vector<IloIntVarArray> ivars_xor;
    vector<IloBoolVarArray> bvars_xor;
    vector<IloConstraintArray> const_xor;

    // Majority
    IloBoolVarArray _var_maj;
    IloConstraint _const_majority;

    // All constraints
    vector<IloConstraint> _const_all_union;  // all constraints including majority


    // Functions
    void readFiles();
    void initializeVariables();
    void genCnfConstraints();
    void genUaiConstraints();
    void genXORConstraints();
    void genMajorityConstraints();

    // utils
    bool makeHashFuncSolvable(vector<vector<bool>> &coeffA);
    void extractXorVarConst(vector<vector<bool>> coeffA, int t);

    //
    void genAllConstraints();
    void prepareModel();

};
#endif //XOR_SMC_SMC_H
