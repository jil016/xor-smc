#include "supply_chain.h"

#define DEBUG_SUPPLYCHAIN
using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


SupplyChain::SupplyChain(){
    // Initialize
    //   timer = IloTimer(env);
    //   model = IloModel(env);
    //   cplex = IloCplex(env);
	env.end();
	
	env = IloEnv();
	model = IloModel(env);
	cplex = IloCplex(model);
    
    sparsify_coeff = true;
    yannakis =true;

    // global parameters, etc.
    timelimit = -1;
}

SupplyChain::~SupplyChain(){
    env.end();
}

// start modification
void SupplyChain::loadParameters(const string& net_folder, int target, int T, char output_dir[]){
    _net_folder = net_folder;
    _prec_prob = 4;
    _prec_cap = 6;
    _prec_cst = 4;   // not used in our application
    _prec_bgt = 4;   // not used

    _network = SupplyNet(net_folder, _prec_cap, _prec_cst, _prec_bgt, _prec_prob);
    _N = _network._N;
    _N_edges = _network._N_edges;
    _N_dedges = _network._N_dedges;
    _N_end = _network._N_end;   // suppose single for now!

    _target = target;   // number of xor constraints
    _end_nodes = _network._end_nodes;

    _edges = _network._edges;
    _edge_map = _network._edge_map;

    _dedges = _network._dedges;
    _dedge_map = _network._dedge_map;

    _T = T;

    _N_inedge = 0;
    for(auto & row : _edge_map){
        if(row[_target] != 0)
            _N_inedge++;
    }

    strcpy(_output_dir, output_dir);
    initializeVariables();
}

void SupplyChain::initializeVariables(){
    // Edge selection related
    _var_select = IloBoolVarArray(env, _N_edges);    // trade plan
    _const_budget = IloConstraintArray(env);


    // Add disaster variables
    for (int t = 0; t < _T; t++){
        _var_disaster.emplace_back(env, _network._N_dedges);  // disaster variables \theta
    }

    // Add variables for probability calculation
    _var_prob_dis.resize(_T);
    _var_prob_add.resize(_T);
    _expr_prob.resize(_T);
    for (int t = 0; t < _T; t++){
        _const_prob_dis.emplace_back(env);  // disaster variables \theta
    }
    _const_prob_add.resize(_T);  // constraints on added variables, T, N_factors, vars


    // add assist variables for node connection
    for (int t = 0; t < _T; t++){
        _var_node_conn.emplace_back(env, _N);    // node connection, indicates the accessibility from primary suppliers
        _var_edge_conn.emplace_back(env, _N_edges);
    }
    for (int t = 0; t < _T; t++){
        _const_connect.emplace_back(env);   // constraints for connectivity
    }

    // add capacity variables
    _var_cap_dis.resize(_T);

    // add majority variables
    _vars_maj = IloBoolVarArray(env, _T);


    /// Begin set variable names
    for (int i = 0; i < _N; i++) {
        for (int j = 0; j < _N; j++) {
            if (_edge_map[i][j] != -1) {
                char *name = new char[32];
                snprintf(name, 32, "s_%d_%d", (int) i, (int) j);
                _var_select[_edge_map[i][j]].setName(name);
            }
        }
    }
    // set node, edge connectivity assist variable names
    // set disaster variable names
    for (int t = 0; t < _T; t++) {
        for (int i = 0; i < _N; i++) {
            char *name = new char[32];
            snprintf(name, 32, "v_t%d_%d", (int) t, (int) i);
            _var_node_conn[t][i].setName(name);

            for (int j = 0; j < _N; j++) {
                if (_edge_map[i][j] != -1) {
                    snprintf(name, 32, "_e_t%d_%d_%d", (int) t, (int) i, (int) j);
                    _var_edge_conn[t][_edge_map[i][j]].setName(name);
                }

                if (_edge_map[i][j] != -1 && _dedge_map[i][j] != -1) {
                    snprintf(name, 32, "de_t%d_%d_%d", (int) t, (int) i, (int) j);
                    _var_disaster[t][_dedge_map[i][j]].setName(name);
                }
            }
        }
    }

    for(int t = 0; t < _T; t++){
        char *name = new char[32];
        snprintf(name, 32, "mj_t%d", (int) t);
        _vars_maj[t].setName(name);
    }
    /// End set names
}


void SupplyChain::genProbConstraints(){
    for (int t = 0; t < _T; t++){
        // read through UAI
        for(int i = 0; i < _network._N_factors; i++){
            IloNumExpr expr_i(env);
            IloBoolVarArray var_add_i(env);
            IloConstraintArray const_i(env);

            // TODO: Now assume variables in different factors
            int f_size = _network._factors[i].size();   // number of variables
            int tab_size = _network._tables[i].size();  // number of terms in this factor

            // then we need tab_size == 2^f_size number of variables
            // sort by low to high binary
            for (int j = 0; j < tab_size; j++){
                char *name = new char[32];
                snprintf(name, 32, "mu_t%d_f%d_%d", (int) t, (int) i, (int) j);

                IloBoolVar mu_i_j (env, 0, 1, name);
                var_add_i.add(mu_i_j);
                expr_i += _network._tables[i][j] * mu_i_j;
            }
            _expr_prob[t].push_back(expr_i);

            // also add constraints for mu_i_j
            for (int k = 0; k < f_size; k++) {
                IloNumExpr add_pos_expr(env);
                IloNumExpr add_neg_expr(env);
                // for each variable
                for (int j = 0; j < tab_size; j++){
                    if((j >> (f_size - k - 1)) % 2 == 1 ){
                        add_pos_expr += var_add_i[j];
                    }
                    else{
                        add_neg_expr += var_add_i[j];
                    }
                }
                const_i.add(add_pos_expr == _var_disaster[t][_network._factors[i][k]]);
                const_i.add(add_neg_expr == 1 - _var_disaster[t][_network._factors[i][k]]);
            }
            _var_prob_add[t].push_back(var_add_i);
            _const_prob_add[t].push_back(const_i);
        }

        _var_prob_dis[t].resize(_network._N_factors);
        for(int i = 0; i < _network._N_factors; i++){
            _var_prob_dis[t][i] = IloBoolVarArray(env, _prec_prob);
            IloNumExpr expr_i(env);

            for(int j = 0; j < _prec_prob; j++){
                expr_i += pow(2, j) * _var_prob_dis[t][i][j];
            }

            _const_prob_dis[t].add(expr_i <= ((int) pow(2, _prec_prob) - 1) * _expr_prob[t][i] + 1);
        }
    }
#ifdef DEBUG_SUPPLYCHAIN
    IloModel model_prob(env);
    IloCplex cplex_prob(model_prob);
    for (int t = 0; t < _T; t++){
        model_prob.add(_var_disaster[t]);
        for (const auto & v : _var_prob_add[t])
            model_prob.add(v);
        model_prob.add(_const_prob_dis[t]);
        for (const auto & c : _const_prob_add[t])
            model_prob.add(c);
    }
    cplex_prob.exportModel("model_prob.lp");
#endif
}

void SupplyChain::genConnectionConstraints(){
    for (int t = 0; t < _T; t++){
        for (int i = 0; i < _N; i++) {
            // Node connection constraints
            int in_degree = 0;
            IloConstraintArray in_const(env);

            for (int j = 0; j < _N; j++) {
                // traverse the whole edge map
                if (_edge_map[i][j] != -1) {
                    // outgoing edge, set edge connectivity
                    IloConstraintArray out_const(env);
                    out_const.add(_var_node_conn[t][i] == 1);    // starting node is connected
                    out_const.add(_var_select[_edge_map[i][j]] == 1);    // is selected
                    if(_dedge_map[i][j] != -1){
                        out_const.add(_var_disaster[t][_dedge_map[i][j]] == 0);    // no disaster
                    }
                    IloAnd out_and(env);
                    out_and.add(out_const);
                    _const_connect[t].add(_var_edge_conn[t][_edge_map[i][j]] == out_and);
                }
                if (_edge_map[j][i] != -1) {
                    // incoming edge, set node connectivity
                    in_degree++;
                    in_const.add(_var_edge_conn[t][_edge_map[j][i]] == 1);
                }
            }
            if (in_degree == 0) {
                // primary supplier, always connected
                _const_connect[t].add((_var_node_conn[t][i] == 1));
            }
            else{
                IloOr in_or(env);
                in_or.add(in_const);
                _const_connect[t].add(_var_node_conn[t][i] == in_or);
                in_degree = 0;
            }
        }
    }
#ifdef DEBUG_SUPPLYCHAIN
    IloModel model_cnt(env);
    IloCplex cplex_cnt(model_cnt);
    model_cnt.add(_var_select);

    for(const auto & vars : _var_node_conn){
        model_cnt.add(vars);
    }
    for(const auto & vars : _var_edge_conn){
        model_cnt.add(vars);
    }
    for(const auto & vars : _var_disaster){
        model_cnt.add(vars);
    }
    for(const auto & cons : _const_connect){
        model_cnt.add(cons);
    }
    cplex_cnt.exportModel("model_connect.lp");
#endif
}

void SupplyChain::genCapacityConstraints(){
    for (int t = 0; t < _T; t++){
        // Capacity, only for those income edges
        vector<int> in_edges;
        vector<int> in_edge_cap;

        IloNumExpr expr_cap(env);

        for (const auto & node : _end_nodes){
            for(int i = 0; i < _N; i++){
                if(_edge_map[i][node] >= 0){    // detected one incoming edge of one end node
                    in_edges.push_back(_edge_map[i][node]); // only for debug purpose
                    in_edge_cap.push_back(_network._dis_capacity[i][node]); // only for debug purpose

                    expr_cap += _var_edge_conn[t][_edge_map[i][node]] * _network._dis_capacity[i][node];
                }
            }
        }
        _expr_cap.push_back(expr_cap);

        // create discretization variables
        _var_cap_dis[t] = IloBoolVarArray(env, _prec_cap);
        IloNumExpr cap_decode_dis(env);
        for(int i = 0; i < _prec_cap; i++){
            cap_decode_dis += pow(2, i) * _var_cap_dis[t][i];
        }
        _const_cap_dis.push_back(cap_decode_dis <= expr_cap);
    }

    #ifdef DEBUG_SUPPLYCHAIN
        IloModel model_cap(env);
        IloCplex cplex_cap(model_cap);
        for(int t = 0; t<_T;t++){
            model_cap.add(_var_edge_conn[t]);
            model_cap.add(_const_cap_dis[t]);
        }
        cplex_cap.exportModel("model_cap.lp");
    #endif
}

void SupplyChain::genBudgetConstraints(){
    for(int i = 0; i < _N; i++){
        // each node i
        IloNumExpr expr_cst(env);
        int in_degree = 0;
        for(int j = 0; j < _N; j++){
            // look at edge (j,i)
            if(_edge_map[j][i] >= 0){
                expr_cst += _var_select[_edge_map[j][i]] * _network._raw_cost[j][i];
                in_degree++;
            }
        }
        if(in_degree > 0) _const_budget.add(expr_cst <= _network._raw_budget[i]);
    }
    #ifdef DEBUG_SUPPLYCHAIN
        IloModel model_bgt(env);
        IloCplex cplex_bgt(model_bgt);
        model_bgt.add(_var_select);
        model_bgt.add(_const_budget);
        cplex_bgt.exportModel("model_bgt.lp");
    #endif
}

void SupplyChain::genXORConstraints(){
    // extract all variables that appear in XOR
    for( int t=0; t < _T; t++){
        _var_all_inxor.emplace_back(env);
        _var_all_inxor[t].add(_var_disaster[t]);
        for (int i = 0; i < _network._N_factors; i++){
            _var_all_inxor[t].add(_var_prob_dis[t][i]);
        }
        _var_all_inxor[t].add(_var_cap_dis[t]);
    }

    for( int t=0; t < _T; t++) {
        int n_vars = _var_all_inxor[t].getSize();
        int n_parity = _network._demand;

        std::cout << ">> Number of variables in XOR is: " << n_vars << endl;
        vector<vector<bool>> coeffA = generate_Toeplitz_matrix(n_parity, n_vars);
//        vector<vector<bool>> coeffA = generate_matrix(n_parity, n_vars);
        while(!makeHashFuncSolvable(coeffA)){
            coeffA = generate_Toeplitz_matrix(n_parity, n_vars);
//            coeffA = generate_matrix(n_parity, n_vars);
        }
        print_matrix(coeffA);
        extractXorVarConst(coeffA, t);
    }
#ifdef DEBUG_SUPPLYCHAIN
    IloModel model_xor(env);
    IloCplex cplex_xor(model_xor);

    for(int t = 0; t<_T;t++){
        model_xor.add(_var_all_inxor[t]);
        model_xor.add(ivars_xor[t]);
        model_xor.add(bvars_xor[t]);
        model_xor.add(const_xor[t]);
    }
    cplex_xor.exportModel("model_xor.lp");
#endif
}

void SupplyChain::genAllConstraints(){
    // generate everything
    genProbConstraints();
    cout << "Disaster constraints generated!" << endl;
    genConnectionConstraints();
    cout << "Connection constraints generated!" << endl;
    genCapacityConstraints();
    cout << "Capacity constraints generated!" << endl;
    genBudgetConstraints();
    cout << "Budget constraints generated!" << endl;
    genXORConstraints();
    cout << "XOR constraints generated!" << endl;
    genMajorityConstraints();
    cout << "Majority constraints generated!" << endl;
}


void SupplyChain::prepareModel(){

    // Objective Variables
    model.add(_var_select);
    model.add(_const_budget);
    cout << "Selection variables and budget constraints added!" << endl;

    // Disaster Probability
    for(int t = 0; t < _T; t++){
        model.add(_var_disaster[t]);

        for (int i = 0; i<_network._N_factors; i++){
            model.add(_var_prob_dis[t][i]);
            model.add(_var_prob_add[t][i]);
            // model.add(_const_prob_add[t][i]);
        }
        // model.add(_const_prob_dis[t]);
    }
    cout << "Disaster variables added!" << endl;


    // Connection
    for(int t = 0; t < _T; t++){
        model.add(_var_node_conn[t]);
        model.add(_var_edge_conn[t]);
        // model.add(_const_connect[t]);
    }
    cout << "Connection variables added!" << endl;

    // Capacity
    for(int t = 0; t < _T; t++){
        model.add(_var_cap_dis[t]);
        // model.add(_const_cap_dis[t]);
    }
    cout << "Capacity variables added!" << endl;

    // XOR Constraints
    for(int t = 0; t < _T; t++){
        model.add(ivars_xor[t]);
        model.add(bvars_xor[t]);
        model.add(const_xor[t]);
    }
    cout << "XOR variables and constraints added!" << endl;

    // Majority and all constraints
    model.add(_vars_maj);
    model.add(_const_majority);
    for(int t = 0; t < _T; t++){
        model.add(_const_all_union[t]);
    }

    cout << "Majority constraints (union with all constraints) added!" << endl;

    // DEBUG
    cplex.exportModel("model_all_w_maj.lp");
}

bool SupplyChain::solveInstance() {
    cout << "\n================== Set Output Path ===================\n" << endl;
//    #ifdef DEBUG_SUPPLYCHAIN
//    time_t now = time(0);
//    char* date_time = ctime(&now);
//    char log_folder_name[100] = "LOG-Shelter_";
//    strcat(log_folder_name, date_time);
//    #else
    filesystem::path log_folder(_output_dir);
//    #endif

    if (!filesystem::is_directory(log_folder) || !filesystem::exists(log_folder)) { // Check if src folder exists
        filesystem::create_directory(log_folder); // create src folder
    }

    filesystem::path result_log_path(log_folder);
    filesystem::path cplex_log_path(log_folder);
    filesystem::path params_log_path(log_folder);

    result_log_path = log_folder / "result.log";
    cplex_log_path = log_folder / "cplex.log";
    params_log_path = log_folder / "params.log";

    cout << "log file path:"<< endl;
    cout << result_log_path << endl;

    ofstream fs_res(result_log_path);
    ofstream fs_cplex(cplex_log_path);
    ofstream fs_params(params_log_path);
    env.setOut(fs_res);
    cplex.setOut(fs_cplex);
    cout << "\n======================================================\n" << endl;

    fs_params << "N: " << _N << "\n"
              << "T: " << _T << "\n"
              << "demand: " << _network._demand << "\n"
              << "Data path: " << _net_folder << endl;

    fs_params.close();
    
    prepareModel();

    // Cplex parameters
    cplex.setParam(IloCplex::Param::WorkMem, 2048);
    cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 2048);
    cplex.setParam(IloCplex::Threads, 4);    // number of parallel threads (automatic by default)

    bool solved = cplex.solve();
    cplex.exportModel("model.lp");
    // env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << cplex.getCplexStatus() << endl;

    if(solved){
        // read the trade selection
        vector<vector<int>> parse_result;
        parse_result.resize(_N);
        for (int n = 0; n < _N; n++){
            parse_result[n].resize(_N);
            fill(parse_result[n].begin(), parse_result[n].end(), 0);
        }

        for (int i = 0; i < _N; i++){
            for (int j = 0; j < _N; j++){
                if(_edge_map[i][j] >= 0){
                    parse_result[i][j] = int(cplex.getValue(_var_select[_edge_map[i][j]]));
                }
            }
        }
        
        for(int i = 0; i<_N; i++){
            for(int j = 0; j<_N; j++){
                env.out() << parse_result[i][j] << " ";
            }
            env.out() << endl;
        }

        for(int t = 0; t < _T; t++){
            env.out() << int(cplex.getValue(_vars_maj[t])) << " ";
        }
    }

    else{
        env.out() << "NO SOLUTION " << endl;
    }

    return true;
}


bool SupplyChain::makeHashFuncSolvable(vector<vector<bool>> &coeffA) {
    if (coeffA.empty()) 
        return true;

    bool solvable = true;            // is the system A x + b =0 solvable?
    
    size_t m = coeffA.size();
    size_t n = coeffA[0].size()-1;
    
    vector <int> indep_columns;
    vector <int> indep_columns_rindex;
    set <int> indep_vars;
        
    // put A in row echelon form
    for (size_t i = 0;i<m;i++) {
        //Find pivot for column k:
        bool empty_row=true;
        int j_max =0;
        for (int s = 0;s<n;s++)
        if (coeffA[i][s]) {
            empty_row=false;
            j_max = s;
            break;
        }
        if (empty_row){        // low rank
        if (coeffA[i][n]) {    //0=1
            solvable = false;
            cout << "Sampled Ax = b is not solvable!" << endl;
            return false;
        }
        }
        else {
        indep_vars.insert(j_max);
        indep_columns.push_back(j_max);          // index of a basis of coeffA
        indep_columns_rindex.push_back(i);        // row index of pivot
            
        for (size_t h=i+1;h<m;h++)
            if (coeffA[h][j_max]) {      // if not already zero
            for (int q=0;q<n+1;q++)      // sum the two rows
                coeffA[h][q] = coeffA[h][q] xor coeffA[i][q];
            }
        }
    }
    
    for (size_t i = 0;i<indep_columns.size();i++) {
        int j_max = indep_columns[i];
        int p = indep_columns_rindex[i];

        for (int h=p-1;h>=0;h--)
        if (coeffA[h][j_max]) {      // if not already zero
            //print_matrix(coeffA);

            for (int q=0;q<n+1;q++)      // sum the two rows
                coeffA[h][q] = coeffA[h][q] xor coeffA[p][q];
            //print_matrix(coeffA);
            }
    }
    
    // produce a solution
    vector <bool>  b;
    vector <bool>  y;
    y.resize(n);
    
    // initialize b to the last column of A
    b.resize(m);
    for (size_t i =0;i<m;i++)
        b[i] = coeffA[i][n];
    
    for (size_t i =0;i<n;i++)
        y[i] = rand()%2;
        
    // sum all the dependent variables that are already fixed
    for (size_t i =0;i<n;i++) {
        if ( (indep_vars.count(i)==0) && (y[i]==1)) {    // dependent variable, and non zero
        // b = b + x[i] * A[] [i]
        for (size_t j =0;j<m;j++)
            b[j] = b[j] xor coeffA[j][i];
        }
    }
        
    // backsubstitute r
    for (int i =indep_columns_rindex.size()-1;i>=0;i--) {
        int c = indep_columns_rindex[i];    // lowest pivot
        if (b[c]==1) {    // we need to add a 1
        y[indep_columns[i]] = 1;
        for (size_t j =0;j<m;j++)
            b[j] = b[j] xor coeffA[j][indep_columns[i]];
        }
        else {
        y[indep_columns[i]] = 0;
        }
    }

    // print_matrix(coeffA);


    #ifdef DEBUG_SUPPLYCHAIN
    cout << "row echelon form:" << endl;
    print_matrix(coeffA);
    #endif
    
    // sparsify
    if (sparsify_coeff && (!coeffA.empty())) {
        cout << "Bits saved: ";
        for (int i=0; i<2; i++)
        cout << sparsify(coeffA) << " ";
        cout << endl;  
    }
    return true;
}

void SupplyChain::extractXorVarConst(vector<vector<bool>> coeffA, int t) {
    // Encode XOR constraints, and extract corresponding variables and constraints
    IloConstraintArray xor_const_array(env);
    IloIntVarArray xor_int_vars_array(env);
    IloBoolVarArray xor_bool_vars_array(env);
    int nbvars = _var_all_inxor[t].getSize();

    std::vector < std::set <size_t> > varAppearancesInXors;
    if (!coeffA.empty())
        varAppearancesInXors.resize(coeffA[0].size());            // dummy parity var

    IloArray<IloArray<IloIntVarArray> > zeta_vars(env);
    IloArray <IloIntVarArray> alpha_vars(env);

    IloBoolVar dummy_parity (env, 0, 1, "dummy");

    xor_bool_vars_array.add(dummy_parity);
    xor_const_array.add((dummy_parity==1));


    vector <size_t> xors_length;

    if (!coeffA.empty()) {
        xors_length.resize(coeffA.size());

        for (size_t j = 0; j < coeffA.size(); j++) {
        size_t f = 0;

        for (size_t l = 0; l < coeffA[j].size(); l++)      // last column is the parity bit b
            if (coeffA[j][l]) {
            f++;
            }

        xors_length[j] = f;    // save length of j-th xor
        }

        for (size_t j = 0; j<coeffA.size();j++) {
        // use yannakis encoding for longer ones
        for (size_t l = 0; l<coeffA[j].size();l++)      // last column is the parity bit b
            if (coeffA[j][l]) {
            varAppearancesInXors[l].insert(j);    // for each var, save list of xors involved
            }
        }

        cout << "XOR minimum length: "
            << *std::min_element(xors_length.begin(),xors_length.end())
            <<" . XOR maximum length: "
            << *std::max_element(xors_length.begin(),xors_length.end())
            << endl;

        if (yannakis) {
            for (size_t j= 0; j<coeffA.size();j++) {
                IloIntVarArray alphas(env);

                // compute xor length
                size_t f =xors_length[j];

                // add alpha_j_k var
                IloNumExpr alpha_sum_to_one(env);
                for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2) {
                char *name = new char[32];
                snprintf(name, 32, "alpha_%d_%d", (int) j, (int) k);
                IloBoolVar alpha_j_k (env, 0, 1, name);          // (18)

                alphas.add(alpha_j_k);
                alpha_sum_to_one = alpha_sum_to_one + alpha_j_k;
                }

                xor_int_vars_array.add(alphas);
                xor_const_array.add((alpha_sum_to_one==1));

                alpha_vars.add(alphas);
            }

            // add zeta_i_j_k var nbvar

            for (size_t i= 0; i < coeffA[0].size();i++) {
                IloArray<IloIntVarArray> zet_jk (env, coeffA.size());

                std::set< size_t > XorsContainingi =  varAppearancesInXors[i];
                for (std::set<size_t >::iterator it=XorsContainingi.begin(); it!=XorsContainingi.end(); ++it) {
                size_t f = xors_length[*it];

                IloIntVarArray zet(env);

                IloNumExpr zeta_sum_to_f(env);

                for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2) {
                    char *name = new char[32];
                    snprintf(name, 32, "zeta_%d_%d_%d", (int) i, (int) *it, (int) k);

                    IloBoolVar zeta_i_j_k (env, 0, 1, name);

                    zet.add(zeta_i_j_k);
                    zeta_sum_to_f = zeta_sum_to_f + zeta_i_j_k;

                    xor_const_array.add((zeta_i_j_k<=alpha_vars[*it][k/2]));
                }
                xor_int_vars_array.add(zet);

                if (i < nbvars){
                    xor_const_array.add((zeta_sum_to_f == _var_all_inxor[t][i]));
                }
                else if (i == nbvars) {
                    xor_const_array.add((zeta_sum_to_f == dummy_parity)); // dummy
                }
                zet_jk[*it]=zet;
                }

                zeta_vars.add(zet_jk);
            }

            for (size_t j= 0; j<coeffA.size();j++){
                // use yannakis encoding for longer ones
                size_t f = xors_length[j];
                for (size_t k = 0; k<= 2* (size_t) floor(f/2);k=k+2) {
                IloNumExpr zeta_sum_to_alpha_k(env);
                for (size_t l = 0; l<coeffA[j].size();l++)  {    // last column is the parity bit b
                    if (coeffA[j][l])
                    zeta_sum_to_alpha_k = zeta_sum_to_alpha_k + zeta_vars[l][j][k/2];
                }

                xor_const_array.add((zeta_sum_to_alpha_k==IloInt(k)*alpha_vars[j][k/2]));
                }
            }
        }
    }
    ivars_xor.push_back(xor_int_vars_array);
    bvars_xor.push_back(xor_bool_vars_array);
    const_xor.push_back(xor_const_array);
}


void SupplyChain::genMajorityConstraints() {
    IloNumExpr majority_sum(env);
    for (int t = 0; t < _T; t++){
        majority_sum = majority_sum + _vars_maj[t];
    }
    _const_majority = (majority_sum > _T/2);

    vector<IloAnd> const_and;
    for (int t = 0; t< _T; t++){
        // load all constraints in t-th repetition
        // Disaster Probability
        IloAnd temp_const_and(env);
        temp_const_and.add(_const_prob_dis[t]);   // constraints on discretization variables
        for(const auto & carray:_const_prob_add[t]){   // constraints on added variables, T, N_factors, vars
            temp_const_and.add(carray);
        }
        // Connection
        temp_const_and.add(_const_connect[t]);
        // Capacity
        temp_const_and.add(_const_cap_dis[t]);

        const_and.push_back(temp_const_and);
    }

    for (int t = 0; t< _T; t++){
        IloConstraint temp_const;
        temp_const = ((const_and[t]) || (_vars_maj[t] == 0));
        _const_all_union.push_back(temp_const);
    }
}
