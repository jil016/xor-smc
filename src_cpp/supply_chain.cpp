#include "supply_chain.h"

// #define DEBUG_SHELTER
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
void SupplyChain::loadParameters(const string& net_folder, int T, char output_dir[]){
    _net_folder = net_folder;
    int prec = 1;   // precision
    _network = SupplyNet(net_folder, prec, prec, prec, prec);
    _N = _network._N;
    _M = _network._M;
    _N_connect = _network._N_connect;
    _N_end = _network._N_end;
    _end_nodes = _network._end_nodes;
    _edges = _network._edges;
    _edge_map = _network._edge_map;
    _T = T;

    strcpy(_output_dir, output_dir);
    initializeVariables();
}

void SupplyChain::initializeVariables(){
    for (int t = 0; t < _T; t++){
        _plan.push_back(IloBoolVarArray(env));
    }

}

//void SupplyChain::initializeVariables(){
//    bvars.resize(_T);
//    for (int t = 0; t < _T; t++){
//        bvars[t].resize(_N);
//        for (int n = 0; n < _N; n++){
//            for (int m = 0; m < _M; m++){
//                bvars[t][n].push_back(IloBoolVarArray(env));
//            }
//        }
//    }
//
//    constraints.resize(_T);
//    for (int t = 0; t < _T; t++){
//        constraints[t].resize(_N);
//        for (int n = 0; n < _N; n++){
//            for (int m = 0; m < _M; m++){
//                constraints[t][n].push_back(IloConstraintArray(env));
//            }
//        }
//    }
//
//    supply_selection.resize(_N);
//    for (int n = 0; n < _N; n++){
//        for (int m = 0; m < _M; m++){
//            if(_network._demand[n][m] != 0){
//                supply_selection[n].push_back(IloBoolVarArray(env, _network._producers[m].size()));
//            }
//            else{
//                supply_selection[n].push_back(IloBoolVarArray(env));
//            }
//        }
//    }
//
//    const_budget = IloConstraintArray(env);
//}

void SupplyChain::genSupplyConstraints(){
    for (int t = 0; t < _T; t++){
        for (int n = 0; n < _N; n++) {
            // each node needs to be below its budget
            // its input must ...
            // its output must be <= its input ...
        }
    }
}

//void SupplyChain::genSupplyConstraints(){
//    //
//    for (int t = 0; t < _T; t++){
//        for (int n = 0; n < _N; n++){
//            for (int m = 0; m < _M; m++){
//                if(_network._demand[n][m] != 0){
//                    // supply
//                    IloBoolVarArray bvars_supply(env);  // represents only one supplier
//                    IloConstraintArray const_supply(env);
//
//                    // capacity
//                    IloBoolVarArray bvars_cap(env, 4);  // dim(capacity) 4 digits
//                    IloNumExpr capacity_encoded = 8*bvars_cap[0] + 4*bvars_cap[1] + 2*bvars_cap[2] + bvars_cap[3];
//
//                    // disasters
//                    vector<IloBoolVarArray> bvars_dis; // _Nd * dim(disaster) 4 digits each
//                    vector<IloNumExpr> disaster_encoded;
//                    for(int d = 0; d < _network._Nd; d++){
//                        IloBoolVarArray bvars_dis_d(env, 4);
//                        IloNumExpr dis_d_encoded(env);
//                        dis_d_encoded = 8*bvars_dis_d[0] + 4*bvars_dis_d[1] + 2*bvars_dis_d[2] + bvars_dis_d[3];
//
//                        bvars_dis.push_back(bvars_dis_d);
//                        disaster_encoded.push_back(dis_d_encoded);
//                    }
//
//                    for (int j = 0; j < _network._producers[m].size(); j ++){
//                        int p = _network._producers[m][j];
//                        // considering (p -> n) supply
//
//                        IloBoolVar bvar_pn(env);
//                        IloConstraint const_supply_pn;
//
//                        // capacity is large enough
//                        int capacity_pn = _network._capacity[p][n];  // precision = 4. this is an integer 0~16
//                        IloConstraint const_cap_pn = (capacity_encoded < capacity_pn); // capacity success
//
//                        // survived in all disasters
//                        IloConstraintArray const_dis_pn(env);
//                        for(int d = 0; d < _network._Nd; d++){
//                            // disaster pass
//                            int prob_pn_d = _network._disaster_map[d][p][n];
//
//                            const_dis_pn.add(disaster_encoded[d] >= prob_pn_d);
//                        }
//
//                        IloAnd const_dis_pn_and(env);
//                        const_dis_pn_and.add(const_dis_pn);
//
//                        // overall
//                        const_supply_pn = (bvar_pn == 0) || (const_cap_pn && const_dis_pn_and);
//
//                        bvars_supply.add(bvar_pn);
//                        const_supply.add(const_supply_pn);
//                    }
//                    // for this specific material,
//
//                    // one edge sampled a time
//                    IloNumExpr expr_one_edge(env);
//                    for (int j = 0; j < _network._producers[m].size(); j++){
//                        expr_one_edge = expr_one_edge + bvars_supply[j];
//                    }
//                    IloConstraint const_one_edge = (expr_one_edge == 1);
//
//
//                    IloConstraintArray const_be_selected(env);
//                    for (int j = 0; j < _network._producers[m].size(); j++){
//                        const_be_selected.add(bvars_supply[j] <= supply_selection[n][m][j]);
//                    }
//
//                    // push_back the following:
//
//                    // const_be_selected
//                    // const_one_edge
//                    // const_supply
//
//                    // bvars_supply
//                    // bvars_cap
//                    // bvars_dis
//
//                    bvars[t][n][m].add(bvars_supply);
//                    bvars[t][n][m].add(bvars_cap);
//                    for(int j = 0; j < bvars_dis.size(); j++){
//                        bvars[t][n][m].add(bvars_dis[j]);
//                    }
//
//                    constraints[t][n][m].add(const_be_selected);
//                    constraints[t][n][m].add(const_one_edge);
//                    constraints[t][n][m].add(const_supply);
//                }
//            }
//        }
//    }
//
//    // limit the budget
//    for (int n = 0; n < _N; n++){
//        IloNumExpr expr_cost_sum(env);
//        int budget_n = _network._budget[n];
//
//        for (int m = 0; m < _M; m++){
//            if(_network._demand[n][m] != 0){
//                for (int j = 0; j < _network._producers[m].size(); j ++){
//                    // cout<< "node :" << n << " material: " << m << " cost: " << _network._cost[j][n]<<endl;
//                    int p = _network._producers[m][j];
//                    expr_cost_sum = expr_cost_sum + supply_selection[n][m][j] * _network._cost[p][n];
//                }
//            }
//        }
//        const_budget.add(expr_cost_sum <= budget_n);
//    }
//    //
//}


void SupplyChain::prepareModel(){

    for (int n = 0; n < _N; n++){
        for (int m = 0; m < _M; m++){
            model.add(supply_selection[n][m]);
        }
    }

    for (int t = 0; t < _T; t++){
        for (int n = 0; n < _N; n++){
            for (int m = 0; m < _M; m++){
                model.add(bvars[t][n][m]);
            }
        }
    }
    cout << "Supply variables added!" << endl;

    for (int t = 0; t < _T; t++){
        for (int n = 0; n < _N; n++){
            for (int m = 0; m < _M; m++){
                model.add(constraints[t][n][m]);
            }
        }
    }
    model.add(const_budget);
    cout << "Supply constraints added!" << endl;

    // XOR constraints ...

    for (int t = 0; t < _T; t++){
        for (int n = 0; n < _N; n++){
            for (int m = 0; m < _M; m++){
                model.add(ivars_xor[t][n][m]);
                model.add(bvars_xor[t][n][m]);
                model.add(const_xor[t][n][m]);
            }
        }
    }
    cout << "XOR constraints added!" << endl;
    // vector<vector<vector<IloConstraintArray>>> const_xor;
    // vector<vector<vector<IloIntVarArray>>> ivars_xor;
    // vector<vector<vector<IloBoolVarArray>>> bvars_xor;
}


bool SupplyChain::solveInstance() {
    cout << "\n================== Set Output Path ===================\n" << endl;
    #ifdef DEBUG_SUPPLY_CHAIN
    time_t now = time(0);
    char* date_time = ctime(&now);
    char log_folder_name[100] = "LOG-Shelter_";
    strcat(log_folder_name, date_time);
    #else
    filesystem::path log_folder(_output_dir);
    #endif

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
              << "M: " << _M << "\n"
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
        // read the supply_selection
        vector<vector<int>> parse_result;
        parse_result.resize(_N);
        for (int n = 0; n < _N; n++){
            parse_result[n].resize(_N);
            fill(parse_result[n].begin(), parse_result[n].end(), 0);
        }

        for (int n = 0; n < _N; n++){
            for (int m = 0; m < _M; m++){
                if(_network._demand[n][m] != 0){
                    for (int j = 0; j < _network._producers[m].size(); j ++){
                        parse_result[_network._producers[m][j]][n] = int(cplex.getValue(supply_selection[n][m][j]));
                    }
                }

            }
        }
        
        for(int i = 0; i<_N; i++){
            for(int j = 0; j<_N; j++){
                env.out() << parse_result[i][j] << " ";
            }
            env.out() << endl;
        }
    }

    return true;
}

void SupplyChain::genXORConstraints(){
    const_xor.resize(_T);
    ivars_xor.resize(_T);
    bvars_xor.resize(_T);

    for( int t=0; t < _T; t++){
        const_xor[t].resize(_N);
        ivars_xor[t].resize(_N);
        bvars_xor[t].resize(_N);

        for(int n=0; n< _N; n++){
            for(int m = 0; m < _M; m++){
                int n_vars = bvars[t][n][m].getSize();
                int n_parity = _network._demand[n][m];  // _demand determines the number of XORs

                std::cout << ">> Number of variables in XOR is: " << n_vars << endl;
                vector<vector<bool>> coeffA = generate_Toeplitz_matrix(n_parity, n_vars);

                // coeffA: n_parity x (n_vars + 1)
                // be sure to generate a set of non-self-conflict XOR constraints
                while(!makeHashFuncSolvable(coeffA)){
                    coeffA = generate_Toeplitz_matrix(n_parity, n_vars);
                }
                extractXorVarConst(coeffA, t, n, m);
            }
        }
    }

}


void SupplyChain::genAllConstraints(){
    // generate everything
    genSupplyConstraints();
    cout << "Supply constraints generated!" << endl;
    genXORConstraints(); 
    cout << "XOR constraints generated!" << endl;
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


    #ifdef DEBUG_SHELTER
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

void SupplyChain::extractXorVarConst(vector<vector<bool>> coeffA, int t, int n, int m) {
    // Encode XOR constraints, and extract corresponding variables and constraints
    IloConstraintArray xor_const_array(env);
    IloIntVarArray xor_int_vars_array(env);
    IloBoolVarArray xor_bool_vars_array(env);
    int nbvars = bvars[t][n][m].getSize();

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
                sprintf(name, "alpha_%d_%d", (int) j, (int) k);
                
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
                    sprintf(name, "zeta_%d_%d_%d", (int) i, (int) *it, (int) k);
                    
                    IloBoolVar zeta_i_j_k (env, 0, 1, name);
                    
                    zet.add(zeta_i_j_k);
                    zeta_sum_to_f = zeta_sum_to_f + zeta_i_j_k;
                    
                    xor_const_array.add((zeta_i_j_k<=alpha_vars[*it][k/2]));
                }
                xor_int_vars_array.add(zet);
                    
                if (i < nbvars){
                    xor_const_array.add((zeta_sum_to_f == bvars[t][n][m][i]));
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
    const_xor[t][n].push_back(xor_const_array);
    ivars_xor[t][n].push_back(xor_int_vars_array);
    bvars_xor[t][n].push_back(xor_bool_vars_array);
}



// void SupplyChain::genMajorityConstraints(){
//     bvars_maj = IloBoolVarArray(env, _T);

//     IloNumExpr majority_sum(env);
//     for (int t = 0; t < _T; t++){
//         majority_sum = majority_sum + bvars_maj[t];
//     }

//     const_majority = (majority_sum > _T/2);

//     const_flow_union.resize(_T);
//     const_is_shelter_union.resize(_T);


//     vector<vector<IloAnd>> const_flow_and;
//     vector<vector<IloAnd>> const_is_shelter_and;

//     const_flow_and.resize(_T);
//     const_is_shelter_and.resize(_T);

//     for (int t = 0; t< _T; t++){
//         for (int s = 0; s< _src.size(); s++){
//             IloAnd temp_cf(env);
//             temp_cf.add(const_flow[t][s]);
//             const_flow_and[t].push_back(temp_cf);

//             IloAnd temp_is(env);
//             temp_is.add(const_is_shelter[t][s]);
//             const_is_shelter_and[t].push_back(temp_is);
//         }
//     }

//     for (int t = 0; t< _T; t++){
//         for (int s = 0; s< _src.size(); s++){
//             IloConstraint temp_const;
//             temp_const = ((const_flow_and[t][s]) || (bvars_maj[t] == 0));
//             const_flow_union[t].push_back(temp_const);
//             temp_const = ((const_is_shelter_and[t][s]) || (bvars_maj[t] == 0));
//             const_is_shelter_union[t].push_back(temp_const);
//         }
//     }
// }




//     bvars.resize(_T);

//     index_supplies.resize(_T);      // _T * _N * _N
//     index_disasters.resize(_T);     // _T * _N * _N * _Nd * dim(disaster)
//     index_capacities.resize(_T);    // _T * _N * _M * dim(capacity)

//     for (int t = 0; t < _T; t++){
//         bvars[t].resize(_N);

//         index_supplies[t].resize(_N);
//         index_disasters[t].resize(_N);
//         index_capacities[t].resize(_N);

//         for (int i = 0; i < _N; i++){
//             index_supplies[t][i].resize(_N);
//             fill(index_supplies[t][i].begin(), index_supplies[t][i].end(), -1); 

//             index_disasters[t][i].resize(_M);
//             for(int j = 0; j < _M; j++){
//                 index_disasters[t][i][j].resize(_network._Nd);
//             }

//             index_capacities[t][i].resize(_M);
//         }


//         for (int i = 0; i < _N; i++){

//             for (int m = 0; m < _M; m++){
//                 bvars[t][i].push_back(IloBoolVarArray(env));

//                 // add supply variables
//                 if(_network._demand[i][m] != 0){
//                     // node i requires _demand[i][m] unit of material m
//                     // will create a bunch of variables for those trade edges
//                     for (int p = 0; p < _network._producers[m].size(); p ++){
//                         // node _network._producers[m][p] produces material m
//                         IloBoolVar var_supply(env);
//                         bvars[t][i][m].add(var_supply);

//                         int idx = _network._producers[m][p];
//                         index_supplies[t][idx][i] = bvars[t][i][m].getSize() - 1;
//                     }
//                 }

//                 // add disaster variables
//                 for (int d = 0; d < _network._Nd; d++){
//                     for (int dd = 0; dd < _network._disaster_precision[d]; dd++){
//                         IloBoolVar var_dis(env);
//                         bvars[t][i][m].add(var_dis);

//                         index_disasters[t][i][m][d].push_back(bvars[t][i][m].getSize() - 1);
//                     }
//                 }

//                 // add capacity variables
//                 for (int c = 0; c < _network._capacity_precision; c++){
//                         IloBoolVar var_cap(env);
//                         bvars[t][i][m].add(var_cap);

//                         index_capacities[t][i][m].push_back(bvars[t][i][m].getSize() - 1);
//                 }
//             }
//         }
//     }

//     supplier_selection = IloBoolVarArray(env);
//     index_selection.resize(_N);
//     for (int i = 0; i< _N; i++){
//         index_selection[i].resize(_N);
//         fill(index_selection[i].begin(), index_selection[i].end(), -1);
//     }

//     for (int i = 0; i< _N; i++){
//         for (int m = 0; m < _M; m++){
//             if(_network._demand[i][m] != 0){
//                 for (int p = 0; p < _network._producers[m].size(); p ++){
//                     IloBoolVar var_select(env);
//                     supplier_selection.add(var_select);

//                     int idx = _network._producers[m][p];
//                     index_selection[idx][i] = supplier_selection.getSize() - 1;
//                 }
//             }
//         }
//     }
// }