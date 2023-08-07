#include "shelter_location.h"

// #define DEBUG_SHELTER
using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


ShelterLocation::ShelterLocation(){
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
    allow_bf = false;

    // global parameters, etc.
    timelimit = -1;
}

ShelterLocation::~ShelterLocation(){
    env.end();
}

void ShelterLocation::loadParameters(char graph_file[], int N, int T, int M, 
                                     vector<int> src, vector<int> q, char output_dir[]){
    
    if(N < 0){
        graph = new Graph(graph_file);
        _N = graph->_N;
    }
    else{
        // fib graph
        _N = N;
        graph = new Graph(N, 1, 2);
    }

    _T = T;
    _M = M;
    _src = src;
    _q = q;
    strcpy(_output_dir, output_dir);
}


void ShelterLocation::genFlowConstraints(){
    // There will be T * n_src complete flows
    bvars.resize(_T);
    flows.resize(_T);
    shelters.resize(_T);

    // One-hot encoding for source nodes
    // Save it to the instance if necessary
    vector<vector<IloBool>> sources(_src.size(), vector<IloBool> (_N, 0));

    for (int s = 0; s< _src.size(); s++){
        sources[s][_src[s]] = 1;
    } 

    // Preprocess
    for (int t = 0; t< _T; t++){
        flows[t].resize(_src.size());
        shelters[t].resize(_src.size());

        for (int s = 0; s< _src.size(); s++){  
            flows[t][s].resize(_N);
            shelters[t][s].resize(_N);
            bvars[t].push_back(IloBoolVarArray(env));

            for (int i = 0; i< _N; i++){
                flows[t][s][i].resize(_N);
                for(int j = 0; j< _N; j++){
                    if(graph->Adj[i][j] == 1){
                        IloBoolVar f_var(env);
                        bvars[t][s].add(f_var);
                        flows[t][s][i][j] = bvars[t][s].getSize() - 1;
                    }
                }
            }            
            for (int i = 0; i< _N; i++){
                IloBoolVar s_var(env);
                bvars[t][s].add(s_var);
                shelters[t][s][i] = bvars[t][s].getSize() - 1;
            }
        }
    }

    // Extract flow constraints
    // vector<vector<IloConstraintArray>> const_flow;
    // this is the main output
    const_flow.resize(_T);

    for (int t = 0; t < _T; t++){
        for (int s = 0; s < _src.size(); s++){
            IloConstraintArray const_flow_t_s(env);
            if (const_flow_t_s.getSize() != 0)
                cout << "Wrong!" << endl;

            // all nodes
            for (int i = 0; i< _N; i++){
                // consider node i
                // in & out
                IloNumExpr degree_abs(env);
                for (int j = 0; j< _N; j++){
                    if(graph->Adj[j][i] == 1){
                        // in 
                        int var_idx = flows[t][s][j][i];
                        degree_abs += bvars[t][s][var_idx];
                    }
                    if(graph->Adj[i][j] == 1){
                        // out 
                        int var_idx = flows[t][s][i][j];
                        degree_abs -= bvars[t][s][var_idx];
                    }
                }
                int var_idx = shelters[t][s][i]; 
                const_flow_t_s.add((degree_abs == (bvars[t][s][var_idx] - sources[s][i])));
            }

            const_flow[t].push_back(const_flow_t_s);
        }
    }

}


void ShelterLocation::genShelterConstraints(){
    shelter_assign = IloBoolVarArray(env, _N);

    // Extract the is_shelter constraints
    ////////////////////////////////////////
    // vector<vector<IloConstraintArray>> const_is_shelter;

    const_is_shelter.resize(_T);
    
    for (int t = 0; t < _T; t++){
        for (int s = 0; s < _src.size(); s++){
            // shelters[t][s] = N dimensional (0,0,0,1,0,0,0,...)
            IloConstraintArray const_number_t_s(env);
            
            for (int i = 0; i< _N; i++){
                int var_idx = shelters[t][s][i];
                const_number_t_s.add(((shelter_assign[i] - bvars[t][s][var_idx]) >= 0));
            }
            const_is_shelter[t].push_back(const_number_t_s);
        }
    }

    // Extract the number constraint
    ////////////////////////////////////////
    // IloConstraint const_max_shelter;
    IloNumExpr count_shelters(env);
    for(int i = 0; i < _N; i++){
        count_shelters = count_shelters + shelter_assign[i];
    }
    const_max_shelter = (count_shelters < _M + 1); 
}


void ShelterLocation::genMajorityConstraints(){
    bvars_maj = IloBoolVarArray(env, _T);

    IloNumExpr majority_sum(env);
    for (int t = 0; t < _T; t++){
        majority_sum = majority_sum + bvars_maj[t];
    }

    const_majority = (majority_sum > _T/2);

    const_flow_union.resize(_T);
    const_is_shelter_union.resize(_T);


    vector<vector<IloAnd>> const_flow_and;
    vector<vector<IloAnd>> const_is_shelter_and;

    const_flow_and.resize(_T);
    const_is_shelter_and.resize(_T);

    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            IloAnd temp_cf(env);
            temp_cf.add(const_flow[t][s]);
            const_flow_and[t].push_back(temp_cf);

            IloAnd temp_is(env);
            temp_is.add(const_is_shelter[t][s]);
            const_is_shelter_and[t].push_back(temp_is);
        }
    }

    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            IloConstraint temp_const;
            temp_const = ((const_flow_and[t][s]) || (bvars_maj[t] == 0));
            const_flow_union[t].push_back(temp_const);
            temp_const = ((const_is_shelter_and[t][s]) || (bvars_maj[t] == 0));
            const_is_shelter_union[t].push_back(temp_const);
        }
    }
}

void ShelterLocation::genNBFConstraints(){
    const_nbf.resize(_T);
    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){  
            IloConstraintArray const_nbf_t_s(env);
            for (int i = 0; i< _N; i++){
                for(int j = i; j< _N; j++){
                    if((graph->Adj[i][j] == 1) && (graph->Adj[j][i] == 1)){
                        // identify bf edges
                        const_nbf_t_s.add((bvars[t][s][flows[t][s][i][j]] + bvars[t][s][flows[t][s][j][i]] < 2));
                    }
                }
            }
            const_nbf[t].push_back(const_nbf_t_s);
        }
    }
}


void ShelterLocation::genXORConstraints(){
    const_xor.resize(_T);
    ivars_xor.resize(_T);
    bvars_xor.resize(_T);

    for( int t=0; t < _T; t++){
        for(int s=0; s<_src.size(); s++){
            int n_vars = bvars[t][s].getSize();
            int n_parity = _q[s];

            std::cout << "Number of variables in XOR is: " << n_vars << endl;

            vector<vector<bool>> coeffA = generate_Toeplitz_matrix(n_parity, n_vars);
            // coeffA: n_parity x (n_vars + 1)
            // be sure to generate a set of non-self-conflict XOR constraints
            while(!makeHashFuncSolvable(coeffA)){
                coeffA = generate_Toeplitz_matrix(n_parity, n_vars);
            }
            extractXorVarConst(coeffA, t, s);
        }
    }

}


void ShelterLocation::genAllConstraints(){
    // generate everything
    genFlowConstraints();
    cout << "Flow constraints generated!" << endl;
    genShelterConstraints();
    cout << "Shelter constraints generated!" << endl;
    genMajorityConstraints();
    cout << "Majority constraints generated!" << endl;
    if(!allow_bf){
        genNBFConstraints();
        cout << "NBF constraints generated!" << endl;
    }
    genXORConstraints(); 
    cout << "XOR constraints generated!" << endl;
}



void ShelterLocation::prepareModel(){

    model.add(shelter_assign);
    model.add(bvars_maj);

    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            model.add(bvars[t][s]);
            model.add(bvars_xor[t][s]);
            model.add(ivars_xor[t][s]);
        }
    }
    cout << "Variables added!" << endl;

    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            model.add(const_flow_union[t][s]);
            model.add(const_is_shelter_union[t][s]);
            model.add(const_xor[t][s]);
        }
    }

    model.add(const_max_shelter);
    model.add(const_majority);

    if(!allow_bf){
        for (int t = 0; t< _T; t++){
            for (int s = 0; s< _src.size(); s++){
                model.add(const_nbf[t][s]);
            }
        }
    }
    cout << "Constraints added!" << endl;
}


bool ShelterLocation::solveInstance() {
    cout << "\n================== Set Output Path ===================\n" << endl;
    #ifdef DEBUG_SHELTER
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
              << "Graph size: " << graph->Adj.size()  << endl;

    ostream_iterator<int> output_iterator(fs_params, ", ");

    fs_params << "\nsources: " << endl;
    std::copy(_src.begin(), _src.end(), output_iterator);    
    fs_params << "\nqlist: " << endl;
    std::copy(_q.begin(), _q.end(), output_iterator);
    
    cplex.exportModel("model.lp");
    prepareModel();
    // cplex.clearModel();  // clear existing model
    // cplex.extract(model);
    // cplex.setParam(IloCplex::Threads, 4);    // number of parallel threads (automatic by default)
    
    bool solved = cplex.solve();
    cplex.exportModel("model.lp");
    // env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << cplex.getCplexStatus() << endl; 

    if(solved){
        vector<int> shelter_idx;

        for(int i = 0; i<shelter_assign.getSize(); i++){
            env.out() << int(cplex.getValue(shelter_assign[i])) << " "; 

            if(cplex.getValue(shelter_assign[i]) == 1){
                shelter_idx.push_back(i);
            }
        }

        #ifdef DEBUG_SHELTER
        int sat_idx;
        for(int i = 0; i<_T; i++){
            if(cplex.getValue(bvars_maj[i]) == 1){
                sat_idx = i;
            }
            env.out() << "\nMajority t = " << cplex.getValue(bvars_maj[i]) << endl; 
        }

        for(int i = 0; i<bvars[sat_idx][0].getSize(); i++){
            env.out() << "flow vars = " << cplex.getValue(bvars[sat_idx][0][i]) << endl; 
        }
        env.out() << "t = " << sat_idx << " is satisfied" << endl; 
        #endif

        env.out() << "\nShelter at:" << endl;
        for (int i = 0; i < shelter_idx.size(); i++){
            env.out() << shelter_idx[i] << " ";
        }
    }
    return true;
}



bool ShelterLocation::makeHashFuncSolvable(vector<vector<bool>> &coeffA) {
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

void ShelterLocation::extractXorVarConst(vector<vector<bool>> coeffA, int t, int s) {
    // Encode XOR constraints, and extract corresponding variables and constraints
    IloConstraintArray xor_const_array(env);
    IloIntVarArray xor_int_vars_array(env);
    IloBoolVarArray xor_bool_vars_array(env);
    int nbvars = bvars[t][s].getSize();

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
                    xor_const_array.add((zeta_sum_to_f == bvars[t][s][i]));
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
    const_xor[t].push_back(xor_const_array);
    ivars_xor[t].push_back(xor_int_vars_array);
    bvars_xor[t].push_back(xor_bool_vars_array);
}