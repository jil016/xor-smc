#include "smc.h"

SMC::SMC(){
    // Initialize
    env.end();

    env = IloEnv();
    model = IloModel(env);
    cplex = IloCplex(model);

    sparsify_coeff = true;
    yannakis =true;

    // global parameters, etc.
    timelimit = -1;
}

SMC::~SMC(){
    env.end();
}

void SMC::loadParameters(const string& cnf_file, const string& uai_file, const string& mar_file,
                    int threshold, int T, char output_dir[]){
    _cnf_file = cnf_file;
    _uai_file = uai_file;
    _mar_file = mar_file;
    _threshold = threshold;   // number of xor constraints
    _T = T;
    _prec_n_digits = 3;
    _prec_n_int = 3;

    strcpy(_output_dir, output_dir);
    readFiles();
    initializeVariables();
}


void SMC::readFiles() {
    ifstream fp;
    fp.open(_cnf_file);

    string line;
    while (getline(fp, line)) {
        if (line.empty()) continue;

        if (line[0] == 'c') {
            continue;
        } else if (line[0] == 'p') {
            std::istringstream iss(line);
            char tmp[1024];
            iss >> tmp;
            iss >> tmp;
            iss >> _N_vars_cnf;
        } else if (isdigit(line[0]) || line[0] == '-') {
            std::vector<int> lits;
            std::istringstream iss(line);
            int l;
            while (iss >> l) {
                if(l != 0) lits.push_back(l);
                else break;
            }
            _clauses.push_back(lits);
        }
    }
    fp.close();


    // Load uai model
    fp.open(_uai_file);

    char pbname[1024];
    fp >> pbname;   // MARKOV or BAYES
    fp >> _N_vars_uai;

    int tmp;
    // read scopes
    for (int j = 0; j < _N_vars_uai; j++){
        fp >> tmp;  // scopes are all 2 by default; otherwise complicated
        if(tmp != 2) {
            cout << "Error: UAI file has non-boolean variable!" << endl;
            exit(0);
        }
    }

    fp >> _N_factors;
    int arity;
    for (int i = 0; i < _N_factors; i++) {
        fp >> arity;
        vector<int> tmp_factor;
        int id;
        for (int j=0; j<arity; j++) {
            fp >> id;
            tmp_factor.push_back(id);
        }
        _factors.push_back(tmp_factor);
    }

    // read in values of CPT tables
    int table_size;
    for (int i=0; i<_N_factors; i++) {
        fp >> table_size;
        double prob;
        vector<double> tmp_table;
        for (int j=0; j<table_size; j++) {
            fp >> prob;
            tmp_table.push_back(prob);
        }
        _tables.push_back(tmp_table);
    }
    cout << "done reading UAI"<< endl;
    fp.close();

    // update scope_evi and scope_mar
    fp.open(_mar_file);
    vector<int> indi_mar(_N_vars_uai, 0);
    int n_mar;
    fp >> n_mar;
    for(int i = 0; i < n_mar; i++){
        int v;
        fp >> v;
        indi_mar[v] = 1;
    }
    for(int i = 0; i < _N_vars_uai; i++){
        if(indi_mar[i] == 0) _scope_evi.push_back(i);
        else _scope_mar.push_back(i);
    }
    fp.close();
}


void SMC::initializeVariables() {
    _var_cnf = IloBoolVarArray(env, _N_vars_cnf);    // cnf
    _const_cnf = IloConstraintArray(env);

    for (int t = 0; t < _T; t++){
        _var_uai.emplace_back(env, _N_vars_uai);  // uai variables
    }

    // Add variables for probability calculation
    _var_prob_dis.resize(_T);
    _var_prob_add.resize(_T);
    _expr_prob.resize(_T);
    for (int t = 0; t < _T; t++){
        _const_prob_dis.emplace_back(env);
    }
    _const_prob_add.resize(_T);  // constraints on added variables, T, N_factors, vars

    // add majority variables
    _var_maj = IloBoolVarArray(env, _T);

    /// Begin set variable names
    for (int i = 0; i < _N_vars_cnf; i++) {
        char *name = new char[32];
        snprintf(name, 32, "x_%d", (int) i);
        _var_cnf[i].setName(name);
    }

    for (int t = 0; t < _T; t++){
        for (int i = 0; i < _N_vars_uai; i++) {
            char *name = new char[32];
            snprintf(name, 32, "y_t%d_%d", (int) t, (int) i);
            _var_uai[t][i].setName(name);
        }
    }

    for(int t = 0; t < _T; t++){
        char *name = new char[32];
        snprintf(name, 32, "mj_t%d", (int) t);
        _var_maj[t].setName(name);
    }
}

void SMC::genCnfConstraints(){
    for(const auto& clause:  _clauses){
        IloNumExpr expr_clause(env);
        for (auto lit: clause){
            int var = abs(lit) - 1;
            if(lit > 0) expr_clause += _var_cnf[var];
            else expr_clause += (1 - _var_cnf[var]);
        }
        _const_cnf.add(expr_clause >= 1);
    }
}

void SMC::genUaiConstraints(){
    for (int t = 0; t < _T; t++){
        // read through UAI
        for(int i = 0; i < _N_factors; i++){
            IloNumExpr expr_i(env);
            IloBoolVarArray var_add_i(env);
            IloConstraintArray const_i(env);

            int f_size = _factors[i].size();   // number of variables
            int tab_size = _tables[i].size();  // number of terms in this factor

            // then we need tab_size == 2^f_size number of variables
            // sort by low to high binary
            for (int j = 0; j < tab_size; j++){
                char *name = new char[32];
                snprintf(name, 32, "mu_t%d_f%d_%d", (int) t, (int) i, (int) j);

                IloBoolVar mu_i_j (env, 0, 1, name);
                var_add_i.add(mu_i_j);
                expr_i += _tables[i][j] * mu_i_j;
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
                const_i.add(add_pos_expr == _var_uai[t][_factors[i][k]]);
                const_i.add(add_neg_expr == 1 - _var_uai[t][_factors[i][k]]);
            }
            _var_prob_add[t].push_back(var_add_i);
            _const_prob_add[t].push_back(const_i);
        }

        _var_prob_dis[t].resize(_N_factors);
        for(int i = 0; i < _N_factors; i++){
            _var_prob_dis[t][i] = IloBoolVarArray(env, _prec_n_digits);
            IloNumExpr expr_i(env);

            for(int j = 0; j < _prec_n_digits; j++){
                expr_i += pow(2, j) * _var_prob_dis[t][i][j];
            }

            _const_prob_dis[t].add(expr_i + 1 <= ((int) pow(2, _prec_n_digits - _prec_n_int)) * _expr_prob[t][i]);
        }
        // uai and cnf share evidence variables
        _const_uai_evi.emplace_back(env);
        for(const auto& v : _scope_evi){
            _const_uai_evi[t].add(_var_cnf[v] == _var_uai[t][v]);
        }
    }
}

void SMC::genXORConstraints(){
    // extract all variables that appear in XOR
    for( int t=0; t < _T; t++){
        _var_all_inxor.emplace_back(env);
        for(const auto& v : _scope_mar){
            _var_all_inxor[t].add(_var_uai[t][v]);
        }
        for (int i = 0; i < _N_factors; i++){
            _var_all_inxor[t].add(_var_prob_dis[t][i]);
        }
    }

    for( int t=0; t < _T; t++) {
        int n_vars = _var_all_inxor[t].getSize();
        int n_parity = 0;
        if(_threshold >= 0)
            n_parity = _threshold;

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
}

void SMC::genMajorityConstraints() {
    IloNumExpr majority_sum(env);
    for (int t = 0; t < _T; t++){
        majority_sum = majority_sum + _var_maj[t];
    }
    _const_majority = (majority_sum > _T/2);

    vector<IloAnd> const_and;
    for (int t = 0; t< _T; t++){
        // load all constraints in t-th repetition
        IloAnd temp_const_and(env);
        temp_const_and.add(_const_prob_dis[t]);   // constraints on discretization variables
        for(const auto & carray:_const_prob_add[t]){   // constraints on added variables, T, N_factors, vars
            temp_const_and.add(carray);
        }
        temp_const_and.add(_const_uai_evi[t]);

        const_and.push_back(temp_const_and);
    }

    for (int t = 0; t< _T; t++){
        IloConstraint temp_const;
        temp_const = ((const_and[t]) || (_var_maj[t] == 0));
        _const_all_union.push_back(temp_const);
    }
}

void SMC::genAllConstraints(){
    // generate everything
    genCnfConstraints();
    cout << "CNF constraints generated!" << endl;
    genUaiConstraints();
    cout << "UAI constraints generated!" << endl;
    genXORConstraints();
    cout << "XOR constraints generated!" << endl;
    genMajorityConstraints();
    cout << "Majority constraints generated!" << endl;
}


void SMC::prepareModel(){
    model.add(_var_cnf);
    model.add(_const_cnf);
    cout << "CNF variables and constraints added!" << endl;

    for(int t = 0; t < _T; t++){
        model.add(_var_uai[t]);

        for (int i = 0; i<_N_factors; i++){
            model.add(_var_prob_dis[t][i]);
            model.add(_var_prob_add[t][i]);
        }
    }
    cout << "UAI variables added!" << endl;

    // XOR Constraints
    for(int t = 0; t < _T; t++){
        model.add(ivars_xor[t]);
        model.add(bvars_xor[t]);
        model.add(const_xor[t]);
    }
    cout << "XOR variables and constraints added!" << endl;

    // Majority and all constraints
    model.add(_var_maj);
    model.add(_const_majority);
    for(int t = 0; t < _T; t++){
        model.add(_const_all_union[t]);
    }
    cout << "Majority constraints (union with all constraints) added!" << endl;
}

bool SMC::solveInstance() {
    cout << "\n================== Set Output Path ===================\n" << endl;

    filesystem::path log_folder(_output_dir);

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

    fs_params << "T: " << _T << "\n"
              << "demand: " << _threshold << "\n"
              << "CNF file: " << _cnf_file << "\n"
              << "UAI file: " << _uai_file << "\n"
              << "MAR file: " << _mar_file << "\n"
              << endl;

    fs_params.close();

    auto start = std::chrono::high_resolution_clock::now();

    genAllConstraints();
    prepareModel();

    // Cplex parameters
    cplex.setParam(IloCplex::Param::WorkMem, 2048);
    cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 2048);
    cplex.setParam(IloCplex::Threads, 4);    // number of parallel threads (automatic by default)

    bool solved = cplex.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cplex.exportModel("model_solved.lp");

    // env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << cplex.getCplexStatus() << endl;

    if(solved){
        // read the trade selection
        vector<int> parse_result;
        parse_result.resize(_N_vars_cnf);

        for (int i = 0; i < _N_vars_cnf; i++){
            parse_result[i] = int(cplex.getValue(_var_cnf[i]));
        }

        for (int i = 0; i < _N_vars_cnf; i++){
            if(parse_result[i] > 0){
                env.out() << (i+1) << " ";
            }
            else{
                env.out() << -(i+1) << " ";
            }
        }
        env.out() << "0\n"
                  << duration.count() * 0.001 << " ms";
    }
    else{
        env.out() << "NO SOLUTION " << endl;
    }

    return true;
}

bool SMC::makeHashFuncSolvable(vector<vector<bool>> &coeffA) {
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
    // sparsify
    if (sparsify_coeff && (!coeffA.empty())) {
        cout << "Bits saved: ";
        for (int i=0; i<2; i++)
            cout << sparsify(coeffA) << " ";
        cout << endl;
    }
    return true;
}

void SMC::extractXorVarConst(vector<vector<bool>> coeffA, int t) {
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