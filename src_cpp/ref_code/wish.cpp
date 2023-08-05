#include "wish.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

//////////////////////////////////////////////////////////

unsigned long get_seed(void) {
  struct timeval tv;
  struct timezone tzp;
  gettimeofday(&tv,&tzp);
  return (( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec;
}

void parseParityArgs(WishInstance & ins, int & argc, char **argv)
{
  // this method eats up all arguments that are relevant for the
  // parity constraint, and returns the rest in argc, argv

  char residualArgv[argc][64];
  strcpy(residualArgv[0], argv[0]);
  int residualArgc = 1;

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-number") ) {
      argIndex++;
      ins.parity_number = (IloInt)atol( argv[argIndex] );
    }
    else if ( !strcmp(argv[argIndex], "-yannakis") ) {
      ins.yannakis = true;
    }
    else {
      // save this option to be returned back
      strcpy(residualArgv[residualArgc++], argv[argIndex]);
    }
  }

  argc = residualArgc;
  for (int i=1; i<argc; ++i) {
    // free(argv[i]);
    argv[i] = new char[strlen(residualArgv[i])+1];
    strcpy(argv[i], residualArgv[i]);
  }
}

void parseArgs(WishInstance & ins, int argc, char **argv) 
{
  // one argument must be the instance filename
  if (argc <= 1) {
    cerr << "ERROR: instance name must be specified" << endl;
    exit(1);
  }

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-timelimit") ) {
      argIndex++;
      ins.timelimit = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-seed") ) {
      argIndex++;
      ins.seed =  atol( argv[argIndex] );
      ins.use_given_seed = true;
    }
    else if ( !strcmp(argv[argIndex], "-verbosity") ) {
      argIndex++;
    }
    else if ( !strcmp(argv[argIndex], "-h") || !strcmp(argv[argIndex], "-help") ) {
      cout << endl
           << "USAGE: iloglue_uai [options] instance.uai" << endl
           << endl
           << "   -timelimit          Timelimit in seconds (default None)" << endl
           << "   -seed               Random seed" << endl
           << endl;
      // print parity constraint options usage
      //printParityUsage(cout);
      exit(0);
    }
    else if (argv[argIndex][0] != '-') {
      // must be the instance name
      strcpy(ins.instance_name, argv[argIndex]);
    }
    else {
      cerr << "ERROR: Unexpected option: " << argv[argIndex] << endl
           << " See usage (iloglue_uai -h)" << endl;
      exit(1);
    }
  }
}

double median(vector<double> &v)
{
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

void print_matrix (vector <vector <bool> > A)
{
for (unsigned int i =0;i<A.size();i++)
  {
    for (unsigned int j =0;j<A[i].size();j++)          // last column is for coefficients b
      cout << A[i][j] << ",";
    cout << endl;
  }
}

vector <vector <bool> > generate_matrix(int m, int n) {
  vector <vector <bool> > A;
  A.resize(m);
  for (int i =0;i<m;i++) {
  A[i].resize(n+1);
  for (int j =0;j<n+1;j++)          // last column is for coefficients b
    //if (rnd_uniform()<0.5)
    if (rand()%2==0)
      A[i][j] = true;
    else
      A[i][j]=false;
  }
  
  return A;  
}

int sparsify(vector <vector <bool> > & A) {
  vector< bitset<1000> > bv;

  bv.resize(A.size());
  size_t m = A.size();
  size_t n = A[0].size();
  int saved_bits = 0;
  size_t initialnbr=0;

  for (size_t i = 0;i<m;i++)
    {
    bitset<1000> row;
    for (size_t s = 0;s<n;s++)
      if (A[i][s]) 
        {
        row.set(s,1);
        initialnbr++;
        }
    bv[i]=row;  
    }

  cout << "Initial # of bits: " <<   initialnbr << endl;


  if (m<100)
  {
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      for (size_t z = 0;z<m;z++)
        for (size_t g = 0;g<m;g++)
          if (i!=l && i!=z && z!=l && i!=g && z!=g && l!=g)
          {
          int cursize = bv[i].count();
          int newsize =  (bv[i]^bv[l]^bv[z]^bv[g]).count();
          if (newsize<cursize)
            {
            saved_bits = saved_bits + cursize - newsize;
            bv[i]^=bv[l]^bv[z]^bv[g];
            for (int q=0;q<n;q++)      // sum the two rows
              A[i][q] = A[i][q] xor A[l][q] xor A[z][q] xor A[g][q];
            }
          }

  }
  if (m<500)
  {
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      for (size_t z = 0;z<m;z++)
      if (i!=l && i!=z && z!=l)
      {
      int cursize = bv[i].count();
      int newsize =  (bv[i]^bv[l]^bv[z]).count();
      if (newsize<cursize)
        {
        saved_bits = saved_bits + cursize - newsize;
        bv[i]^=bv[l]^bv[z];
        for (int q=0;q<n;q++)      // sum the two rows
          A[i][q] = A[i][q] xor A[l][q] xor A[z][q];
        }
      
      }
  }
  if (m<10000)
  {    
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      if (i!=l)
      {
      int cursize = bv[i].count();
      int newsize =  (bv[i]^bv[l]).count();
      if (newsize<cursize)
        {
        saved_bits = saved_bits + cursize - newsize;
        bv[i]^=bv[l];
      //  print_matrix(A);
        
        for (int q=0;q<n;q++)      // sum the two rows
          A[i][q] = A[i][q] xor A[l][q];
        }
      
      }
  }

  cout << "final # of bits: " <<   (int) initialnbr- saved_bits<< endl;    
  return saved_bits;  
}

vector <vector <bool> > generate_Toeplitz_matrix(int m, int n) {
  vector <vector <bool> > A;
  if (m==0)
    return A;
    
  A.resize(m);
  int i;
  for (i =0;i<m;i++) {
    A[i].resize(n+1);
  }
  
  // first column
  for (i =0;i<m;i++) {
    if (rand()%2==0)
      A[i][0] = true;
    else
      A[i][0]=false;
    for (int j =1;j<m-i;j++)
      if (j<n)
        A[i+j][j] = A[i][0];
  }
  
  // last column
  for (i =0;i<m;i++) {
    if (rand()%2==0)
      A[i][n] = true;
    else
      A[i][n]=false;

  }
  
  // first row
  for (int j =1;j<n;j++) {
    if (rand()%2==0)
      A[0][j] = true;
    else
      A[0][j]=false;
      
    for (i =1;i<m;i++)
      if (j+i<n)
        A[i][j+i] = A[0][j];
  }

  return A;  
}

//////////////////////////////////////////////////////////

WishInstance::WishInstance(){
  // Initialize
  timer = new IloTimer(env);
  model = new IloModel(env);
  cost = new IloArray<IloNumArray>(env);
  cplex = new IloCplex(env);
  
  sparsify_coeff = true;
  parity_number = 0;

  yannakis =true;

  // global parameters, etc.
  seed = 10;
  use_given_seed = true;
  timelimit = -1;
}

WishInstance::~WishInstance() {}

void WishInstance::loadInstance() {
  // open the instance file
  ifstream file(instance_name);
  if (!file) {
    cerr << "Could not open file " << instance_name << endl;
    exit(EXIT_FAILURE);
  }

  // stefano mod, read uai file
  // reads uai file to parse domain sizes; creates variables along the way
  cerr << "Creating variables"<< endl;
  file >> pbname;
  file >> nbvar;
  vars = new IloIntVarArray(env, nbvar, 0, 100);

  nbval = 0;
  int tmp;
  for (int i=0; i<nbvar; i++) {
    file >> tmp;
    if (tmp>nbval)
      nbval = tmp;
    (*vars)[i].setBounds(0, tmp-1);
    char *name = new char[16];
    sprintf(name, "x%d", i);
    (*vars)[i].setName(name);
  }
  model->add((*vars));
  file >> nbconstr;
  cerr << "Var:"<< nbvar <<" max dom size:" <<nbval<<" constraints:"<<nbconstr << endl;

  // use a native CP Optimizer representation of the .uai file
  // read in variable scopes of CPT tables
  int arity;
  for (int i=0; i<nbconstr; i++) {
    file >> arity;
    std::vector< int> scope;
    int id;
    for (int j=0; j<arity; j++) {
      file >> id;
      scope.push_back(id);
    }
    scopes.push_back(scope);
  }

  // read in values of CPT tables
  int TableSize;
  for (int i=0; i<nbconstr; i++) {
    file >> TableSize;
    double prod;
    IloNum entry;
    IloNumArray table(env);
    for (int j=0; j<TableSize; j++) {
      file >> prod;        
      // entry = log10(prod);    // log likelihood objective function
      entry = prod;    // normal objective function
      table.add(entry);
    }
    cost->add(table);
  }
  cout << "done reading CPTs"<< endl;
  
  // get random seed
  if (!use_given_seed) seed = get_seed();
  // srand(seed);
}

void WishInstance::readInstance(char *file_path) {
  strcpy(instance_name, file_path);
  loadInstance();
}

void WishInstance::getObjectExpression() {
  // define cost expression
  objexpr = new IloNumExpr(env);

  // my modification
  for (int l = 0; l < nbconstr; l++) {
    int scope_size = scopes[l].size(); // variables in the factor scopes[l][0:num_bvar-1]
    int nbvar_mu = (int)pow(2, scope_size);

    IloBoolVarArray Mu(env, nbvar_mu);

    for (int k = 0; k < nbvar_mu; k++){
      char *name = new char[32];
      sprintf(name, "mu_%d_%d", l, k);
      Mu[k].setName(name);
    }
    model->add(Mu);

    // add constraints
    for (int k = 0; k < scope_size; k++){
      int intvl = (int)pow(2, k);
      // expression
      IloNumExpr constr_expr(env);

      int v_idx = scopes[l][scope_size - k - 1];

      // scopes[l][0], scopes[l][1], scopes[l][2]
      for (int m = 0; m * intvl < nbvar_mu; m+=2){
        for (int n = 0; n < intvl; n++){
          constr_expr += Mu[m * intvl + n];
        }
      }

      model->add((constr_expr == (1 - (*vars)[v_idx])));
    }

    IloNumExpr sum_expr(env);
    for (int k = 0; k < nbvar_mu; k++){
      sum_expr += Mu[k];
    }
    model->add((sum_expr == 1));

    for (int k = 0; k < nbvar_mu; k++){
      if (isfinite((*cost)[l][k]))
        *objexpr += (*cost)[l][k] * Mu[k];
      else {
        if (isinf((*cost)[l][k])) {
          model->add(Mu[k]==0);
        }
        else {
          cout << "Cannot generate ILP"<< endl;
          exit(-1);
        }
      }  
    }
  }
}

void WishInstance::addOptimizationGoal() {
  model->add(IloMaximize(env, *objexpr));
}

void WishInstance::sampleHashCoeffA() {
  // generate matrix of coefficients A x = b. b is the last column
  coeffA.clear();
  coeffA = generate_Toeplitz_matrix(parity_number, nbvar);
}

bool WishInstance::getFeasibleSolution() {
  feasiblesol.clear();
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
  feasiblesol.resize(n);
  y.resize(n);
  
  // initialize b to the last column of A
  b.resize(m);
  for (size_t i =0;i<m;i++)
    b[i] = coeffA[i][n];
  
  for (size_t i =0;i<n;i++)
    y[i] = rand()%2;
    
  // sum all the dependent variables that are already fixed
  for (size_t i =0;i<n;i++) {
    feasiblesol[i] = y[i];
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
      feasiblesol[indep_columns[i]] = 1;
      for (size_t j =0;j<m;j++)
        b[j] = b[j] xor coeffA[j][indep_columns[i]];
    }
    else {
      y[indep_columns[i]] = 0;
      feasiblesol[indep_columns[i]] = 0;
    }
  }

  // print_matrix(coeffA);
  cout << "row echelon form:" << endl;
  print_matrix(coeffA);
  
  // sparsify
  if (sparsify_coeff && (!coeffA.empty())) {
    cout << "Bits saved: ";
    for (int i=0; i<2; i++)
      cout << sparsify(coeffA) << " ";
    cout << endl;  
  }
  return true;
}

void WishInstance::extractXorConstraints() {
  // XOR
  xor_constr_array = new IloConstraintArray(env);
  xor_int_vars_array = new IloIntVarArray(env);
  xor_bool_vars_array = new IloBoolVarArray(env);

  std::vector < std::set <size_t> > varAppearancesInXors;
  if (!coeffA.empty())
    varAppearancesInXors.resize(coeffA[0].size());            // dummy parity var

  IloArray<IloArray<IloIntVarArray> > zeta_vars(env);
  IloArray <IloIntVarArray> alpha_vars(env);

  IloBoolVar dummy_parity (env, 0, 1, "dummy");
  IloConstraint cons = (dummy_parity==1);

  model->add(dummy_parity);
  model->add(cons);

  xor_bool_vars_array->add(dummy_parity);
  xor_constr_array->add(cons);
  


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
        
        cons = (alpha_sum_to_one==1);
        model->add(alphas);
        model->add(cons);              // (15)

        xor_int_vars_array->add(alphas);
        xor_constr_array->add(cons);

        
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
            
            cons = (zeta_i_j_k<=alpha_vars[*it][k/2]);
            model->add(cons);
            xor_constr_array->add(cons);
          }
          model->add(zet);
          xor_int_vars_array->add(zet);
            
          if (i<nbvar){
            cons = (zeta_sum_to_f==(*(vars))[i]);
            model->add(cons); 
            xor_constr_array->add(cons);
          }
          else if (i==nbvar) {       
            cons = (zeta_sum_to_f==dummy_parity);                // dummy
            model->add(cons);
            xor_constr_array->add(cons);
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
          cons = (zeta_sum_to_alpha_k==IloInt(k)*alpha_vars[j][k/2]);
          model->add(cons);
          xor_constr_array->add(cons);
        }
      }
    }
  }
}

bool WishInstance::solveInstance() {
  cplex->clearModel();  // clear existing model
  cplex->extract(*model);
  if (timelimit > 0)
    cplex->setParam(IloCplex::TiLim, timelimit);

  cplex->setParam(IloCplex::Threads, 1);    // number of parallel threads
  
  if (!coeffA.empty()) {
    IloNumArray feasibleinit(env);
    //double [] feasibleinit;
    IloNumVarArray startVar(env);
  
    for (size_t l= 0; l<nbvar;l++) {
      startVar.add((*(vars))[l]);
      feasibleinit.add(feasiblesol[l]);
    }
    cplex->addMIPStart(startVar, feasibleinit);
    // feasibleinit.end();  // https://or.stackexchange.com/questions/3530/no-solution-found-from-n-mip-starts
    // startVar.end();
  }

  if ( !cplex->solve() ) {
    env.out() << "Failed to optimize LP." << endl;
    return false;
  }

  IloNumArray vals(env);
  env.out() << "Solution status = " << cplex->getStatus() << endl;
  env.out() << "Solution value = " << cplex->getObjValue() << endl;

  cplex->getValues(vals, *vars);
  // env.out() << "Values = " << vals << endl;
  return true;
}

void WishInstance::removeXorConstraints() {
  // delete variables and constraints
  //
  // Then free:
  // IloConstraintArray *xor_constr_array;
  // IloIntVarArray *xor_int_vars_array;
  // IloBoolVarArray *xor_bool_vars_array;
  // 
  // Clear vectors:
  // std::vector < std::vector <bool> > coeffA;
  // std::vector <bool> feasiblesol;
  //
	model->remove(*xor_bool_vars_array);
  model->remove(*xor_int_vars_array);
	model->remove(*xor_constr_array);

  xor_bool_vars_array->clear();
	xor_int_vars_array->clear();
	xor_constr_array->clear();

  coeffA.clear();
  feasiblesol.clear();

}

//////////////////////////////////////////////////////////

// Export model
// test_ins.cplex->exportModel("model.lp");