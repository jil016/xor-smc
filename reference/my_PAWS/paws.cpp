#include "wish.h"
#include "paws.h"

// #define DEBUG_PAWS

double total_find_P_time;
double max_find_P_time;
int n_find_P;

SInstance::SInstance() : WishInstance() {
  nbvar_x = 0;
  nbvar_y = 0;
}


void SInstance::addYVariables(int l, int b) {
  //
  // current w(x) doesn't have y_ik or XOR constraint
  // add binary variables: y_i_k, 1<=i<=l-1, 1<=k<=b
  // 
  nbvar_x = nbvar;
  nbvar_y = (l - 1) * b;
  vars_y = new IloIntVarArray(env, nbvar_y, 0, 100);

  for (int i = 0; i < l - 1; i++){
    // constraints
    for (int k = 0; k < b; k++){
      char *name = new char[16];
      sprintf(name, "y_%d_%d", i, k);
      (*vars_y)[i * b + k].setName(name);
      (*vars_y)[i * b + k].setBounds(0, 1);
    }

  }

  vars->add(*vars_y);
  model->add(*vars);
  nbvar += nbvar_y;
}

void SInstance::addSConstraints(IloNum M, int l, int b) {
  // add constraints for y

  IloNum r = pow(2.0, b) / (pow(2.0, b) - 1);

  IloNumExpr new_test(env);

  for (int i = 0; i < l - 1; i++){
    // constraints
    IloNumExpr y_sum_expr(env);

    char *name = new char[16];
    sprintf(name, "g%d", i);
    IloBoolVar g_i(env, 0, 1, name);

    for (int k = 0; k < b; k++){
      y_sum_expr += (*vars_y)[i * b + k];
    }

    // (g_i == 1) <=> (*objexpr) <= (M / pow(r, i))
    // (g_i == 0) <=> (*objexpr) > (M / pow(r, i))
    // (g_i == 1) => y_sum_expr > 1

    model->add(g_i);
    model->add((y_sum_expr >= g_i));

    IloInt bigN = 10 * ceil(M);

    model->add(((*objexpr) >= (ceil(M / pow(r, i + 1)) - bigN * g_i)));
    model->add(((*objexpr) <= (floor(M / pow(r, i + 1)) + bigN * (1 - g_i)))); 
  }

  model->add(((*objexpr) >= ceil(M / pow(r, l)))); // w(x) > M/r^l
}

void SInstance::addNaiveOptimizationGoal() {
  IloNumExpr naive_goal(env);
  naive_goal += (*vars)[0];
  model->add(IloMaximize(env, naive_goal));
}

int SInstance::findPSolutions(int P) {
  IloConstraintArray sol_constraints(env);
  solutions.clear();

  int count = 0;
  while (count < P) {
    if (!solveInstance()) {
      cout << "Cannot solve instance " << count << endl;
      break;
    }

    #ifdef DEBUG_PAWS
    cplex->exportModel("fds_Model1.lp");
    #endif

    count++;
    IloNumArray vals(env);
    cplex->getValues(vals, *vars);
    solutions.push_back(vals);

    cout << "Solution " << count 
         << " : " << vals << endl;

    IloNumExpr constr_expr(env);
    for (int j = 0; j < nbvar; j++) {
      if (vals[j] == 1)
        constr_expr += (1 - (*vars)[j]);
      else
        constr_expr += (*vars)[j];
    }
    IloConstraint cons = (constr_expr > 0);

    sol_constraints.add(cons);
    model->add(cons);
    #ifdef DEBUG_PAWS
    cplex->exportModel("fds_Model2.lp");
    #endif
  }

  model->remove(sol_constraints);  // remove constraints from previous solutions
  #ifdef DEBUG_PAWS
  cplex->exportModel("fds_Model3.lp");
  #endif
  return count;
}

int computeK(SInstance &S_base, int n, double delta, int P) {
  int T = 24 * ceil(log(n / delta));
  int max_count = T / 2;
  int k = -1;
  int count;

  do {
    count = 0;
    k += 1;
    for (int t = 0; t < T; t++) {
      cout << "=====================================\n"
           << "Iteration k: " << k << " t: "<< t << endl;

      cout << "Set parity number to " << k << endl;
      S_base.parity_number = k;
      
      S_base.sampleHashCoeffA();

      int n_solution = 0;

      // Try to find P elements in S 
      cout << "Try to find " << P << " solutions" << endl;
      const clock_t begin_time = clock();

      if(S_base.getFeasibleSolution()){
        cout << "Sample A matrix " << k << "x" << S_base.nbvar  << endl;
        print_matrix(S_base.coeffA);

        S_base.extractXorConstraints();

        n_solution = S_base.findPSolutions(P);

        #ifdef DEBUG_PAWS
        S_base.cplex->exportModel("cpk_Model1.lp");
        #endif

        // remove coefficients A and XOR constraints
        S_base.removeXorConstraints();  

        #ifdef DEBUG_PAWS
        S_base.cplex->exportModel("cpk_Model2.lp");
        #endif
      }
      else {
        cout << "No feasible solution" << endl;
      }

      if(n_solution < P) {
          count++;
          cout << "Didn't find >=P solutions" << endl;
        }
      else
        cout << "Found >=P solutions" << endl;
      
      double time_taken = double( clock () - begin_time ) /  CLOCKS_PER_SEC;
      cout << "Time taken: " << time_taken
            << " sec" << endl;
      
      total_find_P_time += time_taken;
      if(time_taken > max_find_P_time) max_find_P_time = time_taken;
      n_find_P++;
      
      if (count >= max_count)  // early break
        break;
    }
  } while((k < n) && (count < max_count));

  return k;
}


IloNumArray paws(char *file_path, int l, int b, double delta, int P, 
                 int alpha, unsigned long seed, filesystem::path cplex_log_path) {
  // solve　max(w(x))
  WishInstance w_func;

  // set parameters
  std::ofstream cplex_log_fs(cplex_log_path);

  w_func.readInstance(file_path);
  w_func.seed = seed;
  w_func.getObjectExpression(); //
  w_func.addOptimizationGoal(); 

  // solve without Ax=b constraint

  w_func.solveInstance();
  IloNum M = w_func.cplex->getObjValue();
  
  cout << "max{w(x)} is :" << M << endl;

  cout << "\n=====================================\n" << endl;

  // create a paws instance
  try {

    cout << "Preparing S" << endl;

    SInstance S_base;
    S_base.seed = seed;

    // S_base.cplex->setOut(S_base.env.getNullStream()); // mute CPLEX    
    S_base.cplex->setOut(cplex_log_fs);  // redirect to log file

    cout << "Reading w(x)" << endl;

    S_base.readInstance(file_path); // uai_file
    S_base.getObjectExpression(); // get w(x)
    
    cout << "Adding Y variables" << endl;

    S_base.addYVariables(l, b); // add y_l^k

    cout << "nbvar_x: " << S_base.nbvar_x << " nbvar_y: " << S_base.nbvar_y << endl;

    cout << "Adding constraints for S" << endl;

    S_base.addSConstraints(M, l, b); // add constraints for S set

    cout << "Adding naive obj function for S" << endl;

    S_base.addNaiveOptimizationGoal(); // add a naive objective function

    cout << "computeK(S, " 
    << S_base.nbvar << ", "
    << delta << ", "
    << P << ")"
    << endl;

    #ifdef DEBUG_PAWS
    S_base.solveInstance();
    S_base.cplex->exportModel("finalModel0.lp");
    #endif

    int k = computeK(S_base, S_base.nbvar, delta, P) + alpha;
    cout << k << " = computeK(S)" << endl;

    #ifdef DEBUG_PAWS
    S_base.cplex->exportModel("finalModel1.lp");
    #endif

    S_base.parity_number = k;

    S_base.sampleHashCoeffA();

    int count = 0;
    if (S_base.getFeasibleSolution()){
      S_base.extractXorConstraints();
      // S_base is ready, try to find P solutions 
      count = S_base.findPSolutions(P);

      #ifdef DEBUG_PAWS
      S_base.cplex->exportModel("finalModel2.lp");
      #endif
    }

    cout << "Found " << count << " solutions" << endl;
    
    if(count == P){
      printf("No result! Count == P\n");
      return NULL;
    }
    
    if(count == 0) {
      printf("No result! Count == 0\n");
      return NULL;
    }

    int random_p = rand() % P;
    if(random_p <= count) {
      cout << "Sample found:\n" << S_base.solutions[random_p] << endl;
      return S_base.solutions[random_p];
    }
    else {
      cout << "Dropped at random selection. Force output:\n";
      for(int m = 0; m < count; m++){
        cout << S_base.solutions[m] << endl;
      }
      return NULL;
    }
  }
  catch (IloException& ex) {
      cout << "Error: " << ex << endl;
  }
  
  return NULL;
}

IloNumArray pawsWrapper(char *file_path, int n, double epsilon, int b, 
                  double delta, int P, int alpha, unsigned long seed) {
  // epsilon > 0
  // b >= 1
  // P >= 2
  // 0 < delta < 1
  // alpha > gamma
  const clock_t begin_time = clock();
  total_find_P_time = 0;
  max_find_P_time = -1;
  n_find_P = 0;

  cout << "=====================================\n" 
       << "Parameters"<< endl;

  double gamma = log((P + 2 * sqrt(P + 1.0) + 2) / P);
  double c = 1 - ( pow(2, gamma - alpha) / pow(1 - 1.0/P - pow(2, gamma - alpha), 2) );
  double r = pow(2, b) / (pow(2, b) - 1);
  double rho = pow(r, 2) / (1 - epsilon);
  double kappa = 1 / c;

  int l = ceil(log(pow(2, n) / epsilon) / log(r));  // log_r(x) = log10(x) / log10(r)
  // l = 10;
  double cons_coeff = rho * kappa;
  double prob = (1 - delta) * c * pow(2, -(gamma + alpha + 1)) * P / (P - 1);

  cout << "l: " << l << "\n"
       << "b: " << b << "\n"
       << "delta: " << delta << "\n"
       << "P: " << P << "\n"
       << "alpha: " << alpha << "\n"
       << "gamma: " << gamma << "\n"
       << "prob: " << prob << "\n"
       << endl;
  
  if(epsilon <= 0) {
    cerr << "epsilon<=0" <<endl;
    exit(0);
  }
  else if(b < 1) {
    cerr << "b < 1" <<endl;
    exit(0);
  }
  else if(P < 2) {
    cerr << "P < 2" <<endl;
    exit(0);
  }
  else if((delta <= 0) || (delta >= 1)) {
    cerr << "(delta <= 0) || (delta >= 1)" <<endl;
    exit(0);
  }
  else if(alpha <= gamma) {
    cerr << "alpha <= gamma" <<endl;
    exit(0);
  }

  cout << "\n=====================================\n" << endl;

  filesystem::path data_file_path(file_path);

  time_t now = time(0);
  char* date_time = ctime(&now);
  char log_folder_name[100] = "LOG-PAWS_";
  strcat(log_folder_name, date_time);

  filesystem::path log_folder(log_folder_name);
  filesystem::path log_base_path;

  if (!filesystem::is_directory(log_folder) || !filesystem::exists(log_folder)) { // Check if src folder exists
    filesystem::create_directory(log_folder); // create src folder
  }

  log_base_path += log_folder / data_file_path.filename();

  filesystem::path main_log_path(log_base_path);
  filesystem::path cplex_log_path(log_base_path);
  filesystem::path result_log_path(log_base_path);

  main_log_path += ".log";
  cplex_log_path += ".cplex.log";
  result_log_path += ".result.log";

  cout << "log file path:"<< endl;
  cout << main_log_path << endl;
  cout << cplex_log_path << endl;
  cout << result_log_path << endl;

  std::ofstream main_log_fs(main_log_path);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(main_log_fs.rdbuf()); //redirect std::cout to out.txt!

  IloNumArray res = paws(file_path, l, b, delta, P, alpha, seed, cplex_log_path);
  // std::cout.rdbuf(coutbuf);

  std::ofstream fs(result_log_path);

  fs << "Parameters:\n"
     << "epsilon: " << epsilon << "\n"
     << "n: " << n << "\n"
     << "l: " << l << "\n"
     << "b: " << b << "\n"
     << "P: " << P << "\n"
     << "delta: " << delta << "\n"
     << "alpha: " << alpha << "\n"
     << "seed: " << seed << "\n"
     << "gamma: " << gamma << "\n"
     << "prob: " << prob << "\n"
     << "cons_coeff: " << cons_coeff << "\n"
     << endl;

  fs << "Total find P solution time: " << total_find_P_time << "\n";
  fs << "Max find P solution time: " << max_find_P_time << "\n";
  fs << "Number of find P solution: " << n_find_P << "\n";

  fs << "Time taken: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC 
       << " sec\n";
  fs << "Sample result:" << res << endl;
  cout << "Finished! LOG written to " << result_log_path << endl;

  return res;
}
