#include "shelter_location.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


ShelterLocation::ShelterLocation(){
  // Initialize
  //   timer = new IloTimer(env);
  //   model = new IloModel(env);
  //   cplex = new IloCplex(env);
  model = IloModel(env);
  cplex = IloCplex(env);
  
  sparsify_coeff = true;
  n_vars = 0;

  yannakis =true;

  // global parameters, etc.
  timelimit = -1;
}

ShelterLocation::~ShelterLocation(){}

void ShelterLocation::loadParameters(int N, int T, int M, vector<int> src, vector<int> q){
    graph = new Graph(N, 1, 2);
    _N = N;
    _T = T;
    _M = M;
    _src = src;
    _q = q;
}


void ShelterLocation::addFlowConstraints(){
    // There will be T * n_src complete flows
    shelter_assign = IloBoolVarArray(env, _N);

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
    vector<vector<IloConstraintArray>> const_flow;
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
                const_flow_t_s.add(degree_abs == (bvars[t][s][var_idx] - sources[s][i]));
            }

            const_flow[t].push_back(const_flow_t_s);
        }
    }

    // Extract the shelter
    ////////////////////////////////////////
    vector<vector<IloConstraintArray>> const_is_shelter;
    const_is_shelter.resize(_T);
    
    for (int t = 0; t < _T; t++){
        for (int s = 0; s < _src.size(); s++){
            // shelters[t][s] = N dimensional (0,0,0,1,0,0,0,...)
            IloConstraintArray const_number_t_s(env);
            
            for (int i = 0; i< _N; i++){
                int var_idx = shelters[t][s][i];
                const_number_t_s.add((shelter_assign[i] - bvars[t][s][var_idx]) >= 0);
            }
            const_is_shelter[t].push_back(const_number_t_s);
        }
    }

    // Extract the number constraint
    ////////////////////////////////////////
    IloConstraintArray const_max_shelter(env);
    IloNumExpr count_shelters(env);
    for(int i = 0; i < _N; i++){
        count_shelters = count_shelters + shelter_assign[i];
    }
    const_max_shelter.add(count_shelters < _M + 1); 


    ////////////////////////////////////////
    model.add(shelter_assign);

    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            model.add(bvars[t][s]);
        }
    }
    for (int t = 0; t< _T; t++){
        for (int s = 0; s< _src.size(); s++){
            model.add(const_flow[t][s]);
            model.add(const_is_shelter[t][s]);
        }
    }
    model.add(const_max_shelter);
    

    cplex.clearModel();  // clear existing model
    cplex.extract(model);
    cplex.setParam(IloCplex::Threads, 1);    // number of parallel threads
    cplex.solve();
    cplex.exportModel("model.lp");
    // env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << "Solution value = " << cplex.getCplexStatus() << endl; 
    for(int i = 0; i<shelter_assign.getSize(); i++){
        env.out() << "shelter_assign = " << cplex.getValue(shelter_assign[i]) << endl; 
    }
    for(int i = 0; i<bvars[0][0].getSize(); i++){
        env.out() << "bvars = " << cplex.getValue(bvars[0][0][i]) << endl; 
    }

}


void ShelterLocation::addXORConstraints(){

}


void ShelterLocation::parseAllConstraints(){

}


bool ShelterLocation::solveInstance() {
    // cplex->clearModel();  // clear existing model
    // cplex->extract(*model);

    // if (timelimit > 0)
    //     cplex->setParam(IloCplex::TiLim, timelimit);

    // cplex->setParam(IloCplex::Threads, 1);    // number of parallel threads

    // if (!coeffA.empty()) {
    //     IloNumArray feasibleinit(env);
    //     //double [] feasibleinit;
    //     IloNumVarArray startVar(env);
    //     for (size_t l= 0; l<nbvar;l++) {
    //         startVar.add((*(vars))[l]);
    //         feasibleinit.add(feasiblesol[l]);
    //     }
    // cplex->addMIPStart(startVar, feasibleinit);
    // // feasibleinit.end();  // https://or.stackexchange.com/questions/3530/no-solution-found-from-n-mip-starts
    // // startVar.end();
    // }

    // if ( !cplex->solve() ) {
    // env.out() << "Failed to optimize LP." << endl;
    // return false;
    // }

    // IloNumArray vals(env);
    // env.out() << "Solution status = " << cplex->getStatus() << endl;
    // env.out() << "Solution value = " << cplex->getObjValue() << endl;

    // cplex->getValues(vals, *vars);
    // // env.out() << "Values = " << vals << endl;
    return true;
}