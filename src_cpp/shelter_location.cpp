#include "shelter_location.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


ShelterLocation::ShelterLocation(){
  // Initialize
  timer = new IloTimer(env);
  model = new IloModel(env);
  cplex = new IloCplex(env);

  
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

    // T ...
    IloBoolVarArray shelter_assign(env);
    vector<vector<IloBoolVarArray>> bvars;
    
    // The following variables appears in the optimization
    vector<vector<vector<vector<int>>>> flows;
    vector<vector<vector<int>>> shelters;

    IloConstraintArray constraints(env); 

    bvars.resize(_T);
    flows.resize(_T);
    shelters.resize(_T);

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

    vector<vector<IloConstraintArray>> const_flow;
    const_flow.resize(_T);

    

    for (int t = 0; t < _T; t++){
        for (int s = 0; s < _src.size(); s++){
            int src = _src[s];

            IloConstraintArray const_flow_t_s(env);
            // // for the source node
            // // source in == 0
            // IloNumExpr const_src_in(env);
            // for (int i = 0; i< _N; i++){
            //     if(graph->Adj[i][src] == 1){
            //         int var_idx = flows[t][s][i][src];
            //         const_src_in = const_src_in + bvars[t][s][var_idx];
            //     }
            // }
            // const_flow_t_s.add(const_src_in == 0);

            // // source out == 1
            // IloNumExpr const_src_out(env);
            // for (int i = 0; i< _N; i++){
            //     if(graph->Adj[src][i] == 1){
            //         int var_idx = flows[t][s][src][i];
            //         const_src_out = const_src_out + bvars[t][s][var_idx];
            //     }
            // }
            // const_flow_t_s.add(const_src_out == 1);

            // all nodes
            for (int i = 0; i< _N; i++){
                // node i

                // in & out
                IloNumExpr const_mid_in(env);
                IloNumExpr const_mid_out(env);
                for (int j = 0; j< _N; j++){
                    if(graph->Adj[j][i] == 1){
                        // in 
                        int var_idx = flows[t][s][j][i];
                        const_mid_in = const_mid_in + bvars[t][s][var_idx];
                    }
                    if(graph->Adj[i][j] == 1){
                        // out 
                        int var_idx = flows[t][s][i][j];
                        const_mid_out = const_mid_out + bvars[t][s][var_idx];
                    }
                }
                const_flow_t_s.add(const_mid_in == const_mid_out);
            }

        }
    }
    // Add 
    // IloNumExpr constr_expr(env);

    // process graph
}


void addShelterConstraints(){

}


bool ShelterLocation::solveInstance() {
    cplex->clearModel();  // clear existing model
    cplex->extract(*model);

    if (timelimit > 0)
        cplex->setParam(IloCplex::TiLim, timelimit);

    cplex->setParam(IloCplex::Threads, 1);    // number of parallel threads

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