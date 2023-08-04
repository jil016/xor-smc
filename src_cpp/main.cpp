#include <stdio.h>
#include <iostream>
#include <fstream>

#include "graph.h"
#include "shelter_location.h"

using namespace std;


void shelterDesign(){
    
}


int main(int argc, char **argv)
{
    // Need to specify type of graph, 
    // source node
    // params: graph, src, q_list, N, T, M
    int N = 100;
    int T = 1;
    IloEnv env;

   
    vector<int> src{0, 1, 2};
    vector<int> q_vect(N, 2);
    ShelterLocation  sl;

    sl.loadParameters(N, T, 1, src, q_vect);
    sl.genAllConstraints();
    sl.solveInstance();

    // IloModel model;
    // model = IloModel(env);
    // IloCplex cplex;
    // cplex = IloCplex(env);

    // IloBoolVarArray bvars(env, 3);
    // IloConstraintArray cons(env);

    // int mode = 1;
    // int degree = 2;
    // Graph g(N, mode, degree);
    // vector<vector<vector<IloBoolVar>>> flows;

    // flows.resize(3);
    // for (int i = 0; i< T; i++){
    //     flows[i].resize(N);
    //     for (int j = 0; j< N; j++){
    //         flows[i][j].resize(N);
    //         for(int k = 0; k< N; k++){
    //             flows[i][j][k] = IloBoolVar(env);
    //             model.add(flows[i][j][k]);
    //         }
    //     }
    // }

    // cons.add(flows[0][0][0] - flows[0][0][1] - flows[0][1][0] == -1);
    // cons.add(flows[0][0][0] + flows[0][0][1] - flows[0][1][0] == -1);

    // model.add(cons);

    // cplex.clearModel();  // clear existing model
    // cplex.extract(model);
    // cplex.setParam(IloCplex::Threads, 1);    // number of parallel threads
    // cplex.solve();
    // // env.out() << "Solution status = " << cplex.getStatus() << endl;
    // env.out() << "Solution value = " << cplex.getCplexStatus() << endl; 
    // env.out() << "Solution bvars value = " << cplex.getValue(flows[0][0][0]) << endl; 
    // env.out() << "Solution bvars value = " << cplex.getValue(flows[0][0][1]) << endl; 
    // env.out() << "Solution bvars value = " << cplex.getValue(flows[0][1][0]) << endl; 

    return 0;
}