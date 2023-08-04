#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "graph.h"
#include "shelter_location.h"

using namespace std;


void parseArgs(int argc, char **argv, int &N, int &M, int &T, vector<int> &source, vector<int> &qlist) 
{
  // one argument must be the instance filename
  if (argc <= 1) {
    // cerr << "ERROR: instance name must be specified" << endl;
    // exit(1);
    return;
  }

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-N") ) {
        argIndex++;
        N = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-M") ) {
        argIndex++;
        M = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-T") ) {
        argIndex++;
        T = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-source") ) {
        argIndex++;
        stringstream ss(argv[argIndex]);
        while (ss.good()) {
            string substr;
            getline(ss, substr, ',');
            source.push_back(stoi(substr));
        }
    }
    else if ( !strcmp(argv[argIndex], "-qlist") ) {
        argIndex++;
        stringstream ss(argv[argIndex]);
        while (ss.good()) {
            string substr;
            getline(ss, substr, ',');
            qlist.push_back(stoi(substr));
        }
    }
    else if ( !strcmp(argv[argIndex], "-h") || !strcmp(argv[argIndex], "-help") ) {
        cout << endl
            << "USAGE: SMC [options]" << endl
            << endl
            << "   -N          Number of Nodes" << endl
            << "   -M          Maximum Number of Shelters" << endl
            << "   -T          Parameter T" << endl
            << endl;
        // print parity constraint options usage
        //printParityUsage(cout);
        exit(0);
    }
    else {
      cerr << "ERROR: Unexpected option: " << argv[argIndex] << endl;
      exit(1);
    }
  }
}


int main(int argc, char **argv)
{
    // Need to specify type of graph, 
    // source node
    // params: graph, src, q_list, N, T, M

    int N = 10;
    int M = 1;
    int T = 1;
    vector<int> source;
    vector<int> qlist;

    parseArgs(argc, argv, N, M, T, source, qlist);

    if (source.size() == 0){
        source.push_back(0);
        source.push_back(1);
        source.push_back(2);
    }

    if (qlist.size() == 0){
        qlist.resize(source.size());
        fill(qlist.begin(), qlist.end(), 0);
    }

    ShelterLocation  sl;

    // try to calculate T

    sl.loadParameters(N, T, M, source, qlist);
    sl.genAllConstraints();
    sl.solveInstance();

    return 0;
}