#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "graph.h"
#include "shelter_location.h"

using namespace std;


void parseArgs(int argc, char **argv, char graph[], 
               int &N, int &M, int &T, 
               vector<int> &source, vector<int> &qlist,
               char output[]) 
{
  // one argument must be the instance filename
  if (argc <= 1) {
    // cerr << "ERROR: instance name must be specified" << endl;
    // exit(1);
    return;
  }

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-graph") ) {
        argIndex++;
        strcpy(graph, argv[argIndex]);
    }
    else if (argv[argIndex][0] != '-') {
        // must be the graph name
        strcpy(graph, argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-N") ) {
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
    else if ( !strcmp(argv[argIndex], "-output") ) {
        argIndex++;
        strcpy(output, argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-h") || !strcmp(argv[argIndex], "-help") ) {
        cout << endl
            << "USAGE: SMC [options]" << endl
            << endl
            << "   -graph      Graph file name" << endl
            << "   -N          Number of nodes (must skip this if read map from file)" << endl
            << "   -M          Maximum Number of Shelters" << endl
            << "   -T          Parameter T" << endl
            << "   -output     Output directory" << endl
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

    int N = -1;
    int M = 1;
    int T = 1;
    vector<int> source;
    vector<int> qlist;
    char graph_file[1024];
    char output_dir[1024];

    parseArgs(argc, argv, graph_file, N, M, T, source, qlist, output_dir);

    // if (source.size() == 0){
    //     source.push_back(0);
    //     source.push_back(1);
    //     source.push_back(2);
    // }

    // if (qlist.size() == 0){
    //     qlist.resize(source.size());
    //     fill(qlist.begin(), qlist.end(), 0);
    // }

    ShelterLocation  sl;

    // try to calculate T
    sl.loadParameters(graph_file, N, T, M, source, qlist, output_dir);
    sl.genAllConstraints();
    sl.solveInstance();

    return 0;
}