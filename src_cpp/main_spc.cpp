#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "supply_net.h"
#include "supply_chain.h"

using namespace std;


void parseArgs(int argc, char **argv, string &net_folder, 
               int &n_disaster, int &T,
               char output[]) 
{
  // one argument must be the instance filename
  if (argc <= 1) {
    // cerr << "ERROR: instance name must be specified" << endl;
    // exit(1);
    return;
  }

  for (int argIndex=1; argIndex < argc; ++argIndex) {
    if ( !strcmp(argv[argIndex], "-network") ) {
        argIndex++;
        net_folder = argv[argIndex];
    }
    else if (argv[argIndex][0] != '-') {
        // must be the graph name
        net_folder = argv[argIndex];
    }
    else if ( !strcmp(argv[argIndex], "-Nd") ) {
        argIndex++;
        n_disaster = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-T") ) {
        argIndex++;
        T = atol(argv[argIndex]);
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

    int T;
    string net_folder("./networks/net1");
    int n_disaster = 1;
    char output_dir[] = "./LOG-SPC\0";

    parseArgs(argc, argv, net_folder, n_disaster, T, output_dir);

    cout << "T: " <<  T << endl
         << "net_folder: " << net_folder << endl
         << "n_disaster: " << n_disaster << endl
         << "output_dir: " << output_dir << endl;


    SupplyChain sc;
    sc.loadParameters(net_folder, n_disaster, T, output_dir);
    sc.genAllConstraints();
    sc.solveInstance();

    return 0;
}