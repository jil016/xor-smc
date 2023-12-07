#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "supply_net.h"
#include "supply_chain.h"

using namespace std;


void parseArgs(int argc, char **argv, string &net_folder, 
               int &target, int &T, int &seed,
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
    else if ( !strcmp(argv[argIndex], "-target") ) {
        argIndex++;
        target = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-T") ) {
        argIndex++;
        T = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-seed") ) {
        argIndex++;
        seed = atol(argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-output") ) {
        argIndex++;
        strcpy(output, argv[argIndex]);
    }
    else if ( !strcmp(argv[argIndex], "-h") || !strcmp(argv[argIndex], "-help") ) {
        cout << endl
            << "USAGE: Supply Chain [options]" << endl
            << endl
            << "   -network      Network folder name" << endl
            << "   -target       Target node" << endl
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

    int T = 1;
    string net_folder("./test_net_new");
    char output_dir[1024] = "./LOG-SPC\0";
    int target = 5;
    int seed = 20;
    parseArgs(argc, argv, net_folder, target, T, seed, output_dir);
    cout << "T: " <<  T << endl
         << "net_folder: " << net_folder << endl
         << "target: " << target << endl
         << "output_dir: " << output_dir << endl
         << "seed: " << seed << endl;

    srand(seed);
    SupplyChain sc;
    sc.loadParameters(net_folder, target, T, output_dir);
    sc.genAllConstraints();
    sc.prepareModel();
    sc.solveInstance();

    return 0;
}