#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "smc.h"

using namespace std;

void parseArgs(int argc, char **argv, string &cnf_file,
               string &uai_file, string &mar_file,
               int &threshold, int &T, int &seed,
               char output[])
{
    // one argument must be the instance filename
    if (argc <= 1) {
        // cerr << "ERROR: instance name must be specified" << endl;
        // exit(1);
        return;
    }

    for (int argIndex=1; argIndex < argc; ++argIndex) {
        if ( !strcmp(argv[argIndex], "-cnf") ) {
            argIndex++;
            cnf_file = argv[argIndex];
        }
        else if ( !strcmp(argv[argIndex], "-uai") ) {
            argIndex++;
            uai_file = argv[argIndex];
        }
        else if ( !strcmp(argv[argIndex], "-mar") ) {
            argIndex++;
            mar_file = argv[argIndex];
        }
        else if ( !strcmp(argv[argIndex], "-threshold") ) {
            argIndex++;
            threshold = atol(argv[argIndex]);
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
                 << "USAGE: XOR-SMC [options]" << endl
                 << endl
                 << "   -cnf      cnf file path" << endl
                 << "   -uai      uai file path" << endl
                 << "   -mar      mar file path" << endl
                 << "   -threshold       Threshold" << endl
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
    string cnf_file, uai_file, mar_file;

    char output_dir[1024] = "./LOG\0";
    int threshold = 0;
    int seed = 20;
    parseArgs(argc, argv, cnf_file, uai_file, mar_file, threshold, T, seed, output_dir);
    cout << "T: " <<  T << endl
         << "cnf file: " << cnf_file << endl
         << "uai file: " << uai_file << endl
         << "mar file: " << mar_file << endl
         << "threshold: " << threshold << endl
         << "output_dir: " << output_dir << endl
         << "seed: " << seed << endl;

    srand(seed);
    SMC smc_inst;
    smc_inst.loadParameters(cnf_file, uai_file, mar_file, threshold, T, output_dir);
    smc_inst.solveInstance();

    return 0;
}