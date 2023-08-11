#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "supply_net.h"
#include "supply_chain.h"

using namespace std;


int main(int argc, char **argv)
{
    // Need to specify type of graph, 
    // source node
    // params: graph, src, q_list, N, T, M

    int N = -1;
    int M = 1;
    int T = 1;
    vector<int> source;
    vector<int> qlist{0,0,0,0};
    string net_folder("./networks/net1");
    int n_disaster = 1;
    SupplyNet sn(net_folder, n_disaster);

    char output_dir[] = "output";

    SupplyChain sc;
    sc.loadParameters(net_folder, n_disaster, 1, output_dir);
    sc.initializeVariables();

    return 0;
}