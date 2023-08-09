#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "graph.h"
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
    char graph_file[1024] = "./graphs/graph1.txt";
    char output_dir[1024] = "./";

    // if (source.size() == 0){
    //     source.push_back(0);
    //     source.push_back(1);
    //     source.push_back(2);
    // }

    // if (qlist.size() == 0){
    //     qlist.resize(source.size());
    //     fill(qlist.begin(), qlist.end(), 0);
    // }

    SupplyChain  sc;

    // try to calculate T
    sc.loadParameters(graph_file, N, T, M, source, qlist, output_dir);
    sc.initializeVariables();
    sc.genFlowConstraints();
    // sc.genAllConstraints();

    return 0;
}