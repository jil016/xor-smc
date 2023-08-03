#ifndef GRAPH_H
#define GRAPH_H

#include <new>
#include <set>
#include <cmath>
#include <bitset>
#include <vector>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include <stdio.h>
#include <sys/time.h>

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

class Graph {
  public:
    int _N;
    int _M;
    int _mode;
    int _degree;

    // From arguments 
    unsigned long seed;
    bool use_given_seed;

    std::vector < std::vector <int> > Adj;

    char instance_name[1024];

    Graph();
    Graph(int n, int mode, int degree);
    ~Graph();

    void emptyInit();
    void normalInit();
    void loopFreeInit();
    void addEdge(int u, int v);
    void removeEdge(int u, int v);
    void countLoopFreePathsDFS();
    // void forceBreakLoopsBFS();
};


#endif