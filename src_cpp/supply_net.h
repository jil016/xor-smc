#ifndef SUPPLY_NET_H
#define SUPPLY_NET_H

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

class SupplyNet {
  public:

    int _N;  // number of nodes
    int _N_edges; // number of edges
    int _N_end; // number of end nodes
    vector<int> _end_nodes;


    vector<vector<double>> _raw_capacity; // size: _N * _N; capacity[i][j] is the capacity of trading edge from i to j
    vector<vector<int>> _dis_capacity;  // discretized capacity

    vector<vector<double>> _raw_cost;
    vector<double> _raw_budget;
    
    vector<vector<int>> _edges;
    vector<vector<int>> _edge_map;

    int _demand;


    // disaster model -- UAI
    int _N_dedges;  // disaster edges
    vector<vector<int>> _dedges;
    vector<vector<int>> _dedge_map;
    vector<vector<int>> _dedge_watchers;    // i-th dedge appears in factors _dedge_watchers[i]

    // UAI properties
    int _N_factors; //
    vector<vector<int>> _factors;
    vector<vector<double>> _tables;

    // For discretization
    int _prec_cap, _prec_prob;    // precision represented by the number of bits
    double _min_cap, _max_cap, _min_prob, _max_prob;
    void discretizeAll();

    SupplyNet();
    SupplyNet(string net_folder, int prec_cap, int prec_cst, int prec_bgt, int prec_prob);
    ~SupplyNet();
    void loadFromFile(string net_folder);

    /// Unused discretization
    // vector<vector<int>> _dis_cost;   // not used
    // vector<int> _dis_budget;         // not used
    int _prec_cst, _prec_bgt;   // not used
    double _min_cst, _max_cst, _min_bgt, _max_bgt;  // not used
};


#endif