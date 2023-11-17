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
    int _M;  // number of layers (raw materials)
    int _N_connect;
    int _N_end;
    vector<int> _end_nodes;


    vector<vector<double>> _raw_capacity; // size: _N * _N; capacity[i][j] is the capacity of trading edge from i to j
                                          // capacity == 0 means no demand or supply relation.
    vector<vector<double>> _raw_cost;
    vector<double> _raw_budget;

    vector<vector<int>> _edges;
    vector<vector<int>> _edge_map;

    // disaster model -- UAI
    int _N_dedges;  // disaster edges
    vector<vector<int>> _dedges;
    vector<vector<int>> _dedge_map;

    int _N_factors; //
    vector<vector<int>> _factors;
    vector<vector<double>> _tables;

    // discretization
    int _prec_cap, _prec_cst, _prec_bgt, _prec_prob;
    void discretization();

    SupplyNet();
    SupplyNet(string net_folder, int prec_cap, int prec_cst, int prec_bgt, int prec_prob);
    ~SupplyNet();
    void loadFromFile(string net_folder);


    ///////////////////////////////////////////////////////////
    // Read from input
    vector<int> _produce;    // size: _N * 1; node i produces produce[i]
    vector<vector<int>> _demand;  // size: _N * _M; node i requires demand[i][m] unit of material m
    vector<vector<int>> _capacity;
    vector<int> _budget;
    vector<vector<int>> _cost;
    int _capacity_precision; //

    // disasters model ...
    // easy to read
    int _Nd;           // the number of disasters
    vector<int> _disaster_precision; // the number of discretized values, e.g., [0.25, 0.5, 0.75, 1] => _precision=4
    vector<vector<vector<int>>> _disaster_map; // disaster i affects edges in disaster_map[i] with a probabiliy disaster_map[i][j][k]


    // Advanced information
    vector<vector<int>> _producers; // size: _M * UNSURE; material m has produces producers[m]


    SupplyNet(string net_folder, int n_disasters);
};


#endif