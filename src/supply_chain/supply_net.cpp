#include "supply_net.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

//////////////////////////////////////////////////////////
SupplyNet::SupplyNet() {}
SupplyNet::~SupplyNet() {}

void SupplyNet::loadFromFile(string net_folder){
    string capacity_file("/capacity.txt");
    string demand_file("/demand.txt");
    string budget_file("/budget.txt");
    string cost_file("/cost.txt");
    string disaster_file("/disaster.txt");
    string disaster_uai_file("/disaster.uai");
    string disaster_uai_edge_file("/disaster.uai.edges");

    ifstream fp;

    ///// graph connection: edge capacity
    fp.open(net_folder + capacity_file);
    fp >> _N;   // N nodes
    _N_edges = 0;
    _raw_capacity.resize(_N);

    _edges.resize(2);
    _edge_map.resize(_N);
    for(int i = 0; i < _N; i++){ // initialize log_demand
        _edge_map[i].resize(_N);
        fill(_edge_map[i].begin(), _edge_map[i].end(), -1);
    }

    _min_cap = 999999;
    _max_cap = 0;
    for(int i = 0; i < _N; i++){
        _raw_capacity[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _raw_capacity[i][j];
            if(_raw_capacity[i][j] > 0){
                _N_edges += 1;
                _edges[0].push_back(i);
                _edges[1].push_back(j);
                _edge_map[i][j] = _N_edges - 1;

                _min_cap = std::min(_raw_capacity[i][j], _min_cap);
                _max_cap = std::max(_raw_capacity[i][j], _max_cap);;
            }
        }
    }
    fp.close();

    ///// node budget
    fp.open(net_folder + budget_file);
    _min_bgt = 999999;
    _max_bgt = 0;
    _raw_budget.resize(_N);
    for(int i = 0; i < _N; i++){
        fp >> _raw_budget[i];
        _min_bgt = std::min(_raw_budget[i], _min_bgt);
        _max_bgt = std::max(_raw_budget[i], _max_bgt);
    }
    fp.close();

    ///// edge cost
    fp.open(net_folder + cost_file);
    _min_cst = 999999;
    _max_cst = 0;
    _raw_cost.resize(_N);
    for(int i = 0; i < _N; i++){
        _raw_cost[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _raw_cost[i][j];
            if(_raw_capacity[i][j] > 0){    // there is an edge
                _min_cst = std::min(_raw_cost[i][j], _min_cst);
                _max_cst = std::max(_raw_cost[i][j], _max_cst);
            }
        }
    }
    fp.close();


    //// terminal nodes & demand
    fp.open(net_folder + demand_file);
    fp >> _N_end;
    for(int i = 0; i < _N_end; i++){
        int id;
        fp >> id;
        _end_nodes.push_back(id);
    }
    fp >> _demand;
    fp.close();

    // Read disaster edges
    fp.open(net_folder + disaster_uai_edge_file);
    fp >> _N_dedges;
    _dedges.resize(2);
    // read edges
    int tmp;
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < _N_dedges; j++){
            fp >> tmp;
            _dedges[i].push_back(tmp);
        }
    }
    _dedge_map.resize(_N);
    for(int i = 0; i < _N; i++){ // initialize demand
        _dedge_map[i].resize(_N);
        fill(_dedge_map[i].begin(), _dedge_map[i].end(), -1);
    }
    for (int i = 0; i < _N_dedges; i++){
        _dedge_map[_dedges[0][i]][_dedges[1][i]] = i;
    }
    fp.close();

    // Load disaster models
    fp.open(net_folder + disaster_uai_file);

    char pbname[1024];
    fp >> pbname;
    fp >> tmp;
    if(tmp != _N_dedges) exit(0);   // doesn't match the uai.edge file
    // read scopes
    for (int j = 0; j < _N_dedges; j++){
        fp >> tmp;  // scopes are all 2 by default
        if(tmp != 2) exit(0);
    }

    fp >> _N_factors;
    int arity;
    _dedge_watchers.resize(_N_dedges);
    for (int i = 0; i < _N_factors; i++) {
        fp >> arity;
        vector<int> tmp_factor;
        int id;
        for (int j=0; j<arity; j++) {
            fp >> id;
            tmp_factor.push_back(id);
            _dedge_watchers[id].push_back(i);
        }
        _factors.push_back(tmp_factor);
    }

    // read in values of CPT tables
    int table_size;
    for (int i=0; i<_N_factors; i++) {
        fp >> table_size;
        double prob;
        vector<double> tmp_table;
        for (int j=0; j<table_size; j++) {
            fp >> prob;
            tmp_table.push_back(prob);
        }
        _tables.push_back(tmp_table);
    }
    cout << "done reading UAI"<< endl;
    fp.close();
}


SupplyNet::SupplyNet(string net_folder, int prec_cap, int prec_cst, int prec_bgt, int prec_prob):
        _prec_cap(prec_cap), _prec_cst(prec_cst), _prec_bgt(prec_bgt), _prec_prob(prec_prob){
    loadFromFile(net_folder);
    // Don't discretize things NOW!
    // discretizeAll();
}


//void SupplyNet::discretizeAll() {
//    // Discretize capacities.
//    // [min_cap, max_cap] -> [1, ..., 2^prec_cap-1]
//    _dis_capacity.resize(_N);
//    int range = (1 << _prec_cap) - 1;
//    for (int i = 0; i < _N; i++) {
//        _dis_capacity[i].resize(_N);
//        for (int j = 0; j < _N; j++) {
//            if (_raw_capacity[i][j] == 0) {   // no edge
//                _dis_capacity[i][j] = 0;
//            } else if (_min_cap == _max_cap) {
//                _dis_capacity[i][j] = range;
//            } else if (_min_cap < 1 || _max_cap > range) {
//                // exceed effective digits, normalize
//                double clamped_value = std::max(_min_cap, std::min(_max_cap, _raw_capacity[i][j]));
//                double normalized_value = (clamped_value - _min_cap) / (_max_cap - _min_cap);
//                _dis_capacity[i][j] = static_cast<int>(std::round(normalized_value * (range - 1))) + 1;
//            } else {
//                // simply cast
//                _dis_capacity[i][j] = static_cast<int>(std::round(_raw_capacity[i][j]));
//            }
//        }
//    }
//
//    // TODO: Discretize probabilities
//    //  Suppose the discretization is done by special design
//}
