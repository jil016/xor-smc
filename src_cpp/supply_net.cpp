#include "supply_net.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

//////////////////////////////////////////////////////////
SupplyNet::SupplyNet() {}


SupplyNet::SupplyNet(string net_folder, int n_disasters){
    string produce_file("/produce.txt");
    string capacity_file("/capacity.txt");
    string demand_file("/demand.txt");
    string budget_file("/budget.txt");
    string cost_file("/cost.txt");

    ///// produce
    ifstream fp(net_folder + produce_file);
    if (! fp) {
        cout << "Error, file couldn't be opened" << endl;
        exit(1);
    }
    fp >> _N;  // read number of nodes
    fp >> _M;  // read number of raw materials
    _produce.resize(_N);

    // cout <<"_N:" << _N <<endl;
    // cout <<"_M:" << _M <<endl;

    for(int i = 0; i < _N; i++){
        fp >> _produce[i];
    }
    fp.close();

    ///// demand
    fp.open(net_folder + demand_file);
    _demand.resize(_N);
    for(int i = 0; i < _N; i++){ // initialize log_demand
        _demand[i].resize(_M);
        fill(_demand[i].begin(), _demand[i].end(), 0);
    }

    int n_demand;
    fp >> n_demand;
    for(int i = 0; i < n_demand; i++){
        int n_idx;
        int m_idx;
        fp >> n_idx;
        fp >> m_idx;
        fp >> _demand[n_idx][m_idx];
        // cout << "input demand" << endl;
        // cout << n_idx << m_idx << _demand[n_idx][m_idx] << endl;
    }
    fp.close();

    ///// capacity
    fp.open(net_folder + capacity_file);
    _capacity.resize(_N);
    fp >> _capacity_precision;
    for(int i = 0; i < _N; i++){
        _capacity[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _capacity[i][j];
            // if(_capacity[i][j] != 0){
            //     cout << "capacity (i,j)" << i << ", " << j <<"; = " << _capacity[i][j] << endl;
            // }
        }
    }
    fp.close();

    ///// budget
    fp.open(net_folder + budget_file);
    _budget.resize(_N);
    for(int i = 0; i < _N; i++){
        fp >> _budget[i];
    }
    fp.close();

    ///// cost
    fp.open(net_folder + cost_file);
    _cost.resize(_N);
    for(int i = 0; i < _N; i++){
        _cost[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _cost[i][j];
            // if(_cost[i][j] != 0){
            //     cout << "_cost (i,j)" << i << ", " << j <<"; = " << _cost[i][j] << endl;
            // }
        }
    }
    fp.close();

    // load disaster models
    _Nd = n_disasters;
    _disaster_map.resize(_Nd);
    _disaster_precision.resize(_Nd);

    for (int i = 0; i < n_disasters; i++) {
        // read each map
        string prefix("/disaster");
        string suffix(".txt");
        string full_path = net_folder + prefix + std::to_string(i) + suffix;

        fp.open(full_path);
        fp >> _disaster_precision[i];

        _disaster_map[i].resize(_N);
        for (int j = 0; j < _N; j++) {
            _disaster_map[i][j].resize(_N);
            for (int k = 0; k < _N; k++) {
                fp >> _disaster_map[i][j][k];
            }
        }
        fp.close();
    }
}

SupplyNet::~SupplyNet() {}

void SupplyNet::loadFromFile(string net_folder){

    string capacity_file("/capacity.txt");
    string demand_file("/demand.txt");
    string budget_file("/budget.txt");
    string cost_file("/cost.txt");
    string disaster_file("/disaster.txt");

    ifstream fp;

    ///// graph connection: edge capacity
    fp.open(net_folder + capacity_file);
    fp >> _N;   // N nodes
    fp >> _M;   // M layers
    _N_connect = 0;
    _raw_capacity.resize(_N);

    _edges.resize(2);
    _edge_map.resize(_N);
    for(int i = 0; i < _N; i++){ // initialize log_demand
        _edge_map[i].resize(_N);
        fill(_edge_map[i].begin(), _edge_map[i].end(), -1);
    }

    for(int i = 0; i < _N; i++){
        _raw_capacity[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _raw_capacity[i][j];
            if(_raw_capacity[i][j] > 0){
                _N_connect += 1;
                _edges[0].push_back(i);
                _edges[1].push_back(j);
                _edge_map[i][j] = _N_connect - 1;
            }
        }
    }
    fp.close();


    ///// node budget
    fp.open(net_folder + budget_file);
    _raw_budget.resize(_N);
    for(int i = 0; i < _N; i++){
        fp >> _raw_budget[i];
    }
    fp.close();

    ///// edge cost
    fp.open(net_folder + cost_file);
    _raw_cost.resize(_N);
    for(int i = 0; i < _N; i++){
        _raw_cost[i].resize(_N);
        for (int j = 0; j < _N; j++){
            fp >> _raw_cost[i][j];
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
    fp.close();

    // load disaster models
    fp.open(net_folder + disaster_file);

    char pbname[1024];
    fp >> pbname;
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
    for(int i = 0; i < _N; i++){ // initialize log_demand
        _dedge_map[i].resize(_N);
        fill(_dedge_map[i].begin(), _dedge_map[i].end(), -1);
    }

    for (int i = 0; i < _N_dedges; i++){
        _dedge_map[_dedges[0][i]][_dedges[1][i]] = i;
    }

    fp >> _N_factors;

    int arity;
    for (int i = 0; i < _N_factors; i++) {
        fp >> arity;
        vector<int> tmp_factor;
        int id;
        for (int j=0; j<arity; j++) {
            fp >> id;
            tmp_factor.push_back(id);
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
        _prec_cap(prec_cap), _prec_bgt(prec_bgt), _prec_prob(prec_prob){
    loadFromFile(net_folder);
}


void SupplyNet::discretization() {

}