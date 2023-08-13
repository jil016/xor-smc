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

    for (int i = 0; i < n_disasters; i++){
        // read each map
        string prefix("/disaster");
        string suffix(".txt");
        string full_path = net_folder + prefix + std::to_string(i) + suffix;

        fp.open(full_path);
        fp >> _disaster_precision[i];

        _disaster_map[i].resize(_N);
        for(int j = 0; j < _N; j++){
            _disaster_map[i][j].resize(_N);
            for (int k = 0; k < _N; k++){
                fp >> _disaster_map[i][j][k];
            }
        }
        fp.close();
    }

    genAdvancedInfo();
}

SupplyNet::~SupplyNet() {}

void SupplyNet::genAdvancedInfo(){
    _producers.resize(_M);

    for (int i = 0; i < _N; i++){
        _producers[_produce[i]].push_back(i);
    }
}

