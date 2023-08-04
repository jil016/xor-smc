#include "graph.h"

using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

//////////////////////////////////////////////////////////
Graph::Graph() {}


Graph::Graph(int N, int mode, int degree) {
    _N = N;
    _mode = mode;
    _degree = degree;
    
    emptyInit();
    if(mode == 0){ // full graph
        normalInit();
    }
    else if(mode == 1){ // loop free graph
        loopFreeInit();
    }
    else{
        std::cout<<"Didn't Specify Graph Init Mode.\n Empty by Default" << endl;
    }
}

Graph::~Graph() {}

void Graph::emptyInit(){
    Adj.resize(_N);
    for (int i = 0; i < _N; i++) {
        Adj[i].resize(_N);
    }

    for(int i = 0; i < _N; i++) {
        for(int j = 0; j < _N; j++){
                Adj[i][j] = 0;
        }
    }
    _M = 0;
}

void Graph::normalInit(){
    for(int i = 0; i < _N; i++) {
        for(int j = 0; j < _N; j++){
            if(i != j)
                Adj[i][j] = 1;
        }
    }
    _M = _N * _N - _N;
}

void Graph::loopFreeInit(){
    for(int i = 0; i < _N; i++) {
        for(int j = i+1; j < _N; j++){
            if(i != j)
                Adj[i][j] = 1;
        }
    }
    
    if(_degree == -1)
        _M = (_N * _N - _N) / 2;
    else{
        for(int i = 0; i < _N-_degree-1; i++) {
            for (int j = _degree + i + 1; j < _N; j++){
                Adj[i][j] = 0;
                _M++;
            }
        }
    }
}

void Graph::readFromFile(){

}


void Graph::addEdge(int u, int v){

}


void Graph::removeEdge(int u, int v){

}


void Graph::countLoopFreePathsDFS(){}
