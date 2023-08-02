#ifndef MODEL_H
#define MODEL_H

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

class BehaviorModel {
  // stores all variables, and can export a UAI file at anytime
  public:
    char instance_name[1024];
    int n_location;
    double w_reward;
    double w_feature;
    vector <double> reward;
    vector < vector <double> > feature_mat;
    vector < vector <double> > W_mat;
    
    
    BehaviorModel();
    void loadModelFromUai(char *file_path);
    void exportModelToUai(char *file_path);

};

#endif