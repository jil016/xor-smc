#ifndef UTILS_H
#define UTILS_H

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


void print_matrix (vector <vector <bool> > A);
int sparsify(vector <vector <bool> > & A);
vector <vector <bool> > generate_matrix(int m, int n);
vector <vector <bool> > generate_Toeplitz_matrix(int m, int n);


#endif