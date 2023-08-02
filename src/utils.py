import numpy as np


def generate_Toeplitz_matrix(m, n):
    A = np.random.randint(2, size=(m, n), dtype=bool)
    # first column
    for i in range(m):
        for j in range(m - i):
            if j < n:
                A[i+j, j] = A[i, 0]

    
    # first row
    for j in range(1, n):
        for i in range(1, m):
            if(i + j < n):
                A[i, j+i] = A[0, j]

    return A


def generate_XOR_matrix(m, n):
    A = generate_Toeplitz_matrix(m,n)
    b = np.random.randint(2, size=(m, 1), dtype=bool)

    pass



if __name__ == "__main__":
    print(generate_Toeplitz_matrix(4,4))



	

# bool WishInstance::getFeasibleSolution() {
#   feasiblesol.clear();
#   if (coeffA.empty()) 
#     return true;

#   bool solvable = true;            // is the system A x + b =0 solvable?
  
#   size_t m = coeffA.size();
#   size_t n = coeffA[0].size()-1;
  
#   vector <int> indep_columns;
#   vector <int> indep_columns_rindex;
#   set <int> indep_vars;
    
#   // put A in row echelon form
#   for (size_t i = 0;i<m;i++) {
#     //Find pivot for column k:
#     bool empty_row=true;
#     int j_max =0;
#     for (int s = 0;s<n;s++)
#       if (coeffA[i][s]) {
#         empty_row=false;
#         j_max = s;
#         break;
#       }
#     if (empty_row){        // low rank
#       if (coeffA[i][n]) {    //0=1
#         solvable = false;
#         cout << "Sampled Ax = b is not solvable!" << endl;
#         return false;
#       }
#     }
#     else {
#       indep_vars.insert(j_max);
#       indep_columns.push_back(j_max);          // index of a basis of coeffA
#       indep_columns_rindex.push_back(i);        // row index of pivot
        
#       for (size_t h=i+1;h<m;h++)
#         if (coeffA[h][j_max]) {      // if not already zero
#           for (int q=0;q<n+1;q++)      // sum the two rows
#             coeffA[h][q] = coeffA[h][q] xor coeffA[i][q];
#         }
#     }
#   }
  
#   for (size_t i = 0;i<indep_columns.size();i++) {
#     int j_max = indep_columns[i];
#     int p = indep_columns_rindex[i];

#     for (int h=p-1;h>=0;h--)
#       if (coeffA[h][j_max]) {      // if not already zero
#           //print_matrix(coeffA);

#           for (int q=0;q<n+1;q++)      // sum the two rows
#             coeffA[h][q] = coeffA[h][q] xor coeffA[p][q];
#           //print_matrix(coeffA);
#         }
#   }
  
#   // produce a solution
#   vector <bool>  b;
#   vector <bool>  y;
#   feasiblesol.resize(n);
#   y.resize(n);
  
#   // initialize b to the last column of A
#   b.resize(m);
#   for (size_t i =0;i<m;i++)
#     b[i] = coeffA[i][n];
  
#   for (size_t i =0;i<n;i++)
#     y[i] = rand()%2;
    
#   // sum all the dependent variables that are already fixed
#   for (size_t i =0;i<n;i++) {
#     feasiblesol[i] = y[i];
#     if ( (indep_vars.count(i)==0) && (y[i]==1)) {    // dependent variable, and non zero
#       // b = b + x[i] * A[] [i]
#       for (size_t j =0;j<m;j++)
#         b[j] = b[j] xor coeffA[j][i];
#     }
#   }
    
#   // backsubstitute r
#   for (int i =indep_columns_rindex.size()-1;i>=0;i--) {
#     int c = indep_columns_rindex[i];    // lowest pivot
#     if (b[c]==1) {    // we need to add a 1
#       y[indep_columns[i]] = 1;
#       feasiblesol[indep_columns[i]] = 1;
#       for (size_t j =0;j<m;j++)
#         b[j] = b[j] xor coeffA[j][indep_columns[i]];
#     }
#     else {
#       y[indep_columns[i]] = 0;
#       feasiblesol[indep_columns[i]] = 0;
#     }
#   }

#   // print_matrix(coeffA);
#   cout << "row echelon form:" << endl;
#   print_matrix(coeffA);
  
#   // sparsify
#   if (sparsify_coeff && (!coeffA.empty())) {
#     cout << "Bits saved: ";
#     for (int i=0; i<2; i++)
#       cout << sparsify(coeffA) << " ";
#     cout << endl;  
#   }
#   return true;
# }