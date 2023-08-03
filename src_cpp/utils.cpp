#include "utils.h"


using namespace std;

typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;


void print_matrix (vector <vector <bool> > A)
{
for (unsigned int i =0;i<A.size();i++)
  {
    for (unsigned int j =0;j<A[i].size();j++)          // last column is for coefficients b
      cout << A[i][j] << ",";
    cout << endl;
  }
}

vector <vector <bool> > generate_matrix(int m, int n) {
  vector <vector <bool> > A;
  A.resize(m);
  for (int i =0;i<m;i++) {
  A[i].resize(n+1);
  for (int j =0;j<n+1;j++)          // last column is for coefficients b
    //if (rnd_uniform()<0.5)
    if (rand()%2==0)
      A[i][j] = true;
    else
      A[i][j]=false;
  }
  
  return A;  
}

int sparsify(vector <vector <bool> > & A) {
  vector< bitset<1000> > bv;

  bv.resize(A.size());
  size_t m = A.size();
  size_t n = A[0].size();
  int saved_bits = 0;
  size_t initialnbr=0;

  for (size_t i = 0;i<m;i++)
    {
    bitset<1000> row;
    for (size_t s = 0;s<n;s++)
      if (A[i][s]) 
        {
        row.set(s,1);
        initialnbr++;
        }
    bv[i]=row;  
    }

  cout << "Initial # of bits: " <<   initialnbr << endl;


  if (m<100)
  {
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      for (size_t z = 0;z<m;z++)
        for (size_t g = 0;g<m;g++)
          if (i!=l && i!=z && z!=l && i!=g && z!=g && l!=g)
          {
          int cursize = bv[i].count();
          int newsize =  (bv[i]^bv[l]^bv[z]^bv[g]).count();
          if (newsize<cursize)
            {
            saved_bits = saved_bits + cursize - newsize;
            bv[i]^=bv[l]^bv[z]^bv[g];
            for (int q=0;q<n;q++)      // sum the two rows
              A[i][q] = A[i][q] xor A[l][q] xor A[z][q] xor A[g][q];
            }
          }

  }
  if (m<500)
  {
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      for (size_t z = 0;z<m;z++)
      if (i!=l && i!=z && z!=l)
      {
      int cursize = bv[i].count();
      int newsize =  (bv[i]^bv[l]^bv[z]).count();
      if (newsize<cursize)
        {
        saved_bits = saved_bits + cursize - newsize;
        bv[i]^=bv[l]^bv[z];
        for (int q=0;q<n;q++)      // sum the two rows
          A[i][q] = A[i][q] xor A[l][q] xor A[z][q];
        }
      
      }
  }
  if (m<10000)
  {    
  for (size_t i = 0;i<m;i++)
    for (size_t l = 0;l<m;l++)
      if (i!=l)
      {
      int cursize = bv[i].count();
      int newsize =  (bv[i]^bv[l]).count();
      if (newsize<cursize)
        {
        saved_bits = saved_bits + cursize - newsize;
        bv[i]^=bv[l];
      //  print_matrix(A);
        
        for (int q=0;q<n;q++)      // sum the two rows
          A[i][q] = A[i][q] xor A[l][q];
        }
      
      }
  }

  cout << "final # of bits: " <<   (int) initialnbr- saved_bits<< endl;    
  return saved_bits;  
}

vector <vector <bool> > generate_Toeplitz_matrix(int m, int n) {
  vector <vector <bool> > A;
  if (m==0)
    return A;
    
  A.resize(m);
  int i;
  for (i =0;i<m;i++) {
    A[i].resize(n+1);
  }
  
  // first column
  for (i =0;i<m;i++) {
    if (rand()%2==0)
      A[i][0] = true;
    else
      A[i][0]=false;
    for (int j =1;j<m-i;j++)
      if (j<n)
        A[i+j][j] = A[i][0];
  }
  
  // last column
  for (i =0;i<m;i++) {
    if (rand()%2==0)
      A[i][n] = true;
    else
      A[i][n]=false;

  }
  
  // first row
  for (int j =1;j<n;j++) {
    if (rand()%2==0)
      A[0][j] = true;
    else
      A[0][j]=false;
      
    for (i =1;i<m;i++)
      if (j+i<n)
        A[i][j+i] = A[0][j];
  }

  return A;  
}