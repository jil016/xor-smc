## Global Variable

```c++
ILOSTLBEGIN // use ILOG's STL namespace

static const int DEFAULT_FILTER_THRESHOLD = 4;
static const int WHEN_GRAPH_FILTERING_IS_BETTER = 5;

IloInt parity_number          = 0;
IloInt parity_minlength       = -1;
IloInt parity_maxlength       = -1;
unsigned long parity_seed;
IloBool parity_use_given_seed = false;
IloInt parity_filterLevel     = 2; 		// 0: binary representation, individual xors
 										// 1: binary representation, Gaussian elimination
 										// 2: CP variable representation, individual xors
IloInt parity_filterThreshold = DEFAULT_FILTER_THRESHOLD; 

bool PARITY_DONT_HANDLE_RANDOM_SEED = false;

bool yannakis =true;
bool jaroslow = false;
bool wainr = false;
long short_xor_max_length  = 10;
bool use_pairwise_subs = false;

// global parameters, etc.
unsigned long seed;
bool          use_given_seed = false;
IloInt        timelimit      = -1;
bool          use_tb2        = false;
char          instanceName[1024];

// before row_echelon
vector <bool>  feasiblesol;

// before main func
IloCP IlogSolver;
```



Defined

```c++
typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;
```



#### vector <vector <bool> > generate_matrix(int m, int n)

Generate a random binary matrix



#### void row_echelon(vector <vector <bool> > & A)

A[0:m, 0:n] is A

A[0:m,n] is b

First it created:

```c++
bool solvable = true;		// actually, never used
vector <int> indep_columns;
vector <int> indep_columns_rindex;
set <int> indep_vars;
```

For each row **i**, if this row is non-empty, do nothing

if empty, remove. no matter 0=0(no significance) or 0=1(unsolvable row)

else (valid constraints) insert

```c++
indep_vars.insert(j_max);			// j_max is the column index of the first 1-element
indep_columns.push_back(j_max);		// index of a basis of A
indep_columns_rindex.push_back(i);	// row index of pivot, i is the current row index
```

then apply **Gaussian Elimination** to i+1 -- m rows. Thus ensure each row i is linear independent to the rest lines.

After the loop above, **A is a row-echelon matrix**.

Then get solution which is written to the global variable `feasiblesol`. `feasiblesol` is initialized with a **n-length random binary vector**. Thus `feasiblesol = [*,*,1,*,0,*,...,*,1,*]` where `*`s are 0/1/ random numbers.



#### int sparsify(vector <vector <bool> > & A)

Sparsify matrix A



#### void add_linear_combinations(vector <vector <bool> > & A, size_t M)

From row 0 to M, xor them 2 by 2 and insert to A. In the end, there will be M_C_2 new lines added to A.



#### vector <vector <bool> > generate_matrix_maxlength(int m, int n, int k)

Another way to generate random matrix A [m * (n+1)]. For each row, set k random elements to true and others to be random 0/1, which ensures at least k true elements in each row.



#### vector <vector <bool> > generate_Toeplitz_matrix(int m, int n)

Another way to generate random matrix A [m * (n+1)]. A[0:m, 0:n] is a Toeplitz matrix which only contains (m+n-1) different elements. A[0:m, n] is still a normal random binary vector.



#### powerset_type powerset2(set_type const& set)

```c++
typedef std::set<int> set_type;
typedef std::set<set_type> powerset_type;

typedef set_type::const_iterator set_iter;
typedef std::vector<set_iter> vec;
typedef vec::iterator vec_iter;
```

Find all subsets of `set`



#### set <set <int> > powerset(set <int> s)

Find all subsets of `s`

Recursion version (not good).  



#### vector <vector <bool> > substitute_pairwise_vars(vector <vector <bool> > & A, std::vector < std::vector< std::vector < std::vector< IloBoolVar> > > > Mu)

A: m * n 

B: n * (n + n * n)  // n*n is from `Mu.size()^2`

Mu: n * n * 2 * 2 ; record the variable pairs -- this is actually new binary variables introduced by `main()`

Copy all A to B[0:m, 0:n].

For each row of B (which is A if only consider first n columns), if there are two elements A[i, s1] and A[i, s2] in a row i of A that are both 1, and `!Mu[s1][s2].empty()`, which means we have the pairwise variables already, then set

```c++
B[i][s1] = false;
B[i][s2] = false;
B[i][A[i].size() + s1 * Mu.size() + s2] = true;
// add the pairwise var
```

Otherwise, do nothing.

Finally, return B. Won't change A  or Mu.



#### int main(int argc, char **argv)

```c++
char pbname[1024];
int nbvar, nbval, nbconstr;
IloEnv env;
IloTimer timer(env);
```



```c++
// associate the CP solver with the environment
IlogSolver = IloCP(env); // global var
// associate a model with the environment
IloModel model(env);
```



Then read the instance from a `.uai` file.

`test_instance.uai` file:

```
MARKOV
10
2 2 2 2 2 2 2 2 2 2
45
2 1 0
2 2 0
2 2 1
...

4
1
0.506026334357
0.506026334357
1

4
1
0.643189361967
0.643189361967
1

...
```

`pbname` = MARKOV

`nbvar` = 10 binary variables

`nbval` = 2 // each variable can only take at most 2 different values, i.e. **maximum** domain size: 2

`vars[i].setBounds(0, tmp - 1);` set bounds [0, 1] (because vars are general `IloIntVar`, not binary)

`vars[i].setName(name);` set name to x0, x1, ..., x10

Add variables to model by `model.add(vars)`

`nbconstr` = 45 // 45 fatcors/constarints in total

`IloIntVar obj(env, 0, IloIntMax, "objective");`  // define variable that captures the value of the objective function



First read the #fatcors and related variables

In a for loop with `nbconstr` = 45 iterations:

- read `arity` = 2, //  This factor involves 2 variables

- push `1` and `0` to `vector< int> scope` // variable x1 and x0

- push `scope` to `vector<vector< int>> scopes`

- move to the next line `2, 2, 0` // next factor

Then read function tables

In a for loop with `nbconstr` = 45 iterations:

- Read `TableSize`=4

- add those 4 points to `IloNumArray table(env)` by `table.add(log10(x))` -- why log10? integer vs real?

- add `table` to `IloArray<IloNumArray> cost(env)`

**UAI Loaded**

------

**Questions so far:**

​	What are

​	`IloIntVar obj(env, 0, IloIntMax, "objective");`  // define variable that captures the value of the objective function

​	`IloNumExpr objexpr(env);`  // define cost expression

​	`std::vector < std::vector< std::vector < std::vector< IloBoolVar> > > > Mu;`

​	used for?

------

Use `Mu.resize()` to intialize `Mu` to an empty `nbvar` * `nbvar` matrix 

**(later it became  `nbvar` * `nbvar` * 2 * 2)**



A Huge For Loop (l = 0 ; l < `nbconstr`, l++) -- iterate over all factors

- `IloIntExpr pos(env);`

- Add a `IloBoolVar` to `model` by  

  ```c++
  sprintf(name, "mu_%d_%d (0,0)", (int) i, (int) j);		// i, j are 2 variables in a scope (i.e. a factor) from 0-9 (10 bvars in total)
  IloBoolVar mu_i_j_0_0 (env, 0, 1, name);				// same for mu_i_j_0_1, mu_i_j_1_0, mu_i_j_1_1
  model.add(mu_i_j_0_0);									// add
  ```

  Now, binary variables `vars` and `mu_i_j_x_x` are added to model

- Then, for general factor involves 2 bvars

  ```c++
  model.add((mu_i_j_0_0+mu_i_j_1_0 == 1-vars[j]));
  model.add((mu_i_j_0_1+mu_i_j_1_1 == vars[j]));
  model.add((mu_i_j_0_0+mu_i_j_0_1 == 1-vars[i]));
  model.add((mu_i_j_1_0+mu_i_j_1_1 == vars[i]));
  
  model.add((mu_i_j_0_1+mu_i_j_1_1 <= 1));
  model.add((mu_i_j_0_0+mu_i_j_1_0 <= 1));
  model.add((mu_i_j_1_0+mu_i_j_1_1 <= 1));
  model.add((mu_i_j_0_0+mu_i_j_0_1 <= 1));
  
  Mu[i][j][0][0]= mu_i_j_0_0 ;
  Mu[i][j][0][1]= mu_i_j_0_1 ;
  Mu[i][j][1][0]= mu_i_j_1_0 ;
  Mu[i][j][1][1]= mu_i_j_1_1 ;
  ```

  **It seems that it transformed factors to binary variables + proper constraints**. GIven xi, xj, only one of those mu_i_j_x_x could be "1".

  Besides, I don't think

  ```
  model.add((mu_i_j_0_1+mu_i_j_1_1 <= 1));
  model.add((mu_i_j_0_0+mu_i_j_1_0 <= 1));
  model.add((mu_i_j_1_0+mu_i_j_1_1 <= 1));
  model.add((mu_i_j_0_0+mu_i_j_0_1 <= 1));
  ```

  are necessary here.

- `objexpr += cost[l][0]* mu_i_j_0_0`;   // Obvoiously, `objexpr` is definitely the object function that is a combinition of `cost[l][k]* mu_i_j_x_y` where `k in {0,1,2,3}` and `(x,y) in {0,1}^2`

  if `cost[l][0]` (which is one value in range of a factor function l that can be achieved when i=0 and j = 0) is finite, otherwise set `mu_i_j_0_0`==0. Same for other `mu_i_j_x_y`

**So this loop successfully transformed the combination of 2 bvars to 4 separate bvars with a few constaints. And the object function `objexpr` has been set!**

------

**Questions Answered so far:**

​	`std::vector < std::vector< std::vector < std::vector< IloBoolVar> > > > Mu;`   // `n * n * 2 * 2`new binary variables

​	`IloNumExpr objexpr(env);`  // is the object function of new binary variables

------

`model.add(IloMaximize(env, objexpr ));`  maximizing the object function

Generate a Toeplitz martix `A` with `generate_Toeplitz_matrix(parity_number=1, nbvar=10)`  (this cpp only handles a specific parity number. We use python script `WISHCPLEX.py` to increase it incrementally from 1 to nbvars)

In the paper, we should solve for `max(w(x)) s.t. Ax = b mod 2`

`row_echelon(A)`. A becomes a row-echelon matrix. One solution is written to the global variable `feasiblesol`. `feasiblesol` is initialized with a **n-length random binary vector**. Thus `feasiblesol = [*,*,1,*,0,*,...,*,1,*]` where `*`s are 0/1/ random numbers.

------

**Summary**

We got a solution together with a feasible solution of the sampled hash function.

The objective function was loaded with variables & new variables loaded. 

**New questions**

  `IloIntVar obj(env, 0, IloIntMax, "objective");`  // define variable that captures the value of the objective function

  `use_pairwise_subs` -- not used in examples

------

Created `std::vector < std::set <size_t> > varAppearancesInXors` , size: `A[0].size() = 11` (`A is 1*10 and b is 1*1`) (see notes below)

Created `IloArray<IloArray<IloIntVarArray> > zeta_vars(env);` -- used in yannakis encoding

Created `IloArray <IloIntVarArray> alpha_vars(env); ` -- used in yannakis encoding

Created `IloBoolVar dummy_parity (env, 0, 1, "dummy");` new binary variable, also added to `model`. And more importantly, `model.add((dummy_parity==1));` 

Created `vector <size_t> xors_length;` size: `A.size() = 1` (`parity_num == 1`), `xors_length[i]` records #non-zero elements in the i-th row of `A` (actually `A` + `b` annoted in the paper)

Not 100% confident with the following code!

```c++
for (size_t j = 0; j<A.size();j++)
{
    if (!( (wainr || jaroslow) && xors_length[j]<=short_xor_max_length ))	
    // use yannakis encoding for longer ones
    for (size_t l = 0; l<A[j].size();l++)	// last column is the parity bit b
    	if (A[j][l]){
	    	// for each var, save list of xors involved	
    		varAppearancesInXors[l].insert(j);	
    	}
}
```

In test cases, `wainr` and `jaroslow` are both false. So it goes into the first `if` every time.

`varAppearancesInXors[l]` is a vector that records all `j` (row index) that `A[j][l]` (the l-th element in j-th row). Thus `varAppearancesInXors[l]` actually tells us how oftern does the l-th element (corresponding to the l-th original variable (x0-x9) is non-zero. ****



Next, `if (yannakis)` Do yannakis coding for long xors. (Otherwise, skip the yannakis encoding part. In the test cases `yannakis == true` ! )

- First,  if row `A[j]`satisfies `(!( (wainr || jaroslow) && xors_length[j]<=short_xor_max_length ))`. Then introduce new binary variables  `IloBoolVar alpha_j_k (env, 0, 1, "alpha_j_k")` where j is the row index in A and k = {0, 2, 4, ..., 2m}, 2m is the largest ever number that 2m <=  `xors_length[j]` 

  Also, add constraints `IloNumExpr alpha_sum_to_one(env);` that `\sum_{k} alpha_j_k == 1`, i.e. `alpha_sum_to_one == 1` for each row j of A, add 2m new variables and only one of them must be 1 and others are 0.

  `IloIntVarArray alphas(env);` saves all alpha_j_k with same j (a 1-D array). The formerly defined `IloArray <IloIntVarArray> alpha_vars(env)` saves all `alpha` s that have all different j. 

  In summary, `\sum_{j}(m_{j})`new binary variables are added, and `A.size() == parity_number` new constaints are added.

- Then, introduce new varibles -- `zeta`




After yannakis encoding, do jaroslaw encoding and wainwright for short xors

Since `jaroslow` and `wainr` are false in the tests, neither jaroslaw encoding or wainwright were performed.

------

**Why encoding?**

Add the constraint Aσ = b to the model properly. There are several possible encodings for the parity constraints Aσ = b mod 2. 

------



`IloCplex cplex(model);` -- Get the solver ready!

`cplex.setParam(IloCplex::Threads, 1);` -- Number of threads

`IloNumArray feasibleinit(env);`

`IloNumVarArray startVar(env);`

```c++
IloNumArray feasibleinit(env);
//double [] feasibleinit;
IloNumVarArray startVar(env);

for (size_t l= 0; l<nbvar;l++) {
    startVar.add(vars[l]);
    feasibleinit.add(feasiblesol[l]);
}
cplex.addMIPStart(startVar, feasibleinit);
```

Then

```c++
if ( !cplex.solve() ) {
env.out() << "Failed to optimize LP." << endl;
throw(-1);
}
//cout << objexpr;

IloNumArray vals(env);
env.out() << "Solution status = " << cplex.getStatus() << endl;
env.out() << "Solution value log10lik = " << cplex.getObjValue() << endl;

cplex.getValues(vals, vars);
env.out() << "Values = " << vals << endl;

} catch (IloException& ex) {
cout << "Error: " << ex << endl;
```



## Run

1. Install CPLEX Academic & Education. Default path: /opt/ibm/ILOG/CPLEX_Studio201

2. Compile (CPLEX 20.1)

   ```makefile
   ILOGBASE  = /opt/ibm/ILOG/CPLEX_Studio201
   # Add -ldl to the compiling command
   # WH_cplex: WH_cplex.cpp 
   # 	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< $(ILOGLIBS) -L. -lgmp
   WH_cplex: WH_cplex.cpp 
   	$(CC) $(OFLAGS) $(CFLAGS) -o $@ $< $(ILOGLIBS) -L. -lgmp -ldl
   ```

   Then run

   ```shell
   $ make
   ```

   Use python 2.6/2.7, install dependencies

   ```shell
   $ pip install argparse numpy scipy matplotlib
   ```

3. Run with

   ```shell
   $ python WISHCPLEX.py testInstances/clique_attractive_n15_w0.1_fg.uai LOG-test15
   ```

4. Debug cpp

   ```shell
   $ WH_cplex -paritylevel 1 -number 1 -seed 10 testInstances/clique_attractive_n10_w1.0_fg.uai
   ```

5. VSCode tips

   - VSCode No such file or directory when running c++ code
     This question was a result of confusion between the **tasks.json** and the **c_cpp_properties.json** files. I was treating c_cpp_properties.json as though it was used for compilation. **c_cpp_properties.json** is used with Intellisense and in no way deals with compilation. **tasks.json** is used for compilation. If you're unfamiliar with tasks.json as I was you need to specify the include paths here as well.

     In the args section of your tasks.json use "-I" to add an include path, followed by the path you wish to include.

     For my problem that command looked like this:

     > "-I", "C:\Users\Dill\Desktop\temp\header"





## How to read UAI file

Ref: https://www.cs.huji.ac.il/project/PASCAL/fileFormat.php



Example:

```
MARKOV
10
2 2 2 2 2 2 2 2 2 2
45
2 1 0
2 2 0
(42 lines in middle) 
2 9 8 
4
1
0.506026334357
0.506026334357
1

4
1
0.643189361967
0.643189361967
1

(42 blocks in middle)

4
1
0.517966960098
0.517966960098
1
```



- Name: MARKOV

- Number of variables: 10

- Number of elements in the domain of each variable: 2 2 2 2 2 2 2 2 2 2 (all binary, 2 possible values)

- Number of factors in object function: 45

- For each factor, e.g. factor0, 2 1 0

  - 2 variables, variable v1 and v0

  | v1   | v0   | F0             |
  | ---- | ---- | -------------- |
  | 0    | 0    | 1              |
  | 0    | 1    | 0.506026334357 |
  | 1    | 0    | 0.506026334357 |
  | 1    | 1    | 1              |

- Object Function: Sum of all factors
