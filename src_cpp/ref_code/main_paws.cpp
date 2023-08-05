// #include <stdio.h>
// #include "wish.h"
// #include "paws.h"
// #include <iostream>
// #include <fstream>

// using namespace std;


// // Usage: solve one instance
// int main(int argc, char **argv)
// {
//   WishInstance param_ins;

//   parseParityArgs(param_ins, argc, argv);  // first parse and remove parity-related command-line arguments
//   parseArgs(param_ins, argc, argv);  // now parse regular arguments



//   try {
//     param_ins.loadInstance();
//     // param_ins.getObjectExpression();
//     // param_ins.addOptimizationGoal();
//     // param_ins.sampleHashCoeffA();
//     // param_ins.getFeasibleSolution();
//     // param_ins.extractXorConstraints();
//     // param_ins.solveInstance();
//     // param_ins.removeXorConstraints();
//   } 
//   catch (IloException& ex) {
//       cout << "Error: " << ex << endl;
//   }

//   char file_name[100];
//   unsigned long seed = param_ins.seed;

//   srand(seed);
//   strcpy(file_name, param_ins.instance_name);


//   int b = 1;  // b>=1
//   int P = 4;  // P >= 2
//   int alpha = 1; // alpha > gamma0 

//   // assist params
//   int n = param_ins.nbvar;
//   double epsilon = 0.2; // epsilon > 0
//   double delta = 0.5;  // 0 < delta0 < 1
  
//   // paws gridmod_mixed_n8_w3.5_f0.5.uai -paritylevel 1 -samples 100 -nbauxv 15 -b 1 -alpha 1 -pivot 4
//   IloNumArray sol = pawsWrapper(file_name, n, epsilon, b, delta, P, alpha, seed); 

//   return 0;

// }
