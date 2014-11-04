#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;
using std::vector;

#include "include.h"

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);

   int D_aux = atoi(argv[4]);

   //initialize the dimensions of the problem, set the trial
   global::init(D,D_aux,d,L,L);

   Walker walker;

   PEPS<double> peps;
   peps.initialize_jastrow(0.74);

   Environment env(D,D_aux,1);

   env.sU(true,peps,walker);

   env.calc('H',true);
   env.calc('V',true);

   env.test(true);

/*
   double dtau = 0.01;
   int Nw = 1000;

   GFMC gfmc(dtau,Nw);
   gfmc.walk(1);
*/
   return 0;
  
}
