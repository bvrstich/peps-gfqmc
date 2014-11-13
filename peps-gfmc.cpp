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

   global::env[0].calc('H',false,walker);
   global::env[0].calc('H',true,walker);

   global::env[0].calc('V',false,walker);
   global::env[0].calc('V',true,walker);

   global::env[0].test(true);
   global::env[0].test(false);

  // double dtau = 0.01;
   //int Nw = 1000;

 //  GFMC gfmc(dtau,Nw);
//   gfmc.walk(1000000);

   return 0;

}
