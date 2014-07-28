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
   global::init(D,d,L,L);
   
   //and some static objects
   Environment::init(D,D_aux);

   Walker walker;

   PEPS<double> peps;

   char filename[200];
   sprintf(filename,"/home/bright/bestanden/results/peps/output/%dx%d/D=%d",L,L,D);

   peps.load(filename);

   Environment::calc_env('H',peps,walker);

}
