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

   double tau = 0.01;
   int Nw = 1000;

   vector< Walker > walker(Nw);

   ifstream in("debug.walk");

   for(int i = 0;i < Nw;++i){

      bool tmp;
      double weight;

      //read
      for(int j = 0;j < global::Lx*global::Ly;++j){

         in >> tmp;

         walker[i][j] = tmp;

      }

      in >> weight;

      walker[i].sWeight(weight);

      walker[i].calc_EL();

      cout << i << "\t" << weight << "\t" << walker[i].gEL() << endl;

   }
   /*
           GFMC gfmc(tau,Nw);
    */
   return 0;

}
