#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;
using std::complex;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

using namespace global;

/**
 * constructor of the VMC object, takes input parameters that define the QMC walk.
 * @param Nw_in number of Walker states
 */
VMC::VMC(int Nw_in){

   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
VMC::~VMC(){ }

/**
 * initialize the walkers
 */
void VMC::SetupWalkers(){

   walker.resize(Nw);

   walker[0] = Walker(0);
   walker[0].calc_EL();

   walker[1] = Walker(1);
   walker[1].calc_EL();

   for(int i = 2;i < Nw;++i){

      if(i % 2 == 0)
         walker[i] = walker[0];
      else 
         walker[i] = walker[1];

   }

}

void VMC::walk(const int n_steps){

   char filename[200];
   sprintf(filename,"output/VMC/%dx%d/DT=%d/D_aux=%d.txt",Lx,Ly,DT,D_aux);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t" << endl;
   output.close();

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      propagate();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
void VMC::propagate(){

   double sum = 0.0;
   int num = 0;

#pragma omp parallel for reduction(+:sum,num)
   for(unsigned int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(walker[i]);

      //draw new walker
      int pick = dist[myID].metropolis();

      if(pick == 0)
         num++;
      else{

         walker[i] = dist[myID].gwalker(pick);
         walker[i].calc_EL();

      }

      sum += walker[i].gEL();

   }

   EP = sum / (double) Nw;

   num_stable = num;

}

/**
 * write output to file
 */
void VMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/VMC/%dx%d/DT=%d/D_aux=%d.txt",Lx,Ly,DT,D_aux);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t" << num_stable << endl;
   output.close();

}

/**
 * dump the walkers to a single file
 */
void VMC::dump(const char *filename){

   ofstream out(filename);
   out.precision(16);

   for(unsigned int i = 0;i < walker.size();++i){

      for(int row = 0;row < Ly;++row)
         for(int col = 0;col < Ly;++col)
            out << walker[i][row*Lx + col] << " ";
      out << endl;

   }

}
