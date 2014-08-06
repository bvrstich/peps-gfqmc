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
 * constructor of the GFQMC object, takes input parameters that define the QMC walk.
 * @param dtau_in timestep
 * @param Nw_in number of Walker states
 */
GFQMC::GFQMC(double dtau_in,int Nw_in){

   this->dtau = dtau_in;
   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
GFQMC::~GFQMC(){ }

/**
 * initialize the walkers
 */
void GFQMC::SetupWalkers(){

   Walker init_walker;
   init_walker.calc_EL(peps);

   dist[0].construct(init_walker,dtau,init_walker.gEL());
   dist[0].normalize();

   walker.resize(Nw);

#pragma omp parallel for
   for(int i = 0;i < Nw;++i){

      int pick = dist[0].draw();

      walker[i] = dist[0].gwalker(pick);

      walker[i].calc_EL(peps);

   }

}

void GFQMC::walk(const int n_steps){

   //set projected energy
   sEP();

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << "\t" << EP_abs << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      double wsum = propagate();

      double scaling = Nw / wsum;

      double ET = log(scaling)/dtau;

      //calculate the energy
      sEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP << "\t" << EP_abs << endl;
      cout << "         E_T = " << ET << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double max_ov = 0.0;
      double min_ov = 1.0;

      for(int i = 0;i < walker.size();++i){
         
         if(max_ov < std::abs(walker[i].gOverlap()))
            max_ov = std::abs(walker[i].gOverlap());

         if(min_ov > std::abs(walker[i].gOverlap()))
            min_ov = std::abs(walker[i].gOverlap());

      }

#ifdef _DEBUG
      cout << "Minimal Overlap:\t" << min_ov << endl;
      cout << "Maximal Overlap:\t" << max_ov << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double GFQMC::propagate(){

   double sum = 0.0;

#pragma omp parallel for reduction(+: sum)
   for(int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(walker[i],dtau,0.0);
      double nrm = dist[myID].normalize();

      //draw new walker
      int pick = dist[myID].draw();

      walker[i] = dist[myID].gwalker(pick);

      //multiply weight
      walker[i].multWeight(nrm);

      //calculate new properties
      walker[i].calc_EL(peps);

      sum += walker[i].gWeight();

   }

   return sum;

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void GFQMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   double sum = 0.0;

   for(int i = 0;i < walker.size();i++){

      walker[i].multWeight(scaling);

      double weight = walker[i].gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.25){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + rgen_pos<double>());

         if(nCopies == 0){

#ifdef _DEBUG
            cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

            walker.erase(walker.begin() + i);

         }
         else
            walker[i].sWeight(1.0);

      }

      if(weight > 1.5){ //statically energy doesn't change

         int nCopies =(int) ( weight + rgen_pos<double>());
         double new_weight = weight / (double) nCopies;

         walker[i].sWeight(new_weight);

#ifdef _DEBUG
         cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

         for(int n = 1;n < nCopies;++n){

            Walker nw = walker[i];

            walker.push_back(nw);

         }

      }

      sum += weight;

   }

#ifdef _DEBUG
   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight is " << minw << " ." << endl;
   cout << "The max. encountered weight is " << maxw << " ." << endl;
#endif

}

/**
 * set total projected energy of the walkers at a certain timestep
 */
void GFQMC::sEP(){

   double projE_num = 0.0;
   double projE_den = 0.0;

   double projE_abs_num = 0.0;
   double projE_abs_den = 0.0;

   for(int wi = 0;wi < walker.size();wi++){

      double w_loc_en = walker[wi].gEL(); // <Psi_T | H | walk > / <Psi_T | walk >
      double overlap = walker[wi].gOverlap();

      //For the projected energy
      projE_num += walker[wi].gsign() * walker[wi].gWeight() * w_loc_en * overlap;
      projE_den += walker[wi].gsign() * walker[wi].gWeight() * overlap;

      projE_abs_num += std::abs(walker[wi].gWeight() * w_loc_en * overlap);
      projE_abs_den += std::abs(walker[wi].gWeight() * overlap);

   }

   EP = projE_num / projE_den;
   EP_abs = -projE_abs_num / projE_abs_den;

}

/**
 * write output to file
 */
void GFQMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t" << EP_abs << endl;
   output.close();

}
