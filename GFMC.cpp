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
 * constructor of the GFMC object, takes input parameters that define the QMC walk.
 * @param dtau_in timestep
 * @param Nw_in number of Walker states
 */
GFMC::GFMC(double dtau_in,int Nw_in){

   this->dtau = dtau_in;
   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * destructor...
 */
GFMC::~GFMC(){ 

   //remove the walkers
   for(int i = 0;i < walker.size();++i)
      delete walker[i];
   
}

/**
 * initialize the walkers: read in distribution from walkers dir
 */
void GFMC::SetupWalkers(){

   walker.resize(Nw);

   walker[0] = new Walker(0);
   walker[0]->update_env(-1);
   walker[0]->calc_EL();

   walker[1] = new Walker(1);
   walker[1]->update_env(-1);
   walker[1]->calc_EL();

   for(int i = 2;i < Nw;++i){

      if(i % 2 == 0)
         walker[i] = new Walker(*walker[0]);
      else 
         walker[i] = new Walker(*walker[1]);

   }

}

void GFMC::walk(const int n_steps){

   //set projected energy
   sEP();

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << "\t" << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      double wsum = propagate();

      double scaling = Nw / wsum;

      ET = 1.0/dtau * ( 1 - 1.0/scaling);

      //calculate the energy
      sEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP << endl;
      cout << "         E_T = " << ET << endl;
      cout << "    # stable = " << num_stable << endl;
      cout << "    avg sign = " << avs << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(step,scaling);

      double max_en = -100.0;
      double min_en = 100.0;

      for(unsigned int i = 0;i < walker.size();++i){

         if(max_en < walker[i]->gEL())
            max_en = walker[i]->gEL();

         if(min_en > walker[i]->gEL())
            min_en = walker[i]->gEL();

      }

#ifdef _DEBUG
      cout << endl;
      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Maximal Energy:\t" << max_en << endl;
      cout << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double GFMC::propagate(){

   double sum = 0.0;
   int num = 0;

#pragma omp parallel for reduction(+: sum,num)
   for(unsigned int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(*walker[i],dtau,0.0);
      dist[myID].check_negative();

      double nrm = dist[myID].normalize();

      //draw new walker
      int pick = dist[myID].draw();

      if(pick == 0)//nothing happens
         num++;
      else{

         walker[i]->sconf( dist[myID].gconf(pick) );

         if(dist[myID].gsign_flip(pick))
            walker[i]->sign_flip();

         //calculate new properties
         walker[i]->update_env(-1);//dist[myID].gind(pick));

         walker[i]->calc_EL();

      }

      //multiply weight
      walker[i]->multWeight(nrm);

      sum += walker[i]->gWeight();

   }

   num_stable = num;

   return sum;

}

/**
 * @param step when step is certain number, do the population control
 * @param scaling factor with which to rescale the walkers
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void GFMC::PopulationControl(int step,double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   for(unsigned int i = 0;i < walker.size();i++)
      walker[i]->multWeight(scaling);

   //if(step % 100 == 0){

      //remove those too small
      for(unsigned int i = 0;i < walker.size();i++){

         double weight = walker[i]->gWeight();

         if(weight < minw)
            minw = weight;

         if(weight > maxw)
            maxw = weight;

         if (weight < 0.10){ //Energy doesn't change statistically

            int nCopies = (int) ( weight + RN());

            if(nCopies == 0){

#ifdef _DEBUG
               cout << "Walker with weight " << weight << " will be deleted." << endl;
#endif

               delete walker[i];

               walker.erase(walker.begin() + i);
               --i;

            }
            else
               walker[i]->sWeight(1.0);

         }

      }

      //multiply those too large
      int pop = walker.size();

      for(int i = 0;i < pop;i++){

         double weight = walker[i]->gWeight();

         if(weight > 2.0){ //statically energy doesn't change

            int nCopies =(int) ( weight + RN());
            double new_weight = weight / (double) nCopies;

            walker[i]->sWeight(new_weight);

#ifdef _DEBUG
            cout << "Walker with weight " << weight << " will be " << nCopies << " times copied with weight " << new_weight << "." << endl;
#endif

            for(int n = 1;n < nCopies;++n){

               Walker *nw = new Walker(*walker[i]);

               walker.push_back(nw);

            }

         }

      }

      //get the new weight
      double sum = 0.0;

      for(unsigned int i = 0;i < walker.size();++i)
         sum += walker[i]->gWeight();

      //rescale the weights to unity for correct ET estimate in next iteration
      for(unsigned int i = 0;i < walker.size();++i)
         walker[i]->multWeight((double)Nw/sum);

#ifdef _DEBUG
      cout << endl;
      cout << "total weight:\t" << sum << endl;
      cout << endl;

      cout << "The min. encountered weight is " << minw << " ." << endl;
      cout << "The max. encountered weight is " << maxw << " ." << endl;
      cout << endl;
#endif

//   }

}

/**
 * set total projected energy of the walkers at a certain timestep
 */
void GFMC::sEP(){

   double projE_num = 0.0;
   double projE_den = 0.0;

   double sign_num = 0.0;
   double sign_den = 0.0;

   for(unsigned int wi = 0;wi < walker.size();wi++){

      double w_loc_en = walker[wi]->gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num += walker[wi]->gsign() * walker[wi]->gWeight() * w_loc_en;
      projE_den += walker[wi]->gsign() * walker[wi]->gWeight();

      sign_num += walker[wi]->gsign() * walker[wi]->gWeight();
      sign_den += walker[wi]->gWeight();

   }

   EP = projE_num / projE_den;
   avs = sign_num / sign_den;

}

/**
 * write output to file
 */
void GFMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/%dx%d/D=%d.txt",Lx,Ly,DT);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t" << ET << "\t" << num_stable << "\t" << avs << endl;
   output.close();

}

/**
 * dump the walkers to a single file
 */
void GFMC::dump(const char *filename){

   ofstream out(filename);
   out.precision(16);

   for(unsigned int i = 0;i < walker.size();++i){

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c)
            out << (*walker[i])[r*Lx + c] << " ";
      out << walker[i]->gWeight() << endl;

   }

}
