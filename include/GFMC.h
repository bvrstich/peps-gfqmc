#ifndef GFMC_H
#define GFMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

class Walker;
class Distribution;

class GFMC {

   public:
   
      //constructor with input trialwavefunction
      GFMC(double,int);
      
      //Destructor
      virtual ~GFMC();
      
      //Let the walkers propagate for n_steps steps
      void walk(int);

      void test();

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double propagate();
      
      //Control the population of walkers based on scaling * weight
      void PopulationControl(double);

      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      void sEP();

      //Write the projected energy, target energy
      void write(int);

      //Setup the walkers
      void SetupWalkers();

      void dump(const char *);

   private:
      
      //!The total desired number of walkers
      int Nw;

      //!timestep
      double dtau;

      //!projected energy at current timestep
      double EP;

      //!target energy, growth estimate
      double ET;

      //!number of walkers that stay the same in a given generation
      int num_stable;

      //!The walkers
      std::vector<Walker> walker;

      //!Distribution of possible final states given a walker state
      std::vector< Distribution > dist;
      
};

#endif
