#ifndef GFQMC_H
#define GFQMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

class GFQMC {

   public:
   
      //constructor with input trialwavefunction
      GFQMC(double,int);
      
      //Destructor
      virtual ~GFQMC();
      
      //Let the walkers propagate for n_steps steps
      void walk(int);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      double propagate();
      
      //Control the population of walkers based on scaling * weight
      void PopulationControl(double);

      //Calculate the single walker projected energies, update the energy history, calculate the fluctuation metric, and the total projected energy
      void sEP();

      //Write the projected energy, target energy
      void write(int,double);

      //Setup the walkers
      void SetupWalkers();

   private:
      
      //!The total desired number of walkers
      int Nw;

      //!timestep
      double dtau;

      //!projected energy at current timestep
      double EP;
      
      //!The walkers
      std::vector<Walker> walker;

      //!Distribution of possible final states given a walker state
      Distribution dist;
      
};

#endif
