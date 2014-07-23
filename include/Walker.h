#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

using namespace btas;

#include "PEPS.h"

/**
 * class definition of Walker, made to describe the product state walkers. An array of L*L booleans vector. Each site represents a spin up (0) or spin down (1)
 */
class Walker : public vector< bool > {

   public:

      //empty contstructor
      Walker();
   
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();

      double gWeight() const;

      void sWeight(double);

      void multWeight(double);

      double gOverlap() const;

      double gEL() const;

  private:

      //!The walker weight
      double weight;

      //!The walker overlap with the trial wfn
      double overlap;

      //!local energy
      double EL;

};

#endif
