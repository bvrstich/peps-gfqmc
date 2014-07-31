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

   friend ostream &operator<<(ostream &output,const Walker &walker_p);

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

      double exp_en(const Walker &);

      void calc_EL(const PEPS<double> &);

      void save(const char *filename);

      int gsign() const;

      void sign_flip();

  private:

      //!The walker weight
      double weight;

      //!sign of the walker
      int sign;

      //!The walker overlap with the trial wfn
      double overlap;

      //!local energy
      double EL;

};

#endif
