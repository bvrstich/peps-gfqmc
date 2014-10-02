#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

namespace global{

   int DT;

   int D_aux;

   int Lx;
   int Ly;

   int d;

   PEPS<double> peps;

   int omp_num_threads;

   Random RN;

   /**
    * @param DT_in virtual dimension of the trial
    * @param D_aux_in auxiliary dimension for peps contraction
    * @param d_in physical dimension
    * @param Lx_in x dimension of the square lattice
    * @param Ly_in y dimension of the square lattice
    */
   void init(int DT_in,int D_aux_in,int d_in,int Lx_in,int Ly_in){

      Lx = Lx_in;
      Ly = Ly_in;

      d = d_in;

      DT = DT_in;
      D_aux = D_aux_in;

#ifdef _OPENMP
      omp_num_threads = omp_get_max_threads();
#else
      omp_num_threads = 1;
#endif

      peps.resize(Lx*Ly);
      peps.set_jastrow(0.74);

   }

   //!function which generates random complex numbers uniformly on a square of side 2 [(-1,1):(-1,1)]
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

   //!function which generates uniform random numbers between [0:1]
   template<>
      double rgen_pos(){ 

         return RN();

      }

}
