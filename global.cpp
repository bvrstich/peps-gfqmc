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
using std::ifstream;

#include "include.h"

namespace global{

   int DT;

   int D_aux;

   int Lx;
   int Ly;

   int d;

   std::vector< PEPS<double> > peps;

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

      int mult = D_aux/DT;

      char filename_in[200];
      sprintf(filename_in,"/home/bright/bestanden/results/peps/%dx%d/D=%d/D_aux=%d/peps",Lx,Ly,DT,mult*DT*DT);

      ifstream in(filename_in);

      peps.resize(4);

      peps[0] = PEPS<double>(filename_in);

      //now construct permuted peps'
      peps[1] = peps[0];
      peps[2] = peps[0];
      peps[3] = peps[0];

      for(int r = 0;r < Ly;++r)
         for(int c = 0;c < Lx;++c)
            for(int s = 0;s < d;++s){

               Permute(peps[0](r,c,s),shape(2,0,3,1),peps[1](r,c,s));
               Permute(peps[0](r,c,s),shape(0,2,1,3),peps[2](r,c,s));
               Permute(peps[0](r,c,s),shape(3,1,2,0),peps[3](r,c,s));

            }

      Walker walker;

      double tmp = walker.overlap();

      peps[0].scal(1.0/tmp);
      peps[1].scal(1.0/tmp);
      peps[2].scal(1.0/tmp);
      peps[3].scal(1.0/tmp);

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
