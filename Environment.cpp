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

using namespace global;

//statics
vector< vector< MPS > > Environment::l;
vector< vector< MPS > > Environment::r;
vector< vector< MPS > > Environment::t;
vector< vector< MPS > > Environment::b;

vector< SL_PEPS > Environment::U;
vector< SL_PEPS > Environment::I;

/** 
 * initialize all the static variables
 * @param D bond dimension of the trial peps
 * @param D_aux_in auxiliary bond dimension for the contractions
 */
void Environment::init(){

   t.resize(omp_num_threads);
   b.resize(omp_num_threads);

   for(int thr = 0;thr < omp_num_threads;++thr){

      t[thr].resize(Ly - 1);
      b[thr].resize(Ly - 1);

      b[thr][0] = MPS(DT);
      t[thr][Ly-2] = MPS(DT);

      int dim = DT;

      for(int i = 1;i < Ly - 2;++i){

         dim *= DT;

         if(dim < D_aux){

            b[thr][i] = MPS(dim);
            t[thr][Ly-2-i] = MPS(dim);

         }
         else{

            b[thr][i] = MPS(D_aux);
            t[thr][Ly-2-i] = MPS(D_aux);

         }

      }

      b[thr][Ly-2] = MPS(D_aux);
      t[thr][0] = MPS(D_aux);

   }

   r.resize(omp_num_threads);
   l.resize(omp_num_threads);

   for(int thr = 0;thr < omp_num_threads;++thr){

      r[thr].resize(Lx - 1);
      l[thr].resize(Lx - 1);

      l[thr][0] = MPS(DT);
      r[thr][Lx-2] = MPS(DT);

      int dim = DT;

      for(int i = 1;i < Lx - 2;++i){

         dim *= DT;

         if(dim < D_aux){

            l[thr][i] = MPS(dim);
            r[thr][Lx-2-i] = MPS(dim);

         }
         else{

            l[thr][i] = MPS(D_aux);
            r[thr][Lx-2-i] = MPS(D_aux);

         }

      }

      l[thr][Lx-2] = MPS(D_aux);
      r[thr][0] = MPS(D_aux);

   }

   U.resize(omp_num_threads);
   I.resize(omp_num_threads);

   for(int thr = 0;thr < omp_num_threads;++thr){

      U[thr] = SL_PEPS(DT);
      I[thr] = SL_PEPS(DT);

   }

}

/**
 * construct the environment mps's for the input PEPS:
 * @param option if 'H' construct top and bottom environment
 *               if 'V' construct left and right environment
 * @param peps input PEPS<double>
 */
void Environment::calc_env(char option,const PEPS< double > &peps,const Walker &walker){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   if(option == 'H'){

      //construct bottom layer
      b[myID][0].fill('b',U[myID]);

      for(int r = 1;r < Ly - 1;++r){

         MPS tmp(b[myID][r - 1]);

         //apply to form MPS with bond dimension D^2
         tmp.gemv('L','H',r,U[myID]);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         b[myID][r].compress(D_aux,tmp,1);

     }

      //then construct top layer
      t[myID][Ly - 2].fill('t',U[myID]);

      for(int r = Ly - 2;r > 0;--r){

         //apply to form MPS with bond dimension D^4
         MPS tmp(t[myID][r]);

         tmp.gemv('U','H',r,U[myID]);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         t[myID][r - 1].compress(D_aux,tmp,1);

      }

   }
   else{//Vertical

      //then left layer
      l[myID][0].fill('l',U[myID]);

      for(int c = 1;c < Lx - 1;++c){

         //i'th col as MPO
         MPS tmp(l[myID][c - 1]);

         //apply to form MPS with bond dimension D^4
         tmp.gemv('L','V',c,U[myID]);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         l[myID][c].compress(D_aux,tmp,1);

      }

      //finally construct right layer
      r[myID][Lx - 2].fill('r',U[myID]);

      for(int c = Lx - 2;c > 0;--c){

         //apply to form MPS with bond dimension D^4
         MPS tmp(r[myID][c]);

         tmp.gemv('U','V',c,U[myID]);

         //reduce the dimensions of the edge states using thin svd
         tmp.cut_edges();

         //compress in sweeping fashion
         r[myID][c - 1].compress(D_aux,tmp,1);

      }

   }

}

/**
 * test if the enviroment is correctly contracted
 */
void Environment::test_env(char option){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   if(option == 'H'){

      cout << endl;
      cout << "FROM BOTTOM TO TOP" << endl;
      cout << endl;
      for(int i = 0;i < Ly - 1;++i)
         cout << i << "\t" << b[myID][i].dot(t[myID][i]) << endl;
   }
   else{

      cout << endl;
      cout << "FROM LEFT TO RIGHT" << endl;
      cout << endl;
      for(int i = 0;i < Lx - 1;++i)
         cout << i << "\t" << r[myID][i].dot(l[myID][i]) << endl;
      cout << endl;

   }

}
