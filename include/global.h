#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

using namespace btas;

#include "Walker.h"
#include "PEPS.h"

namespace global {

   extern Random RN;

   //!x dimension of the lattice, nr of cols
   extern int Lx;

   //!y dimension of the lattice, nr of rows
   extern int Ly;

   //!physical dimension of sites
   extern int d;

   //!virtual dimension of the trial
   extern int DT;

   //!backup walker state for stability in GFQMC algorithm
   extern vector<Walker> backup_walker;

   //!number of omp threads
   extern int omp_num_threads;

   //!trial wavefunction
   extern PEPS<double> peps;

   void init(int,int,int,int);

   template<typename T>
      T rgen();

   template<typename T>
      T rgen_pos();
   
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
