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

/** 
 * empty constructor: just sets the length of the vector
 */
SL_PEPS::SL_PEPS() : vector< TArray<double,4> >( Lx*Ly ) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
SL_PEPS::SL_PEPS(int D_in) : vector< TArray<double,4> >( Lx*Ly ) {

   D = D_in;

   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(1,D,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx - 1 ].resize(D,D,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(1,1,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx - 1 ].resize(D,1,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(D,D,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(D,1,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(1,D,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(D,D,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(D,D,D,D);

}

/**
 * copy constructor
 */
SL_PEPS::SL_PEPS(const SL_PEPS &slp_copy) : vector< TArray<double,4> >(slp_copy) {

   this->D = slp_copy.gD();

}

/**
 * empty destructor
 */
SL_PEPS::~SL_PEPS(){ }

/**
 * @return virtual dimension of the SL_PEPS
 */
int SL_PEPS::gD() const {

   return D;

}

/**
 * fill the SL_PEPS object by contracting a peps with a walker
 * @param inverse if true take the inverse walker
 * @param peps input PEPS<> object
 * @param walker the Walker object
 */
void SL_PEPS::fill(bool inverse,const PEPS< double > &peps,const Walker &walker){

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         int dim = peps(r,c).size()/2;

         if(inverse)
            blas::copy(dim,peps(r,c).data() + (!walker[r*Lx + c])*dim,1,(*this)[r*Lx + c].data(),1);
         else
            blas::copy(dim,peps(r,c).data() + walker[r*Lx + c]*dim,1,(*this)[r*Lx + c].data(),1);

      }

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
const TArray<double,4> &SL_PEPS::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
TArray<double,4> &SL_PEPS::operator()(int r,int c) {

   return (*this)[r*Lx + c];

}
