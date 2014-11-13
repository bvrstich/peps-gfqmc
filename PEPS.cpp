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
 * construct an empty PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 */
template<typename T>
PEPS<T>::PEPS() : vector< TArray<T,4> >(2 * Lx * Ly) { }

/**
 * construct from file
 * @param filename name of the file
 */
template<typename T>
PEPS<T>::PEPS(const char *filename) : vector< TArray<T,4> >(2 * Lx * Ly){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ifstream fin(name);

         int Da,Db,Dc,Dd,De;

         fin >> Da >> Db >> Dc >> Dd >> De;

         for(int s = 0;s < d;++s)
            (*this)[s*Lx*Ly + row*Lx + col].resize(Da,Db,Dd,De);

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fin >> a_ >> b_ >> c_ >> d_ >> e_ >> (*this)[c_*Lx*Ly + row*Lx + col](a_,b_,d_,e_);

      }

   this->D = (*this)[0].shape(3);

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,4> >(peps_copy) {

   D = peps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
PEPS<T>::~PEPS(){ }

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @param s spin index
 * @return the tensor on site (r,c) with spin s
 */
template<typename T>
const TArray<T,4> &PEPS<T>::operator()(int r,int c,int s) const {

   return (*this)[s*Lx*Ly + r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,4> &PEPS<T>::operator()(int r,int c,int s) {

   return (*this)[s*Lx*Ly + r*Lx + c];

}

/**
 * @return the cutoff virtual dimension
 */
template<typename T>
int PEPS<T>::gD() const {

   return D;

}

/**
 * @param D_in value to the D to
 */
template<typename T>
void PEPS<T>::sD(int D_in) {

   this->D = D_in;

}

/**
 * scale the peps with a number
 * @param val scalar to be multiplied with the peps
 */
template<typename T>
void PEPS<T>::scal(T val){

   val = pow(val,(T)1.0/(T)this->size());

   //now scale every tensor in the array
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         for(int s = 0;s < d;++s)
            Scal(val,(*this)[ s*Lx*Ly + r*Lx + c ]);

}

//forward declarations for types to be used!
template PEPS<double>::PEPS();
template PEPS< complex<double> >::PEPS();

template PEPS<double>::PEPS(const char *);
template PEPS< complex<double> >::PEPS(const char *);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,4> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,4> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,4> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,4> &PEPS< complex<double> >::operator()(int r,int c) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::scal(double val);
template void PEPS< complex<double> >::scal(complex<double> val);
