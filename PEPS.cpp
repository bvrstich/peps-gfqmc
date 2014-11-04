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
PEPS<T>::PEPS() : vector< TArray<T,5> >(Lx * Ly) { }

/**
 * construct constructs a standard PEPS object, note: be sure to initialize the Lattice object before calling the constructor
 * @param D_in cutoff virtual dimension
 */
template<typename T>
PEPS<T>::PEPS(int D_in) : vector< TArray<T,5> >(Lx * Ly) {

   D = D_in;
   
   //corners first

   //r == 0 : c == 0
   (*this)[ 0 ].resize(1,D,d,1,D);

   //r == 0 : c == L - 1
   (*this)[ Lx-1 ].resize(D,D,d,1,1);

   //r == L - 1 : c == 0
   (*this)[ (Ly-1)*Lx ].resize(1,1,d,D,D);

   //r == L - 1 : c == L - 1
   (*this)[ (Ly-1)*Lx + Lx-1 ].resize(D,1,d,D,1);

   //sides:

   //r == 0
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ c ].resize(D,D,d,1,D);

   //r == Ly - 1
   for(int c = 1;c < Lx - 1;++c)
      (*this)[ (Ly-1)*Lx + c ].resize(D,1,d,D,D);

   //c == 0
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx ].resize(1,D,d,D,D);

   //c == Lx - 1
   for(int r = 1;r < Ly - 1;++r)
      (*this)[ r*Lx + Lx - 1 ].resize(D,D,d,D,1);

   //the rest is full
   for(int r = 1;r < Ly - 1;++r)
      for(int c = 1;c < Lx - 1;++c)
         (*this)[ r*Lx + c ].resize(D,D,d,D,D);

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         (*this)[ r*Lx + c ].generate(rgen<T>);

         Normalize((*this)[ r*Lx + c ]);
         Scal((T)D,(*this)[ r*Lx + c ]);

      }

}

/**
 * copy constructor
 */
template<typename T>
PEPS<T>::PEPS(const PEPS<T> &peps_copy) : vector< TArray<T,5> >(peps_copy) {

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
 * @return the tensor on site (r,c)
 */
template<typename T>
const TArray<T,5> &PEPS<T>::operator()(int r,int c) const {

   return (*this)[r*Lx + c];

}

/**
 * access to the individual tensors: const version
 * @param r row index
 * @param c col index
 * @return the tensor on site (r,c)
 */
template<typename T>
TArray<T,5> &PEPS<T>::operator()(int r,int c) {

   return (*this)[r*Lx + c];

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
 * initialize the peps to the direct sum of two antiferromagnetic D=1 structures
 * @param f jastrow factor
 */
template<>
void PEPS<double>::initialize_jastrow(double f) {

   D = 2;

   //bottom row, first site
   (*this)[0].resize(d,1,D,1,D);
   (*this)[0] = 0.0;

   (*this)[0](0,0,0,0,0) = 1.0;
   (*this)[0](1,0,1,0,1) = 1.0;

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[col].resize(d,D,D,1,D);
      (*this)[col] = 0.0;

      (*this)[col](0,0,0,0,0) = f;
      (*this)[col](0,1,0,0,0) = 1.0;
      (*this)[col](1,0,1,0,1) = 1.0;
      (*this)[col](1,1,1,0,1) = f;

   }

   //bottom row, last site
   (*this)[Lx-1].resize(d,D,D,1,1);
   (*this)[Lx-1] = 0.0;

   (*this)[Lx-1](0,0,0,0,0) = f;
   (*this)[Lx-1](0,1,0,0,0) = 1.0;

   (*this)[Lx-1](1,0,1,0,0) = 1.0;
   (*this)[Lx-1](1,1,1,0,0) = f;

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      (*this)[row*Lx].resize(d,1,D,D,D);
      (*this)[row*Lx] = 0.0;

      (*this)[row*Lx](0,0,0,0,0) = f;
      (*this)[row*Lx](0,0,0,1,0) = 1.0;
      (*this)[row*Lx](1,0,1,0,1) = 1.0;
      (*this)[row*Lx](1,0,1,1,1) = f;

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         (*this)[row*Lx + col].resize(d,D,D,D,D);
         (*this)[row*Lx + col] = 0.0;

         (*this)[row*Lx + col](0,0,0,0,0) = f*f;
         (*this)[row*Lx + col](0,0,0,1,0) = f;
         (*this)[row*Lx + col](0,1,0,0,0) = f;
         (*this)[row*Lx + col](0,1,0,1,0) = 1.0;

         (*this)[row*Lx + col](1,0,1,0,1) = 1.0;
         (*this)[row*Lx + col](1,0,1,1,1) = f;
         (*this)[row*Lx + col](1,1,1,0,1) = f;
         (*this)[row*Lx + col](1,1,1,1,1) = f*f;

      }

      //rightmost site on row 'row'
      (*this)[row*Lx + Lx - 1].resize(d,D,D,D,1);
      (*this)[row*Lx + Lx - 1] = 0.0;

      (*this)[row*Lx + Lx - 1](0,0,0,0,0) = f*f;
      (*this)[row*Lx + Lx - 1](0,0,0,1,0) = f;
      (*this)[row*Lx + Lx - 1](0,1,0,0,0) = f;
      (*this)[row*Lx + Lx - 1](0,1,0,1,0) = 1.0;

      (*this)[row*Lx + Lx - 1](1,0,1,0,0) = 1.0;
      (*this)[row*Lx + Lx - 1](1,1,1,0,0) = f;
      (*this)[row*Lx + Lx - 1](1,0,1,1,0) = f;
      (*this)[row*Lx + Lx - 1](1,1,1,1,0) = f*f;

   }

   //top row
   //leftmost site
   (*this)[(Ly - 1)*Lx].resize(d,1,1,D,D);
   (*this)[(Ly - 1)*Lx] = 0.0;

   (*this)[(Ly - 1)*Lx](0,0,0,0,0) = f;
   (*this)[(Ly - 1)*Lx](0,0,0,1,0) = 1.0;
   (*this)[(Ly - 1)*Lx](1,0,0,0,1) = 1.0;
   (*this)[(Ly - 1)*Lx](1,0,0,1,1) = f;

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      (*this)[(Ly - 1)*Lx + col].resize(d,D,1,D,D);
      (*this)[(Ly - 1)*Lx + col] = 0.0;

      (*this)[(Ly - 1)*Lx + col](0,0,0,0,0) = f*f;
      (*this)[(Ly - 1)*Lx + col](0,0,0,1,0) = f;
      (*this)[(Ly - 1)*Lx + col](0,1,0,0,0) = f;
      (*this)[(Ly - 1)*Lx + col](0,1,0,1,0) = 1.0;

      (*this)[(Ly - 1)*Lx + col](1,0,0,0,1) = 1.0;
      (*this)[(Ly - 1)*Lx + col](1,0,0,1,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,1,0,0,1) = f;
      (*this)[(Ly - 1)*Lx + col](1,1,0,1,1) = f*f;

   }

   //top row rightmost site
   (*this)[(Ly - 1)*Lx + Lx - 1].resize(d,D,1,D,1);
   (*this)[(Ly - 1)*Lx + Lx - 1] = 0.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,0,0) = f*f;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,0,0,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,1,0,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](0,1,0,1,0) = 1.0;

   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,0,0) = 1.0;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,0,0,1,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,1,0,0,0) = f;
   (*this)[(Ly - 1)*Lx + Lx - 1](1,1,0,1,0) = f*f;

}

/**
 * increase the bond dimension by one
 * @param D_in bond dimension to grow to
 * @param noise level of noise added to the initial state
 */
template<>
void PEPS<double>::grow_bond_dimension(int D_in,double noise) {

   D = D_in;

   DArray<5> tmp;

   //bottom row, first site
   tmp.resize(1,D,d,1,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[0].shape(0);++i)
      for(int j = 0;j < (*this)[0].shape(1);++j)
         for(int k = 0;k < (*this)[0].shape(2);++k)
            for(int l = 0;l < (*this)[0].shape(3);++l)
               for(int m = 0;m < (*this)[0].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[0](i,j,k,l,m);

   (*this)[0] = std::move(tmp);

   //bottom row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,D,d,1,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[col].shape(0);++i)
         for(int j = 0;j < (*this)[col].shape(1);++j)
            for(int k = 0;k < (*this)[col].shape(2);++k)
               for(int l = 0;l < (*this)[col].shape(3);++l)
                  for(int m = 0;m < (*this)[col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[col](i,j,k,l,m);

      (*this)[col] = std::move(tmp);

   }

   //bottom row, last site
   tmp.clear();
   tmp.resize(D,D,d,1,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[Lx-1].shape(0);++i)
      for(int j = 0;j < (*this)[Lx-1].shape(1);++j)
         for(int k = 0;k < (*this)[Lx-1].shape(2);++k)
            for(int l = 0;l < (*this)[Lx-1].shape(3);++l)
               for(int m = 0;m < (*this)[Lx-1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[Lx - 1](i,j,k,l,m);

   (*this)[Lx-1] = std::move(tmp);

   //middle sites
   for(int row = 1;row < Ly - 1;++row){

      //leftmost middle site: col == 0
      tmp.clear();
      tmp.resize(1,D,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx](i,j,k,l,m);

      (*this)[row*Lx] = std::move(tmp);

      //middle sites on row 'row'
      for(int col = 1;col < Lx - 1;++col){

         tmp.clear();
         tmp.resize(D,D,d,D,D);

         tmp.generate(rgen<double>);
         Scal(noise,tmp);

         for(int i = 0;i < (*this)[row*Lx + col].shape(0);++i)
            for(int j = 0;j < (*this)[row*Lx + col].shape(1);++j)
               for(int k = 0;k < (*this)[row*Lx + col].shape(2);++k)
                  for(int l = 0;l < (*this)[row*Lx + col].shape(3);++l)
                     for(int m = 0;m < (*this)[row*Lx + col].shape(4);++m)
                        tmp(i,j,k,l,m) += (*this)[row*Lx + col](i,j,k,l,m);

         (*this)[row*Lx + col] = std::move(tmp);

      }

      //rightmost site on row 'row'
      tmp.clear();
      tmp.resize(D,D,d,D,1);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[row*Lx + Lx - 1].shape(0);++i)
         for(int j = 0;j < (*this)[row*Lx + Lx - 1].shape(1);++j)
            for(int k = 0;k < (*this)[row*Lx + Lx - 1].shape(2);++k)
               for(int l = 0;l < (*this)[row*Lx + Lx - 1].shape(3);++l)
                  for(int m = 0;m < (*this)[row*Lx + Lx - 1].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[row*Lx + Lx - 1](i,j,k,l,m);

      (*this)[row*Lx + Lx - 1] = std::move(tmp);

   }

   //top row
   //leftmost site
   tmp.clear();
   tmp.resize(1,1,d,D,D);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx] = std::move(tmp);

   //top row, middle sites
   for(int col = 1;col < Lx - 1;++col){

      tmp.clear();
      tmp.resize(D,1,d,D,D);

      tmp.generate(rgen<double>);
      Scal(noise,tmp);

      for(int i = 0;i < (*this)[(Ly - 1)*Lx + col].shape(0);++i)
         for(int j = 0;j < (*this)[(Ly - 1)*Lx + col].shape(1);++j)
            for(int k = 0;k < (*this)[(Ly - 1)*Lx + col].shape(2);++k)
               for(int l = 0;l < (*this)[(Ly - 1)*Lx + col].shape(3);++l)
                  for(int m = 0;m < (*this)[(Ly - 1)*Lx + col].shape(4);++m)
                     tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + col](i,j,k,l,m);

      (*this)[(Ly - 1)*Lx + col] = std::move(tmp);

   }

   //top row rightmost site
   tmp.clear();
   tmp.resize(D,1,d,D,1);

   tmp.generate(rgen<double>);
   Scal(noise,tmp);

   for(int i = 0;i < (*this)[(Ly - 1)*Lx + Lx - 1].shape(0);++i)
      for(int j = 0;j < (*this)[(Ly - 1)*Lx + Lx - 1].shape(1);++j)
         for(int k = 0;k < (*this)[(Ly - 1)*Lx + Lx - 1].shape(2);++k)
            for(int l = 0;l < (*this)[(Ly - 1)*Lx + Lx - 1].shape(3);++l)
               for(int m = 0;m < (*this)[(Ly - 1)*Lx + Lx - 1].shape(4);++m)
                  tmp(i,j,k,l,m) += (*this)[(Ly - 1)*Lx + Lx - 1](i,j,k,l,m);

   (*this)[(Ly - 1)*Lx + Lx - 1] = std::move(tmp);

}

/**
 * scale the peps with a number
 * @param val scalar to be multiplied with the peps
 */
template<typename T>
void PEPS<T>::scal(T val){

   val = pow(val,(T)1.0/(T)this->size());

   //now initialize with random numbers
   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c)
         Scal(val,(*this)[ r*Lx + c ]);

}

/**
 * @param mpx will be written to file
 * @param filename name of the file
 * save the MPX object to a file in binary format.
 */

template<typename T>
void PEPS<T>::save(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ofstream fout(name);
         fout.precision(16);

         int Da = (*this)(row,col).shape(0);
         int Db = (*this)(row,col).shape(1);
         int Dc = (*this)(row,col).shape(2);
         int Dd = (*this)(row,col).shape(3);
         int De = (*this)(row,col).shape(4);

         fout << Da << "\t" << Db << "\t" << Dc << "\t" << Dd << "\t" << De << endl;

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fout << a_ << "\t" << b_ << "\t" << c_ << "\t" << d_ << "\t" << e_ << "\t" << (*this)(row,col)(a_,b_,c_,d_,e_) << endl;

      }

}

/**
 * @param mpx will be constructed from file
 * @param filename name of the file
 * load the MPX object from a file in binary format.
 */
template<typename T>
void PEPS<T>::load(const char *filename){

   for(int row = 0;row < Ly;++row)
      for(int col = 0;col < Lx;++col){

         char name[200];

         sprintf(name,"%s/site_(%d,%d).peps",filename,row,col);

         std::ifstream fin(name);

         int Da,Db,Dc,Dd,De;

         fin >> Da >> Db >> Dc >> Dd >> De;

         (*this)(row,col).resize(Da,Db,Dc,Dd,De);

         for(int a_ = 0;a_ < Da;++a_)
            for(int b_ = 0;b_ < Db;++b_)
               for(int c_ = 0;c_ < Dc;++c_)
                  for(int d_ = 0;d_ < Dd;++d_)
                     for(int e_ = 0;e_ < De;++e_)
                        fin >> a_ >> b_ >> c_ >> d_ >> e_ >> (*this)(row,col)(a_,b_,c_,d_,e_);

      }

}

//forward declarations for types to be used!
template PEPS<double>::PEPS();
template PEPS< complex<double> >::PEPS();

template PEPS<double>::PEPS(int);
template PEPS< complex<double> >::PEPS(int);

template PEPS<double>::PEPS(const PEPS<double> &);
template PEPS< complex<double> >::PEPS(const PEPS< complex<double> > &);

template PEPS<double>::~PEPS();
template PEPS< complex<double> >::~PEPS();

template TArray<double,5> &PEPS<double>::operator()(int r,int c);
template TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c);

template const TArray<double,5> &PEPS<double>::operator()(int r,int c) const;
template const TArray<complex<double>,5> &PEPS< complex<double> >::operator()(int r,int c) const;

template int PEPS<double>::gD() const;
template int PEPS< complex<double> >::gD() const;

template void PEPS<double>::sD(int);
template void PEPS< complex<double> >::sD(int);

template void PEPS<double>::scal(double val);
template void PEPS< complex<double> >::scal(complex<double> val);

template void PEPS<double>::load(const char *filename);
template void PEPS< complex<double> >::load(const char *filename);

template void PEPS<double>::save(const char *filename);
template void PEPS< complex<double> >::save(const char *filename);
