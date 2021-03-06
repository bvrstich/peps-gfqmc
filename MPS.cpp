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
 * empty constructor
 */
MPS::MPS() : vector< TArray<double,3> >( ) { }

/** 
 * standard constructor: just takes in, sets the lenght of the tensors
 * @param L_in length of the chain
 */
MPS::MPS(int L_in) : vector< TArray<double,3> >( L_in ) { }

/** 
 * standard constructor: just takes in, sets the lenght of the tensors
 * @param L_in length of the chain
 * @param d_in physical dimension
 * @param D_in virtual dimension
 */
MPS::MPS(int L_in,int d_in,int D_in) : vector< TArray<double,3> >( L_in ) {

   D = D_in;
   d = d_in;

   vector<int> vdim(L_in + 1);

   vdim[0] = 1;

   for(int i = 1;i < L_in;++i){

      int tmp = vdim[i - 1] * d;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L_in] = 1;

   for(int i = L_in - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i)
      (*this)[i].resize(vdim[i],d,vdim[i+1]);

}

/**
 * copy constructor
 */
MPS::MPS(const MPS &mps_copy) : vector< TArray<double,3> >(mps_copy) {

   this->D = mps_copy.gD();
   this->d = mps_copy.gd();

}

/**
 * empty destructor
 */
MPS::~MPS(){ }

/**
 * fill with random entries
 */
void MPS::fill_Random(){

   DArray<2> L;

   for(int i = this->size() - 1;i >= 0;--i){

      (*this)[i].generate(rgen<double>);
      Gelqf(L,(*this)[i]);

   }

}

/**
 * @return virtual dimension of the MPS
 */
int MPS::gD() const {

   return D;

}

/**
 * @return physical dimension of the MPS
 */
int MPS::gd() const {

   return d;

}

/**
 * canonicalize the mps
 * @param dir Left or Right canonicalization
 * @param norm if true: normalize, else not
 */
void MPS::canonicalize(const BTAS_SIDE &dir,bool norm){

   int length = this->size();

   if(dir == Left){//QR

      TArray<double,2> R;
      TArray<double,3> tmp;

      for(int i = 0;i < length - 1;++i){

         R.clear();

         //do QR
         Geqrf((*this)[i],R);

         //paste to next matrix
         tmp.clear();

         Contract(1.0,R,shape(1),(*this)[i + 1],shape(0),0.0,tmp);

         (*this)[i + 1] = std::move(tmp);

      }

      if(norm){

         double nrm = sqrt(Dotc((*this)[length-1],(*this)[length-1]));
         Scal(1.0/nrm,(*this)[length-1]);

      }

   }
   else{//LQ

      TArray<double,2> L;
      TArray<double,3> tmp;

      for(int i = length - 1;i > 0;--i){

         L.clear();

         //do QR
         Gelqf(L,(*this)[i]);

         //paste to previous matrix
         tmp.clear();

         Contract(1.0,(*this)[i - 1],shape(2),L,shape(0),0.0,tmp);

         (*this)[i - 1] = std::move(tmp);

      }

      if(norm){

         double nrm = sqrt(Dotc((*this)[0],(*this)[0]));
         Scal(1.0/nrm,(*this)[0]);

      }

   }

}

/**
 * @param bra the bra of the inner product
 * @return the inner product of two MPS's, with *this being the ket
 */
double MPS::dot(const MPS &bra) const {

   TArray<double,2> E;

   Contract(1.0,bra[0],shape(0,1),(*this)[0],shape(0,1),0.0,E);

   TArray<double,3> I;

   for(int i = 1;i < this->size();++i){

      I.clear();

      Contract(1.0,bra[i],shape(0),E,shape(0),0.0,I);

      E.clear();

      Contract(1.0,I,shape(2,0),(*this)[i],shape(0,1),0.0,E);

   }

   return E(0,0);

}

/**
 * normalize the MPS
 */
void MPS::normalize(){

   double nrm2 = this->dot(*this);

   nrm2 = pow(nrm2,(double)1.0/(double)(2.0*this->size()));

   for(int i = 0;i < this->size();++i)
      Scal(1.0/nrm2,(*this)[i]);

}

/**
 * fill a MPS object, by creating a single layer from contracting a peps with a physical vector, i.e. a walker
 * @param option 't'op, 'b'ottom, 'l'eft or 'r'ight
 * @param inverse if true use inverse walker
 * @param walker the input Walker object
 */
void MPS::fill(char option,bool inverse,const Walker &walker) {

   if(option == 'b'){

      //just a shallow copy of bottom row
      for(int col = 0;col < Lx;++col){

         bool s = inverse ^ walker[col];
         int dim = (*this)[col].size();

         blas::copy(dim, peps[0](0,col,s).data(), 1, (*this)[col].data(), 1);

      }

   }
   else if(option == 't'){

      //just a shallow copy of the top row
      for(int col = 0;col < Lx;++col){

         bool s = inverse ^ walker[(Ly-1)*Lx + col];
         int dim = (*this)[col].size();

         blas::copy(dim, peps[0](Ly-1,col,s).data(), 1, (*this)[col].data(), 1);

      }

   }
   else if(option == 'l'){

      //just a shallow copy of the left col
      for(int row = 0;row < Ly;++row){

         bool s = inverse ^ walker[row*Lx];
         int dim = (*this)[row].size();

         blas::copy(dim, peps[1](row,0,s).data(), 1, (*this)[row].data(), 1);

      }

   }
   else{//right

      //just a shallow copy of the right col
      for(int row = 0;row < Ly;++row){

         bool s = inverse ^ walker[row*Lx + Lx - 1];
         int dim = (*this)[row].size();

         blas::copy(dim, peps[1](row,Lx - 1,s).data(), 1, (*this)[row].data(), 1);

      }

   }

}

ostream &operator<<(ostream &output,const MPS &mps_p){

   for(int i = 0;i < mps_p.size();++i){

      output << endl;
      output << "Tensor on site\t" << i << endl;
      output << endl;
      output << i << "\t" << mps_p[i] << endl;

   }

   return output;

}

/**
 * scale the MPS with a constant factor
 * @param alpha scalingfactor
 */
void MPS::scal(double alpha){

   int sign;

   if(alpha > 0)
      sign = 1;
   else
      sign = -1;

   alpha = pow(fabs(alpha),1.0/(double)this->size());

   Scal(sign * alpha,(*this)[0]);

   for(int i = 1;i < this->size();++i)
      Scal(alpha,(*this)[i]);

}
