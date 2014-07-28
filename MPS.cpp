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
MPS::MPS() : vector< TArray<double,3> >( Lx ) { }

/** 
 * standard constructor: just takes in
 * @param L_in length of the chain
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
MPS::MPS(int D_in) : vector< TArray<double,3> >( Lx ) {

   this->D = D_in;

   (*this)[0].resize(1,D,D);

   for(int c = 1;c < Lx-1;++c)
      (*this)[c].resize(D,D,D);

   (*this)[Lx-1].resize(D,D,1);

}

/**
 * copy constructor
 */
MPS::MPS(const MPS &mps_copy) : vector< TArray<double,3> >(mps_copy) {

   this->D = mps_copy.gD();

}

/**
 * empty destructor
 */
MPS::~MPS(){ }

/**
 * fill with random entries
 */
void MPS::fill_Random(){

   for(int i = 0;i < this->size();++i)
      (*this)[i].generate(rgen< double >);

}

/**
 * @return virtual dimension of the MPS
 */
int MPS::gD() const {

   return D;

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
 * find an approximate form of the state 'mps' compressed to a bond dimension 'Dc' by performing an SVD on an non-canonical state.
 * @param dir Left or Right - going compression
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
void MPS::guess(const BTAS_SIDE &dir,int Dc,const MPS &mps){

   int L = mps.size();

   if(dir == Left){

      TArray<double,3> U;
      TArray<double,2> V;
      TArray<double,1> S;

      Gesvd('S','S',mps[0],S,U,V,Dc);

      (*this)[0] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      //and contract V with mps on next site
      (*this)[1].clear();

      Contract(1.0,V,shape(1),mps[1],shape(0),0.0,(*this)[1]);

      for(int i = 1;i < L - 1;++i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(U);

         //multiply S to V
         Dimm(S,V);

         //and contract V with mps on next site
         (*this)[i + 1].clear();

         Contract(1.0,V,shape(1),mps[i + 1],shape(0),0.0,(*this)[i + 1]);

      }

   }
   else{

      TArray<double,2> U;
      TArray<double,3> V;
      TArray<double,1> S;

      Gesvd('S','S',mps[L - 1],S,U,V,Dc);

      (*this)[L - 1] = std::move(V);

      //multiply U and S
      Dimm(U,S);

      //and contract U with mps on previous site
      (*this)[L - 2].clear();

      Contract(1.0,mps[L - 2],shape(2),U,shape(0),0.0,(*this)[L - 2]);

      for(int i = L - 2;i > 0;--i){

         Gesvd('S','S',(*this)[i],S,U,V,Dc);

         (*this)[i] = std::move(V);

         //multiply S to V
         Dimm(U,S);

         //and contract V with mps on next site
         (*this)[i - 1].clear();

         Contract(1.0,mps[i - 1],shape(2),U,shape(0),0.0,(*this)[i - 1]);

      }

   }

   if(Dc < mps.gD())
      this->D = Dc;
   else
      Dc = mps.gD();

}

/**
 * find the best compression of the state 'mps' a bond dimension 'Dc' by optimizing the tensor in a sweeping fashion
 * @param Dc the compressed dimension
 * @param mps state to be compressed
 */
void MPS::compress(int Dc,const MPS &mps,int n_iter){

   int L = mps.size();

   //initial guess by performing svd compression of uncanonicalized state: output is right-canonicalized state
   guess(Right,Dc,mps);

   //construct renormalized operators
   std::vector< TArray<double,2> > RO(L - 1);
   std::vector< TArray<double,2> > LO(L - 1);

   compress::init_ro(Right,RO,mps,*this);

   int iter = 0;

   while(iter < n_iter){

      //first site
      int M = mps[0].shape(1);
      int N = RO[0].shape(0);
      int K = mps[0].shape(2);

      blas::gemm(CblasRowMajor,CblasNoTrans,CblasConjTrans, M, N, K, 1.0, mps[0].data(),K,RO[0].data(),K,0.0,(*this)[0].data(),N);

      //QR
      Geqrf((*this)[0],RO[0]);

      //paste to next matrix
      TArray<double,3> tmp;

      M = RO[0].shape(0);
      N = (*this)[1].shape(1) * (*this)[1].shape(2);
      K = RO[0].shape(1);

      tmp.resize(shape(M,(*this)[1].shape(1),(*this)[1].shape(2)));

      blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, RO[0].data(),K,(*this)[1].data(),N,0.0,tmp.data(),N);

      (*this)[1] = std::move(tmp);

      compress::update_L(0,LO,mps,*this);

      for(int i = 1;i < L - 1;++i){

         M = mps[i].shape(0) * mps[i].shape(1);
         N = RO[i].shape(0);
         K = RO[i].shape(1);

         tmp.resize(shape(mps[i].shape(0),mps[i].shape(1),N));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasConjTrans, M, N, K, 1.0, mps[i].data(),K,RO[i].data(),K,0.0,tmp.data(),N);

         M = LO[i-1].shape(0);
         N = tmp.shape(1)*tmp.shape(2);
         K = tmp.shape(0);

         (*this)[i].resize(shape(LO[i-1].shape(0),mps[i].shape(1),RO[i].shape(0)));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, LO[i-1].data(),K,tmp.data(),N,0.0,(*this)[i].data(),N);

         Geqrf((*this)[i],RO[i]);

         //paste to next matrix
         M = RO[i].shape(0);
         N = (*this)[i+1].shape(1) * (*this)[i+1].shape(2);
         K = RO[i].shape(1);

         tmp.resize(shape(M,(*this)[i+1].shape(1),(*this)[i+1].shape(2)));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, RO[i].data(),K,(*this)[i+1].data(),N,0.0,tmp.data(),N);

         (*this)[i+1] = std::move(tmp);

         compress::update_L(i,LO,mps,*this);

      }

      //and backward!
      M = LO[L-2].shape(0);
      N = (*this)[L-1].shape(1);
      K = LO[L-2].shape(1);

      blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, LO[L-2].data(),K,mps[L-1].data(),N,0.0,(*this)[L-1].data(),N);

      //LQ
      Gelqf(LO[L - 2],(*this)[L - 1]);

      //paste to previous matrix
      M = (*this)[L-2].shape(0)*(*this)[L-2].shape(1);
      N = LO[L-2].shape(1);
      K = LO[L-2].shape(0);

      tmp.resize(shape((*this)[L-2].shape(0),(*this)[L-2].shape(1),N));

      blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, (*this)[L-2].data(),K,LO[L-2].data(),N,0.0,tmp.data(),N);

      (*this)[L - 2] = std::move(tmp);

      compress::update_R(L-1,RO,mps,*this);

      for(int i = L - 2;i > 0;--i){

         M = mps[i].shape(0) * mps[i].shape(1);
         N = RO[i].shape(0);
         K = RO[i].shape(1);

         tmp.resize(shape(mps[i].shape(0),mps[i].shape(1),N));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasConjTrans, M, N, K, 1.0, mps[i].data(),K,RO[i].data(),K,0.0,tmp.data(),N);

         M = LO[i-1].shape(0);
         N = tmp.shape(1)*tmp.shape(2);
         K = tmp.shape(0);

         (*this)[i].resize(shape(LO[i-1].shape(0),mps[i].shape(1),RO[i].shape(0)));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, LO[i-1].data(),K,tmp.data(),N,0.0,(*this)[i].data(),N);

         Gelqf(LO[i],(*this)[i]);

         //paste to previous matrix
         M = (*this)[i-1].shape(0)*(*this)[i-1].shape(1);
         N = LO[i].shape(1);
         K = LO[i].shape(0);

         tmp.resize(shape((*this)[i-1].shape(0),(*this)[i-1].shape(1),N));

         blas::gemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K, 1.0, (*this)[i-1].data(),K,LO[i].data(),N,0.0,tmp.data(),N);

         (*this)[i-1] = std::move(tmp);

         compress::update_R(i,RO,mps,*this);

      }

      ++iter;

   }

   this->D = Dc;

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
 * reduce the dimension of the edge states after MPO action using thin svd.
 */
void MPS::cut_edges() {

   int L = this->size();

   //Left
   TArray<double,3> U;
   TArray<double,2> V;
   TArray<double,1> S;

   int i = 0;

   //easy compression
   while( (*this)[i].shape(0)*(*this)[i].shape(1) < (*this)[i].shape(2) ){

      U.clear();
      S.clear();
      V.clear();

      Gesvd('S','S',(*this)[i],S,U,V);

      (*this)[i] = std::move(U);

      //multiply S to V
      Dimm(S,V);

      U.clear();

      //and contract V with mps on next site
      Contract(1.0,V,shape(1),(*this)[i+1],shape(0),0.0,U);

      (*this)[i+1] = std::move(U);

      ++i;

   }

   i = L - 1;

   while( (*this)[i].shape(0) > (*this)[L - 1].shape(1)*(*this)[i].shape(2) ){

      //Right
      U.clear();
      V.clear();
      S.clear();

      Gesvd('S','S',(*this)[i],S,V,U);

      (*this)[i] = std::move(U);

      //multiply U and S
      Dimm(V,S);

      //and contract U with mps on previous site
      U.clear();

      Contract(1.0,(*this)[i-1],shape(2),V,shape(0),0.0,U);
      (*this)[i-1] = std::move(U);

      --i;

   }

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
 * fill the MPS by taking the correct part of the PEPS
 * @param option 'b'ottom , 't'op, 'l'eft or 'r'ight
 */
void MPS::fill(char option,const PEPS<double> &peps,const Walker &walker){

   if(option == 'b'){

   }

}
