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
Environment::Environment(){ }

/** 
 * constructor with allocation
 * @param D_in bond dimension of peps state
 * @param D_aux_in contraction bond dimension
 * @param comp_sweeps_in sets the number of sweeps done for MPS compression
 */
Environment::Environment(int D_in,int D_aux_in,int comp_sweeps_in){

   t.resize( 2*(Ly - 1) );
   b.resize( 2*(Ly - 1) );

   r.resize( 2*(Lx - 1) );
   l.resize( 2*(Lx - 1) );

   D = D_in;
   D_aux = D_aux_in;
   comp_sweeps = comp_sweeps_in;

   //allocate the memory
   
   //bottom
   int tmp = D;

   for(int i = 0;i < Ly - 1;++i){

      if(tmp < D_aux){

         b[i] = MPS(Lx,D,tmp);
         b[i + Ly - 1] = MPS(Lx,D,tmp);
         tmp *= D;

      }
      else{

         b[i] = MPS(Lx,D,D_aux);
         b[i + Ly-1] = MPS(Lx,D,D_aux);

      }

   }
   
   //top
   tmp = D;

   for(int i = Ly - 2;i >= 0;--i){

      if(tmp < D_aux){

         t[i] = MPS(Lx,D,tmp);
         t[i + Ly - 1] = MPS(Lx,D,tmp);
         tmp *= D;

      }
      else{

         t[i] = MPS(Lx,D,D_aux);
         t[i + Ly-1] = MPS(Lx,D,D_aux);

      }

   }

   //left
   tmp = D;

   for(int i = 0;i < Lx - 1;++i){

      if(tmp < D_aux){

         l[i] = MPS(Ly,D,tmp);
         l[i + Lx - 1] = MPS(Ly,D,tmp);
         tmp *= D;

      }
      else{

         l[i] = MPS(Ly,D,D_aux);
         l[i + Lx - 1] = MPS(Ly,D,D_aux);

      }

   }
   
   //finally right
   tmp = D;

   for(int i = Lx - 2;i >= 0;--i){

      if(tmp < D_aux){

         r[i] = MPS(Ly,D,tmp);
         r[i + Lx - 1] = MPS(Ly,D,tmp);
         tmp *= D;

      }
      else{

         r[i] = MPS(Ly,D,D_aux);
         r[i + Lx - 1] = MPS(Ly,D,D_aux);

      }

   }

   flag_b = true;
   flag_t = true;
   flag_l = true;
   flag_r = true;

}

/** 
 * copy constructor with allocation
 */
Environment::Environment(const Environment &env_copy){

   t = env_copy.gt();
   b = env_copy.gb();

   r = env_copy.gr();
   l = env_copy.gl();

   D = env_copy.gD();
   D_aux = env_copy.gD_aux();

   comp_sweeps = env_copy.gcomp_sweeps();

   flag_b = env_copy.gflag_b();
   flag_t = env_copy.gflag_t();
   flag_l = env_copy.gflag_l();
   flag_r = env_copy.gflag_r();

}

/**
 * empty destructor
 */
Environment::~Environment(){ }

/**
 * construct the enviroment mps's for the input PEPS and Walker
 * @param dir    if 'T' construct full top  environment
 *               if 'B' construct full bottom environment
 *               if 'L' construct full left  environment
 *               if 'R' construct full right environment
 *               if 'A' construct all environments, t,l,r and b
 * @param inverse if true inverse walker
 * @param walker input Walker object
 */
void Environment::calc(const char dir,bool inverse,const Walker &walker){

   if(dir == 'B' || dir == 'A'){//bottom

      //bottom 
      b[inverse * (Ly - 1)].fill('b',inverse,walker);

      for(int i = 1;i < Ly - 1;++i)
         this->add_layer('b',i,inverse,walker);

      flag_b = false;

   }

   if(dir == 'T' || dir == 'A'){//top

      //top
      t[Ly - 2 + inverse * (Ly - 1)].fill('t',inverse,walker);

      for(int i = Ly - 3;i >= 0;--i)
         this->add_layer('t',i,inverse,walker);

      flag_t = false;

   }

   if(dir == 'R' || dir == 'A'){//right

      //right
      r[Lx - 2 + inverse * (Lx - 1)].fill('r',inverse,walker);

      for(int i = Lx - 3;i >= 0;--i)
         this->add_layer('r',i,inverse,walker);

      flag_r = false;

   }

   if(dir == 'L' || dir == 'A'){//left

      //and left
      l[inverse * (Lx - 1)].fill('l',inverse,walker);

      for(int i = 1;i < Lx - 1;++i)
         this->add_layer('l',i,inverse,walker);

      flag_l = false;

   }

}

/**
 * test if the enviroment is correctly contracted
 * @param inverse if true inverse walker is used
 */
void Environment::test(bool inverse){

   cout << endl;
   cout << "FROM BOTTOM TO TOP" << endl;
   cout << endl;
   for(int i = 0;i < Ly - 1;++i)
      cout << i << "\t" << b[i + inverse * (Ly - 1)].dot(t[i + inverse * (Ly - 1) ]) << endl;

   cout << endl;
   cout << "FROM LEFT TO RIGHT" << endl;
   cout << endl;
   for(int i = 0;i < Lx - 1;++i)
      cout << i << "\t" << l[i + inverse * (Lx - 1)].dot(r[i + inverse * (Lx - 1)]) << endl;
   cout << endl;

}

/**
 * @return the flag saying if there is a previous guess for the bottom environment or not
 */
bool Environment::gflag_b() const {

   return flag_b;

}

/**
 * @return the flag saying if there is a previous guess for the top environment or not
 */
bool Environment::gflag_t() const {

   return flag_t;

}

/**
 * @return the flag saying if there is a previous guess for the left environment or not
 */
bool Environment::gflag_l() const {

   return flag_l;

}

/**
 * @return the flag saying if there is a previous guess for the right environment or not
 */
bool Environment::gflag_r() const {

   return flag_r;

}

/**
 * const version
 * @param option if true, regular walker, if false inverse walker
 * @param col the column index
 * @return the right boundary 'MPS' environment on column col
 */
const MPS &Environment::gr(bool option,int col) const {

   return r[col + option * (Lx - 1)];

}

/**
 * @param col the column index: access version
 * @param option if true, regular walker, if false inverse walker
 * @return the right boundary 'MPS' environment on column col
 */
MPS &Environment::gr(bool option,int col) {

   return r[col + option * (Lx - 1)];

}

/**
 * const version
 * @param option if true, regular walker, if false inverse walker
 * @param col the column index
 * @return the left boundary 'MPS' environment on column col
 */
const MPS &Environment::gl(bool option,int col) const {

   return l[col + option * (Lx - 1)];

}

/**
 * @param option if true, regular walker, if false inverse walker
 * @param col the column index: access version
 * @return the left boundary 'MPS' environment on column col
 */
MPS &Environment::gl(bool option,int col) {

   return l[col + option * (Lx - 1)];

}

/**
 * const version
 * @param option if true, regular walker, if false inverse walker
 * @param row the row index
 * @return the top boundary 'MPS' environment on row 'row'
 */
const MPS &Environment::gt(bool option,int row) const {

   return t[row + option * (Ly - 1)];

}

/**
 * access version
 * @param option if true, regular walker, if false inverse walker
 * @param row the row index
 * @return the top boundary 'MPS' environment on row 'row'
 */
MPS &Environment::gt(bool option,int row) {

   return t[row + option * (Ly - 1)];

}

/**
 * const version
 * @param option if true, regular walker, if false inverse walker
 * @param row the row index
 * @return the bottom boundary 'MPS' environment on row 'row'
 */
const MPS &Environment::gb(bool option,int row) const {

   return b[row + option * (Ly - 1)];

}

/**
 * access version
 * @param option if true, regular walker, if false inverse walker
 * @param row the row index
 * @return the bottom boundary 'MPS' environment on row 'row'
 */
MPS &Environment::gb(bool option,int row) {

   return b[row + option * (Ly - 1)];

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
int Environment::gD_aux() const {

   return D_aux;

}

/**
 * @return the auxiliary bond dimension for the contraction
 **/
int Environment::gD() const {

   return D;

}

/**
 * @return the number of sweeps performed during compression
 **/
int Environment::gcomp_sweeps() const {

   return comp_sweeps;

}


/**
 * set a new bond dimension
 */
void Environment::sD(int D_in) {

   D = D_in;

}

/**
 * set a new auxiliary bond dimension
 */
void Environment::sD_aux(int D_aux_in) {

   D_aux = D_aux_in;

}

/**
 * @return the full bottom boundary 'MPS'
 */
const vector< MPS > &Environment::gb() const {

   return b;

}

/**
 * @return the full top boundary 'MPS'
 */
const vector< MPS > &Environment::gt() const {

   return t;

}

/**
 * @return the full left boundary 'MPS'
 */
const vector< MPS > &Environment::gl() const {

   return l;

}

/**
 * @return the full right boundary 'MPS'
 */
const vector< MPS > &Environment::gr() const {

   return r;

}

/**
 * construct the (t,b,l or r) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPS
 * @param dir 't'op, 'b'ottom, 'l'eft or 'r'ight environment
 * @param rc row or column index
 * @param inverse if true, inverse of walker is taken
 */
void Environment::add_layer(const char dir,int rc,bool inverse,const Walker &walker){

   if(dir == 'b'){

      DArray<4> tmp4;
      DArray<4> tmp4bis;
      DArray<3> tmp3bis;
      DArray<3> tmp3;

      if(flag_b)
         b[rc + inverse * (Ly - 1)].fill_Random();

      vector< DArray<3> > R(Lx - 1);

      //first construct rightmost operator
      bool s = inverse ^ walker[rc*Lx + Lx - 1];

      Gemm(CblasNoTrans,CblasTrans,1.0,b[rc - 1 + inverse * (Ly - 1)][Lx - 1],peps[0](rc,Lx - 1,s),0.0,tmp3);

      int M = tmp3.shape(0) * tmp3.shape(1);
      int N = b[rc + inverse * (Ly - 1)][Lx - 1].shape(0);
      int K = tmp3.shape(2);

      R[Lx - 2].resize(shape(tmp3.shape(0),tmp3.shape(1),b[rc + inverse * (Ly - 1)][Lx - 1].shape(0)));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,b[rc + inverse * (Ly - 1)][Lx - 1].data(),K,0.0,R[Lx-2].data(),N);

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,b[rc - 1 + inverse * (Ly - 1)][i],R[i],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,3,1,2),tmp4bis);

         //get spin
         s = inverse ^ walker[rc*Lx + i];

         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[0](rc,i,s),0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,2,3,1),tmp4bis);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,b[rc + inverse * (Ly - 1)][i],0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //spin
         s = inverse ^ walker[rc*Lx];

         //now start sweeping to get the compressed boundary MPS
         tmp3.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,b[rc - 1 + inverse * (Ly - 1)][0],peps[2](rc,0,s),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(0,2,1),tmp3bis);

         M = tmp3bis.shape(2);
         N = R[0].shape(2);
         K = R[0].shape(0)*R[0].shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3bis.data(),M,R[0].data(),N,0.0,b[rc + inverse * (Ly - 1)][0].data(),N);

         //QR
         DArray<2> tmp2;
         Geqrf(b[rc + inverse * (Ly - 1)][0],tmp2);

         //construct new left operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = b[rc + inverse * (Ly - 1)][0].shape(2);
         K = tmp3bis.shape(2);

         R[0].resize(shape(tmp3bis.shape(0),tmp3bis.shape(1),N));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, tmp3bis.data(),K,b[rc + inverse * (Ly - 1)][0].data(),N,0.0,R[0].data(),N);

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            //spin
            s = inverse ^ walker[rc*Lx + i];

            tmp4.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],b[rc - 1 + inverse * (Ly - 1)][i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[2](rc,i,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,1,3),tmp4bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,R[i],0.0,b[rc + inverse * (Ly - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(b[rc + inverse * (Ly - 1)][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4bis,b[rc + inverse * (Ly - 1)][i],0.0,R[i]);

         }

         //rightmost site
         s = inverse ^ walker[rc*Lx + Lx - 1];

         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,b[rc - 1 + inverse * (Ly - 1)][Lx - 1],peps[0](rc,Lx - 1,s),0.0,tmp3);

         //construct new right operator
         M = R[Lx - 2].shape(2);
         N = tmp3.shape(2);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, R[Lx - 2].data(),M,tmp3.data(),N,0.0,b[rc + inverse * (Ly - 1)][Lx - 1].data(),N);

         //LQ
         tmp2.clear();
         Gelqf(tmp2,b[rc + inverse * (Ly - 1)][Lx - 1]);

         //construct new right operator
         M = tmp3.shape(0) * tmp3.shape(1);
         N = b[rc + inverse * (Ly - 1)][Lx - 1].shape(0);
         K = tmp3.shape(2);

         R[Lx - 2].resize(shape(tmp3.shape(0),tmp3.shape(1),b[rc + inverse * (Ly - 1)][Lx - 1].shape(0)));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,b[rc + inverse * (Ly - 1)][Lx - 1].data(),K,0.0,R[Lx-2].data(),N);

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            //get spin
            s = inverse ^ walker[rc*Lx + i];

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,b[rc - 1 + inverse * (Ly - 1)][i],R[i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[0](rc,i,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4bis,0.0,b[rc + inverse*(Ly - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,b[rc + inverse * (Ly - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,b[rc + inverse * (Ly - 1)][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<3> tmp3;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,b[rc + inverse * (Ly - 1)][0],tmp2,0.0,tmp3);

         b[rc + inverse*(Ly - 1)][0] = std::move(tmp3);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(b[rc + inverse*(Ly - 1)][0]);

      //rescale the first site
      Scal((1.0/nrm), b[rc + inverse*(Ly - 1)][0]);

      //then multiply the norm over the whole chain
      b[rc + inverse*(Ly - 1)].scal(nrm);

   }
   else if(dir == 't'){

      if(flag_t)
         t[rc + inverse*(Ly - 1)].fill_Random();

      DArray<4> tmp4;
      DArray<4> tmp4bis;
      DArray<3> tmp3;
      DArray<3> tmp3bis;

      vector< DArray<3> > R(Lx - 1);

      //first construct rightmost operator
      bool s = inverse ^ walker[(rc + 1)*Lx + Lx - 1];

      Gemm(CblasNoTrans,CblasTrans,1.0,t[rc + 1 + inverse * (Ly - 1)][Lx - 1],peps[2](rc + 1,Lx - 1,s),0.0,tmp3);

      int M = tmp3.shape(0) * tmp3.shape(1);
      int N = t[rc + inverse * (Ly - 1)][Lx - 1].shape(0);
      int K = tmp3.shape(2);

      R[Lx - 2].resize(shape(tmp3.shape(0),tmp3.shape(1),t[rc + inverse * (Ly - 1)][Lx - 1].shape(0)));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,t[rc + inverse * (Ly - 1)][Lx - 1].data(),K,0.0,R[Lx-2].data(),N);

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc + 1 + inverse * (Ly - 1)][i],R[i],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,3,1,2),tmp4bis);

         //get spin
         s = inverse ^ walker[(rc + 1)*Lx + i];

         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[2](rc + 1,i,s),0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,2,3,1),tmp4bis);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,t[rc + inverse * (Ly - 1)][i],0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //spin
         s = inverse ^ walker[(rc + 1)*Lx];

         //now start sweeping to get the compressed boundary MPS
         tmp3.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,t[rc + 1 + inverse * (Ly - 1)][0],peps[0](rc + 1,0,s),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(0,2,1),tmp3bis);

         M = tmp3bis.shape(2);
         N = R[0].shape(2);
         K = R[0].shape(0)*R[0].shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3bis.data(),M,R[0].data(),N,0.0,t[rc + inverse * (Ly - 1)][0].data(),N);

         //QR
         DArray<2> tmp2;
         Geqrf(t[rc + inverse * (Ly - 1)][0],tmp2);

         //construct new left operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = t[rc + inverse * (Ly - 1)][0].shape(2);
         K = tmp3bis.shape(2);

         R[0].resize(shape(tmp3bis.shape(0),tmp3bis.shape(1),N));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, tmp3bis.data(),K,t[rc + inverse * (Ly - 1)][0].data(),N,0.0,R[0].data(),N);

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            //spin
            s = inverse ^ walker[(rc + 1)*Lx + i];

            tmp4.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],t[rc + 1 + inverse * (Ly - 1)][i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](rc + 1,i,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,1,3),tmp4bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,R[i],0.0,t[rc + inverse * (Ly - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(t[rc + inverse * (Ly - 1)][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4bis,t[rc + inverse * (Ly - 1)][i],0.0,R[i]);

         }

         //rightmost site
         s = inverse ^ walker[(rc + 1)*Lx + Lx - 1];

         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,t[rc + 1 + inverse * (Ly - 1)][Lx - 1],peps[2](rc + 1,Lx - 1,s),0.0,tmp3);

         //construct new right operator
         M = R[Lx - 2].shape(2);
         N = tmp3.shape(2);
         K = tmp3.shape(0) * tmp3.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, R[Lx - 2].data(),M,tmp3.data(),N,0.0,t[rc + inverse * (Ly - 1)][Lx - 1].data(),N);

         //LQ
         tmp2.clear();
         Gelqf(tmp2,t[rc + inverse * (Ly - 1)][Lx - 1]);

         //construct new right operator
         M = tmp3.shape(0) * tmp3.shape(1);
         N = t[rc + inverse * (Ly - 1)][Lx - 1].shape(0);
         K = tmp3.shape(2);

         R[Lx - 2].resize(shape(tmp3.shape(0),tmp3.shape(1),t[rc + inverse * (Ly - 1)][Lx - 1].shape(0)));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,t[rc + inverse * (Ly - 1)][Lx - 1].data(),K,0.0,R[Lx-2].data(),N);

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            //get spin
            s = inverse ^ walker[(rc + 1)*Lx + i];

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc + 1 + inverse * (Ly - 1)][i],R[i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[2](rc + 1,i,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4bis,0.0,t[rc + inverse*(Ly - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,t[rc + inverse * (Ly - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,t[rc + inverse * (Ly - 1)][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<3> tmp3;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc + inverse * (Ly - 1)][0],tmp2,0.0,tmp3);

         t[rc + inverse*(Ly - 1)][0] = std::move(tmp3);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(t[rc + inverse*(Ly - 1)][0]);

      //rescale the first site
      Scal((1.0/nrm), t[rc + inverse*(Ly - 1)][0]);

      //then multiply the norm over the whole chain
      t[rc + inverse*(Ly - 1)].scal(nrm);

   }
   else if(dir == 'r'){

      if(flag_r)
         r[rc + inverse*(Lx - 1)].fill_Random();

      vector< DArray<3> > R(Ly - 1);

      DArray<4> tmp4;
      DArray<4> tmp4bis;

      //first construct rightmost operator
      bool s = inverse ^ walker[(Ly - 1)*Lx + rc + 1];

      DArray<3> tmp3;
      Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc + inverse*(Lx - 1)][Ly - 1],peps[0](Ly - 1,rc+1,s),0.0,tmp3);

      DArray<3> tmp3bis;
      Permute(tmp3,shape(2,1,0),tmp3bis);

      int M = r[rc + 1 + inverse*(Lx - 1)][Ly - 1].shape(0);
      int N = tmp3bis.shape(1) * tmp3bis.shape(2);
      int K = tmp3bis.shape(0);

      R[Lx - 2].resize(shape(M,tmp3bis.shape(1),tmp3bis.shape(2)));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,r[rc + 1 + inverse*(Lx - 1)][Ly - 1].data(),K,tmp3bis.data(),N,0.0,R[Lx-2].data(),N);

      //now move from right to left to construct the rest
      for(int i = Ly - 2;i > 0;--i){

         s = inverse ^ walker[i*Lx + rc + 1];

         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,r[rc + inverse*(Lx - 1)][i],R[i],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,peps[0](i,rc+1,s),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc + 1 + inverse*(Lx - 1)][i],tmp4bis,0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPS
         s = inverse ^ walker[rc + 1];

         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](0,rc+1,s),r[rc + 1 + inverse*(Lx - 1)][0],0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         M = tmp3bis.shape(2);
         N = R[0].shape(2);
         K = tmp3bis.shape(0) * tmp3bis.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),M,R[0].data(),N,0.0,r[rc + inverse*(Lx - 1)][0].data(),N);

         //QR
         DArray<2> tmp2;
         Geqrf(r[rc + inverse*(Lx - 1)][0],tmp2);

         //construct new left operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = r[rc + inverse*(Lx - 1)][0].shape(2);
         K = tmp3bis.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,r[rc + inverse*(Lx - 1)][0].data(),N,0.0,R[0].data(),N);

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Ly-1;++i){

            s = inverse ^ walker[i*Lx + rc + 1];

            tmp4.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],r[rc + 1 + inverse*(Lx - 1)][i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(3,1,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[0](i,rc+1,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,2,0,3),tmp4bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,R[i],0.0,r[rc + inverse*(Lx - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(r[rc + inverse*(Lx - 1)][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4bis,r[rc + inverse*(Lx - 1)][i],0.0,R[i]);

         }

         //rightmost site
         s = inverse ^ walker[(Ly - 1)*Lx + rc + 1];

         tmp4.clear();
         Contract(1.0,r[rc + 1 + inverse*(Lx - 1)][Ly - 1],shape(1,2),peps[0](Ly - 1,rc+1,s),shape(3,1),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(0,2,1),tmp3bis);

         //construct new left operator
         M = R[Ly - 2].shape(2);
         N = tmp3bis.shape(2);
         K = tmp3bis.shape(0) * tmp3bis.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,R[Ly - 2].data(),M,tmp3bis.data(),N,0.0,r[rc  + inverse*(Lx - 1)][Ly - 1].data(),N);

         //LQ
         tmp2.clear();
         Gelqf(tmp2,r[rc + inverse*(Lx - 1)][Ly - 1]);

         //construct new right operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = r[rc  + inverse*(Lx - 1)][Ly - 1].shape(0);
         K = tmp3bis.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3bis.data(),K,r[rc  + inverse*(Lx - 1)][Ly - 1].data(),K,0.0,R[Ly - 2].data(),N);

         //back to the beginning with a leftgoing sweep
         for(int i = Ly-2;i > 0;--i){

            s = inverse ^ walker[i*Lx + rc + 1];

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc + 1 + inverse*(Lx - 1)][i],R[i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[1](i,rc+1,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4bis,0.0,r[rc + inverse*(Lx - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,r[rc + inverse*(Lx - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,r[rc + inverse*(Lx - 1)][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<3> tmp3;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc + inverse*(Lx - 1)][0],tmp2,0.0,tmp3);

         r[rc + inverse*(Lx - 1)][0] = std::move(tmp3);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(r[rc + inverse*(Lx - 1)][0]);

      //rescale the first site
      Scal((1.0/nrm),r[rc + inverse*(Lx - 1)][0]);

      //then multiply the norm over the whole chain
      r[rc + inverse*(Lx - 1)].scal(nrm);

   }
   else{//left

      if(flag_l)
         l[rc + inverse*(Lx - 1)].fill_Random();

      vector< DArray<3> > R(Ly - 1);

      DArray<4> tmp4;
      DArray<4> tmp4bis;

      //first construct rightmost operator
      bool s = inverse ^ walker[(Ly - 1)*Lx + rc];

      DArray<3> tmp3;
      Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc + inverse*(Lx - 1)][Ly - 1],peps[3](Ly - 1,rc,s),0.0,tmp3);

      DArray<3> tmp3bis;
      Permute(tmp3,shape(2,1,0),tmp3bis);

      int M = l[rc - 1 + inverse*(Lx - 1)][Ly - 1].shape(0);
      int N = tmp3bis.shape(1) * tmp3bis.shape(2);
      int K = tmp3bis.shape(0);

      R[Lx - 2].resize(shape(M,tmp3bis.shape(1),tmp3bis.shape(2)));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,l[rc - 1 + inverse*(Lx - 1)][Ly - 1].data(),K,tmp3bis.data(),N,0.0,R[Lx-2].data(),N);

      //now move from right to left to construct the rest
      for(int i = Ly - 2;i > 0;--i){

         s = inverse ^ walker[i*Lx + rc];

         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,l[rc + inverse*(Lx - 1)][i],R[i],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,peps[3](i,rc,s),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc - 1 + inverse*(Lx - 1)][i],tmp4bis,0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPS
         s = inverse ^ walker[rc];

         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[3](0,rc,s),l[rc - 1 + inverse*(Lx - 1)][0],0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         M = tmp3bis.shape(2);
         N = R[0].shape(2);
         K = tmp3bis.shape(0) * tmp3bis.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),M,R[0].data(),N,0.0,l[rc + inverse*(Lx - 1)][0].data(),N);

         //QR
         DArray<2> tmp2;
         Geqrf(l[rc + inverse*(Lx - 1)][0],tmp2);

         //construct new left operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = l[rc + inverse*(Lx - 1)][0].shape(2);
         K = tmp3bis.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,l[rc + inverse*(Lx - 1)][0].data(),N,0.0,R[0].data(),N);

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Ly-1;++i){

            s = inverse ^ walker[i*Lx + rc];

            tmp4.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],l[rc - 1 + inverse*(Lx - 1)][i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(3,1,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,peps[3](i,rc,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,2,0,3),tmp4bis);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,R[i],0.0,l[rc + inverse*(Lx - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(l[rc + inverse*(Lx - 1)][i],tmp2);

            //construct left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4bis,l[rc + inverse*(Lx - 1)][i],0.0,R[i]);

         }

         //rightmost site
         s = inverse ^ walker[(Ly - 1)*Lx + rc];

         tmp4.clear();
         Contract(1.0,l[rc - 1 + inverse*(Lx - 1)][Ly - 1],shape(1,2),peps[3](Ly - 1,rc,s),shape(3,1),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(0,2,1),tmp3bis);

         //construct new left operator
         M = R[Ly - 2].shape(2);
         N = tmp3bis.shape(2);
         K = tmp3bis.shape(0) * tmp3bis.shape(1);

         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,R[Ly - 2].data(),M,tmp3bis.data(),N,0.0,l[rc  + inverse*(Lx - 1)][Ly - 1].data(),N);

         //LQ
         tmp2.clear();
         Gelqf(tmp2,l[rc + inverse*(Lx - 1)][Ly - 1]);

         //construct new right operator
         M = tmp3bis.shape(0) * tmp3bis.shape(1);
         N = l[rc  + inverse*(Lx - 1)][Ly - 1].shape(0);
         K = tmp3bis.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3bis.data(),K,l[rc  + inverse*(Lx - 1)][Ly - 1].data(),K,0.0,R[Ly - 2].data(),N);

         //back to the beginning with a leftgoing sweep
         for(int i = Ly-2;i > 0;--i){

            s = inverse ^ walker[i*Lx + rc];

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc - 1 + inverse*(Lx - 1)][i],R[i],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](i,rc,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4bis,0.0,l[rc + inverse*(Lx - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,l[rc + inverse*(Lx - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,l[rc + inverse*(Lx - 1)][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<3> tmp3;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc + inverse*(Lx - 1)][0],tmp2,0.0,tmp3);

         l[rc + inverse*(Lx - 1)][0] = std::move(tmp3);

         ++iter;

      }

      //redistribute the norm over the chain
      double nrm =  Nrm2(l[rc + inverse*(Lx - 1)][0]);

      //rescale the first site
      Scal((1.0/nrm), l[rc + inverse*(Lx - 1)][0]);

      //then multiply the norm over the whole chain
      l[rc + inverse*(Lx - 1)].scal(nrm);

   }

}
