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

   U.resize(2);

   for(int i = 0;i < 2;++i)
      U[i] = SL_PEPS(D);

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

}

/** 
 * copy constructor with allocation
 */
Environment::Environment(const Environment &env_copy){

   t = env_copy.gt();
   b = env_copy.gb();

   r = env_copy.gr();
   l = env_copy.gl();

   U = env_copy.gU();

   D = env_copy.gD();
   D_aux = env_copy.gD_aux();

   comp_sweeps = env_copy.gcomp_sweeps();

}

/**
 * empty destructor
 */
Environment::~Environment(){ }

/**
 * construct the enviroment mps's for the input PEPS and Walker: make sure the correct SL_PEPS object is filled first
 * @param dir    if 'L' construct full left environment
 *               if 'R' construct full right environment
 *               if 'T' construct full top environment
 *               if 'B' construct full bottom environment
 *               if 'A' construct all environments
 * @param inverse if true inverse walker
 * @param peps input PEPS<double>
 * @param D_aux dimension to which environment will be compressed
 */
void Environment::calc(const char dir,bool inverse){

   if(dir == 'H'){//only bottom and top

      //bottom 
      b[inverse * (Ly - 1)].fill('b',U[inverse]);

      for(int i = 1;i < Ly - 1;++i)
         this->add_layer('b',i,inverse);

      //top
      t[Ly - 2 + inverse * (Ly - 1)].fill('t',U[inverse]);

      for(int i = Ly - 3;i >= 0;--i)
         this->add_layer('t',i,inverse);

   }
   else{//Vertical

      //right
      r[Lx - 2 + inverse * (Lx - 1)].fill('r',U[inverse]);

      for(int i = Lx - 3;i >= 0;--i)
         this->add_layer('r',i,inverse);

      //and left
      l[inverse * (Lx - 1)].fill('l',U[inverse]);

      for(int i = 1;i < Lx - 1;++i)
         this->add_layer('l',i,inverse);

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
 * @param inverse if true return inverse
 * @return the overlap between the regular walker and the PEPS
 */
const SL_PEPS &Environment::gU(bool inverse) const {

   return U[inverse];

}

/**
 * @param inverse if true return inverse
 * @return the overlap between the regular walker and the PEPS
 */
const std::vector< SL_PEPS > &Environment::gU() const {

   return U;

}

/**
 * constract the physical indices of the PEPS with the walker, put it in U object
 * @param inverse if true, take the inverse of the Walker
 * @param peps the PEPS<double> object
 * @param Walker the regular walker
 */
void Environment::sU(bool inverse,const PEPS<double> &peps,const Walker &walker){

   U[inverse].fill(inverse,peps,walker);

}

/**
 * construct the (t,b,l or r) environment on row/col 'rc' by adding a the appropriate peps row/col and compressing the boundary MPS
 * @param dir 't'op, 'b'ottom, 'l'eft or 'r'ight environment
 * @param rc row or column index
 * @param inverse if true, inverse of walker is taken
 */
void Environment::add_layer(const char dir,int rc,bool inverse){

   if(dir == 'b'){

      b[rc + inverse * (Ly - 1)].fill_Random();

      vector< DArray<3> > R(Lx - 1);

      //first construct rightmost operator
      DArray<5> tmp5;
      Contract(1.0,b[rc - 1 + inverse * (Ly - 1)][Lx - 1],shape(1),U[inverse](rc,Lx - 1),shape(2),0.0,tmp5);

      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),b[rc + inverse * (Ly - 1)][Lx - 1],shape(1),0.0,tmp6);

      R[Lx - 2] = tmp6.reshape_clear(shape(tmp6.shape(0),tmp6.shape(2),tmp6.shape(4)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<4> tmp4;
         Contract(1.0,b[rc - 1 + inverse * (Ly - 1)][i],shape(2),R[i],shape(0),0.0,tmp4);

         DArray<4> tmp4bis;
         Contract(1.0,tmp4,shape(1,2),U[inverse](rc,i),shape(2,3),0.0,tmp4bis);

         Contract(1.0,tmp4bis,shape(3,1),b[rc + inverse * (Ly - 1)][i],shape(1,2),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPS
         DArray<4> tmp4;
         Contract(1.0,b[rc - 1 + inverse * (Ly - 1)][0],shape(2),R[0],shape(0),0.0,tmp4);

         DArray<4> tmp4bis;
         Contract(1.0,U[inverse](rc,0),shape(2,3),tmp4,shape(1,2),0.0,tmp4bis);

         b[rc + inverse * (Ly - 1)][0] = tmp4bis.reshape_clear(shape(1,D,tmp4bis.shape(3)));

         //QR
         DArray<2> tmp2;
         Geqrf(b[rc + inverse * (Ly - 1)][0],tmp2);

         //construct new left operator
         tmp5.clear();
         Contract(1.0,b[rc-1 + inverse * (Ly - 1)][0],shape(1),U[inverse](rc,0),shape(2),0.0,tmp5);

         tmp6.clear();
         Contract(1.0,tmp5,shape(3),b[rc + inverse * (Ly - 1)][0],shape(1),0.0,tmp6);

         R[0] = tmp6.reshape_clear(shape(tmp6.shape(1),tmp6.shape(3),tmp6.shape(5)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp4.clear();
            Contract(1.0,R[i - 1],shape(0),b[rc - 1 + inverse * (Ly - 1)][i],shape(0),0.0,tmp4);

            tmp4bis.clear();
            Contract(1.0,tmp4,shape(0,2),U[inverse](rc,i),shape(0,2),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,2,1,3),tmp4);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,R[i],0.0,b[rc + inverse * (Ly - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(b[rc + inverse * (Ly - 1)][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4,b[rc + inverse * (Ly - 1)][i],0.0,R[i]);

         }

         //rightmost site
         tmp4.clear();
         Contract(1.0,R[Lx - 2],shape(0),b[rc - 1 + inverse * (Ly - 1)][Lx - 1],shape(0),0.0,tmp4);

         tmp4bis.clear();
         Contract(1.0,tmp4,shape(0,2),U[inverse](rc,Lx - 1),shape(0,2),0.0,tmp4bis);

         b[rc + inverse * (Ly - 1)][Lx - 1] = tmp4bis.reshape_clear(shape(tmp4bis.shape(0),D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,b[rc + inverse * (Ly - 1)][Lx - 1]);

         //construct new right operator
         tmp5.clear();
         Contract(1.0,b[rc - 1 + inverse * (Ly - 1)][Lx - 1],shape(1),U[inverse](rc,Lx - 1),shape(2),0.0,tmp5);

         tmp6.clear();
         Contract(1.0,tmp5,shape(3),b[rc + inverse * (Ly - 1)][Lx - 1],shape(1),0.0,tmp6);

         R[Lx - 2] = tmp6.reshape_clear(shape(tmp6.shape(0),tmp6.shape(2),tmp6.shape(4)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp4.clear();
            Contract(1.0,b[rc - 1 + inverse * (Ly - 1)][i],shape(2),R[i],shape(0),0.0,tmp4);

            tmp4bis.clear();
            Contract(1.0,tmp4,shape(1,2),U[inverse](rc,i),shape(2,3),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,2,3,1),tmp4);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4,0.0,b[rc + inverse*(Ly - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,b[rc + inverse * (Ly - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4,b[rc + inverse * (Ly - 1)][i],0.0,R[i - 1]);

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

      t[rc + inverse*(Ly - 1)].fill_Random();

      vector< DArray<3> > R(Lx - 1);

      //first construct rightmost operator
      DArray<5> tmp5;
      Contract(1.0,t[rc + 1 + inverse*(Ly - 1)][Lx - 1],shape(1),U[inverse](rc + 1,Lx - 1),shape(1),0.0,tmp5);

      DArray<6> tmp6;
      Contract(1.0,tmp5,shape(3),t[rc + inverse*(Ly - 1)][Lx - 1],shape(1),0.0,tmp6);

      R[Lx - 2] = tmp6.reshape_clear(shape(tmp6.shape(0),tmp6.shape(2),tmp6.shape(4)));

      //now move from right to left to construct the rest
      for(int i = Lx - 2;i > 0;--i){

         DArray<4> tmp4;
         Contract(1.0,t[rc + 1 + inverse*(Ly - 1)][i],shape(2),R[i],shape(0),0.0,tmp4);

         DArray<4> tmp4bis;
         Contract(1.0,tmp4,shape(1,2),U[inverse](rc + 1,i),shape(1,3),0.0,tmp4bis);

         Contract(1.0,tmp4bis,shape(3,1),t[rc + inverse*(Ly - 1)][i],shape(1,2),0.0,R[i - 1]);

      }

      int iter = 0;

      while(iter < comp_sweeps){

         //now start sweeping to get the compressed boundary MPS
         DArray<4> tmp4;
         Contract(1.0,t[rc + 1 + inverse*(Ly - 1)][0],shape(2),R[0],shape(0),0.0,tmp4);

         DArray<4> tmp4bis;
         Contract(1.0,U[inverse](rc + 1,0),shape(1,3),tmp4,shape(1,2),0.0,tmp4bis);

         t[rc + inverse*(Ly - 1)][0] = tmp4bis.reshape_clear(shape(1,D,tmp4bis.shape(3)));

         //QR
         DArray<2> tmp2;
         Geqrf(t[rc + inverse*(Ly - 1)][0],tmp2);

         //construct new left operator
         tmp5.clear();
         Contract(1.0,t[rc+1 + inverse*(Ly - 1)][0],shape(1),U[inverse](rc+1,0),shape(1),0.0,tmp5);

         tmp6.clear();
         Contract(1.0,tmp5,shape(3),t[rc + inverse*(Ly - 1)][0],shape(1),0.0,tmp6);

         R[0] = tmp6.reshape_clear(shape(tmp6.shape(1),tmp6.shape(3),tmp6.shape(5)));

         //now for the rest of the rightgoing sweep.
         for(int i = 1;i < Lx-1;++i){

            tmp4.clear();
            Contract(1.0,R[i - 1],shape(0),t[rc + 1 + inverse*(Ly - 1)][i],shape(0),0.0,tmp4);

            tmp4bis.clear();
            Contract(1.0,tmp4,shape(0,2),U[inverse](rc + 1,i),shape(0,1),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,2,1,3),tmp4);

            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,R[i],0.0,t[rc + inverse*(Ly - 1)][i]);

            //QR
            tmp2.clear();
            Geqrf(t[rc + inverse*(Ly - 1)][i],tmp2);

            //construct new left operator
            Gemm(CblasTrans,CblasNoTrans,1.0,tmp4,t[rc + inverse*(Ly - 1)][i],0.0,R[i]);

         }

         //rightmost site
         tmp4.clear();
         Contract(1.0,R[Lx - 2],shape(0),t[rc + 1 + inverse*(Ly - 1)][Lx - 1],shape(0),0.0,tmp4);

         tmp4bis.clear();
         Contract(1.0,tmp4,shape(0,2),U[inverse](rc + 1,Lx - 1),shape(0,1),0.0,tmp4bis);

         t[rc + inverse*(Ly - 1)][Lx - 1] = tmp4bis.reshape_clear(shape(tmp4bis.shape(0),D,1));

         //LQ
         tmp2.clear();
         Gelqf(tmp2,t[rc + inverse*(Ly - 1)][Lx - 1]);

         //construct new right operator
         tmp5.clear();
         Contract(1.0,t[rc + 1 + inverse*(Ly - 1)][Lx - 1],shape(1),U[inverse](rc + 1,Lx - 1),shape(1),0.0,tmp5);

         tmp6.clear();
         Contract(1.0,tmp5,shape(3),t[rc + inverse*(Ly - 1)][Lx - 1],shape(1),0.0,tmp6);

         R[Lx - 2] = tmp6.reshape_clear(shape(tmp6.shape(0),tmp6.shape(2),tmp6.shape(4)));

         //back to the beginning with a leftgoing sweep
         for(int i = Lx-2;i > 0;--i){

            tmp4.clear();
            Contract(1.0,t[rc + 1 + inverse*(Ly - 1)][i],shape(2),R[i],shape(0),0.0,tmp4);

            tmp4bis.clear();
            Contract(1.0,tmp4,shape(1,2),U[inverse](rc+1,i),shape(1,3),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,2,3,1),tmp4);

            Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp4,0.0,t[rc + inverse*(Ly - 1)][i]);

            //LQ
            tmp2.clear();
            Gelqf(tmp2,t[rc + inverse*(Ly - 1)][i]);

            //construct right operator
            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4,t[rc + inverse*(Ly - 1)][i],0.0,R[i - 1]);

         }

         //multiply the last L matrix with the first matrix:
         DArray<3> tmp3;
         Gemm(CblasNoTrans,CblasNoTrans,1.0,t[rc + inverse*(Ly - 1)][0],tmp2,0.0,tmp3);

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
   /*
      else if(option == 'r'){

      if(!flag_r)
      r[rc].fill_Random();

      vector< DArray<4> > R(Ly - 1);

//first construct rightmost operator
DArray<7> tmp7;
Contract(1.0,r[rc + 1][Ly - 1],shape(1),peps(Ly - 1,rc+1),shape(4),0.0,tmp7);

DArray<8> tmp8;
Contract(1.0,tmp7,shape(1,5),peps(Ly-1,rc+1),shape(4,2),0.0,tmp8);

DArray<8> tmp8bis;
Contract(1.0,tmp8,shape(2,5),r[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

//now move from right to left to construct the rest
for(int i = Ly - 2;i > 0;--i){

DArray<6> tmp6;
Contract(1.0,r[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

tmp7.clear();
Contract(1.0,tmp6,shape(1,3),peps(i,rc+1),shape(4,1),0.0,tmp7);

tmp6.clear();
Contract(1.0,tmp7,shape(1,2,5),peps(i,rc+1),shape(4,1,2),0.0,tmp6);

Contract(1.0,tmp6,shape(2,4,1),r[rc][i],shape(1,2,3),0.0,R[i - 1]);

}

int iter = 0;

while(iter < comp_sweeps){

//now start sweeping to get the compressed boundary MPS
DArray<6> tmp6;
Contract(1.0,r[rc + 1][0],shape(3),R[0],shape(0),0.0,tmp6);

tmp7.clear();
Contract(1.0,peps(0,rc+1),shape(1,4),tmp6,shape(4,2),0.0,tmp7);

tmp6.clear();
Contract(1.0,peps(0,rc+1),shape(2,4,1),tmp7,shape(1,4,5),0.0,tmp6);

r[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

//QR
DArray<2> tmp2;
Geqrf(r[rc][0],tmp2);

//construct new left operator
tmp7.clear();
Contract(1.0,r[rc+1][0],shape(1),peps(0,rc+1),shape(4),0.0,tmp7);

tmp8.clear();
Contract(1.0,tmp7,shape(1,5),peps(0,rc+1),shape(4,2),0.0,tmp8);

tmp8bis.clear();
Contract(1.0,tmp8,shape(2,5),r[rc][0],shape(1,2),0.0,tmp8bis);

R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(7)));

//now for the rest of the rightgoing sweep.
for(int i = 1;i < Ly-1;++i){

tmp6.clear();
Contract(1.0,R[i - 1],shape(0),r[rc + 1][i],shape(0),0.0,tmp6);

tmp7.clear();
Contract(1.0,tmp6,shape(0,3),peps(i,rc+1),shape(3,4),0.0,tmp7);

tmp6.clear();
Contract(1.0,tmp7,shape(0,2,6),peps(i,rc+1),shape(3,4,2),0.0,tmp6);

DArray<6> tmp6bis;
Permute(tmp6,shape(0,2,4,1,3,5),tmp6bis);

Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,r[rc][i]);

//QR
tmp2.clear();
Geqrf(r[rc][i],tmp2);

//construct left operator
Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,r[rc][i],0.0,R[i]);

}

//rightmost site
tmp6.clear();
Contract(1.0,R[Ly - 2],shape(0),r[rc + 1][Ly - 1],shape(0),0.0,tmp6);

tmp7.clear();
Contract(1.0,tmp6,shape(0,3),peps(Ly - 1,rc+1),shape(3,4),0.0,tmp7);

tmp6.clear();
Contract(1.0,tmp7,shape(0,2,6),peps(Ly - 1,rc + 1),shape(3,4,2),0.0,tmp6);

r[rc][Ly - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

//LQ
tmp2.clear();
Gelqf(tmp2,r[rc][Ly - 1]);

//construct new right operator
tmp7.clear();
Contract(1.0,r[rc + 1][Ly - 1],shape(1),peps(Ly - 1,rc+1),shape(4),0.0,tmp7);

tmp8.clear();
Contract(1.0,tmp7,shape(1,5),peps(Ly - 1,rc+1),shape(4,2),0.0,tmp8);

tmp8bis.clear();
Contract(1.0,tmp8,shape(2,5),r[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

//back to the beginning with a leftgoing sweep
for(int i = Ly-2;i > 0;--i){

   tmp6.clear();
   Contract(1.0,r[rc + 1][i],shape(3),R[i],shape(0),0.0,tmp6);

   tmp7.clear();
   Contract(1.0,tmp6,shape(1,3),peps(i,rc+1),shape(4,1),0.0,tmp7);

   tmp6.clear();
   Contract(1.0,tmp7,shape(1,2,5),peps(i,rc+1),shape(4,1,2),0.0,tmp6);

   DArray<6> tmp6bis;
   Permute(tmp6,shape(0,3,5,2,4,1),tmp6bis);

   Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,r[rc][i]);

   //LQ
   tmp2.clear();
   Gelqf(tmp2,r[rc][i]);

   //construct right operator
   Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,r[rc][i],0.0,R[i - 1]);

}

//multiply the last L matrix with the first matrix:
DArray<4> tmp4;
Gemm(CblasNoTrans,CblasNoTrans,1.0,r[rc][0],tmp2,0.0,tmp4);

r[rc][0] = std::move(tmp4);

++iter;

}

//redistribute the norm over the chain
double nrm =  Nrm2(r[rc][0]);

//rescale the first site
Scal((1.0/nrm), r[rc][0]);

//then multiply the norm over the whole chain
r[rc].scal(nrm);

}
else{//left

   if(!flag_l)
      l[rc].fill_Random();

   vector< DArray<4> > R(Ly - 1);

   //first construct rightmost operator
   DArray<7> tmp7;
   Contract(1.0,l[rc - 1][Ly - 1],shape(1),peps(Ly-1,rc),shape(0),0.0,tmp7);

   DArray<8> tmp8;
   Contract(1.0,tmp7,shape(1,4),peps(Ly - 1,rc),shape(0,2),0.0,tmp8);

   DArray<8> tmp8bis;
   Contract(1.0,tmp8,shape(4,7),l[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

   R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

   //now move from right to left to construct the rest
   for(int i = Ly - 2;i > 0;--i){

      DArray<6> tmp6;
      Contract(1.0,l[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

      tmp7.clear();
      Contract(1.0,tmp6,shape(1,3),peps(i,rc),shape(0,1),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(1,2,4),peps(i,rc),shape(0,1,2),0.0,tmp6);

      Contract(1.0,tmp6,shape(3,5,1),l[rc][i],shape(1,2,3),0.0,R[i - 1]);

   }

   int iter = 0;

   while(iter < comp_sweeps){

      //now start sweeping to get the compressed boundary MPS
      DArray<6> tmp6;
      Contract(1.0,l[rc - 1][0],shape(3),R[0],shape(0),0.0,tmp6);

      tmp7.clear();
      Contract(1.0,peps(0,rc),shape(0,1),tmp6,shape(2,4),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,peps(0,rc),shape(0,1,2),tmp7,shape(4,5,0),0.0,tmp6);

      l[rc][0] = tmp6.reshape_clear(shape(1,D,D,tmp6.shape(5)));

      //QR
      DArray<2> tmp2;
      Geqrf(l[rc][0],tmp2);

      //construct new left operator
      tmp7.clear();
      Contract(1.0,l[rc-1][0],shape(1),peps(0,rc),shape(0),0.0,tmp7);

      tmp8.clear();
      Contract(1.0,tmp7,shape(1,4),peps(0,rc),shape(0,2),0.0,tmp8);

      tmp8bis.clear();
      Contract(1.0,tmp8,shape(4,7),l[rc][0],shape(1,2),0.0,tmp8bis);

      R[0] = tmp8bis.reshape_clear(shape(tmp8bis.shape(1),tmp8bis.shape(2),tmp8bis.shape(4),tmp8bis.shape(7)));

      //now for the rest of the rightgoing sweep.
      for(int i = 1;i < Ly-1;++i){

         tmp6.clear();
         Contract(1.0,R[i - 1],shape(0),l[rc - 1][i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(0,3),peps(i,rc),shape(3,0),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(0,2,5),peps(i,rc),shape(3,0,2),0.0,tmp6);

         DArray<6> tmp6bis;
         Permute(tmp6,shape(0,3,5,1,2,4),tmp6bis);

         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp6bis,R[i],0.0,l[rc][i]);

         //QR
         tmp2.clear();
         Geqrf(l[rc][i],tmp2);

         //construct left operator
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp6bis,l[rc][i],0.0,R[i]);

      }

      //rightmost site
      tmp6.clear();
      Contract(1.0,R[Ly - 2],shape(0),l[rc - 1][Ly - 1],shape(0),0.0,tmp6);

      tmp7.clear();
      Contract(1.0,tmp6,shape(0,3),peps(Ly - 1,rc),shape(3,0),0.0,tmp7);

      tmp6.clear();
      Contract(1.0,tmp7,shape(0,2,5),peps(Ly - 1,rc),shape(3,0,2),0.0,tmp6);

      l[rc][Ly - 1] = tmp6.reshape_clear(shape(tmp6.shape(0),D,D,1));

      //LQ
      tmp2.clear();
      Gelqf(tmp2,l[rc][Ly - 1]);

      //construct new right operator
      tmp7.clear();
      Contract(1.0,l[rc - 1][Ly - 1],shape(1),peps(Ly - 1,rc),shape(0),0.0,tmp7);

      tmp8.clear();
      Contract(1.0,tmp7,shape(1,4),peps(Ly - 1,rc),shape(0,2),0.0,tmp8);

      tmp8bis.clear();
      Contract(1.0,tmp8,shape(4,7),l[rc][Ly - 1],shape(1,2),0.0,tmp8bis);

      R[Ly - 2] = tmp8bis.reshape_clear(shape(tmp8bis.shape(0),tmp8bis.shape(3),tmp8bis.shape(5),tmp8bis.shape(6)));

      //back to the beginning with a leftgoing sweep
      for(int i = Ly-2;i > 0;--i){

         tmp6.clear();
         Contract(1.0,l[rc - 1][i],shape(3),R[i],shape(0),0.0,tmp6);

         tmp7.clear();
         Contract(1.0,tmp6,shape(1,3),peps(i,rc),shape(0,1),0.0,tmp7);

         tmp6.clear();
         Contract(1.0,tmp7,shape(1,2,4),peps(i,rc),shape(0,1,2),0.0,tmp6);

         DArray<6> tmp6bis;
         Permute(tmp6,shape(0,2,4,3,5,1),tmp6bis);

         Gemm(CblasTrans,CblasNoTrans,1.0,R[i - 1],tmp6bis,0.0,l[rc][i]);

         //LQ
         tmp2.clear();
         Gelqf(tmp2,l[rc][i]);

         //construct right operator
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp6bis,l[rc][i],0.0,R[i - 1]);

      }

      //multiply the last L matrix with the first matrix:
      DArray<4> tmp4;
      Gemm(CblasNoTrans,CblasNoTrans,1.0,l[rc][0],tmp2,0.0,tmp4);

      l[rc][0] = std::move(tmp4);

      ++iter;

   }

   //redistribute the norm over the chain
   double nrm =  Nrm2(l[rc][0]);

   //rescale the first site
   Scal((1.0/nrm), l[rc][0]);

   //then multiply the norm over the whole chain
   l[rc].scal(nrm);

}
*/
}
