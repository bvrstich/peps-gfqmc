#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

#include "include.h"

using namespace global;

/**
 * construct a Walker object: initialize on AF state
 * @param n_trot_in number of trotter terms
 */
Walker::Walker() : std::vector< bool >( Lx * Ly ){

   weight = 1.0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         if( (r + c)%2 == 0)
            (*this)[ r*Lx + c ] = true;
         else
            (*this)[ r*Lx + c ] = false;

      }

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< bool >(walker) {

   this->weight = walker.gWeight();
   this->overlap = walker.gOverlap();
   this->EL = walker.gEL();

}

/**
 * destructor
 */
Walker::~Walker(){ }

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}

/** 
 * @return the overlap of the walker with the Trial
 */
double Walker::gOverlap() const{

   return overlap; 

}

/** 
 * @return the local energy
 */
double Walker::gEL() const{

   return EL; 

}

ostream &operator<<(ostream &output,const Walker &walker_p){

   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx;++c)
         output << walker_p[r*Lx + c] << " ";

      output << endl;

   }

   return output;

}

/**
 * get the energy between *this and walker_i <*this|H| walker_i>
 * @param dtau timestep
 * @param walker_i input Walker
 */
double Walker::exp_en(const Walker &walker_i){

   double tmp = 0.0;

   //first horizontal
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx - 1;++c){

         //Sz Sz
         if( ((*this)[r*Lx + c] == walker_i[r*Lx + c]) && ((*this)[r*Lx + (c + 1)] == walker_i[r*Lx + (c + 1)]) ){

            if( (*this)[r*Lx + c] == (*this)[r*Lx + (c + 1)] )//up up or down down
               tmp += 0.25;
            else //up down or down up
               tmp -= 0.25;

         }

         //S+ S- (extra minus sign for sign problem)
         if( ((*this)[r*Lx + c] != walker_i[r*Lx + c]) && ((*this)[r*Lx + (c + 1)] != walker_i[r*Lx + (c + 1)]) ){

            if( (*this)[r*Lx + c] != (*this)[r*Lx + (c + 1)] )
               tmp -= 0.5;

         }

      }

   }

   //then vertical
   for(int c = 0;c < Lx;++c){

      for(int r = 0;r < Ly - 1;++r){

         //Sz Sz
         if( ((*this)[r*Lx + c] == walker_i[r*Lx + c]) && ((*this)[(r + 1)*Lx + c] == walker_i[(r + 1)*Lx + c]) ){

            if( (*this)[r*Lx + c] == (*this)[(r + 1)*Lx + c] )//up up or down down
               tmp += 0.25;
            else //up down or down up
               tmp -= 0.25;

         }

         //S+ S-
         if( ((*this)[r*Lx + c] != walker_i[r*Lx + c]) && ((*this)[(r + 1)*Lx + c] != walker_i[(r + 1)*Lx + c]) ){

            if( (*this)[r*Lx + c] != (*this)[(r + 1)*Lx + c] )
               tmp -= 0.5;

         }

      }

   }

   return tmp;

}

/**
 * calculate the local energy expectation value
 * @param peps trial wave function represented as peps
 */
void Walker::calc_EL(const PEPS< double > &peps){

   // ---- || evaluate the expectation values in an MPO/MPS manner, first from bottom to top, then left to right || ----

   double energy = 0.0;

   int M,N,K;

   //calculate the single layer contractions first:
   Environment::U.fill('H',false,peps,*this);
   Environment::I.fill('H',true,peps,*this);

   //first construct the top and bottom (horizontal) environment layers
   Environment::calc_env('H',peps,*this);

   // #################################################################
   // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || bottom row: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<2> > R(Lx - 1);

   //first the rightmost operator
   DArray<4> tmp4;
   DArray<3> tmp3;

   //tmp comes out index (t,b)
   Contract(1.0,Environment::t[0][Lx - 1],shape(1),Environment::b[0][Lx - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp4.reshape_clear(shape(Environment::t[0][Lx - 1].shape(0),Environment::b[0][Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();
      Contract(1.0,Environment::t[0][col],shape(2),R[col],shape(0),0.0,tmp3);

      Contract(1.0,tmp3,shape(1,2),Environment::b[0][col],shape(1,2),0.0,R[col-1]);

   }

   //2 left going operators: S+/- and 1
   DArray<2> LI;
   DArray<2> LU;

   TArray<double,5> tmp5;

   //tmp comes out index (t,b)
   Contract(1.0,Environment::t[0][0],shape(1),Environment::U(0,0),shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear( shape(Environment::t[0][0].shape(2),Environment::U(0,0).shape(3)) );

   //only calculate LI if it contributes
   if( (*this)[0] != (*this)[1] ){

      //tmp comes out index (t,b)
      Contract(1.0,Environment::t[0][0],shape(1),Environment::I(0,0),shape(1),0.0,tmp5);

      LI = tmp5.reshape_clear( shape(Environment::t[0][0].shape(2),Environment::I(0,0).shape(3)) );

   }

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //first close down the S+/- action:

      //only calculate if it contributes
      if( (*this)[col - 1] != (*this)[col] ){

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(1.0,Environment::t[0][col],shape(2),R[col],shape(0),0.0,tmp3);

         // 1) paste I to the right
         M = tmp3.shape(0);
         N = Environment::I(0,col).shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,Environment::I(0,col).data(),K,0.0,R[col-1].data(),N);

         //contract with left LI 
         energy += 0.5 * Dot(LI,R[col - 1]);

      }

      //construct left renormalized operators for next site: first paste top to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(0),Environment::t[0][col],shape(0),0.0,tmp3);

      //1) construct new unity on the left
      LU.resize(Environment::t[0][col].shape(2),Environment::U(0,col).shape(3));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::U(0,col).data(),N,0.0,LU.data(),N);

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         LI.resize(Environment::t[0][col].shape(2),Environment::I(0,col).shape(3));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::I(0,col).data(),N,0.0,LI.data(),N);

      }

   }

   //last site of bottom row: close down LI
   if((*this)[Lx - 2] != (*this)[Lx - 1]){

      Contract(1.0,Environment::t[0][Lx-1],shape(1),Environment::I(0,Lx-1),shape(1),0.0,tmp5);

      R[Lx-2] = tmp5.reshape_clear(shape(Environment::t[0][Lx-1].shape(0),Environment::I(0,Lx-1).shape(0)));

      energy += 0.5 * Dot(LI,R[Lx-2]);

   }
   
   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< TArray<double,3> > RO(Lx - 1);

   //2 left renormalized operators needed
   TArray<double,3> LOI;
   TArray<double,3> LOU;

   for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operator

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::t[row][Lx - 1],shape(1),Environment::U(row,Lx-1),shape(1),0.0,tmp5);

      //then bottom enviroment
      TArray<double,6> tmp6;
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),Environment::U(row,Lx-1).shape(0),Environment::b[row-1][Lx - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int col = Lx-2;col > 0;--col){

         I4.clear();
         Contract(1.0,Environment::t[row][col],shape(2),RO[col],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::U(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-1]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) construct left renormalized operator with unity

      //paste top environment on local Sz
      tmp5.clear();
      Contract(1.0,Environment::t[row][0],shape(1),Environment::U(row,0),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOU = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),Environment::U(row,0).shape(3),Environment::b[row-1][0].shape(2)));

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[row*Lx] != (*this)[row*Lx + 1]){

         //paste top environment on local inverted Sz
         tmp5.clear();
         Contract(1.0,Environment::t[row][0],shape(1),Environment::I(row,0),shape(1),0.0,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(1.0,tmp5,shape(3),Environment::b[row-1][0],shape(1),0.0,tmp6);

         //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
         LOI = tmp6.reshape_clear(shape(Environment::t[row][0].shape(2),Environment::I(row,0).shape(3),Environment::b[row-1][0].shape(2)));

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int col = 1;col < Lx - 1;++col){

         enum {i,j,k,o,m,n};

         //1) close down LOI with I if it contributes
         if((*this)[row*Lx + col - 1] != (*this)[row*Lx + col]){

            //first add top to the right side, put it in I4
            I4.clear();
            Contract(1.0,Environment::t[row][col],shape(2),RO[col],shape(0),0.0,I4);

            I4bis.clear();
            Contract(1.0,I4,shape(i,j,k,o),Environment::I(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

            Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(1,2),0.0,RO[col-1]);

            //expectation value:
            energy += 0.5 * Dot(LOI,RO[col-1]);

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity
         I4.clear();
         Contract(1.0,Environment::t[row][col],shape(0),LOU,shape(0),0.0,I4);

         // 1) construct new left unity
         Contract(1.0,I4,shape(i,j,k,o),Environment::U(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOU.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOU);

         // 2) if it contributes, construct new left inverted 
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Contract(1.0,I4,shape(i,j,k,o),Environment::I(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

            LOI.clear();
            Contract(1.0,I4bis,shape(2,3),Environment::b[row-1][col],shape(0,1),0.0,LOI);

         }

      }

      //last site on the right: close down LOI if it contributes
      if((*this)[row*Lx + Lx - 2] != (*this)[row*Lx + Lx - 1]){

         //paste top environment on
         tmp5.clear();
         Contract(1.0,Environment::t[row][Lx - 1],shape(1),Environment::I(row,Lx-1),shape(1),0.0,tmp5);

         //then bottom enviroment
         Contract(1.0,tmp5,shape(3),Environment::b[row-1][Lx-1],shape(1),0.0,tmp6);

         //move to a DArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(Environment::t[row][Lx - 1].shape(0),Environment::I(row,Lx-1).shape(0),Environment::b[row-1][Lx - 1].shape(0)));

         //add to energy
         energy += 0.5 * Dot(LOI,RO[Lx - 2]);

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (t,b)
   tmp4.clear();
   Contract(1.0,Environment::t[Ly-2][Lx - 1],shape(1),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp4.reshape_clear(shape(Environment::t[Ly-2][Lx - 1].shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();
      Contract(1.0,Environment::t[Ly-2][col],shape(2),R[col],shape(0),0.0,tmp3);

      R[col - 1].clear();
      Contract(1.0,tmp3,shape(1,2),Environment::b[Ly-2][col],shape(1,2),0.0,R[col-1]);

   }
   
   //construct the left going operators on the first top site

   //unity
   tmp5.clear();
   Contract(1.0,Environment::U(Ly-1,0),shape(2),Environment::b[Ly-2][0],shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear(shape(Environment::U(Ly-1,0).shape(3),Environment::b[Ly-2][0].shape(2)));

   //inverse if it contributes
   if( (*this)[(Ly - 1)*Lx] != (*this)[(Ly - 1)*Lx + 1]){

      tmp5.clear();
      Contract(1.0,Environment::I(Ly-1,0),shape(2),Environment::b[Ly-2][0],shape(1),0.0,tmp5);

      LI = tmp5.reshape_clear(shape(Environment::I(Ly-1,0).shape(3),Environment::b[Ly-2][0].shape(2)));

   }

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the I term from the previous site for the energy
      if( (*this)[(Ly - 1)*Lx + col - 1] != (*this)[(Ly - 1)*Lx + col]){

         //construct the right intermediate contraction (paste bottom to right)
         tmp3.clear();
         Contract(1.0,Environment::b[Ly-2][col],shape(2),R[col],shape(1),0.0,tmp3);

         // 1) paste Sx to the right
         M = Environment::I(Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, Environment::I(Ly-1,col).data(),K,tmp3.data(),K,0.0,R[col-1].data(),N);

         //contract with left LI 
         energy += 0.5 * Dot(LI,R[col - 1]);

      }

      //construct left renormalized operators for next site: first paste bottom to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(1),Environment::b[Ly-2][col],shape(0),0.0,tmp3);

      // 1) construct new left unity operator
      LU.resize(Environment::U(Ly-1,col).shape(3),Environment::b[Ly-2][col].shape(2));

      M = Environment::U(Ly-1,col).shape(3);
      N = tmp3.shape(2);
      K = tmp3.shape(0) * tmp3.shape(1);

      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::U(Ly-1,col).data(),M,tmp3.data(),N,0.0,LU.data(),N);

      // 2) if it contributes, construct new left LI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         LI.resize(Environment::I(Ly-1,col).shape(3),Environment::b[Ly-2][col].shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::I(Ly-1,col).data(),M,tmp3.data(),N,0.0,LI.data(),N);

      }

   }
   
   //finally close down on last top site

   //first calculate overlap
   overlap = Dot(LU,R[Lx-2]);

   // close down last LI
   if( (*this)[(Ly - 1)*Lx + Lx - 2] != (*this)[(Ly - 1)*Lx + Lx - 1]){

      //tmp comes out index (t,b)
      tmp5.clear();
      Contract(1.0,Environment::I(Ly-1,Lx-1),shape(2),Environment::b[Ly-2][Lx - 1],shape(1),0.0,tmp5);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp5.reshape_clear(shape(Environment::I(Ly-1,Lx-1).shape(0),Environment::b[Ly-2][Lx - 1].shape(0)));

      //energy
      energy += 0.5 * Dot(LI,R[Lx-2]);

   }

   // #################################################################
   // ### ----      Horizontal Sz contribution is easy         ---- ### 
   // #################################################################

   int cnt = 0;

   for(int row = 0;row < Ly;++row){

      for(int col = 0;col < Lx - 1;++col){

         if( (*this)[row*Lx + col] != (*this)[row*Lx + col + 1] )
            cnt -= 1;
         else
            cnt += 1;

      }

   }

   energy += 0.25 * cnt * overlap;

/*

   // #################################################################
   // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
   // #################################################################

   //calculate the single layer contractions first:
   Environment::U[myID].fill('V',peps,*this);

   Environment::Sx[myID].fill('V',peps,Sx,*this);
   Environment::Sy[myID].fill('V',peps,Sy,*this);
   Environment::Sz[myID].fill('V',peps,Sz,*this);

   //first construct the top and bottom (horizontal) environment layers
   Environment::calc_env('V',peps,*this);

   // -- (1) -- || left column: similar to overlap calculation

   //first construct the right renormalized operators
   vector< DArray<2> > R(Ly - 1);

   //first the rightmost operator
   DArray<4> tmp4;
   DArray<3> tmp3;

   //tmp comes out index (t,b)
   Contract(1.0,Environment::r[myID][0][Ly - 1],shape(1),Environment::l[myID][0][Ly - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp4.reshape_clear(shape(Environment::r[myID][0][Ly - 1].shape(0),Environment::l[myID][0][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();
      Contract(1.0,Environment::r[myID][0][row],shape(2),R[row],shape(0),0.0,tmp3);

      Contract(1.0,tmp3,shape(1,2),Environment::l[myID][0][row],shape(1,2),0.0,R[row-1]);

   }

   //4 left going operators: Sx, Sy, Sz, and 1
   DArray<2> LSx;
   DArray<2> LSy;
   DArray<2> LSz;
   DArray<2> LU;

   TArray<double,5> tmp5;

   //tmp comes out index (r,l)
   Contract(1.0,Environment::r[myID][0][0],shape(1),Environment::Sx[myID](0,0),shape(1),0.0,tmp5);

   LSx = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sx[myID](0,0).shape(3)));

   //then Sy
   Contract(1.0,Environment::r[myID][0][0],shape(1),Environment::Sy[myID](0,0),shape(1),0.0,tmp5);

   LSy = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sy[myID](0,0).shape(3)));

   //then Sz
   Contract(1.0,Environment::r[myID][0][0],shape(1),Environment::Sz[myID](0,0),shape(1),0.0,tmp5);

   LSz = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

   //finally unity
   Contract(1.0,Environment::r[myID][0][0],shape(1),Environment::U[myID](0,0),shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear(shape(Environment::r[myID][0][0].shape(2),Environment::Sz[myID](0,0).shape(3)));

   int dim = R[0].size();

   //now contract x,y and z with R for local expectation values:
   auxvec[myID][0][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
   auxvec[myID][0][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
   auxvec[myID][0][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //first close down the x,y and z terms from the previous site for the energy

      //construct the right intermediate contraction (paste top to right)
      tmp3.clear();
      Contract(1.0,Environment::r[myID][0][row],shape(2),R[row],shape(0),0.0,tmp3);

      // 1) paste Sx to the right
      int M = tmp3.shape(0);
      int N = Environment::Sx[myID](row,0).shape(0);
      int K = tmp3.shape(1) * tmp3.shape(2);

      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,Environment::Sx[myID](row,0).data(),K,0.0,R[row-1].data(),N);

      //contract with left Sx
      energy += Dot(LSx,R[row - 1]);

      // 2) then paste Sy to the right
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,Environment::Sy[myID](row,0).data(),K,0.0,R[row-1].data(),N);

      //contract with left Sy
      energy += Dot(LSy,R[row - 1]);

      // 3) finally Sz
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,Environment::Sz[myID](row,0).data(),K,0.0,R[row-1].data(),N);

      //contract with left Sz
      energy += Dot(LSz,R[row - 1]);

      //construct left renormalized operators for next site: first paste top to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(0),Environment::r[myID][0][row],shape(0),0.0,tmp3);

      // 1) construct new Sx left operator
      LSx.resize(Environment::r[myID][0][row].shape(2),Environment::Sx[myID](row,0).shape(3));

      M = tmp3.shape(2);
      N = Environment::Sx[myID](row,0).shape(3);
      K = tmp3.shape(0) * tmp3.shape(1);

      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::Sx[myID](row,0).data(),N,0.0,LSx.data(),N);

      // 2) construct new Sy left operator
      LSy.resize(Environment::r[myID][0][row].shape(2),Environment::Sy[myID](row,0).shape(3));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::Sy[myID](row,0).data(),N,0.0,LSy.data(),N);

      // 3) construct new Sz left operator
      LSz.resize(Environment::r[myID][0][row].shape(2),Environment::Sz[myID](row,0).shape(3));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::Sz[myID](row,0).data(),N,0.0,LSz.data(),N);

      //now contract x,y and z with R for local expectation values:
      dim = R[row].size();

      auxvec[myID][row*Lx][0] = blas::dot(dim,LSx.data(),1,R[row].data(),1);
      auxvec[myID][row*Lx][1] = blas::dot(dim,LSy.data(),1,R[row].data(),1);
      auxvec[myID][row*Lx][2] = blas::dot(dim,LSz.data(),1,R[row].data(),1);

      // 4) finally construct new unity on the left
      LU.resize(Environment::r[myID][0][row].shape(2),Environment::U[myID](row,0).shape(3));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),M,Environment::U[myID](row,0).data(),N,0.0,LU.data(),N);

   }

   //last site of left column: close down the left x,y and z

   //1) Sx to close down LSx
   Contract(1.0,Environment::r[myID][0][Ly-1],shape(1),Environment::Sx[myID](Ly-1,0),shape(1),0.0,tmp5);

   R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sx[myID](Ly-1,0).shape(0)));

   energy += Dot(LSx,R[Ly-2]);

   dim = R[Ly-2].size();
   auxvec[myID][(Ly-1)*Lx][0] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   //2) Sy to close down LSy
   Contract(1.0,Environment::r[myID][0][Ly-1],shape(1),Environment::Sy[myID](Ly-1,0),shape(1),0.0,tmp5);

   R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sy[myID](Ly-1,0).shape(0)));

   energy += Dot(LSy,R[Ly-2]);

   auxvec[myID][(Ly-1)*Lx][1] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   //3) Sz to close down Lz
   Contract(1.0,Environment::r[myID][0][Ly-1],shape(1),Environment::Sz[myID](Ly-1,0),shape(1),0.0,tmp5);

   R[Ly-2] = tmp5.reshape_clear(shape(Environment::r[myID][0][Ly-1].shape(0),Environment::Sz[myID](Ly-1,0).shape(0)));

   energy += Dot(LSz,R[Ly-2]);

   auxvec[myID][(Ly-1)*Lx][2] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   // -- (2) -- now move from left to right calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< TArray<double,3> > RO(Ly - 1);

   //4 left renormalized operators needed
   TArray<double,3> LOSx;
   TArray<double,3> LOSy;
   TArray<double,3> LOSz;
   TArray<double,3> LOU;

   for(int col = 1;col < Lx - 1;++col){

      //first create right renormalized operator

      //paste right environment on
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][Ly - 1],shape(1),Environment::U[myID](Ly-1,col),shape(1),0.0,tmp5);

      //then left enviroment
      TArray<double,6> tmp6;
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::U[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int row = Ly-2;row > 0;--row){

         I4.clear();
         Contract(1.0,Environment::r[myID][col][row],shape(2),RO[row],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-1].clear();
         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),0.0,RO[row-1]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) Sx

      //paste top environment on local Sx
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][0],shape(1),Environment::Sx[myID](0,col),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOSx = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sx[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

      // 2) Sy

      //paste top environment on local Sy
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][0],shape(1),Environment::Sy[myID](0,col),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOSy = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sy[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

      // 3) Sz

      //paste top environment on local Sz
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][0],shape(1),Environment::Sz[myID](0,col),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOSz = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::Sz[myID](0,col).shape(3),Environment::r[myID][col-1][0].shape(2)));

      // 4) 1 -- finally construct left renormalized operator with unity

      //paste top environment on local Sz
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][0],shape(1),Environment::U[myID](0,col),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOU = tmp6.reshape_clear(shape(Environment::r[myID][col][0].shape(2),Environment::U[myID](0,col).shape(3),Environment::l[myID][col-1][0].shape(2)));

      //now contract x,y and z with R for local expectation values:
      dim = RO[0].size();

      auxvec[myID][col][0] = blas::dot(dim,LOSx.data(),1,RO[0].data(),1);
      auxvec[myID][col][1] = blas::dot(dim,LOSy.data(),1,RO[0].data(),1);
      auxvec[myID][col][2] = blas::dot(dim,LOSz.data(),1,RO[0].data(),1);

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int row = 1;row < Ly - 1;++row){

         //first add top to the right side, put it in I4
         I4.clear();
         Contract(1.0,Environment::r[myID][col][row],shape(2),RO[row],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         //1) close down LOSx with Sx
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         energy += Dot(LOSx,RO[row-1]);

         //2) close down LOSy with Sy
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         energy += Dot(LOSy,RO[row-1]);

         //3) finally close down LOSz with Sz
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(1,2),0.0,RO[row-1]);

         //expectation value:
         energy += Dot(LOSz,RO[row-1]);

         // now construct the new left going renormalized operators

         //first attach top to left unity
         I4.clear();
         Contract(1.0,Environment::r[myID][col][row],shape(0),LOU,shape(0),0.0,I4);

         // 1) construct new left Sx operator
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sx[myID](row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOSx.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),0.0,LOSx);

         // 2) construct new left Sy operator
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sy[myID](row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOSy.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),0.0,LOSy);

         // 3) construct new left Sz operator
         Contract(1.0,I4,shape(i,j,k,o),Environment::Sz[myID](row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOSz.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),0.0,LOSz);

         // 4) finally construct new left unity
         Contract(1.0,I4,shape(i,j,k,o),Environment::U[myID](row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOU.clear();
         Contract(1.0,I4bis,shape(2,3),Environment::l[myID][col-1][row],shape(0,1),0.0,LOU);

         //now contract x,y and z with R for local expectation values:
         dim = RO[row].size();

         auxvec[myID][row*Lx + col][0] = blas::dot(dim,LOSx.data(),1,RO[row].data(),1);
         auxvec[myID][row*Lx + col][1] = blas::dot(dim,LOSy.data(),1,RO[row].data(),1);
         auxvec[myID][row*Lx + col][2] = blas::dot(dim,LOSz.data(),1,RO[row].data(),1);

      }

      //last site on the right: close down on the incomings

      //1) first LSx with Sx

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sx[myID](Ly-1,col),shape(1),0.0,tmp5);

      //then bottom enviroment
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sx[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

      //add to energy
      energy += Dot(LOSx,RO[Ly - 2]);

      //local expectation value
      dim = RO[Ly-2].size();
      auxvec[myID][(Ly-1)*Lx + col][0] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

      //2) then LSy with Sy

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sy[myID](Ly-1,col),shape(1),0.0,tmp5);

      //then bottom enviroment
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sy[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

      //add to value
      energy += Dot(LOSy,RO[Ly - 2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + col][1] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

      //3) then LSz with Sz

      //paste top environment on
      tmp5.clear();
      Contract(1.0,Environment::r[myID][col][Ly - 1],shape(1),Environment::Sz[myID](Ly-1,col),shape(1),0.0,tmp5);

      //then bottom enviroment
      Contract(1.0,tmp5,shape(3),Environment::l[myID][col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 2] = tmp6.reshape_clear(shape(Environment::r[myID][col][Ly - 1].shape(0),Environment::Sz[myID](Ly-1,col).shape(0),Environment::l[myID][col-1][Ly - 1].shape(0)));

      //add to value
      energy += Dot(LOSz,RO[Ly - 2]);

      //local expectation value
      auxvec[myID][(Ly-1)*Lx + col][2] = blas::dot(dim,LOU.data(),1,RO[Ly-2].data(),1);

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (r,l)
   tmp4.clear();
   Contract(1.0,Environment::r[myID][Lx-2][Ly - 1],shape(1),Environment::l[myID][Lx-2][Ly - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp4.reshape_clear(shape(Environment::r[myID][Lx-2][Ly - 1].shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();
      Contract(1.0,Environment::r[myID][Lx-2][row],shape(2),R[row],shape(0),0.0,tmp3);

      R[row - 1].clear();
      Contract(1.0,tmp3,shape(1,2),Environment::l[myID][Lx-2][row],shape(1,2),0.0,R[row-1]);

   }

   //construct the left going operators on the first top site

   //first Sx

   //tmp comes out index (r,l)
   tmp5.clear();
   Contract(1.0,Environment::Sx[myID](0,Lx-1),shape(2),Environment::l[myID][Lx-2][0],shape(1),0.0,tmp5);

   LSx = tmp5.reshape_clear(shape(Environment::Sx[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

   //then Sy

   //tmp5 comes out index (r,l)
   Contract(1.0,Environment::Sy[myID](0,Lx-1),shape(2),Environment::l[myID][Lx-2][0],shape(1),0.0,tmp5);

   LSy = tmp5.reshape_clear(shape(Environment::Sy[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

   //then Sz 

   //tmp comes out index (r,l)
   Contract(1.0,Environment::Sz[myID](0,Ly-1),shape(2),Environment::l[myID][Ly-2][0],shape(1),0.0,tmp5);

   LSz = tmp5.reshape_clear(shape(Environment::Sz[myID](0,Lx-1).shape(3),Environment::l[myID][Lx-2][0].shape(2)));

   //and finally unity
   Contract(1.0,Environment::r[myID][Lx-2][0],shape(1),Environment::l[myID][Lx-2][0],shape(1),0.0,tmp4);

   LU = tmp4.reshape_clear(shape(Environment::r[myID][Lx-2][0].shape(2),Environment::l[myID][Lx-2][0].shape(2)));

   dim = R[0].size();

   //now contract x,y and z with R for local expectation values:
   auxvec[myID][Lx-1][0] = blas::dot(dim,LSx.data(),1,R[0].data(),1);
   auxvec[myID][Lx-1][1] = blas::dot(dim,LSy.data(),1,R[0].data(),1);
   auxvec[myID][Lx-1][2] = blas::dot(dim,LSz.data(),1,R[0].data(),1);

   //middle of the chain:
   for(int row = 1;row < Ly-1;++row){

      //first close down the x,y and z terms from the previous site for the energy

      //construct the right intermediate contraction (paste bottom to right)
      tmp3.clear();
      Contract(1.0,Environment::l[myID][Lx-2][row],shape(2),R[row],shape(1),0.0,tmp3);

      // 1) paste Sx to the right
      int M = Environment::Sx[myID](row,Lx-1).shape(0);
      int N = tmp3.shape(0);
      int K = tmp3.shape(1) * tmp3.shape(2);

      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, Environment::Sx[myID](row,Lx-1).data(),K,tmp3.data(),K,0.0,R[row-1].data(),N);

      //contract with left Sx
      energy += Dot(LSx,R[row - 1]);

      // 2) paste Sy to the right
      M = Environment::Sy[myID](row,Lx-1).shape(0);
      N = tmp3.shape(0);
      K = tmp3.shape(1) * tmp3.shape(2);

      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, Environment::Sy[myID](row,Lx-1).data(),K,tmp3.data(),K,0.0,R[row-1].data(),N);

      //contract with left Sy
      energy += Dot(LSy,R[row - 1]);

      // 3) paste Sz to the right
      M = Environment::Sz[myID](row,Lx-1).shape(0);
      N = tmp3.shape(0);
      K = tmp3.shape(1) * tmp3.shape(2);

      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, Environment::Sz[myID](row,Lx-1).data(),K,tmp3.data(),K,0.0,R[row-1].data(),N);

      //contract with left Sy
      energy += Dot(LSz,R[row - 1]);

      //construct left renormalized operators for next site: first paste bottom to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(1),Environment::l[myID][Lx-2][row],shape(0),0.0,tmp3);

      // 1) construct new Sx left operator
      LSx.resize(Environment::Sx[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));

      M = Environment::Sx[myID](row,Lx-1).shape(3);
      N = tmp3.shape(2);
      K = tmp3.shape(0) * tmp3.shape(1);

      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::Sx[myID](row,Lx-1).data(),M,tmp3.data(),N,0.0,LSx.data(),N);

      // 2) construct new Sy left operator
      LSy.resize(Environment::Sy[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::Sy[myID](row,Lx-1).data(),M,tmp3.data(),N,0.0,LSy.data(),N);

      // 3) construct new Sz left operator
      LSz.resize(Environment::Sz[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::Sz[myID](row,Lx-1).data(),M,tmp3.data(),N,0.0,LSz.data(),N);

      // 4) finally construct new unity on the left
      LU.resize(Environment::U[myID](row,Lx-1).shape(3),Environment::l[myID][Lx-2][row].shape(2));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0, Environment::r[myID][Lx-2][row].data(),M,tmp3.data(),N,0.0,LU.data(),N);

      //now contract x,y and z with R for local expectation values:
      dim = R[row].size();

      auxvec[myID][row*Lx + Lx-1][0] = blas::dot(dim,LSx.data(),1,R[row].data(),1);
      auxvec[myID][row*Lx + Lx-1][1] = blas::dot(dim,LSy.data(),1,R[row].data(),1);
      auxvec[myID][row*Lx + Lx-1][2] = blas::dot(dim,LSz.data(),1,R[row].data(),1);

   }

   //finally close down on last top site

   //first calculate overlap
   overlap = Dot(LU,R[Ly-2]);

   //1) Sx to close down LSx

   //tmp comes out index (r,l)
   tmp5.clear();
   Contract(1.0,Environment::Sx[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),0.0,tmp5);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sx[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

   //energy
   energy += Dot(LSx,R[Ly-2]);

   //local expectation value
   auxvec[myID][(Ly-1)*Lx + Lx-1][0] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   //2) Sy to close down Ly
   Contract(1.0,Environment::Sy[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),0.0,tmp5);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sy[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

   //energy
   energy += Dot(LSy,R[Ly-2]);

   //local expectation value
   auxvec[myID][(Ly-1)*Lx + Lx-1][1] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   //3) Sz to close down Lz

   //tmp comes out index (t,b)
   Contract(1.0,Environment::Sz[myID](Ly-1,Lx-1),shape(2),Environment::l[myID][Lx-2][Ly - 1],shape(1),0.0,tmp5);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp5.reshape_clear(shape(Environment::Sz[myID](Ly-1,Lx-1).shape(0),Environment::l[myID][Lx-2][Ly - 1].shape(0)));

   //energy
   energy += Dot(LSz,R[Ly-2]);

   //local expectation value
   auxvec[myID][(Ly-1)*Lx + Lx-1][2] = blas::dot(dim,LU.data(),1,R[Ly-2].data(),1);

   // #################################################################
   // ### ----               SET THE PROPERTIES                ---- ### 
   // #################################################################

   EL = energy/overlap;

   //set VL: Sx, Sy and Sz on every site: ready to construct auxiliary expectation values quickly!
   for(int k = 0;k < Trotter::n_trot;++k)
      for(int r = 0;r < 3;++r){

         VL[r*Trotter::n_trot + k] = auxvec[myID][0][r] * Trotter::gV()(k,0);


         for(int row = 0;row < Ly;++row)
            for(int col = 0;col < Lx;++col)
               VL[r*Trotter::n_trot + k] += auxvec[myID][row*Lx + col][r] * Trotter::gV()(k,row*Lx + col);

         VL[r*Trotter::n_trot + k] /= overlap;

      }
   */
}
