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
 * @param optoin start with up or down spin
 */
Walker::Walker(int option) : std::vector< bool >( Lx * Ly ){

   weight = 1.0;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         if( (r + c + option)%2 == 0)
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
   this->nn_over = walker.gnn_over();
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
 * @return the vector containing the overlaps of all the neighbouring walker states with the trial
 */
const vector<double> &Walker::gnn_over() const {

   return nn_over; 

}

/** 
 * @param index of the neighbour walker
 * @return the overlap of the walker with index 'index' with the trial
 */
double Walker::gnn_over(int index) const {

   return nn_over[index]; 

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
 * get the "potential" energy of the walker, which is < \sum_i Sz_i Sz_{i+1} >
 */
double Walker::pot_en() const {

   double tmp = 0.0;

   //first horizontal
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx - 1;++c){

         //Sz Sz
         if( (*this)[r*Lx + c] == (*this)[r*Lx + (c + 1)] )//up up or down down
            tmp += 0.25;
         else //up down or down up
            tmp -= 0.25;

      }

   }

   //then vertical
   for(int c = 0;c < Lx;++c){

      for(int r = 0;r < Ly - 1;++r){

         //Sz Sz
         if( (*this)[r*Lx + c] == (*this)[(r + 1)*Lx + c] )//up up or down down
            tmp += 0.25;
         else //up down or down up
            tmp -= 0.25;

      }

   }

   return tmp;

}

/**
 * calculate the local energy expectation value and overlap with the accesible states
 * @param peps trial wave function represented as peps
 */
void Walker::calc_EL(){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   // ---- || evaluate the expectation values in an MPO/MPS manner, first from bottom to top, then left to right || ----
   double ward;

   EL = 0.0;

   int M,N,K;

   nn_over.clear();

   //calculate the single layer contractions first:
   env[myID].sU(false,peps,*this);
   env[myID].sU(true,peps,*this);

   //first construct the top and bottom (horizontal) environment layers
   env[myID].calc('H',false);
   env[myID].calc('H',true);

   // #################################################################
   // ### ---- from bottom to top: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || bottom row: similar to overlap calculation

   //first construct the right renormalized operators: direct and inverse
   vector< DArray<2> > RU(Lx - 1);
   vector< DArray<2> > RI(Lx - 1);

   //first the rightmost operator
   DArray<4> tmp4;
   DArray<3> tmp3;

   //A: regular
   Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(false,0)[Lx - 1],env[myID].gb(false,0)[Lx - 1],0.0,RU[Lx - 2]);

   //B: inverse
   Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(true,0)[Lx - 1],env[myID].gb(true,0)[Lx - 1],0.0,RI[Lx - 2]);

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();

      //U
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,0)[col],RU[col],0.0,tmp3);
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(false,0)[col],0.0,RU[col - 1]);

      //I
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,0)[col],RI[col],0.0,tmp3);
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(true,0)[col],0.0,RI[col - 1]);

   }

   //4 left going operators: S+/- and 1
   DArray<2> LUU;
   DArray<2> LUI;

   DArray<2> LIU;
   DArray<2> LII;

   //U overlap
   Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,0)[0],env[myID].gb(false,0)[0],0.0,LUU);

   //I overlap
   Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,0)[0],env[myID].gb(true,0)[0],0.0,LII);

   //calculate the overlap with this state
   double tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   nn_over.push_back(1.0);

   //only calculate LUI and LIU if it contributes
   if( (*this)[0] != (*this)[1] ){

      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,0)[0],env[myID].gb(true,0)[0],0.0,LUI); //regular
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,0)[0],env[myID].gb(false,0)[0],0.0,LIU); //inverse

   }

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //only calculate if it contributes
      if( (*this)[col - 1] != (*this)[col] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,0)[col],RU[col],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(true,0)[col],0.0,RU[col - 1]);

         ward = Dot(LUI,RU[col - 1]);

         // B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,0)[col],RI[col],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(false,0)[col],0.0,RI[col - 1]);

         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward /tmp_over;

      }

      //construct left renormalized operators for next site:

      //A: regular 
      tmp3.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,LUU,env[myID].gt(false,0)[col],0.0,tmp3);

      //1) construct new unity on the left
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env[myID].gb(false,0)[col],0.0,LUU);

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env[myID].gb(true,0)[col],0.0,LUI);

      }

      //B: inverse 
      Gemm(CblasTrans,CblasNoTrans,1.0,LII,env[myID].gt(true,0)[col],0.0,tmp3);

      //1) construct new unity on the left
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env[myID].gb(true,0)[col],0.0,LII);

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env[myID].gb(false,0)[col],0.0,LIU);

      }

   }

   //last site of bottom row: close down LUI and LIU
   if((*this)[Lx - 2] != (*this)[Lx - 1]){

      //A: regular LUI
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(false,0)[Lx-1],env[myID].gb(true,0)[Lx-1],0.0,RU[Lx - 2]);

      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse LIU
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(true,0)[Lx-1],env[myID].gb(false,0)[Lx-1],0.0,RI[Lx - 2]);

      ward += Dot(LIU,RI[Lx-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward/tmp_over;

   }

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< TArray<double,3> > ROU(Lx - 1);
   vector< TArray<double,3> > ROI(Lx - 1);

   //4 left renormalized operators needed
   TArray<double,3> LOUU;
   TArray<double,3> LOUI;

   TArray<double,3> LOIU;
   TArray<double,3> LOII;

   DArray<3> tmp3bis;
   DArray<4> tmp4bis;

   for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operators

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(false,row - 1)[Lx-1],env[myID].gU(false)(row,Lx-1),0.0,tmp3);

      tmp3bis.clear();
      Permute(tmp3,shape(2,1,0),tmp3bis);

      M = env[myID].gt(false,row)[Lx - 1].shape(0);
      N = tmp3bis.shape(1) * tmp3bis.shape(2);
      K = tmp3bis.shape(0);

      ROU[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env[myID].gt(false,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROU[Lx - 2].data(),N);

      //B: inverse
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(true,row - 1)[Lx-1],env[myID].gU(true)(row,Lx-1),0.0,tmp3);

      Permute(tmp3,shape(2,1,0),tmp3bis);

      ROI[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env[myID].gt(true,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROI[Lx - 2].data(),N);

      //now construct the middle operators
      for(int col = Lx-2;col > 0;--col){

         //A: regular
         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(false,row - 1)[col],ROU[col],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gU(false)(row,col),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         ROU[col - 1].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,row)[col],tmp4bis,0.0,ROU[col - 1]);

         //B: inverse
         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(true,row - 1)[col],ROI[col],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gU(true)(row,col),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         ROI[col - 1].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,row)[col],tmp4bis,0.0,ROI[col - 1]);

      }

      // --- now move from left to right to get the expectation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) construct left renormalized operator with unity

      //A: regular
      tmp3.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,row)[0],env[myID].gU(false)(row,0),0.0,tmp3);

      tmp3bis.clear();
      Permute(tmp3,shape(0,2,1),tmp3bis);

      M = tmp3bis.shape(0) * tmp3bis.shape(1);
      N = env[myID].gb(false,row-1)[0].shape(2);
      K = tmp3bis.shape(2);

      LOUU.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env[myID].gb(false,row-1)[0].data(),N,0.0,LOUU.data(),N);

      //overlap A
      tmp_over = Dot(LOUU,ROU[0]);

      //B: inverse
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,row)[0],env[myID].gU(true)(row,0),0.0,tmp3);

      Permute(tmp3,shape(0,2,1),tmp3bis);

      //M,N and K are the same
      LOII.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env[myID].gb(true,row-1)[0].data(),N,0.0,LOII.data(),N);

      //overlap B
      tmp_over += Dot(LOII,ROI[0]);

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[row*Lx] != (*this)[row*Lx + 1]){

         //A: regular
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,row)[0],env[myID].gU(true)(row,0),0.0,tmp3);

         Permute(tmp3,shape(0,2,1),tmp3bis);

         //M,N and K are same as before
         LOUI.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env[myID].gb(false,row-1)[0].data(),N,0.0,LOUI.data(),N);

         //B: inverse
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,row)[0],env[myID].gU(false)(row,0),0.0,tmp3);

         Permute(tmp3,shape(0,2,1),tmp3bis);

         //M,N and K are the same
         LOIU.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env[myID].gb(true,row-1)[0].data(),N,0.0,LOIU.data(),N);

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int col = 1;col < Lx - 1;++col){

         //1) close down LO(U/I)I if it contributes
         if((*this)[row*Lx + col - 1] != (*this)[row*Lx + col]){

            //A: regular
            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(false,row - 1)[col],ROU[col],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gU(true)(row,col),tmp4bis,0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,row)[col],tmp4bis,0.0,ROU[col - 1]);

            //first part of expectation value
            ward = Dot(LOUI,ROU[col-1]);

            //B: inverse
            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(true,row - 1)[col],ROI[col],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gU(false)(row,col),tmp4bis,0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,row)[col],tmp4bis,0.0,ROI[col - 1]);

            //B expectation value
            ward += Dot(LOIU,ROI[col-1]);

            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new left going renormalized operators
         DArray<4> perm4;

         //A: regular
         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,row)[col],LOUU,0.0,tmp4);

         Permute(tmp4,shape(1,3,2,0),perm4);

         // 1) construct new LOUU
         tmp4bis.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,env[myID].gU(false)(row,col),0.0,tmp4bis);

         tmp4.clear();
         Permute(tmp4bis,shape(0,3,1,2),tmp4);

         LOUU.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env[myID].gb(false,row-1)[col],0.0,LOUU);

         // 2) if it contributes, construct new left inverted : LOUI
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,env[myID].gU(true)(row,col),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,3,1,2),tmp4);

            LOUI.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env[myID].gb(false,row-1)[col],0.0,LOUI);

         }

         //B: inverse
         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,row)[col],LOII,0.0,tmp4);

         Permute(tmp4,shape(1,3,2,0),perm4);

         //1) construct new LOII
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,env[myID].gU(true)(row,col),0.0,tmp4bis);

         Permute(tmp4bis,shape(0,3,1,2),tmp4);

         LOII.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env[myID].gb(true,row-1)[col],0.0,LOII);

         // 2) if it contributes, construct new left inverted : LOIU
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,env[myID].gU(false)(row,col),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,3,1,2),tmp4);

            LOIU.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env[myID].gb(true,row-1)[col],0.0,LOIU);

         }

      }

      //last site on the right: close down LOUI and LOIU if it contributes
      if((*this)[row*Lx + Lx - 2] != (*this)[row*Lx + Lx - 1]){

         //A: regular
         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(false,row - 1)[Lx-1],env[myID].gU(true)(row,Lx-1),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         M = env[myID].gt(false,row)[Lx - 1].shape(0);
         N = tmp3bis.shape(1) * tmp3bis.shape(2);
         K = tmp3bis.shape(0);

         ROU[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env[myID].gt(false,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROU[Lx - 2].data(),N);

         //expectation value A
         ward = Dot(LOUI,ROU[Lx - 2]);

         //B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gb(true,row - 1)[Lx-1],env[myID].gU(false)(row,Lx-1),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         ROI[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env[myID].gt(true,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROI[Lx - 2].data(),N);

         //expectation value A
         ward += Dot(LOIU,ROI[Lx - 2]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //A: regular
   Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(false,Ly-2)[Lx - 1],env[myID].gb(false,Ly-2)[Lx - 1],0.0,RU[Lx - 2]);

   //B: inverse
   Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(true,Ly-2)[Lx - 1],env[myID].gb(true,Ly-2)[Lx - 1],0.0,RI[Lx - 2]);

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[col],RU[col],0.0,tmp3);

      RU[col - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(false,Ly-2)[col],0.0,RU[col - 1]);

      //B: inverse
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[col],RI[col],0.0,tmp3);

      RI[col - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(true,Ly-2)[col],0.0,RI[col - 1]);

   }

   //construct the right going operators on the first top site

   //A: regular
   LUU.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[0],env[myID].gb(false,Ly-2)[0],0.0,LUU);

   //overlap part A
   tmp_over = Dot(LUU,RU[0]);

   //B: inverse
   LII.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[0],env[myID].gb(true,Ly-2)[0],0.0,LII);

   //overlap part B
   tmp_over += Dot(LII,RI[0]);

   //LUI and LIU if they contribute
   if( (*this)[(Ly - 1)*Lx] != (*this)[(Ly - 1)*Lx + 1]){

      //A: regular
      LUI.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[0],env[myID].gb(false,Ly-2)[0],0.0,LUI);

      //B: inverse
      LIU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[0],env[myID].gb(true,Ly-2)[0],0.0,LIU);

   }

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the I term from the previous site for the energy
      if( (*this)[(Ly - 1)*Lx + col - 1] != (*this)[(Ly - 1)*Lx + col]){

         //A: regular
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[col],RU[col],0.0,tmp3);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(false,Ly-2)[col],0.0,RU[col - 1]);

         //expectation A
         ward = Dot(LUI,RU[col - 1]);

         //B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[col],RI[col],0.0,tmp3);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env[myID].gb(true,Ly-2)[col],0.0,RI[col - 1]);

         //expectation B
         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

      //construct left renormalized operators for next site

      //A regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LUU,env[myID].gb(false,Ly-2)[col],0.0,tmp3);

      // 1) construct new LUU
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[col],tmp3,0.0,LUU);

      // 2) if it contributes, construct new left LUI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[col],tmp3,0.0,LUI);

      }

      //B inverse
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LII,env[myID].gb(true,Ly-2)[col],0.0,tmp3);

      // 1) construct new LII
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(true,Ly-2)[col],tmp3,0.0,LII);

      // 2) if it contributes, construct new left LUI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env[myID].gt(false,Ly-2)[col],tmp3,0.0,LIU);

      }

   }

   //finally close down on last top site

   // close down last LUI and LIU
   if( (*this)[(Ly - 1)*Lx + Lx - 2] != (*this)[(Ly - 1)*Lx + Lx - 1]){

      //A: regular
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(true,Ly-2)[Lx - 1],env[myID].gb(false,Ly-2)[Lx - 1],0.0,RU[Lx - 2]);

      //energy from A part
      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse
      Gemm(CblasNoTrans,CblasTrans,1.0,env[myID].gt(false,Ly-2)[Lx - 1],env[myID].gb(true,Ly-2)[Lx - 1],0.0,RI[Lx - 2]);

      //energy from B part
      ward += Dot(LIU,RI[Lx-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward/tmp_over;

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

   EL += 0.25 * cnt;

   cout << EL << endl;
/*
   // #################################################################
   // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
   // #################################################################

   //construct the left and right (vertical) environment layers
   env[myID].calc('V',false);
   env[myID].calc('V',true);

   // -- (1) -- || right column: similar to overlap calculation

   //A: regular
   Contract(1.0,env[myID].gl(false,Lx - 2)[Ly - 1],shape(1),env[myID].gr(false,Lx - 2)[Ly - 1],shape(1),0.0,tmp4);

   RU[Ly - 2] = tmp4.reshape_clear(shape(env[myID].gl(false,Lx - 2)[Ly - 1].shape(0),env[myID].gr(false,Lx - 2)[Ly - 1].shape(0)));

   //B: inverse
   Contract(1.0,env[myID].gl(true,Lx - 2)[Ly - 1],shape(1),env[myID].gr(true,Lx - 2)[Ly - 1],shape(1),0.0,tmp4);

   RI[Ly - 2] = tmp4.reshape_clear(shape(env[myID].gl(true,Lx - 2)[Ly - 1].shape(0),env[myID].gr(true,Lx - 2)[Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();

      //U
      Contract(1.0,env[myID].gl(false,Lx - 2)[row],shape(2),RU[row],shape(0),0.0,tmp3);
      Contract(1.0,tmp3,shape(1,2),env[myID].gr(false,Lx - 2)[row],shape(1,2),0.0,RU[row-1]);

      //I
      Contract(1.0,env[myID].gl(true,Lx - 2)[row],shape(2),RI[row],shape(0),0.0,tmp3);
      Contract(1.0,tmp3,shape(1,2),env[myID].gr(true,Lx - 2)[row],shape(1,2),0.0,RI[row-1]);

   }

   //4 left going operators: S+/- and 1

   //U overlap
   Contract(1.0,env[myID].gl(false,Lx - 2)[0],shape(1),env[myID].gU(false)(0,Lx - 1),shape(0),0.0,tmp5);
   LUU = tmp5.reshape_clear( shape(env[myID].gt(false,Lx - 2)[0].shape(2),env[myID].gU(false)(0,Lx - 1).shape(1)) );

   //I overlap
   Contract(1.0,env[myID].gl(true,Lx - 2)[0],shape(1),env[myID].gU(true)(0,Lx - 1),shape(0),0.0,tmp5);
   LII = tmp5.reshape_clear( shape(env[myID].gt(true,Lx - 2)[0].shape(2),env[myID].gU(true)(0,Lx - 1).shape(1)) );

   //calculate the overlap with this state
   tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   nn_over.push_back(1.0);

   //only calculate LUI and LIU if it contributes
   if( (*this)[Lx - 1] != (*this)[2*Lx - 1] ){

      //regular
      Contract(1.0,env[myID].gl(false,Lx - 2)[0],shape(1),env[myID].gU(true)(0,Lx - 1),shape(0),0.0,tmp5);
      LUI = tmp5.reshape_clear( shape(env[myID].gt(false,Lx - 2)[0].shape(2),env[myID].gU(true)(0,Lx - 1).shape(1)) );

      //inverse
      Contract(1.0,env[myID].gl(true,Lx - 2)[0],shape(1),env[myID].gU(false)(0,Lx - 1),shape(0),0.0,tmp5);
      LIU = tmp5.reshape_clear( shape(env[myID].gt(true,Lx - 2)[0].shape(2),env[myID].gU(false)(0,Lx - 1).shape(1)) );

   }

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //only calculate if it contributes
      if( (*this)[(row - 1)*Lx + Lx - 1] != (*this)[row*Lx + Lx - 1] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(1.0,env[myID].gl(false,Lx - 2)[col],shape(2),RU[row],shape(0),0.0,tmp3);

         // 1) paste I to the right
         M = tmp3.shape(0);
         N = env[myID].gU(true)(row,Lx - 1).shape(2);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),K,env[myID].gU(true)(row,Lx - 1).data(),N,0.0,RU[row-1].data(),N);

         ward = Dot(LUI,RU[row - 1]);

         // B: inverse

         //construct the right intermediate contraction (paste top to right)
         Contract(1.0,env[myID].gt(true,0)[col],shape(2),RI[row],shape(0),0.0,tmp3);

         // 1) paste U to the right
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, tmp3.data(),K,env[myID].gU(false)(row,Lx - 1).data(),N,0.0,RI[row-1].data(),N);

         ward += Dot(LIU,RI[row - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward / tmp_over;

      }

      //construct left renormalized operators for next site:

      //A: regular 
      tmp3.clear();
      Contract(1.0,LUU,shape(0),env[myID].gt(false,0)[col],shape(0),0.0,tmp3);

      //1) construct new unity on the left
      tmp3bis.clear();
      Contract(1.0,tmp3,shape(0,1),env[myID].gU(false)(0,col),shape(0,1),0.0,tmp3bis);

      LUU = tmp3bis.reshape_clear(shape(env[myID].gt(false,0)[col].shape(2),env[myID].gU(false)(0,col).shape(3)));

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         Contract(1.0,tmp3,shape(0,1),env[myID].gU(true)(0,col),shape(0,1),0.0,tmp3bis);

         LUI = tmp3bis.reshape_clear(shape(env[myID].gt(false,0)[col].shape(2),env[myID].gU(true)(0,col).shape(3)));

      }

      //B: inverse 
      Contract(1.0,LII,shape(0),env[myID].gt(true,0)[col],shape(0),0.0,tmp3);

      //1) construct new inverse on the left
      Contract(1.0,tmp3,shape(0,1),env[myID].gU(true)(0,col),shape(0,1),0.0,tmp3bis);

      LII = tmp3bis.reshape_clear(shape(env[myID].gt(true,0)[col].shape(2),env[myID].gU(true)(0,col).shape(3)));

      //2) if it contributes, calculate LIU
      if((*this)[col] != (*this)[col + 1]){

         Contract(1.0,tmp3,shape(0,1),env[myID].gU(false)(0,col),shape(0,1),0.0,tmp3bis);

         LIU = tmp3bis.reshape_clear(shape(env[myID].gt(true,0)[col].shape(2),env[myID].gU(false)(0,col).shape(3)));

      }

   }

   //last site of bottom row: close down LUI and LIU
   if((*this)[Lx - 2] != (*this)[Lx - 1]){

      //A: regular LUI
      Contract(1.0,env[myID].gt(false,0)[Lx-1],shape(1),env[myID].gU(true)(0,Lx-1),shape(1),0.0,tmp5);

      RU[Lx-2] = tmp5.reshape_clear(shape(env[myID].gt(false,0)[Lx-1].shape(0),env[myID].gU(true)(0,Lx-1).shape(0)));

      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse LIU
      Contract(1.0,env[myID].gt(true,0)[Lx-1],shape(1),env[myID].gU(false)(0,Lx-1),shape(1),0.0,tmp5);

      RI[Lx-2] = tmp5.reshape_clear(shape(env[myID].gt(true,0)[Lx-1].shape(0),env[myID].gU(false)(0,Lx-1).shape(0)));

      ward += Dot(LIU,RI[Lx-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward/tmp_over;

   }

   // -- (2) -- now move from bottom to top calculating everything like an MPO/MPS expectation value

   //Right renormalized operators
   vector< TArray<double,3> > ROU(Lx - 1);
   vector< TArray<double,3> > ROI(Lx - 1);

   //4 left renormalized operators needed
   TArray<double,3> LOUU;
   TArray<double,3> LOUI;

   TArray<double,3> LOIU;
   TArray<double,3> LOII;

   for(int row = 1;row < Ly - 1;++row){

      //first create right renormalized operators

      //A: regular
      tmp5.clear();
      Contract(1.0,env[myID].gt(false,row)[Lx - 1],shape(1),env[myID].gU(false)(row,Lx-1),shape(1),0.0,tmp5);

      TArray<double,6> tmp6;
      Contract(1.0,tmp5,shape(3),env[myID].gb(false,row-1)[Lx-1],shape(1),0.0,tmp6);

      ROU[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt(false,row)[Lx - 1].shape(0),env[myID].gU(false)(row,Lx-1).shape(0),env[myID].gb(false,row-1)[Lx - 1].shape(0)));

      //B: inverse
      Contract(1.0,env[myID].gt(true,row)[Lx - 1],shape(1),env[myID].gU(true)(row,Lx-1),shape(1),0.0,tmp5);

      Contract(1.0,tmp5,shape(3),env[myID].gb(true,row-1)[Lx-1],shape(1),0.0,tmp6);

      ROI[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt(true,row)[Lx - 1].shape(0),env[myID].gU(true)(row,Lx-1).shape(0),env[myID].gb(true,row-1)[Lx - 1].shape(0)));

      DArray<4> tmp4;
      DArray<4> tmp4bis;

      //now construct the middle operators
      for(int col = Lx-2;col > 0;--col){

         //A: regular
         tmp4.clear();
         Contract(1.0,env[myID].gt(false,row)[col],shape(2),ROU[col],shape(0),0.0,tmp4);

         tmp4bis.clear();
         Contract(1.0,tmp4,shape(1,2),env[myID].gU(false)(row,col),shape(1,3),0.0,tmp4bis);

         Contract(1.0,tmp4bis,shape(3,1),env[myID].gb(false,row-1)[col],shape(1,2),0.0,ROU[col-1]);

         //B: inverse
         Contract(1.0,env[myID].gt(true,row)[col],shape(2),ROI[col],shape(0),0.0,tmp4);

         Contract(1.0,tmp4,shape(1,2),env[myID].gU(true)(row,col),shape(1,3),0.0,tmp4bis);

         Contract(1.0,tmp4bis,shape(3,1),env[myID].gb(true,row-1)[col],shape(1,2),0.0,ROI[col-1]);

      }

      // --- now move from left to right to get the expectation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) construct left renormalized operator with unity

      //A: regular
      tmp5.clear();
      Contract(1.0,env[myID].gt(false,row)[0],shape(1),env[myID].gU(false)(row,0),shape(1),0.0,tmp5);

      tmp6.clear();
      Contract(1.0,tmp5,shape(3),env[myID].gb(false,row-1)[0],shape(1),0.0,tmp6);

      LOUU = tmp6.reshape_clear(shape(env[myID].gt(false,row)[0].shape(2),env[myID].gU(false)(row,0).shape(3),env[myID].gb(false,row-1)[0].shape(2)));

      //overlap A
      tmp_over = Dot(LOUU,ROU[0]);

      //B: inverse
      Contract(1.0,env[myID].gt(true,row)[0],shape(1),env[myID].gU(true)(row,0),shape(1),0.0,tmp5);

      Contract(1.0,tmp5,shape(3),env[myID].gb(true,row-1)[0],shape(1),0.0,tmp6);

      LOII = tmp6.reshape_clear(shape(env[myID].gt(true,row)[0].shape(2),env[myID].gU(true)(row,0).shape(3),env[myID].gb(true,row-1)[0].shape(2)));

      //overlap B
      tmp_over += Dot(LOII,ROI[0]);

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[row*Lx] != (*this)[row*Lx + 1]){

         //A: regular
         Contract(1.0,env[myID].gt(false,row)[0],shape(1),env[myID].gU(true)(row,0),shape(1),0.0,tmp5);

         Contract(1.0,tmp5,shape(3),env[myID].gb(false,row-1)[0],shape(1),0.0,tmp6);

         LOUI = tmp6.reshape_clear(shape(env[myID].gt(false,row)[0].shape(2),env[myID].gU(true)(row,0).shape(3),env[myID].gb(false,row-1)[0].shape(2)));

         //B: inverse
         Contract(1.0,env[myID].gt(true,row)[0],shape(1),env[myID].gU(false)(row,0),shape(1),0.0,tmp5);

         Contract(1.0,tmp5,shape(3),env[myID].gb(true,row-1)[0],shape(1),0.0,tmp6);

         LOIU = tmp6.reshape_clear(shape(env[myID].gt(true,row)[0].shape(2),env[myID].gU(false)(row,0).shape(3),env[myID].gb(true,row-1)[0].shape(2)));

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int col = 1;col < Lx - 1;++col){

         enum {i,j,k,o,m,n};

         //1) close down LO(U/I)I if it contributes
         if((*this)[row*Lx + col - 1] != (*this)[row*Lx + col]){

            //A: regular
            tmp4.clear();
            Contract(1.0,env[myID].gt(false,row)[col],shape(2),ROU[col],shape(0),0.0,tmp4);

            tmp4bis.clear();
            Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(m,j,n,k),0.0,tmp4bis,shape(i,m,n,o));

            Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(false,row-1)[col],shape(1,2),0.0,ROU[col-1]);

            //first part of expectation value
            ward = Dot(LOUI,ROU[col-1]);

            //B: inverse
            Contract(1.0,env[myID].gt(true,row)[col],shape(2),ROI[col],shape(0),0.0,tmp4);

            Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(m,j,n,k),0.0,tmp4bis,shape(i,m,n,o));

            Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(true,row-1)[col],shape(1,2),0.0,ROI[col-1]);

            //B expectation value
            ward += Dot(LOIU,ROI[col-1]);

            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new left going renormalized operators

         //A: regular
         tmp4.clear();
         Contract(1.0,env[myID].gt(false,row)[col],shape(0),LOUU,shape(0),0.0,tmp4);

         // 1) construct new LOUU
         tmp4bis.clear();
         Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(k,i,m,n),0.0,tmp4bis,shape(j,n,o,m));

         LOUU.clear();
         Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(false,row-1)[col],shape(0,1),0.0,LOUU);

         // 2) if it contributes, construct new left inverted : LOUI
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(k,i,m,n),0.0,tmp4bis,shape(j,n,o,m));

            LOUI.clear();
            Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(false,row-1)[col],shape(0,1),0.0,LOUI);

         }

         //B: inverse
         Contract(1.0,env[myID].gt(true,row)[col],shape(0),LOII,shape(0),0.0,tmp4);

         //1) construct new LOII
         Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(k,i,m,n),0.0,tmp4bis,shape(j,n,o,m));

         LOII.clear();
         Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(true,row-1)[col],shape(0,1),0.0,LOII);

         // 2) if it contributes, construct new left inverted : LOIU
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Contract(1.0,tmp4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(k,i,m,n),0.0,tmp4bis,shape(j,n,o,m));

            LOIU.clear();
            Contract(1.0,tmp4bis,shape(2,3),env[myID].gb(true,row-1)[col],shape(0,1),0.0,LOIU);

         }

      }

      //last site on the right: close down LOUI and LOIU if it contributes
      if((*this)[row*Lx + Lx - 2] != (*this)[row*Lx + Lx - 1]){

         //A: regular
         tmp5.clear();
         Contract(1.0,env[myID].gt(false,row)[Lx - 1],shape(1),env[myID].gU(true)(row,Lx-1),shape(1),0.0,tmp5);

         //then bottom enviroment
         Contract(1.0,tmp5,shape(3),env[myID].gb(false,row-1)[Lx-1],shape(1),0.0,tmp6);

         //move to a DArray<3> object
         ROU[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt(false,row)[Lx - 1].shape(0),env[myID].gU(true)(row,Lx-1).shape(0),env[myID].gb(false,row-1)[Lx - 1].shape(0)));

         //expectation value A
         ward = Dot(LOUI,ROU[Lx - 2]);

         //B: inverse
         Contract(1.0,env[myID].gt(true,row)[Lx - 1],shape(1),env[myID].gU(false)(row,Lx-1),shape(1),0.0,tmp5);

         //then bottom enviroment
         Contract(1.0,tmp5,shape(3),env[myID].gb(true,row-1)[Lx-1],shape(1),0.0,tmp6);

         //move to a DArray<3> object
         ROI[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt(true,row)[Lx - 1].shape(0),env[myID].gU(false)(row,Lx-1).shape(0),env[myID].gb(true,row-1)[Lx - 1].shape(0)));

         //expectation value A
         ward += Dot(LOIU,ROI[Lx - 2]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //A: regular
   tmp4.clear();
   Contract(1.0,env[myID].gt(false,Ly-2)[Lx - 1],shape(1),env[myID].gb(false,Ly-2)[Lx - 1],shape(1),0.0,tmp4);

   RU[Lx - 2] = tmp4.reshape_clear(shape(env[myID].gt(false,Ly-2)[Lx - 1].shape(0),env[myID].gb(false,Ly-2)[Lx - 1].shape(0)));

   //B: inverse
   Contract(1.0,env[myID].gt(true,Ly-2)[Lx - 1],shape(1),env[myID].gb(true,Ly-2)[Lx - 1],shape(1),0.0,tmp4);

   RI[Lx - 2] = tmp4.reshape_clear(shape(env[myID].gt(true,Ly-2)[Lx - 1].shape(0),env[myID].gb(true,Ly-2)[Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      //A: regular
      tmp3.clear();
      Contract(1.0,env[myID].gt(false,Ly-2)[col],shape(2),RU[col],shape(0),0.0,tmp3);

      RU[col - 1].clear();
      Contract(1.0,tmp3,shape(1,2),env[myID].gb(false,Ly-2)[col],shape(1,2),0.0,RU[col-1]);

      //B: regular
      Contract(1.0,env[myID].gt(true,Ly-2)[col],shape(2),RI[col],shape(0),0.0,tmp3);

      RI[col - 1].clear();
      Contract(1.0,tmp3,shape(1,2),env[myID].gb(true,Ly-2)[col],shape(1,2),0.0,RI[col-1]);

   }

   //construct the left going operators on the first top site

   //A: regular
   tmp5.clear();
   Contract(1.0,env[myID].gU(false)(Ly-1,0),shape(2),env[myID].gb(false,Ly-2)[0],shape(1),0.0,tmp5);

   LUU = tmp5.reshape_clear(shape(env[myID].gU(false)(Ly-1,0).shape(3),env[myID].gb(false,Ly-2)[0].shape(2)));

   //overlap part A
   tmp_over = Dot(LUU,RU[0]);

   //B: inverse
   tmp5.clear();
   Contract(1.0,env[myID].gU(true)(Ly-1,0),shape(2),env[myID].gb(true,Ly-2)[0],shape(1),0.0,tmp5);

   LII = tmp5.reshape_clear(shape(env[myID].gU(false)(Ly-1,0).shape(3),env[myID].gb(false,Ly-2)[0].shape(2)));

   //overlap part B
   tmp_over += Dot(LII,RI[0]);

   //LUI and LIU if they contribute
   if( (*this)[(Ly - 1)*Lx] != (*this)[(Ly - 1)*Lx + 1]){

      //A: regular
      Contract(1.0,env[myID].gU(true)(Ly-1,0),shape(2),env[myID].gb(false,Ly-2)[0],shape(1),0.0,tmp5);

      LUI = tmp5.reshape_clear(shape(env[myID].gU(true)(Ly-1,0).shape(3),env[myID].gb(false,Ly-2)[0].shape(2)));

      //B: inverse
      Contract(1.0,env[myID].gU(false)(Ly-1,0),shape(2),env[myID].gb(true,Ly-2)[0],shape(1),0.0,tmp5);

      LIU = tmp5.reshape_clear(shape(env[myID].gU(false)(Ly-1,0).shape(3),env[myID].gb(true,Ly-2)[0].shape(2)));

   }

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the I term from the previous site for the energy
      if( (*this)[(Ly - 1)*Lx + col - 1] != (*this)[(Ly - 1)*Lx + col]){

         //A: regular
         tmp3.clear();
         Contract(1.0,env[myID].gb(false,Ly-2)[col],shape(2),RU[col],shape(1),0.0,tmp3);

         // 1) paste I to the right
         M = env[myID].gU(true)(Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, env[myID].gU(true)(Ly-1,col).data(),K,tmp3.data(),K,0.0,RU[col-1].data(),N);

         //expectation A
         ward = Dot(LUI,RU[col - 1]);

         //B: inverse
         Contract(1.0,env[myID].gb(true,Ly-2)[col],shape(2),RI[col],shape(1),0.0,tmp3);

         // 1) paste I to the right
         M = env[myID].gU(false)(Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, env[myID].gU(false)(Ly-1,col).data(),K,tmp3.data(),K,0.0,RI[col-1].data(),N);

         //expectation B
         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

      //construct left renormalized operators for next site

      //A regular
      tmp3.clear();
      Contract(1.0,LUU,shape(1),env[myID].gb(false,Ly-2)[col],shape(0),0.0,tmp3);

      // 1) construct new LII
      tmp3bis.clear();
      Contract(1.0,env[myID].gU(false)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

      LUU = tmp3bis.reshape_clear(shape(env[myID].gU(false)(Ly-1,col).shape(3),env[myID].gb(false,Ly-2)[col].shape(2)));

      // 2) if it contributes, construct new left LUI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         Contract(1.0,env[myID].gU(true)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

         LUI = tmp3bis.reshape_clear(shape(env[myID].gU(true)(Ly-1,col).shape(3),env[myID].gb(false,Ly-2)[col].shape(2)));

      }

      //B inverse
      tmp3.clear();
      Contract(1.0,LII,shape(1),env[myID].gb(true,Ly-2)[col],shape(0),0.0,tmp3);

      // 1) construct new LII
      Contract(1.0,env[myID].gU(true)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

      LII = tmp3bis.reshape_clear(shape(env[myID].gU(true)(Ly-1,col).shape(3),env[myID].gb(true,Ly-2)[col].shape(2)));

      // 2) if it contributes, construct new left LIU
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         Contract(1.0,env[myID].gU(false)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

         LIU = tmp3bis.reshape_clear(shape(env[myID].gU(false)(Ly-1,col).shape(3),env[myID].gb(true,Ly-2)[col].shape(2)));

      }
   }

   //finally close down on last top site

   // close down last LUI and LIU
   if( (*this)[(Ly - 1)*Lx + Lx - 2] != (*this)[(Ly - 1)*Lx + Lx - 1]){

      //A: regular
      tmp5.clear();
      Contract(1.0,env[myID].gU(true)(Ly-1,Lx-1),shape(2),env[myID].gb(false,Ly-2)[Lx - 1],shape(1),0.0,tmp5);

      //reshape tmp to a 2-index array
      RU[Lx - 2] = tmp5.reshape_clear(shape(env[myID].gU(true)(Ly-1,Lx-1).shape(0),env[myID].gb(false,Ly-2)[Lx - 1].shape(0)));

      //energy from A part
      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse
      tmp5.clear();
      Contract(1.0,env[myID].gU(false)(Ly-1,Lx-1),shape(2),env[myID].gb(true,Ly-2)[Lx - 1],shape(1),0.0,tmp5);

      //reshape tmp to a 2-index array
      RI[Lx - 2] = tmp5.reshape_clear(shape(env[myID].gU(false)(Ly-1,Lx-1).shape(0),env[myID].gb(true,Ly-2)[Lx - 1].shape(0)));

      //energy from B part
      ward += Dot(LIU,RI[Lx-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward/tmp_over;

   }


   // #################################################################
   // ### ----          Vertical Sz contribution is easy       ---- ### 
   // #################################################################

   cnt = 0;

   for(int col = 0;col < Lx;++col){

      for(int row = 0;row < Ly - 1;++row){

         if( (*this)[row*Lx + col] != (*this)[(row + 1)*Lx + col] )
            cnt -= 1;
         else
            cnt += 1;

      }

   }

   EL += 0.25 * cnt;
   */
}
