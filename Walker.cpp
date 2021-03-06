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
 * @param option start with up or down spin
 */
Walker::Walker(int option) : std::vector< bool >( Lx * Ly ){

   weight = 1.0;
   sign = 1;

   for(int r = 0;r < Ly;++r)
      for(int c = 0;c < Lx;++c){

         if( (r + c + option)%2 == 0)
            (*this)[ r*Lx + c ] = true;
         else
            (*this)[ r*Lx + c ] = false;

      }

   env = Environment(DT,D_aux,1);

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< bool >(walker) {

   this->weight = walker.gWeight();
   this->nn_over = walker.gnn_over();
   this->EL = walker.gEL();
   this->sign = walker.gsign();
   this->env = walker.genv();

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
 * @return the weight corresponding to the walker
 */
const Environment &Walker::genv() const {

   return env; 

}

/**
 * flip the sign of the walker
 */
void Walker::sign_flip(){

   sign *= -1;

}

/** 
 * @return the sign of the walker
 */
int Walker::gsign() const{

   return sign; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set the walker configuration, bool vector
 */
void Walker::sconf(const vector<bool> &conf){

   this->std::vector<bool>::operator=(conf);

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
 * make sure the environment is calculated
 */
void Walker::calc_EL(){

   // ---- || evaluate the expectation values in an MPO/MPS manner, first from bottom to top, then left to right || ----
   double ward;

   EL = 0.0;

   int M,N,K;

   nn_over.clear();

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

   DArray<4> perm4;

   //A: regular
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(false,0)[Lx - 1],env.gb(false,0)[Lx - 1],0.0,RU[Lx - 2]);

   //B: inverse
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(true,0)[Lx - 1],env.gb(true,0)[Lx - 1],0.0,RI[Lx - 2]);

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();

      //U
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,0)[col],RU[col],0.0,tmp3);
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(false,0)[col],0.0,RU[col - 1]);

      //I
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,0)[col],RI[col],0.0,tmp3);
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(true,0)[col],0.0,RI[col - 1]);

   }

   //4 left going operators: S+/- and 1
   DArray<2> LUU;
   DArray<2> LUI;

   DArray<2> LIU;
   DArray<2> LII;

   //U overlap
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,0)[0],env.gb(false,0)[0],0.0,LUU);

   //I overlap
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,0)[0],env.gb(true,0)[0],0.0,LII);

   //calculate the overlap with this state
   double tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   nn_over.push_back(1.0);

   //only calculate LUI and LIU if it contributes
   if( (*this)[0] != (*this)[1] ){

      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,0)[0],env.gb(true,0)[0],0.0,LUI); //regular
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,0)[0],env.gb(false,0)[0],0.0,LIU); //inverse

   }

   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //only calculate if it contributes
      if( (*this)[col - 1] != (*this)[col] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,0)[col],RU[col],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(true,0)[col],0.0,RU[col - 1]);

         ward = Dot(LUI,RU[col - 1]);

         // B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,0)[col],RI[col],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(false,0)[col],0.0,RI[col - 1]);

         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward /tmp_over;

      }

      //construct left renormalized operators for next site:

      //A: regular 
      tmp3.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,LUU,env.gt(false,0)[col],0.0,tmp3);

      //1) construct new unity on the left
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gb(false,0)[col],0.0,LUU);

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gb(true,0)[col],0.0,LUI);

      }

      //B: inverse 
      Gemm(CblasTrans,CblasNoTrans,1.0,LII,env.gt(true,0)[col],0.0,tmp3);

      //1) construct new unity on the left
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gb(true,0)[col],0.0,LII);

      //2) if it contributes, calculate inverse on the left
      if((*this)[col] != (*this)[col + 1]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gb(false,0)[col],0.0,LIU);

      }

   }

   //last site of bottom row: close down LUI and LIU
   if((*this)[Lx - 2] != (*this)[Lx - 1]){

      //A: regular LUI
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(false,0)[Lx-1],env.gb(true,0)[Lx-1],0.0,RU[Lx - 2]);

      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse LIU
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(true,0)[Lx-1],env.gb(false,0)[Lx-1],0.0,RI[Lx - 2]);

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
      bool s = (*this)[row*Lx + Lx - 1];

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(false,row - 1)[Lx-1],peps[0](row,Lx-1,s),0.0,tmp3);

      tmp3bis.clear();
      Permute(tmp3,shape(2,1,0),tmp3bis);

      M = env.gt(false,row)[Lx - 1].shape(0);
      N = tmp3bis.shape(1) * tmp3bis.shape(2);
      K = tmp3bis.shape(0);

      ROU[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(false,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROU[Lx - 2].data(),N);

      //B: inverse
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(true,row - 1)[Lx-1],peps[0](row,Lx-1,!s),0.0,tmp3);

      Permute(tmp3,shape(2,1,0),tmp3bis);

      ROI[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(true,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROI[Lx - 2].data(),N);

      //now construct the middle operators
      for(int col = Lx-2;col > 0;--col){

         s = (*this)[row*Lx + col];

         //A: regular
         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(false,row - 1)[col],ROU[col],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,s),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         ROU[col - 1].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,row)[col],tmp4bis,0.0,ROU[col - 1]);

         //B: inverse
         tmp4.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(true,row - 1)[col],ROI[col],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,!s),tmp4bis,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(1,3,0,2),tmp4bis);

         ROI[col - 1].clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,row)[col],tmp4bis,0.0,ROI[col - 1]);

      }

      // --- now move from left to right to get the expectation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) construct left renormalized operator with unity

      s = (*this)[row*Lx];

      //A: regular
      tmp3.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,row)[0],peps[0](row,0,s),0.0,tmp3);

      tmp3bis.clear();
      Permute(tmp3,shape(0,2,1),tmp3bis);

      M = tmp3bis.shape(0) * tmp3bis.shape(1);
      N = env.gb(false,row-1)[0].shape(2);
      K = tmp3bis.shape(2);

      LOUU.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env.gb(false,row-1)[0].data(),N,0.0,LOUU.data(),N);

      //overlap A
      tmp_over = Dot(LOUU,ROU[0]);

      //B: inverse
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,row)[0],peps[0](row,0,!s),0.0,tmp3);

      Permute(tmp3,shape(0,2,1),tmp3bis);

      //M,N and K are the same
      LOII.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env.gb(true,row-1)[0].data(),N,0.0,LOII.data(),N);

      //overlap B
      tmp_over += Dot(LOII,ROI[0]);

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[row*Lx] != (*this)[row*Lx + 1]){

         //A: regular
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,row)[0],peps[0](row,0,!s),0.0,tmp3);

         Permute(tmp3,shape(0,2,1),tmp3bis);

         //M,N and K are same as before
         LOUI.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env.gb(false,row-1)[0].data(),N,0.0,LOUI.data(),N);

         //B: inverse
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,row)[0],peps[0](row,0,s),0.0,tmp3);

         Permute(tmp3,shape(0,2,1),tmp3bis);

         //M,N and K are the same
         LOIU.resize( shape(tmp3bis.shape(0),tmp3bis.shape(1),N) );
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,tmp3bis.data(),K,env.gb(true,row-1)[0].data(),N,0.0,LOIU.data(),N);

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int col = 1;col < Lx - 1;++col){

         s = (*this)[row*Lx + col];

         //1) close down LO(U/I)I if it contributes
         if((*this)[row*Lx + col - 1] != (*this)[row*Lx + col]){

            //A: regular
            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(false,row - 1)[col],ROU[col],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,!s),tmp4bis,0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,row)[col],tmp4bis,0.0,ROU[col - 1]);

            //first part of expectation value
            ward = Dot(LOUI,ROU[col-1]);

            //B: inverse
            tmp4.clear();
            Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(true,row - 1)[col],ROI[col],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,s),tmp4bis,0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(1,3,0,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,row)[col],tmp4bis,0.0,ROI[col - 1]);

            //B expectation value
            ward += Dot(LOIU,ROI[col-1]);

            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new left going renormalized operators

         //A: regular
         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,row)[col],LOUU,0.0,tmp4);

         perm4.clear();
         Permute(tmp4,shape(1,3,2,0),perm4);

         // 1) construct new LOUU
         tmp4bis.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,peps[0](row,col,s),0.0,tmp4bis);

         tmp4.clear();
         Permute(tmp4bis,shape(0,3,1,2),tmp4);

         LOUU.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env.gb(false,row-1)[col],0.0,LOUU);

         // 2) if it contributes, construct new left inverted : LOUI
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,peps[0](row,col,!s),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,3,1,2),tmp4);

            LOUI.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env.gb(false,row-1)[col],0.0,LOUI);

         }

         //B: inverse
         tmp4.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,row)[col],LOII,0.0,tmp4);

         Permute(tmp4,shape(1,3,2,0),perm4);

         //1) construct new LOII
         Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,peps[0](row,col,!s),0.0,tmp4bis);

         Permute(tmp4bis,shape(0,3,1,2),tmp4);

         LOII.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env.gb(true,row-1)[col],0.0,LOII);

         // 2) if it contributes, construct new left inverted : LOIU
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,perm4,peps[0](row,col,s),0.0,tmp4bis);

            Permute(tmp4bis,shape(0,3,1,2),tmp4);

            LOIU.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4,env.gb(true,row-1)[col],0.0,LOIU);

         }

      }

      s = (*this)[row*Lx + Lx - 1];

      //last site on the right: close down LOUI and LOIU if it contributes
      if((*this)[row*Lx + Lx - 2] != (*this)[row*Lx + Lx - 1]){

         //A: regular
         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(false,row - 1)[Lx-1],peps[0](row,Lx-1,!s),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         M = env.gt(false,row)[Lx - 1].shape(0);
         N = tmp3bis.shape(1) * tmp3bis.shape(2);
         K = tmp3bis.shape(0);

         ROU[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(false,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROU[Lx - 2].data(),N);

         //expectation value A
         ward = Dot(LOUI,ROU[Lx - 2]);

         //B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,env.gb(true,row - 1)[Lx-1],peps[0](row,Lx-1,s),0.0,tmp3);

         tmp3bis.clear();
         Permute(tmp3,shape(2,1,0),tmp3bis);

         ROI[Lx - 2].resize(M,tmp3bis.shape(1),tmp3bis.shape(2));
         blas::gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0,env.gt(true,row)[Lx - 1].data(),K,tmp3bis.data(),N,0.0,ROI[Lx - 2].data(),N);

         //expectation value A
         ward += Dot(LOIU,ROI[Lx - 2]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //A: regular
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(false,Ly-2)[Lx - 1],env.gb(false,Ly-2)[Lx - 1],0.0,RU[Lx - 2]);

   //B: inverse
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(true,Ly-2)[Lx - 1],env.gb(true,Ly-2)[Lx - 1],0.0,RI[Lx - 2]);

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[col],RU[col],0.0,tmp3);

      RU[col - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(false,Ly-2)[col],0.0,RU[col - 1]);

      //B: inverse
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[col],RI[col],0.0,tmp3);

      RI[col - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(true,Ly-2)[col],0.0,RI[col - 1]);

   }

   //construct the right going operators on the first top site

   //A: regular
   LUU.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[0],env.gb(false,Ly-2)[0],0.0,LUU);

   //overlap part A
   tmp_over = Dot(LUU,RU[0]);

   //B: inverse
   LII.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[0],env.gb(true,Ly-2)[0],0.0,LII);

   //overlap part B
   tmp_over += Dot(LII,RI[0]);

   //LUI and LIU if they contribute
   if( (*this)[(Ly - 1)*Lx] != (*this)[(Ly - 1)*Lx + 1]){

      //A: regular
      LUI.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[0],env.gb(false,Ly-2)[0],0.0,LUI);

      //B: inverse
      LIU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[0],env.gb(true,Ly-2)[0],0.0,LIU);

   }

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the I term from the previous site for the energy
      if( (*this)[(Ly - 1)*Lx + col - 1] != (*this)[(Ly - 1)*Lx + col]){

         //A: regular
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[col],RU[col],0.0,tmp3);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(false,Ly-2)[col],0.0,RU[col - 1]);

         //expectation A
         ward = Dot(LUI,RU[col - 1]);

         //B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[col],RI[col],0.0,tmp3);

         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gb(true,Ly-2)[col],0.0,RI[col - 1]);

         //expectation B
         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

      //construct left renormalized operators for next site

      //A regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LUU,env.gb(false,Ly-2)[col],0.0,tmp3);

      // 1) construct new LUU
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[col],tmp3,0.0,LUU);

      // 2) if it contributes, construct new left LUI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[col],tmp3,0.0,LUI);

      }

      //B inverse
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LII,env.gb(true,Ly-2)[col],0.0,tmp3);

      // 1) construct new LII
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(true,Ly-2)[col],tmp3,0.0,LII);

      // 2) if it contributes, construct new left LUI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gt(false,Ly-2)[col],tmp3,0.0,LIU);

      }

   }

   //finally close down on last top site

   // close down last LUI and LIU
   if( (*this)[(Ly - 1)*Lx + Lx - 2] != (*this)[(Ly - 1)*Lx + Lx - 1]){

      //A: regular
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(true,Ly-2)[Lx - 1],env.gb(false,Ly-2)[Lx - 1],0.0,RU[Lx - 2]);

      //energy from A part
      ward = Dot(LUI,RU[Lx-2]);

      //B: inverse
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gt(false,Ly-2)[Lx - 1],env.gb(true,Ly-2)[Lx - 1],0.0,RI[Lx - 2]);

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

   // #################################################################
   // ### ---- from left to right: contract in mps/mpo fashion ---- ### 
   // #################################################################

   // -- (1) -- || right column: similar to overlap calculation

   //first the rightmost operator

   //A: regular
   RU[Ly - 2].clear();
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(false,Lx - 2)[Ly - 1],env.gr(false,Lx - 2)[Ly - 1],0.0,RU[Ly - 2]);

   //B: inverse
   RI[Ly - 2].clear();
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(true,Lx - 2)[Ly - 1],env.gr(true,Lx - 2)[Ly - 1],0.0,RI[Ly - 2]);

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();

      //U
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,Lx - 2)[row],RU[row],0.0,tmp3);

      RU[row  - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(false,Lx - 2)[row],0.0,RU[row - 1]);

      //I
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,Lx - 2)[row],RI[row],0.0,tmp3);

      RI[row  - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(true,Lx - 2)[row],0.0,RI[row - 1]);

   }

   //construct 4 right going operators

   //U overlap
   LUU.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,Lx - 2)[0],env.gr(false,Lx - 2)[0],0.0,LUU);

   //I overlap
   LII.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,Lx - 2)[0],env.gr(true,Lx - 2)[0],0.0,LII);

   //calculate the overlap with this state
   tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   //only calculate LUI and LIU if it contributes
   if( (*this)[Lx - 1] != (*this)[2*Lx - 1] ){

      LUI.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,Lx - 2)[0],env.gr(true,Lx - 2)[0],0.0,LUI); //regular

      LIU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,Lx - 2)[0],env.gr(false,Lx - 2)[0],0.0,LIU); //inverse

   }

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //only calculate if it contributes
      if( (*this)[(row - 1)*Lx + Lx - 1] != (*this)[row*Lx + Lx - 1] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,Lx - 2)[row],RU[row],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(true,Lx - 2)[row],0.0,RU[row - 1]);

         ward = Dot(LUI,RU[row - 1]);

         // B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,Lx - 2)[row],RI[row],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(false,Lx - 2)[row],0.0,RI[row - 1]);

         ward += Dot(LIU,RI[row - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward /tmp_over;

      }

      //construct left renormalized operators for next site:

      //A: regular 
      tmp3.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,LUU,env.gl(false,Lx - 2)[row],0.0,tmp3);

      //1) construct new unity on the left
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gr(false,Lx - 2)[row],0.0,LUU);

      //2) if it contributes, calculate inverse on the left
      if((*this)[row*Lx + Lx - 1] != (*this)[(row + 1)*Lx + Lx - 1]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gr(true,Lx - 2)[row],0.0,LUI);

      }

      //B: inverse 
      Gemm(CblasTrans,CblasNoTrans,1.0,LII,env.gl(true,Lx - 2)[row],0.0,tmp3);

      //1) construct new unity on the left
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gr(true,Lx - 2)[row],0.0,LII);

      //2) if it contributes, calculate inverse on the left
      if((*this)[row*Lx + Lx - 1] != (*this)[(row + 1)*Lx + Lx - 1]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,tmp3,env.gr(false,Lx - 2)[row],0.0,LIU);

      }

   }

   //last site of right col: close down LUI and LIU
   if((*this)[(Ly - 2)*Lx + Lx - 1] != (*this)[Lx*Ly - 1]){

      //A: regular LUI
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(false,Lx - 2)[Ly - 1],env.gr(true,Lx - 2)[Ly - 1],0.0,RU[Ly - 2]);

      ward = Dot(LUI,RU[Ly-2]);

      //B: inverse LIU
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(true,Lx - 2)[Ly - 1],env.gr(false,Lx - 2)[Ly - 1],0.0,RI[Ly - 2]);

      ward += Dot(LIU,RI[Ly-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward/tmp_over;

   }

   // -- (2) -- now move from right to left calculating everything like an MPO/MPS expectation value

   for(int col = Lx - 2;col > 0;--col){

      //first create right renormalized operators
      bool s = (*this)[(Ly - 1)*Lx + col];

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[Ly - 1],peps[0](Ly-1,col,s),0.0,tmp3);

      M = tmp3.shape(0) * tmp3.shape(1);
      N = env.gr(false,col)[Ly - 1].shape(0);
      K = tmp3.shape(2);

      ROU[Ly - 2].resize(tmp3.shape(0),tmp3.shape(1),env.gr(false,col - 1)[Ly - 1].shape(0));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3.data(),K,env.gr(false,col)[Ly - 1].data(),K,0.0,ROU[Ly - 2].data(),N);

      //B: inverse
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[Ly - 1],peps[0](Ly-1,col,!s),0.0,tmp3);

      M = tmp3.shape(0) * tmp3.shape(1);
      N = env.gr(true,col)[Ly - 1].shape(0);
      K = tmp3.shape(2);

      ROI[Ly - 2].resize(tmp3.shape(0),tmp3.shape(1),env.gr(true,col)[Ly - 1].shape(0));
      blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3.data(),K,env.gr(true,col)[Ly - 1].data(),K,0.0,ROI[Ly - 2].data(),N);

      //now construct the middle operators
      for(int row = Ly-2;row > 0;--row){

         s = (*this)[row*Lx + col];

         //A: regular
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[row],ROU[row],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,3,1,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](row,col,s),0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,2,3,1),tmp4bis);

         ROU[row - 1].clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,env.gr(false,col)[row],0.0,ROU[row - 1]);

         //B: inverse
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[row],ROI[row],0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,3,1,2),tmp4bis);

         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](row,col,!s),0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(0,2,3,1),tmp4bis);

         ROI[row - 1].clear();
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,env.gr(true,col)[row],0.0,ROI[row - 1]);

      }

      // --- now move from left to right to get the expectation value of the interactions ---

      // --- First construct the right going operators for the first site -----

      // 1) construct left renormalized operator with unity

      s = (*this)[col];

      //A: regular
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](0,col,s),env.gr(false,col)[0],0.0,tmp3);

      M = env.gl(false,col - 1)[0].shape(2);
      N = tmp3.shape(1) * tmp3.shape(2);
      K = tmp3.shape(0);

      LOUU.resize(M,tmp3.shape(1),tmp3.shape(2));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,env.gl(false,col - 1)[0].data(),M,tmp3.data(),N,0.0,LOUU.data(),N);

      //overlap A
      tmp_over = Dot(LOUU,ROU[0]);

      //B: inverse
      Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](0,col,!s),env.gr(true,col)[0],0.0,tmp3);

      LOII.resize(M,tmp3.shape(1),tmp3.shape(2));
      blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,env.gl(true,col - 1)[0].data(),M,tmp3.data(),N,0.0,LOII.data(),N);

      //overlap B
      tmp_over += Dot(LOII,ROI[0]);

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[col] != (*this)[Lx + col]){

         //A: regular
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](0,col,!s),env.gr(false,col)[0],0.0,tmp3);

         LOUI.resize(M,tmp3.shape(1),tmp3.shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,env.gl(false,col - 1)[0].data(),M,tmp3.data(),N,0.0,LOUI.data(),N);

         //B: inverse
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](0,col,s),env.gr(true,col)[0],0.0,tmp3);

         LOIU.resize(M,tmp3.shape(1),tmp3.shape(2));
         blas::gemm(CblasRowMajor, CblasTrans, CblasNoTrans, M, N, K, 1.0,env.gl(true,col - 1)[0].data(),M,tmp3.data(),N,0.0,LOIU.data(),N);

      }

      // --- now for the middle sites, close down the operators on the left and construct new ones --- 
      for(int row = 1;row < Ly - 1;++row){

         s = (*this)[row*Lx + col];

         //1) close down LO(U/I)I if it contributes
         if((*this)[(row - 1)*Lx + col] != (*this)[row*Lx + col]){

            //A: regular
            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[row],ROU[row],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](row,col,!s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,env.gr(false,col)[row],0.0,ROU[row - 1]);

            //first part of expectation value
            ward = Dot(LOUI,ROU[row-1]);

            //B: inverse
            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[row],ROI[row],0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,3,1,2),tmp4bis);

            tmp4.clear();
            Gemm(CblasNoTrans,CblasNoTrans,1.0,tmp4bis,peps[0](row,col,s),0.0,tmp4);

            tmp4bis.clear();
            Permute(tmp4,shape(0,2,3,1),tmp4bis);

            Gemm(CblasNoTrans,CblasTrans,1.0,tmp4bis,env.gr(true,col)[row],0.0,ROI[row - 1]);

            //B expectation value
            ward += Dot(LOIU,ROI[row-1]);

            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new right going renormalized operators

         //A: regular
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,LOUU,env.gr(false,col)[row],0.0,tmp4);

         perm4.clear();
         Permute(tmp4,shape(1,2,0,3),perm4);

         // 1) construct new LOUU
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,s),perm4,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(2,0,1,3),tmp4bis);

         LOUU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[row],tmp4bis,0.0,LOUU);

         // 2) if it contributes, construct new left inverted : LOUI
         if((*this)[row*Lx + col] != (*this)[(row + 1)*Lx + col]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,!s),perm4,0.0,tmp4);

            Permute(tmp4,shape(2,0,1,3),tmp4bis);

            LOUI.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[row],tmp4bis,0.0,LOUI);

         }

         //B: inverse
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,LOII,env.gr(true,col)[row],0.0,tmp4);

         perm4.clear();
         Permute(tmp4,shape(1,2,0,3),perm4);

         // 1) construct new LOII
         tmp4.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,!s),perm4,0.0,tmp4);

         tmp4bis.clear();
         Permute(tmp4,shape(2,0,1,3),tmp4bis);

         LOII.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[row],tmp4bis,0.0,LOII);

         // 2) if it contributes, construct new left inverted : LOIU
         if((*this)[row*Lx + col] != (*this)[(row + 1)*Lx + col]){

            Gemm(CblasNoTrans,CblasNoTrans,1.0,peps[0](row,col,s),perm4,0.0,tmp4);

            Permute(tmp4,shape(2,0,1,3),tmp4bis);

            LOIU.clear();
            Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[row],tmp4bis,0.0,LOIU);

         }

      }

      s = (*this)[(Ly - 1)*Lx + col];

      //last site on the right: close down LOUI and LOIU if it contributes
      if((*this)[(Ly - 2)*Lx + col] != (*this)[(Ly - 1)*Lx + col]){

         //A: regular
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,col - 1)[Ly - 1],peps[0](Ly-1,col,!s),0.0,tmp3);

         M = tmp3.shape(0) * tmp3.shape(1);
         N = env.gr(false,col)[Ly - 1].shape(0);
         K = tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3.data(),K,env.gr(false,col)[Ly - 1].data(),K,0.0,ROU[Ly - 2].data(),N);

         //expectation value A
         ward = Dot(LOUI,ROU[Lx - 2]);

         //B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,col - 1)[Ly - 1],peps[0](Ly-1,col,s),0.0,tmp3);

         M = tmp3.shape(0) * tmp3.shape(1);
         N = env.gr(true,col)[Ly - 1].shape(0);
         K = tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0,tmp3.data(),K,env.gr(true,col)[Ly - 1].data(),K,0.0,ROI[Ly - 2].data(),N);

         //expectation value A
         ward += Dot(LOIU,ROI[Lx - 2]);

         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

   }

   // -- (3) -- || left col = 0: again similar to overlap calculation

   //A: regular
   RU[Ly - 2].clear();
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(false,0)[Ly - 1],env.gr(false,0)[Ly - 1],0.0,RU[Ly - 2]);

   //B: inverse
   RI[Ly - 2].clear();
   Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(true,0)[Ly - 1],env.gr(true,0)[Ly - 1],0.0,RI[Ly - 2]);

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();

      //U
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,0)[row],RU[row],0.0,tmp3);

      RU[row  - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(false,0)[row],0.0,RU[row - 1]);

      //I
      Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,0)[row],RI[row],0.0,tmp3);

      RI[row  - 1].clear();
      Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(true,0)[row],0.0,RI[row - 1]);

   }

   //construct 4 right going operators

   //U overlap
   LUU.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,0)[0],env.gr(false,0)[0],0.0,LUU);

   //I overlap
   LII.clear();
   Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,0)[0],env.gr(true,0)[0],0.0,LII);

   //calculate the overlap with this state
   tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   //only calculate LUI and LIU if it contributes
   if( (*this)[0] != (*this)[Lx] ){

      LUI.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,0)[0],env.gr(false,0)[0],0.0,LUI); //regular

      LIU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,0)[0],env.gr(true,0)[0],0.0,LIU); //inverse

   }

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      //only calculate if it contributes
      if( (*this)[(row - 1)*Lx] != (*this)[row*Lx] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(true,0)[row],RU[row],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(false,0)[row],0.0,RU[row - 1]);

         ward = Dot(LUI,RU[row - 1]);

         // B: inverse
         tmp3.clear();
         Gemm(CblasNoTrans,CblasNoTrans,1.0,env.gl(false,0)[row],RI[row],0.0,tmp3);

         //paste I to the right
         Gemm(CblasNoTrans,CblasTrans,1.0,tmp3,env.gr(true,0)[row],0.0,RI[row - 1]);

         ward += Dot(LIU,RI[row - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward /tmp_over;

      }

      //construct left renormalized operators for next site:

      //A: regular 
      tmp3.clear();
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LUU,env.gr(false,0)[row],0.0,tmp3);

      //1) construct new unity on the left
      LUU.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,0)[row],tmp3,0.0,LUU);

      //2) if it contributes, calculate inverse on the left
      if((*this)[row*Lx] != (*this)[(row + 1)*Lx]){

         LUI.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,0)[row],tmp3,0.0,LUI);

      }

      //B: inverse 
      Gemm(CblasNoTrans,CblasNoTrans,1.0,LII,env.gr(true,0)[row],0.0,tmp3);

      //1) construct new unity on the left
      LII.clear();
      Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(true,0)[row],tmp3,0.0,LII);

      //2) if it contributes, calculate inverse on the left
      if((*this)[row*Lx] != (*this)[(row + 1)*Lx]){

         LIU.clear();
         Gemm(CblasTrans,CblasNoTrans,1.0,env.gl(false,0)[row],tmp3,0.0,LIU);

      }

   }

   //last site of left col: close down LUI and LIU
   if((*this)[(Ly - 2)*Lx] != (*this)[(Ly - 1)*Lx]){

      //A: regular LUI
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(true,0)[Ly - 1],env.gr(false,0)[Ly - 1],0.0,RU[Ly - 2]);

      ward = Dot(LUI,RU[Ly-2]);

      //B: inverse LIU
      Gemm(CblasNoTrans,CblasTrans,1.0,env.gl(false,0)[Ly - 1],env.gr(true,0)[Ly - 1],0.0,RI[Ly - 2]);

      ward += Dot(LIU,RI[Ly-2]);

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

}


/**
 * @return the overlap between the trial peps (reg + inv) and the walker
 */
double Walker::overlap() {

   int half = Ly/2;

   //A: regular

   //bottom 
   env.gb(false,0).fill('b',false,*this);

   for(int i = 1;i <= half;++i)
      env.add_layer('b',i,false,*this);

   //top
   env.gt(false,Ly - 2).fill('t',false,*this);

   for(int i = Ly - 3;i >= half;--i)
      env.add_layer('t',i,false,*this);

   double tmp = env.gb(false,half).dot(env.gt(false,half));

   //A: inverse

   //bottom 
   env.gb(true,0).fill('b',true,*this);

   for(int i = 1;i <= half;++i)
      env.add_layer('b',i,true,*this);

   //top
   env.gt(true,Ly - 2).fill('t',true,*this);

   for(int i = Ly - 3;i >= half;--i)
      env.add_layer('t',i,true,*this);

   tmp += env.gb(true,half).dot(env.gt(true,half));

   return tmp;

}

/**
 * update the environment for the new walker, taking into account that only number on index have changed
 * @param index number containing info about which sites changed
 */
void Walker::update_env(int index) {

   if(index < 0){

      env.calc('A',false,*this);
      env.calc('A',true,*this);

   }
   else{

      if(index >= Lx*Ly){//vertical gate

         index -= Lx*Ly;

         int col = index % Lx;
         int row = (index - col)/Lx;

         //first bottom and top
         if(row == 0){

            env.calc('B',false,*this);
            env.calc('B',true,*this);

            env.add_layer('t',0,false,*this);
            env.add_layer('t',0,true,*this);

         }
         else if(row == Ly - 2){

            env.calc('T',false,*this);
            env.calc('T',true,*this);

            env.add_layer('b',Ly - 2,false,*this);
            env.add_layer('b',Ly - 2,true,*this);

         }
         else{

            for(int i = row;i < Ly - 1;++i){

               env.add_layer('b',i,false,*this);
               env.add_layer('b',i,true,*this);

            }

            for(int i = row;i >= 0;--i){

               env.add_layer('t',i,false,*this);
               env.add_layer('t',i,true,*this);

            }

         }

         //then left and right
         if(col == 0){

            env.calc('L',false,*this);
            env.calc('L',true,*this);

         }
         else if(col == Lx - 1){

            env.calc('R',false,*this);
            env.calc('R',true,*this);

         }
         else{

            for(int i = col;i < Ly - 1;++i){

               env.add_layer('l',i,false,*this);
               env.add_layer('l',i,true,*this);

            }

            for(int i = col - 1;i >= 0;--i){

               env.add_layer('r',i,false,*this);
               env.add_layer('r',i,true,*this);

            }

         }

      }
      else{//horizontal gate

         int col = index % Lx;
         int row = (index - col)/Lx;

         //first bottom and top
         if(row == 0){

            env.calc('B',false,*this);
            env.calc('B',true,*this);

         }
         else if(row == Ly - 1){

            env.calc('T',false,*this);
            env.calc('T',true,*this);

         }
         else{

            for(int i = row;i < Ly - 1;++i){

               env.add_layer('b',i,false,*this);
               env.add_layer('b',i,true,*this);

            }

            for(int i = row - 1;i >= 0;--i){

               env.add_layer('t',i,false,*this);
               env.add_layer('t',i,true,*this);

            }

         }

         //then left and right
         if(col == 0){

            env.calc('L',false,*this);
            env.calc('L',true,*this);

            env.add_layer('r',0,false,*this);
            env.add_layer('r',0,true,*this);

         }
         else if(col == Lx - 2){

            env.calc('R',false,*this);
            env.calc('R',true,*this);

            env.add_layer('l',Lx - 2,false,*this);
            env.add_layer('l',Lx - 2,true,*this);

         }
         else{

            for(int i = col;i < Ly - 1;++i){

               env.add_layer('l',i,false,*this);
               env.add_layer('l',i,true,*this);

            }

            for(int i = col;i >= 0;--i){

               env.add_layer('r',i,false,*this);
               env.add_layer('r',i,true,*this);

            }

         }

      }

   }

}
