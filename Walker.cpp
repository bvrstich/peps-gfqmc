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
   DArray<3> tmp3bis;

   //tmp comes out index (t,b): U
   Contract(1.0,env[myID].gt(false,0)[Lx - 1],shape(1),env[myID].gb(false,0)[Lx - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   RU[Lx - 2] = tmp4.reshape_clear(shape(env[myID].gt(false,0)[Lx - 1].shape(0),env[myID].gb(false,0)[Lx - 1].shape(0)));

   //inverse
   Contract(1.0,env[myID].gt(true,0)[Lx - 1],shape(1),env[myID].gb(true,0)[Lx - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   RI[Lx - 2] = tmp4.reshape_clear(shape(env[myID].gt(true,0)[Lx - 1].shape(0),env[myID].gb(true,0)[Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();

      //U
      Contract(1.0,env[myID].gt(false,0)[col],shape(2),RU[col],shape(0),0.0,tmp3);
      Contract(1.0,tmp3,shape(1,2),env[myID].gb(false,0)[col],shape(1,2),0.0,RU[col-1]);

      //I
      Contract(1.0,env[myID].gt(true,0)[col],shape(2),RI[col],shape(0),0.0,tmp3);
      Contract(1.0,tmp3,shape(1,2),env[myID].gb(true,0)[col],shape(1,2),0.0,RI[col-1]);

   }

   //4 left going operators: S+/- and 1
   DArray<2> LUU;
   DArray<2> LUI;

   DArray<2> LIU;
   DArray<2> LII;

   TArray<double,5> tmp5;

   //U overlap
   Contract(1.0,env[myID].gt(false,0)[0],shape(1),env[myID].gU(false)(0,0),shape(1),0.0,tmp5);
   LUU = tmp5.reshape_clear( shape(env[myID].gt(false,0)[0].shape(2),env[myID].gU(false)(0,0).shape(3)) );

   //I overlap
   Contract(1.0,env[myID].gt(true,0)[0],shape(1),env[myID].gU(true)(0,0),shape(1),0.0,tmp5);
   LII = tmp5.reshape_clear( shape(env[myID].gt(true,0)[0].shape(2),env[myID].gU(true)(0,0).shape(3)) );

   //calculate the overlap with this state
   double tmp_over = Dot(RU[0],LUU) + Dot(RI[0],LII);

   nn_over.push_back(1.0);

   //only calculate LUI and LIU if it contributes
   if( (*this)[0] != (*this)[1] ){

      //regular
      Contract(1.0,env[myID].gt(false,0)[0],shape(1),env[myID].gU(true)(0,0),shape(1),0.0,tmp5);
      LUI = tmp5.reshape_clear( shape(env[myID].gt(false,0)[0].shape(2),env[myID].gU(true)(0,0).shape(3)) );

      //inverse
      Contract(1.0,env[myID].gt(true,0)[0],shape(1),env[myID].gU(false)(0,0),shape(1),0.0,tmp5);
      LIU = tmp5.reshape_clear( shape(env[myID].gt(true,0)[0].shape(2),env[myID].gU(false)(0,0).shape(3)) );

   }
  
   //now for the middle terms
   for(int col = 1;col < Lx - 1;++col){

      //only calculate if it contributes
      if( (*this)[col - 1] != (*this)[col] ){

         // A: regular

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(1.0,env[myID].gt(false,0)[col],shape(2),RU[col],shape(0),0.0,tmp3);

         // 1) paste I to the right
         M = tmp3.shape(0);
         N = env[myID].gU(true)(0,col).shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,env[myID].gU(true)(0,col).data(),K,0.0,RU[col-1].data(),N);

         ward = Dot(LUI,RU[col - 1]);

         // B: inverse

         //construct the right intermediate contraction (paste top to right)
         Contract(1.0,env[myID].gt(true,0)[col],shape(2),RI[col],shape(0),0.0,tmp3);

         // 1) paste U to the right
         M = tmp3.shape(0);
         N = env[myID].gU(false)(0,col).shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,env[myID].gU(false)(0,col).data(),K,0.0,RI[col-1].data(),N);

         ward += Dot(LIU,RI[col - 1]);

         nn_over.push_back(ward/tmp_over);

         //contract with left LI 
         EL -= 0.5 * ward /tmp_over;

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

   cout << EL << endl;
/*
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
      Contract(1.0,env[myID].gt[row][Lx - 1],shape(1),env[myID].gU(false)(row,Lx-1),shape(1),0.0,tmp5);

      //then bottom enviroment
      TArray<double,6> tmp6;
      Contract(1.0,tmp5,shape(3),env[myID].gb[row-1][Lx-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt[row][Lx - 1].shape(0),env[myID].gU(false)(row,Lx-1).shape(0),env[myID].gb[row-1][Lx - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int col = Lx-2;col > 0;--col){

         I4.clear();
         Contract(1.0,env[myID].gt[row][col],shape(2),RO[col],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[col-1].clear();
         Contract(1.0,I4bis,shape(2,3),env[myID].gb[row-1][col],shape(1,2),0.0,RO[col-1]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // 1) construct left renormalized operator with unity

      //paste top environment on local Sz
      tmp5.clear();
      Contract(1.0,env[myID].gt[row][0],shape(1),env[myID].gU(false)(row,0),shape(1),0.0,tmp5);

      //then bottom enviroment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),env[myID].gb[row-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
      LOU = tmp6.reshape_clear(shape(env[myID].gt[row][0].shape(2),env[myID].gU(false)(row,0).shape(3),env[myID].gb[row-1][0].shape(2)));

      tmp_over = Dot(LOU,RO[0]);

      // 2) construct left operator with inverted spin if it contributes
      if((*this)[row*Lx] != (*this)[row*Lx + 1]){

         //paste top environment on local inverted Sz
         tmp5.clear();
         Contract(1.0,env[myID].gt[row][0],shape(1),env[myID].gU(true)(row,0),shape(1),0.0,tmp5);

         //then bottom enviroment on that
         tmp6.clear();
         Contract(1.0,tmp5,shape(3),env[myID].gb[row-1][0],shape(1),0.0,tmp6);

         //move to a DArray<3> object: order (top-env,peps-row,bottom-env)
         LOI = tmp6.reshape_clear(shape(env[myID].gt[row][0].shape(2),env[myID].gU(true)(row,0).shape(3),env[myID].gb[row-1][0].shape(2)));

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int col = 1;col < Lx - 1;++col){

         enum {i,j,k,o,m,n};

         //1) close down LOI with I if it contributes
         if((*this)[row*Lx + col - 1] != (*this)[row*Lx + col]){

            //first add top to the right side, put it in I4
            I4.clear();
            Contract(1.0,env[myID].gt[row][col],shape(2),RO[col],shape(0),0.0,I4);

            I4bis.clear();
            Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

            Contract(1.0,I4bis,shape(2,3),env[myID].gb[row-1][col],shape(1,2),0.0,RO[col-1]);

            //expectation value:
            ward = Dot(LOI,RO[col-1]);
            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity
         I4.clear();
         Contract(1.0,env[myID].gt[row][col],shape(0),LOU,shape(0),0.0,I4);

         // 1) construct new left unity
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOU.clear();
         Contract(1.0,I4bis,shape(2,3),env[myID].gb[row-1][col],shape(0,1),0.0,LOU);

         // 2) if it contributes, construct new left inverted 
         if((*this)[row*Lx + col] != (*this)[row*Lx + col + 1]){

            Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

            LOI.clear();
            Contract(1.0,I4bis,shape(2,3),env[myID].gb[row-1][col],shape(0,1),0.0,LOI);

         }

      }

      //last site on the right: close down LOI if it contributes
      if((*this)[row*Lx + Lx - 2] != (*this)[row*Lx + Lx - 1]){

         //paste top environment on
         tmp5.clear();
         Contract(1.0,env[myID].gt[row][Lx - 1],shape(1),env[myID].gU(true)(row,Lx-1),shape(1),0.0,tmp5);

         //then bottom enviroment
         Contract(1.0,tmp5,shape(3),env[myID].gb[row-1][Lx-1],shape(1),0.0,tmp6);

         //move to a DArray<3> object
         RO[Lx - 2] = tmp6.reshape_clear(shape(env[myID].gt[row][Lx - 1].shape(0),env[myID].gU(true)(row,Lx-1).shape(0),env[myID].gb[row-1][Lx - 1].shape(0)));

         //add to energy
         ward = Dot(LOI,RO[Lx - 2]);
         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (t,b)
   tmp4.clear();
   Contract(1.0,env[myID].gt[Ly-2][Lx - 1],shape(1),env[myID].gb[Ly-2][Lx - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Lx - 2] = tmp4.reshape_clear(shape(env[myID].gt[Ly-2][Lx - 1].shape(0),env[myID].gb[Ly-2][Lx - 1].shape(0)));

   //now construct the rest
   for(int col = Lx - 2;col > 0;--col){

      tmp3.clear();
      Contract(1.0,env[myID].gt[Ly-2][col],shape(2),R[col],shape(0),0.0,tmp3);

      R[col - 1].clear();
      Contract(1.0,tmp3,shape(1,2),env[myID].gb[Ly-2][col],shape(1,2),0.0,R[col-1]);

   }

   //construct the left going operators on the first top site

   //unity
   tmp5.clear();
   Contract(1.0,env[myID].gU(false)(Ly-1,0),shape(2),env[myID].gb[Ly-2][0],shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear(shape(env[myID].gU(false)(Ly-1,0).shape(3),env[myID].gb[Ly-2][0].shape(2)));

   tmp_over = Dot(LU,R[0]);

   //inverse if it contributes
   if( (*this)[(Ly - 1)*Lx] != (*this)[(Ly - 1)*Lx + 1]){

      tmp5.clear();
      Contract(1.0,env[myID].gU(true)(Ly-1,0),shape(2),env[myID].gb[Ly-2][0],shape(1),0.0,tmp5);

      LI = tmp5.reshape_clear(shape(env[myID].gU(true)(Ly-1,0).shape(3),env[myID].gb[Ly-2][0].shape(2)));

   }

   //middle of the chain:
   for(int col = 1;col < Lx-1;++col){

      //first close down the I term from the previous site for the energy
      if( (*this)[(Ly - 1)*Lx + col - 1] != (*this)[(Ly - 1)*Lx + col]){

         //construct the right intermediate contraction (paste bottom to right)
         tmp3.clear();
         Contract(1.0,env[myID].gb[Ly-2][col],shape(2),R[col],shape(1),0.0,tmp3);

         // 1) paste Sx to the right
         M = env[myID].gU(true)(Ly-1,col).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, env[myID].gU(true)(Ly-1,col).data(),K,tmp3.data(),K,0.0,R[col-1].data(),N);

         //contract with left LI 
         ward = Dot(LI,R[col - 1]);
         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward/tmp_over;

      }

      //construct left renormalized operators for next site: first paste bottom to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(1),env[myID].gb[Ly-2][col],shape(0),0.0,tmp3);

      // 1) construct new left unity operator
      tmp3bis.clear();
      Contract(1.0,env[myID].gU(false)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

      LU = tmp3bis.reshape_clear(shape(env[myID].gU(false)(Ly-1,col).shape(3),env[myID].gb[Ly-2][col].shape(2)));

      // 2) if it contributes, construct new left LI
      if( (*this)[(Ly - 1)*Lx + col] != (*this)[(Ly - 1)*Lx + col + 1]){

         tmp3bis.clear();
         Contract(1.0,env[myID].gU(true)(Ly-1,col),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

         LI = tmp3bis.reshape_clear(shape(env[myID].gU(true)(Ly-1,col).shape(3),env[myID].gb[Ly-2][col].shape(2)));

      }

   }
   
   //finally close down on last top site

   // close down last LI
   if( (*this)[(Ly - 1)*Lx + Lx - 2] != (*this)[(Ly - 1)*Lx + Lx - 1]){

      //tmp comes out index (t,b)
      tmp5.clear();
      Contract(1.0,env[myID].gU(true)(Ly-1,Lx-1),shape(2),env[myID].gb[Ly-2][Lx - 1],shape(1),0.0,tmp5);

      //reshape tmp to a 2-index array
      R[Lx - 2] = tmp5.reshape_clear(shape(env[myID].gU(true)(Ly-1,Lx-1).shape(0),env[myID].gb[Ly-2][Lx - 1].shape(0)));

      //energy
      ward =  Dot(LI,R[Lx-2]);
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

   //calculate the single layer contractions first:
   env[myID].gU(false).fill('V',false,peps,*this);
   env[myID].gU(true).fill('V',true,peps,*this);

   //first construct the top and bottom (horizontal) environment layers
   Environment::calc_env('V',peps,*this);

   // -- (1) -- || left column: similar to overlap calculation

   //tmp comes out index (r,l)
   tmp4.clear();
   Contract(1.0,env[myID].gr[0][Ly - 1],shape(1),env[myID].gl[0][Ly - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp4.reshape_clear(shape(env[myID].gr[0][Ly - 1].shape(0),env[myID].gl[0][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();
      Contract(1.0,env[myID].gr[0][row],shape(2),R[row],shape(0),0.0,tmp3);

      R[row - 1].clear();
      Contract(1.0,tmp3,shape(1,2),env[myID].gl[0][row],shape(1,2),0.0,R[row-1]);

   }

   //left going operator:
   tmp5.clear();

   //unity
   Contract(1.0,env[myID].gr[0][0],shape(1),env[myID].gU(false)(0,0),shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear(shape(env[myID].gr[0][0].shape(2),env[myID].gU(false)(0,0).shape(3)));

   tmp_over = Dot(LU,R[0]);

   //inverse if it contributes:
   if( (*this)[0] != (*this)[Lx]){

      tmp5.clear();

      Contract(1.0,env[myID].gr[0][0],shape(1),env[myID].gU(true)(0,0),shape(1),0.0,tmp5);

      LI = tmp5.reshape_clear(shape(env[myID].gr[0][0].shape(2),env[myID].gU(true)(0,0).shape(3)));

   }

   //now for the middle terms
   for(int row = 1;row < Ly - 1;++row){

      if( (*this)[(row - 1)*Lx] != (*this)[row*Lx] ){

         //first close down the LI from the previous site for the energy if necessary

         //construct the right intermediate contraction (paste top to right)
         tmp3.clear();
         Contract(1.0,env[myID].gr[0][row],shape(2),R[row],shape(0),0.0,tmp3);

         // 1) paste Sx to the right
         M = tmp3.shape(0);
         N = env[myID].gU(true)(row,0).shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, tmp3.data(),K,env[myID].gU(true)(row,0).data(),K,0.0,R[row-1].data(),N);

         //contract with left Sx
         ward = Dot(LI,R[row - 1]);
         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward / tmp_over;

      }

      //construct left renormalized operators for next site: first paste top to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(0),env[myID].gr[0][row],shape(0),0.0,tmp3);

      //1) construct new unity on the left
      tmp3bis.clear();
      Contract(1.0,tmp3,shape(0,1),env[myID].gU(false)(row,0),shape(0,1),0.0,tmp3bis);

      LU = tmp3bis.reshape_clear(shape(env[myID].gr[0][row].shape(2),env[myID].gU(false)(row,0).shape(3)));

      // 2) if contribution, construct new inverse on the left
      if( (*this)[row*Lx] != (*this)[(row + 1)*Lx] ){

         tmp3bis.clear();
         Contract(1.0,tmp3,shape(0,1),env[myID].gU(true)(row,0),shape(0,1),0.0,tmp3bis);

         LI = tmp3bis.reshape_clear(shape(env[myID].gr[0][row].shape(2),env[myID].gU(true)(row,0).shape(3)));

      }

   }

   //last site of left column: close down the left LI if necessary
   if( (*this)[(Ly - 2)*Lx] != (*this)[(Ly - 1)*Lx] ){

      Contract(1.0,env[myID].gr[0][Ly-1],shape(1),env[myID].gU(true)(Ly-1,0),shape(1),0.0,tmp5);

      R[Ly-2] = tmp5.reshape_clear(shape(env[myID].gr[0][Ly-1].shape(0),env[myID].gU(true)(Ly-1,0).shape(0)));

      ward = Dot(LI,R[Ly-2]);
      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward / tmp_over;

   }

   // -- (2) -- now move from left to right calculating everything like an MPO/MPS expectation value

   for(int col = 1;col < Lx - 1;++col){

      //first create right renormalized operator

      //paste right environment on
      tmp5.clear();
      Contract(1.0,env[myID].gr[col][Ly - 1],shape(1),env[myID].gU(false)(Ly-1,col),shape(1),0.0,tmp5);

      //then left enviroment
      TArray<double,6> tmp6;
      Contract(1.0,tmp5,shape(3),env[myID].gl[col-1][Ly-1],shape(1),0.0,tmp6);

      //move to a DArray<3> object
      RO[Ly - 2] = tmp6.reshape_clear(shape(env[myID].gr[col][Ly - 1].shape(0),env[myID].gU(false)(Ly-1,col).shape(0),env[myID].gl[col-1][Ly - 1].shape(0)));

      DArray<4> I4;
      DArray<4> I4bis;

      //now construct the middle operators
      for(int row = Ly-2;row > 0;--row){

         I4.clear();
         Contract(1.0,env[myID].gr[col][row],shape(2),RO[row],shape(0),0.0,I4);

         enum {i,j,k,o,m,n};

         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

         RO[row-1].clear();
         Contract(1.0,I4bis,shape(2,3),env[myID].gl[col-1][row],shape(1,2),0.0,RO[row-1]);

      }

      // --- now move from left to right to get the expecation value of the interactions ---
      // --- First construct the left going operators for the first site -----

      // construct left renormalized operator with unity

      //paste left environment on local Sz
      tmp5.clear();
      Contract(1.0,env[myID].gr[col][0],shape(1),env[myID].gU(false)(0,col),shape(1),0.0,tmp5);

      //then right environment on that
      tmp6.clear();
      Contract(1.0,tmp5,shape(3),env[myID].gl[col-1][0],shape(1),0.0,tmp6);

      //move to a DArray<3> object: order (left-env,peps-row,right-env)
      LOU = tmp6.reshape_clear(shape(env[myID].gr[col][0].shape(2),env[myID].gU(false)(0,col).shape(3),env[myID].gl[col-1][0].shape(2)));

      tmp_over =  Dot(LOU,RO[0]);

      //construct left inverse if it contributes
      if( (*this)[col] != (*this)[Lx + col]){

         tmp5.clear();
         Contract(1.0,env[myID].gr[col][0],shape(1),env[myID].gU(true)(0,col),shape(1),0.0,tmp5);

         //then right environment on that
         tmp6.clear();
         Contract(1.0,tmp5,shape(3),env[myID].gl[col-1][0],shape(1),0.0,tmp6);

         //move to a DArray<3> object: order (left-env,peps-row,right-env)
         LOI = tmp6.reshape_clear(shape(env[myID].gr[col][0].shape(2),env[myID].gU(true)(0,col).shape(3),env[myID].gl[col-1][0].shape(2)));

      }

      // --- now for the middle sites, close down the operators on the left and construct new 1.0s --- 
      for(int row = 1;row < Ly - 1;++row){

         enum {i,j,k,o,m,n};

         if( (*this)[ (row - 1)*Lx + col] != (*this)[row*Lx + col]){

            //close down LOI with I if it contributes

            //first add top to the right side, put it in I4
            I4.clear();
            Contract(1.0,env[myID].gr[col][row],shape(2),RO[row],shape(0),0.0,I4);

            I4bis.clear();
            Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(m,j,n,k),0.0,I4bis,shape(i,m,n,o));

            Contract(1.0,I4bis,shape(2,3),env[myID].gl[col-1][row],shape(1,2),0.0,RO[row-1]);

            //expectation value:
            ward = Dot(LOI,RO[row-1]);
            nn_over.push_back(ward/tmp_over);

            EL -= 0.5 * ward / tmp_over;

         }

         // now construct the new left going renormalized operators

         //first attach top to left unity
         I4.clear();
         Contract(1.0,env[myID].gr[col][row],shape(0),LOU,shape(0),0.0,I4);

         // and construct new left unity
         I4bis.clear();
         Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(false)(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

         LOU.clear();
         Contract(1.0,I4bis,shape(2,3),env[myID].gl[col-1][row],shape(0,1),0.0,LOU);

         //if it contributes, construct new left inverse
         if( (*this)[ row*Lx + col] != (*this)[ (row + 1)*Lx + col]){

            Contract(1.0,I4,shape(i,j,k,o),env[myID].gU(true)(row,col),shape(k,i,m,n),0.0,I4bis,shape(j,n,o,m));

            LOI.clear();
            Contract(1.0,I4bis,shape(2,3),env[myID].gl[col-1][row],shape(0,1),0.0,LOI);

         }

      }

      //last site on the right: close down on the incomings if possible

      //first LI with I
      if( (*this)[ (Ly - 2)*Lx + col] != (*this)[ (Ly - 1)*Lx + col]){

         //paste top environment on
         tmp5.clear();
         Contract(1.0,env[myID].gr[col][Ly - 1],shape(1),env[myID].gU(true)(Ly-1,col),shape(1),0.0,tmp5);

         //then bottom enviroment
         Contract(1.0,tmp5,shape(3),env[myID].gl[col-1][Ly-1],shape(1),0.0,tmp6);

         //move to a DArray<3> object
         RO[Ly - 2] = tmp6.reshape_clear(shape(env[myID].gr[col][Ly - 1].shape(0),env[myID].gU(true)(Ly-1,col).shape(0),env[myID].gl[col-1][Ly - 1].shape(0)));

         //add to energy
         ward = Dot(LOI,RO[Ly - 2]);
         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward / tmp_over;

      }

   }

   // -- (3) -- || top row = Ly-1: again similar to overlap calculation

   //first construct the right renormalized operators

   //tmp comes out index (r,l)
   tmp4.clear();
   Contract(1.0,env[myID].gr[Lx-2][Ly - 1],shape(1),env[myID].gl[Lx-2][Ly - 1],shape(1),0.0,tmp4);

   //reshape tmp to a 2-index array
   R[Ly - 2] = tmp4.reshape_clear(shape(env[myID].gr[Lx-2][Ly - 1].shape(0),env[myID].gl[Lx-2][Ly - 1].shape(0)));

   //now construct the rest
   for(int row = Ly - 2;row > 0;--row){

      tmp3.clear();
      Contract(1.0,env[myID].gr[Lx-2][row],shape(2),R[row],shape(0),0.0,tmp3);

      R[row - 1].clear();
      Contract(1.0,tmp3,shape(1,2),env[myID].gl[Lx-2][row],shape(1,2),0.0,R[row-1]);

   }

   //construct the left going operators on the first top site

   //first unity
   tmp5.clear();
   Contract(1.0,env[myID].gU(false)(0,Lx-1),shape(2),env[myID].gl[Lx-2][0],shape(1),0.0,tmp5);

   LU = tmp5.reshape_clear(shape(env[myID].gU(false)(0,Lx-1).shape(3),env[myID].gl[Lx-2][0].shape(2)));

   tmp_over = Dot(LU,R[0]);

   //then inverse if necessary
   if( (*this)[Lx - 1] != (*this)[2*Lx - 1]){

      tmp5.clear();
      Contract(1.0,env[myID].gU(true)(0,Lx-1),shape(2),env[myID].gl[Lx-2][0],shape(1),0.0,tmp5);

      LI = tmp5.reshape_clear(shape(env[myID].gU(true)(0,Lx-1).shape(3),env[myID].gl[Lx-2][0].shape(2)));

   }

   //middle of the chain:
   for(int row = 1;row < Ly-1;++row){

      //first close down the inverse term from the previous site for the energy

      if( (*this)[ (row - 1)*Lx + Lx - 1] != (*this)[row*Lx + Lx - 1] ) {

         //construct the right intermediate contraction (paste bottom to right)
         tmp3.clear();
         Contract(1.0,env[myID].gl[Lx-2][row],shape(2),R[row],shape(1),0.0,tmp3);

         // 1) paste Sx to the right
         M = env[myID].gU(true)(row,Lx-1).shape(0);
         N = tmp3.shape(0);
         K = tmp3.shape(1) * tmp3.shape(2);

         blas::gemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1.0, env[myID].gU(true)(row,Lx-1).data(),K,tmp3.data(),K,0.0,R[row-1].data(),N);

         //contract with left inverse
         ward = Dot(LI,R[row - 1]);
         nn_over.push_back(ward/tmp_over);

         EL -= 0.5 * ward / tmp_over;

      }

      //construct left renormalized operators for next site: first paste bottom to Left unity
      tmp3.clear();
      Contract(1.0,LU,shape(1),env[myID].gl[Lx-2][row],shape(0),0.0,tmp3);

      // ly construct new unity on the left
      tmp3bis.clear();
      Contract(1.0,env[myID].gU(false)(row,Lx-1),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

      LU = tmp3bis.reshape_clear(shape(env[myID].gU(false)(row,Lx-1).shape(3),env[myID].gl[Lx-2][row].shape(2)));

      // construct new left inverse if it contributes
      if( (*this)[ row*Lx + Lx - 1] != (*this)[ (row + 1)*Lx + Lx - 1] ) {

         tmp3bis.clear();
         Contract(1.0,env[myID].gU(true)(row,Lx-1),shape(0,2),tmp3,shape(0,1),0.0,tmp3bis);

         LI = tmp3bis.reshape_clear(shape(env[myID].gU(true)(row,Lx-1).shape(3),env[myID].gl[Lx-2][row].shape(2)));

      }

   }

   //finally close down on last top site

   // I to close down LI
   if( (*this)[ (Ly - 2)*Lx + Lx - 1] != (*this)[ (Ly - 1)*Lx + Lx - 1] ) {

      //tmp comes out index (r,l)
      tmp5.clear();
      Contract(1.0,env[myID].gU(true)(Ly-1,Lx-1),shape(2),env[myID].gl[Lx-2][Ly - 1],shape(1),0.0,tmp5);

      //reshape tmp to a 2-index array
      R[Ly - 2] = tmp5.reshape_clear(shape(env[myID].gU(true)(Ly-1,Lx-1).shape(0),env[myID].gl[Lx-2][Ly - 1].shape(0)));

      //energy
      ward = Dot(LI,R[Ly-2]);

      nn_over.push_back(ward/tmp_over);

      EL -= 0.5 * ward / tmp_over;

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
