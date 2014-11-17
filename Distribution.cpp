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
Distribution::Distribution() : vector<double> () { }

/** 
 * copy constructor
 */
Distribution::Distribution(const Distribution &dist_copy) : vector<double>(dist_copy){

   this->conf = dist_copy.gconf();
   this->ind = dist_copy.gind();
   this->sign_flip = dist_copy.gsign_flip();

}

/** 
 * empty destructor
 */
Distribution::~Distribution(){ }

/** 
 * @return number of grid points
 */
int Distribution::gn() const {

   return this->size();

}

/** 
 * @return configurations of final walker states
 */
const vector< vector<bool> > &Distribution::gconf() const {

   return conf;

}

/** 
 * @return list of indicies of final walker configuration (row,col and dir)
 */
const vector<int> &Distribution::gind() const {

   return ind;

}

/** 
 * @return list of booleans, for sign_flip is true sign of final walker is flipped
 */
const vector<bool> &Distribution::gsign_flip() const {

   return sign_flip;

}

ostream &operator<<(ostream &output,const Distribution &dist_p){

   for(int i = 0;i < dist_p.gn();++i)
      output << i << "\t" << dist_p[i] << endl;

   return output;

}

/** 
 * normalize the distribution
 * @return the total weight/norm of the distrubition
 */
double Distribution::normalize(){

   double nrm = 0.0;

   for(int i = 0;i < this->size();++i)
      nrm += (*this)[i];

   for(int i = 0;i < this->size();++i)
      (*this)[i] /= nrm;

   return nrm;

}

/**
 * draw a number from the distribution
 */
int Distribution::draw() const {

   //Get what you should do
   int trial = RN()*this->size();
   double x = RN();

   while((*this)[trial] < x){

      trial = RN()*this->size();
      x = RN();

   }

   return trial;

}

/** 
 * @return the configuration of the final walker state Walker corresponding to index pick
 * const version
 */
const vector<bool> &Distribution::gconf(int pick) const {

   return conf[pick];

}

/** 
 * @return the index of the final walker state that has been picked
 * const version
 */
int Distribution::gind(int pick) const {

   return ind[pick];

}

/** 
 * @return the sign flip flag of the picked walker
 * const version
 */
bool Distribution::gsign_flip(int pick) const {

   return sign_flip[pick];

}

/**
 * construct and fill the distribution by calculating the matrix elements <0|1-dtau * H|i> for all i = 0,...,n
 * @param walker_i input walker, the list and distribution is constructed from this
 * @param dtau timestep
 * @param ET estimator for ground state energy
 */
void Distribution::construct(const Walker &walker_i,double dtau,double ET){

   //first reset the lists
   conf.clear();
   ind.clear();
   sign_flip.clear();

   this->clear();

   //f == i first 
   conf.push_back(walker_i);
   ind.push_back(-1);
   sign_flip.push_back(false);

   vector<bool> tmp;

   //first horizontal 'final' states
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx - 1;++c){

         if(walker_i[r*Lx + c] != walker_i[r*Lx + (c + 1)]){

            tmp = conf[0];

            tmp[r*Lx + c] = !(conf[0][r*Lx + c]);
            tmp[r*Lx + (c + 1)] = !(conf[0][r*Lx + (c + 1)]);

            conf.push_back(tmp);
            ind.push_back(r*Lx + c);

         }

      }

   }

   //then vertical 'final' states
   for(int c = Lx - 1;c >= 0;--c){

      for(int r = 0;r < Ly - 1;++r){

         if(walker_i[r*Lx + c] != walker_i[(r + 1)*Lx + c]){

            tmp = conf[0];

            tmp[r*Lx + c] = !(conf[0][r*Lx + c]);
            tmp[(r + 1)*Lx + c] = !(conf[0][(r + 1)*Lx + c]);

            conf.push_back(tmp);
            ind.push_back(Lx*Ly + r*Lx + c);

         }

      }

   }

   this->resize(conf.size());

   (*this)[0] = 1.0 - dtau * (walker_i.pot_en() - ET);

   for(int i = 1;i < conf.size();++i){

      double ward = 0.5 * dtau * ( walker_i.gnn_over(i) );

      if(ward < 0.0){

         (*this)[i] = -ward;
         sign_flip.push_back(true);

      }
      else{

         (*this)[i] = ward;
         sign_flip.push_back(false);

      }

   }

}

/**
 * check for negative entries
 */
void Distribution::check_negative() const {

   for(int i = 0;i < this->size();++i)
      if( (*this)[i] < 0.0 )
         cout << "ERROR\t" << (*this)[i] << endl;

}
