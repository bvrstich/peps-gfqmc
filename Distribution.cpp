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
 * standard constructor
 * @param walker_i input walker, the list and distribution is constructed from this
 */
Distribution::Distribution(const Walker &walker_i) : vector<double> (){

   //f == i first 
   list.push_back(walker_i);

   //first horizontal 'final' states
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx - 1;++c){

         if(walker_i[r*Lx + c] != walker_i[r*Lx + (c + 1)]){

            Walker walker_f(walker_i);

            walker_f[r*Lx + c] = !(walker_i[r*Lx + c]);
            walker_f[r*Lx + (c + 1)] = !(walker_i[r*Lx + (c + 1)]);

            list.push_back(walker_f);

         }

      }

   }

   //then vertical 'final' states
   for(int c = 0;c < Lx;++c){

      for(int r = 0;r < Ly - 1;++r){

         if(walker_i[r*Lx + c] != walker_i[(r + 1)*Lx + c]){

            Walker walker_f(walker_i);

            walker_f[r*Lx + c] = !(walker_i[r*Lx + c]);
            walker_f[(r + 1)*Lx + c] = !(walker_i[(r + 1)*Lx + c]);

            list.push_back(walker_f);

         }

      }

   }

   this->resize(list.size());

}

/** 
 * copy constructor
 */
Distribution::Distribution(const Distribution &dist_copy) : vector<double>(dist_copy){

   list = dist_copy.glist();

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
   double x = RN();

   int cnt = 0;

   double sum = (*this)[0];

   while ( (sum <= x) && (cnt < this->size() - 1) ){

      cnt += 1;
      sum += (*this)[cnt];

   }

   return cnt;

}

/** 
 * @return the list of final Walker states
 */
const vector< Walker > &Distribution::glist() const {

   return list;

}

/** 
 * @return the final state Walker with index i
 */
const Walker &Distribution::gwalker(int index) const {

   return list[index];

}

/**
 * fill the distribution by calculating the matrix elements <0|1-dtau * H|i> for all i = 0,...,n
 * @param dtau timestep
 */
void Distribution::fill(double dtau,double ET){

   (*this)[0] = 1.0 - dtau * (list[0].exp_en(list[0]) - ET);

   for(int i = 1;i < list.size();++i)
      (*this)[i] = - dtau * (list[0].exp_en(list[i]));

}
