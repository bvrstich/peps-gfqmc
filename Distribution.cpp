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
 * @return the list of final Walker states
 */
const vector< Walker > &Distribution::glist() const {

   return list;

}

/** 
 * @return the final state Walker with index i
 * const version
 */
const Walker &Distribution::gwalker(int index) const {

   return list[index];

}

/** 
 * @return the final state Walker with index i
 */
Walker &Distribution::gwalker(int index) {

   return list[index];

}

/**
 * construct and fill the distribution by calculating the matrix elements <0|1-dtau * H|i> for all i = 0,...,n
 * @param walker_i input walker, the list and distribution is constructed from this
 */
void Distribution::construct(const Walker &walker_i){

   //first reset the lists
   list.clear();
   this->clear();

   //f == i first 
   list.push_back(walker_i);

   //first horizontal 'final' states
   for(int r = 0;r < Ly;++r){

      for(int c = 0;c < Lx - 1;++c){

         if(walker_i[r*Lx + c] != walker_i[r*Lx + c + 1]){

            Walker walker_f(walker_i);

            walker_f[r*Lx + c] = !(walker_i[r*Lx + c]);
            walker_f[r*Lx + c + 1] = !(walker_i[r*Lx + c + 1]);

            list.push_back(walker_f);

         }

      }

   }

   //then vertical 'final' states
   for(int c = Lx - 1;c >= 0;--c){

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

   (*this)[0] = 1.0;

   for(unsigned int i = 1;i < this->size();++i)
      (*this)[i] = walker_i.gnn_over(i) * walker_i.gnn_over(i);

}

/**
 * draw a new walker using metropolis algorithm
 */
int Distribution::metropolis() const {

   //draw uniform move
   int trial = (RN()*(list.size() - 1) + 1);

   double x = RN();

   if((*this)[trial] > x)
      return trial;
   else
      return 0;

}
