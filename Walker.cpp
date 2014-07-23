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
