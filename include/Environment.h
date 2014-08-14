#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      static void init();

      static void calc_env(char,const PEPS< double > &,const Walker &walker);

      static void calc_overlap_env(const PEPS< double > &,const Walker &);

      static void test_env();

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      static vector< vector< MPS > > l;
      static vector< vector< MPS > > r;
      static vector< vector< MPS > > t;
      static vector< vector< MPS > > b;

      //!contraction between current walker state and peps
      static vector< SL_PEPS > U;

      //!contraction between inverted walker state and peps
      static vector< SL_PEPS > I;

   private:

};

#endif
