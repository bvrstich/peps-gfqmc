#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

template<typename T>
class PEPS;

class MPS;
class Walker;
class SL_PEPS;

/**
 * @author Brecht Verstichel
 * @data 02-05-2014\n\n
 * Class used to calculate the enviroment of a peps. Needed for the calculation of expectation values and the update of tensors.
 */
class Environment {

   public:

      Environment();

      Environment(int,int,int);

      //copy constructor
      Environment(const Environment &);

      //destructor
      virtual ~Environment();

      const MPS &gl(int) const;
      MPS &gl(int);

      const MPS &gr(int) const;
      MPS &gr(int);

      const MPS &gt(int) const;
      MPS &gt(int);

      const MPS &gb(int) const;
      MPS &gb(int);

      const vector< MPS > &gl() const;
      const vector< MPS > &gr() const;
      const vector< MPS > &gt() const;
      const vector< MPS > &gb() const;

      int gD() const;
      int gD_aux() const;

      int gcomp_sweeps() const;

      void sD(int);
      void sD_aux(int);

      void calc(char,const PEPS< double > &,const Walker &walker);

      void test(char);

      void add_layer(const char,int,const PEPS<double> &);

   private:

      //!stores an array environment MPS's for l(eft) , r(ight), t(op) and b(ottom)
      vector< MPS > l;
      vector< MPS > r;
      vector< MPS > t;
      vector< MPS > b;

      //!contraction between current walker state and peps
      SL_PEPS U;

      //!contraction between inverted walker state and peps
      SL_PEPS I;

      //!regular bond dimension of peps
      int D;

      //!Auxiliary dimension, for the contraction
      int D_aux;

      //!nr of sweeps in compression
      int comp_sweeps;

};

#endif
