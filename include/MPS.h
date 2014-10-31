#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

/**
 * @author Brecht Verstichel
 * @date 12-06-2014\n\n
 * This class MPS is a class written for the construction of matrix product states without symmetry.
 * More specifically it will be used for the contraction of PEPS networks. Where the reduction to a MPS-like form is done.
 */
class MPS : public vector< TArray<double,3> > {

   public:

      MPS();

      MPS(int);

      MPS(int,int,int);

      //copy constructor
      MPS(const MPS &);

      //destructor
      virtual ~MPS();

      int gD() const;

      int gd() const;

      void fill_Random();

      void normalize();

      void canonicalize(const BTAS_SIDE &,bool);

      void cut_edges();

      void guess(const BTAS_SIDE &,int ,const MPS &mps);

      void compress(int ,const MPS &mps,int);

      double dot(const MPS &bra) const;

      void fill(char,const SL_PEPS &);

      void gemv(char,char,int,const SL_PEPS &);

   private:

      //!dimension of the bonds
      int D;

      //!physical dimension
      int d;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
