#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Walker.h"

/**
 * @author Brecht Verstichel
 * @date 22-07-2014\n\n
 * The class is written to easily draw from non-smooth/analytical distributions.
 * It contains a list of Walker states which encode the information of the distribution.
 */
class Distribution : public vector<double> {

   friend ostream &operator<<(ostream &output,const Distribution &dist_p);

   public:

      Distribution();

      Distribution(const Distribution &);

      virtual ~Distribution();

      void construct(const Walker &);

      int gn() const;

      const vector< Walker > &glist() const;

      const Walker &gwalker(int) const;

      Walker &gwalker(int);

      int metropolis() const;

      double energy() const;

   private:

      //!list of final walker states
      vector< Walker > list;

            
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
