#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

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

      void construct(const Walker &,double,double);

      int gn() const;

      double normalize();

      int draw() const;

      const vector< vector<bool> > &gconf() const;
      const vector<int> &gind() const;
      const vector<bool> &gsign_flip() const;

      const vector<bool> &gconf(int) const;
      int gind(int) const;
      bool gsign_flip(int) const;

      void check_negative() const;

   private:

      //!list of final walker configurations
      vector< vector<bool> > conf;

      //!list of final walker indices (row,col,dir)
      vector<int> ind;

      //!list of sign flips, false is no sign flip
      vector<bool> sign_flip;
      
};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
