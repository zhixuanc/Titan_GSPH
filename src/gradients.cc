/*
 * =====================================================================================
 *
 *       Filename:  gradients.cc
 *
 *    Description:  
 *
 *        Created:  04/06/2010 01:22:38 PM EDT
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *        License:  GNU General Public License
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * =====================================================================================
 * $Id:$
 */


#include <vector>
#include <cassert>
#include <cstdio>
using namespace std;

#include <hashtab.h>
#include <particle.h>

#include "constants.h"
#include "sph_header.h"

int calc_gradients(HashTable *P_table)
{
  int i,j,k;
  const int SIZE = DIMENSION;
  double dwdx[DIMENSION];
  double dx[DIMENSION], s[DIMENSION];
  double gradU[DIMENSION*DIMENSION], d2U[DIMENSION*DIMENSION];
  double AA[SIZE][SIZE], dd[SIZE][DIMENSION];
  double ff[SIZE][DIMENSION], V1[SIZE], V2[SIZE];
  double t1[DIMENSION], t2[DIMENSION];
  vector<Key> neighs;
  static int icount = 0;

  HTIterator *itr = new HTIterator(P_table);
  Particle *pi;
  while ( pi = (Particle *) itr->next() ) 
    if ( pi->is_real() )
    {
      const double *xi = pi->get_coords();
      const double *vi= pi->get_vel();
      double h = pi->get_smlen();
      double supp = 3*h;
      neighs = pi->get_neighs();      
      int num_neigh = neighs.size();

      // Initialize summation variables to ZERO
      for (i=0; i<DIMSQRD; i++)
        gradU[i]=0.;
     
      for (i=0; i<DIMENSION; i++)
      {
        t1[i] = 0.;
        t2[i] = 0.;
      }
      // iterate over all the (real) neighbors
      for (j=0; j<num_neigh; j++)
      {
        Particle *pj = (Particle *) P_table->lookup(neighs[j]);
        if ( !pj )
        {
          printf("something wrong at:\n");
          printf("pi: (%f, %f), pj: (%u, %u) \n",
                  xi[0], xi[1], neighs[j].key[0], neighs[j].key[1]);
          printf("neigh size: %d\n", num_neigh);
          return (-2);
        }

        if ( *pi == *pj ) continue;
        if ( !pj->is_ghost() )
        {
          const double *xj = pj->get_coords();
          for (i=0; i<DIMENSION; i++)
          {
            dx[i]  = xj[i] - xi[i];
            s[i]   = -dx[i]/h;
          }
          if ( in_support(dx, supp) )
          {
            double rj = pj->get_density();
            double tmp = pj->get_mass()/rj;
            const double *vj = pj->get_vel();
            for (i=0; i<DIMENSION; i++)
            {
              dwdx[i] = d_weight(s,h,i);
              t1[i] += (vj[i]-vi[i])*tmp*dwdx[i];
              t2[i] += dx[i]*tmp*dwdx[i];
            }
            for (i=0; i<DIMENSION; i++)
              for (k=0; k<DIMENSION; k++)
              {
                AA[i][k] += dx[k]*dwdx[i]*tmp;
                dd[i][k] += (vj[k]-vi[k])*tmp;
              }
          }
        }
      }
  
      linsolve(&AA[0][0], DIMENSION, &dd[0][0], DIMENSION, &ff[0][0]);
      for (i=0; i<DIMENSION; i++)
        for (k=0; k<DIMENSION; k++)
          gradU[i*DIMENSION + k] = ff[i][k];

      for (i=0; i<DIMENSION; i++)
        gradU[i*DIMENSION+i] = t1[i]/t2[i];

      for (i=0; i<DIMSQRD; i++)
        if ( isnan(gradU[i]) )
        {
          fprintf(stderr,"FATAL ERROR: calc_gradients() failed\n");
          fprintf(stderr,".. at (%f, %f), no of neighbors: %d \n", 
                          xi[0], xi[1], num_neigh);
          return -1;
        }
      // update values of slopes
      pi->put_d_vel(gradU);
    }
  
  // clean up stuff
  delete itr;
  return 0;
}
