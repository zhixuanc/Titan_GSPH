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
  double dwdx[DIMENSION];
  double dx[DIMENSION], s[DIMENSION];
  double gradU[NOEQxDIM];
  double AA[DIMENSION][DIMENSION], dd[DIMENSION][NO_OF_EQNS];
  double delu[NO_OF_EQNS], ff[DIMENSION][NO_OF_EQNS];
  double dia1[DIMENSION], dia2[DIMENSION];
  vector<Key> neighs;
  static int icount = 0;

  HTIterator *itr = new HTIterator(P_table);
  Particle *pi;
  while ((pi=static_cast<Particle *>(itr->next()))) 
    if ( pi->is_real() )
    {
      const double *xi = pi->get_coords();
      const double *ui= pi->get_state_vars();
      double h = pi->get_smlen();
      double supp = 3*h;
      int nreal=0;
      neighs = pi->get_neighs();      

      // Initialize summation variables to ZERO
      for (i = 0; i < DIMENSION; i++)
      {
        for (j = 0; j < DIMENSION; j++)
          AA[i][j] = 0;
        for (j = 0; j < NO_OF_EQNS; j++)
        {
          dd[i][j] = 0;
          gradU[i * NO_OF_EQNS + j] = 0;
        }
        dia1[i] = 0.;
        dia2[i] = 0.;
      }
     
      // iterate over all the (real) neighbors
      vector<Key>::iterator p_itr;
      for (p_itr=neighs.begin(); p_itr!=neighs.end(); p_itr++)
      {
        Particle *pj=static_cast<Particle *>(P_table->lookup(*p_itr));
        if (!pj)
        {
          fprintf(stderr,"something wrong at:\n");
          fprintf(stderr,"pi: (%f, %f), pj: (%u, %u) \n",
                  xi[0], xi[1], neighs[j].key[0], neighs[j].key[1]);
          fprintf(stderr,"neigh size: %d\n", neighs.size());
          return (-2);
        }

        if ( *pi == *pj ) continue;

        if ( !pj->is_ghost() )
        {
          nreal++;
          const double *xj = pj->get_coords();
          for (i=0; i<DIMENSION; i++)
          {
            dx[i]  = xj[i] - xi[i];
            s[i]   = -dx[i]/h;
          }
          if ( in_support(dx, supp) )
          {
            const double *uj = pj->get_state_vars();
            double volj = pj->get_mass()/uj[0];

            // pre-compute weight-gradients
            for (i=0; i<DIMENSION; i++)
              dwdx[i] = d_weight(s,h,i);

           // pre-compute Uj - Ui
            for (i = 0; i < NO_OF_EQNS; i++)
              delu[i] = uj[i] - ui[i];

            // form the DX matrix
            for (i = 0; i < DIMENSION; i++)
              for (j = 0; j < DIMENSION; j++)
                AA[i][j] += volj * dwdx[i] * dx[j];

            // form Df matrix
            for (i = 0; i < DIMENSION; i++)
              for (j = 0; j < NO_OF_EQNS; j++)
                dd[i][j] += volj * dwdx[i] * delu[j];

            for (i =0; i < DIMENSION; i++)
            {
              dia1[i] += volj * delu[i+1] * dwdx[i];
              dia2[i] += volj * dx[i] * dwdx[i];
            }
          }
        }
      }
  
      linsolve(&AA[0][0], &dd[0][0], NO_OF_EQNS, &ff[0][0]);

      // check and save derivatives
      for (i = 0; i < NO_OF_EQNS; i++)
        for (j = 0; j < DIMENSION; j++)
          gradU[i * DIMENSION + j] = ff[j][i];

      for (i = 0; i < DIMENSION; i++)
        gradU[i * NO_OF_EQNS + DIMENSION] = dia1[i] / dia2[i];

      //  if there is a problem ... then die a violent death
      for (i = 0; i < NOEQxDIM; i++)
        if (isnan(gradU[i]) && (nreal > 5))
        {
          fprintf(stderr,"FATAL ERROR: calc_gradients() failed\n");
          fprintf(stderr,".. at (%f, %f), no of neighbors: %d \n", 
                          xi[0], xi[1], nreal);
#ifdef DEBUG
          vector<Key> neighs = pi->get_neighs();
          vector<Key>::iterator pitr;
          for (pitr=neighs.begin(); pitr!=neighs.end(); pitr++)
          {
            Particle *pj=static_cast<Particle *>(P_table->lookup(*pitr));
            if (*pi == *pj) continue;
            if (!pj->is_ghost())
              fprintf(stderr,"%e, %e, %e, %e\n",
                      (xi[0]-*(pj->get_coords())),
                      (xi[1]-*(pj->get_coords()+1)),
                      (ui[1]-*(pj->get_vel())),
                      (ui[2]-*(pj->get_vel()+1)));
          }
#endif
          return -1;
        }

      // if particle doesn't have many neighbors, 
      // it doesn't have gradients
      if (nreal <= 5)
        for (i = 0; i < NOEQxDIM; i++)
          gradU[i] = 0.;

      // update values of slopes
      pi->put_d_state_vars (gradU);
    }
  
  // clean up stuff
  delete itr;
  return 0;
}
