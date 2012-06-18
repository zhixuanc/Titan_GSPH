
/*
 * =====================================================================================
 *
 *       Filename:  mom_update.cc
 *
 *    Description:  
 *
 *        Created:  04/09/2009 05:13:09 PM EDT
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
#include <cmath>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <bucket.h>

#include "particle.h"
#include "constants.h"
#include "sph_header.h"
#include "riemann_solve.h"

const double sqrt2 = 1.41421356237310;

int
mom_update(int myid, HashTable * P_table, HashTable * BG_mesh,
           MatProps * matprops, TimeProps * timeprops)
{
  int i, j, k;

  vector < Key > neighs;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION], sj[DIMENSION];
  double dx2[DIMENSION];
  double ustar, rhs[DIMENSION], unew[NO_OF_EQNS];
  double dwdx[DIMENSION], uvecj[DIMENSION];
  double bndcrd[DIMENSION];
  double f_coef[PARTICLE_DENSITY][PARTICLE_DENSITY][DIMENSION];
  double bedfrict[DIMENSION], dudx[DIMSQRD];
  double tauxy[DIMENSION][2];
  Particle *pi;

  // three point gauss-quadrature points
  double gqp[3] = { -0.7746, 0, 0.7746 };
  double gqw[3] = { 0.5556, 0.8889, 0.5556 };
  double gravity[3] = { 0., 0., -1.};

  // friction data
  double tanintfrict = matprops->tanintfrict;
  double sinintfrict = matprops->sinintfrict;

  // time-step
  double dt = timeprops->dtime;
  HTIterator *itr = new HTIterator(BG_mesh);
  Bucket *buck = NULL;

  // go through bucket table
  while ((buck = (Bucket *) itr->next()))
    if (buck->is_active() && !(buck->is_guest()))
    {
      Key *neighbors = buck->get_neighbors();

      vector < Key > plist = buck->get_plist();
      vector < Key >::iterator ip;
      for (ip = plist.begin(); ip != plist.end(); ip++)
      {
        pi = (Particle *) P_table->lookup(*ip);
        if (pi->is_real())
        {
          for (i = 0; i < DIMENSION; i++)
            xi[i] = *(pi->get_coords() + i);

          // expanded smoothing length for Momentum equation
          double h = pi->get_smlen();
          double hi = sqrt2 * h;
          double h3 = 3 * h;
          double supp = 3 * hi;

          const double *uvec = pi->get_state_vars();

          // density must always be positive
          assert(uvec[0] > 0);
          double Vi = 1.0 / uvec[0];

          // reset rhs to zero, for current particle
          double wnorm = 0.;

          for (i = 0; i < DIMENSION; i++)
          {
            rhs[i] = 0;
            bedfrict[i] = 0;
            for (j = 0; j < DIMENSION; j++)
              tauxy[i][j] = 0.;
          }

          // compute bed-friction
          if (buck->get_bucket_type() == MIXED)
          {
            for (i = 0; i < PARTICLE_DENSITY; i++)
              for (j = 0; j < PARTICLE_DENSITY; j++)
              {
                bndcrd[0] = *(buck->get_bnd_xcrd() + i);
                bndcrd[1] = *(buck->get_bnd_ycrd() + j);
                bndcrd[2] = buck->get_bndZ(bndcrd);

                for (k = 0; k < DIMENSION; k++)
                  dx[k] = xi[k] - bndcrd[k];
                if (in_support(dx, h3))
                {
                  for (k = 0; k < DIMENSION; k++)
                    si[k] = dx[k] / h;
                  double w = weight(si, h);

                  for (k = 0; k < DIMENSION; k++)
                    bedfrict[k] += buck->get_f_coef(i, j, k) * w;
                  wnorm += w;
                }
              }
          }
          for (int ineigh = 0; ineigh < NEIGH_SIZE; ineigh++)
            if (*(buck->get_neigh_proc() + ineigh) > -1)
            {
              Bucket *buck_neigh =
                (Bucket *) BG_mesh->lookup(neighbors[ineigh]);
              // If *guest*-neighbor is not there, it is not needed
              if ((!buck_neigh) && (*(buck->get_neigh_proc() + ineigh) != myid))
                continue;
              assert(buck_neigh);
              // Only if is is a boundary bucket
              if (buck_neigh->get_bucket_type() == MIXED)
              {
                // get dimensions of boundary particles
                for (i = 0; i < PARTICLE_DENSITY; i++)
                  for (j = 0; j < PARTICLE_DENSITY; j++)
                  {
                    bndcrd[0] = *(buck_neigh->get_bnd_xcrd() + i);
                    bndcrd[1] = *(buck_neigh->get_bnd_ycrd() + j);
                    bndcrd[2] = buck_neigh->get_bndZ(bndcrd);

                    for (k = 0; k < DIMENSION; k++)
                      dx[k] = xi[k] - bndcrd[k];
                    if (in_support (dx, h3))
                    {
                      for (k = 0; k < DIMENSION; k++)
                        si[k] = dx[k] / h;
                      double w = weight(si, h);

                      for (k = 0; k < DIMENSION; k++)
                        bedfrict[k] += buck->get_f_coef (i, j, k) * w;
                      wnorm += w;
                    }
                  }
              }
            }
          if ( wnorm > 0 )
          {
            for (i = 0; i < DIMENSION; i++)
              bedfrict[i] /= wnorm;
            if (isnan(bedfrict[0]) || isnan(bedfrict[1]) || isnan(bedfrict[2]))
            {
              fprintf(stderr, "ERROR: Getting NaN's for bedfriction values\n");
              return -2;
            }
          }
          pi->put_bedfrict(bedfrict);

          // list of neighbors
          vector < Key > neighs = pi->get_neighs();
          vector < Key >::iterator p_itr;
          for (p_itr = neighs.begin(); p_itr != neighs.end(); p_itr++)
          {
            Particle *pj = (Particle *) P_table->lookup(*p_itr);
            assert (pj);

            // self contribution is zero as dw(0)=0
            if (*pi == *pj)
              continue;
            double dist = 0;

            for (i = 0; i < DIMENSION; i++)
            {
              dx[i] = xi[i] - *(pj->get_coords() + i);
              si[i] = dx[i] / hi;
              dist += dx[i] * dx[i];
            }
            dist = sqrt(dist);
            // if dx < 3 * sqrt(2) * h
            if (in_support(dx, supp))
            {
              double mj = pj->get_mass();
              double Vj = 1.0 / pj->get_density();

              // solve the Riemann problem between to particles 
              ustar = riemann_solve(pi, pj, matprops, dt);
              if (ustar <= 0)
              {
                fprintf(stderr, " FATAL ERROR: rho* = %f < 0 \n", ustar);
                return -1;
              }
              double pstar = matprops->pressure (ustar);


              // Cij and Dij linear interpolation constants 
              // variable names used in Inutsuka 2003 paper
              double Cij = (Vi - Vj) / dist;
              double Dij = (Vi + Vj) / 2.0;
              double Vij = (0.5 * (h *h * Cij * Cij)) + (Dij * Dij);

              // pre-compute weight function derivatives
              for (k = 0; k < DIMENSION; k++)
                dwdx[k] = d_weight (si, hi, k);

              // Velocity update 
              rhs[0] += -2 * mj * pstar * Vij * dwdx[0];
              rhs[1] += -2 * mj * pstar * Vij * dwdx[1];
              rhs[2] += -2 * mj * pstar * Vij * dwdx[2];

            }
          } // end loop over neighs

          //  calculate divergence 
          for (k = 0; k < DIMSQRD; k++)
            dudx[k] = *(pi->get_d_vel() + k);

          double divergU = 0.;
          for (k = 0; k < DIMENSION; k++)
            divergU += dudx[k * DIMENSION + k];

          // update state variables
          // density
          unew[0] = uvec[0] - dt * uvec[0] * divergU;
          if (unew[0] <= 0 || isnan(unew[0]))
          {
            fprintf(stderr, "FATAL ERROR: Negative/Zero density\n");
            fprintf(stderr, "at (%f, %f) = %f\n", xi[0], xi[1], unew[0]);
            return (-1);
          }

          // stress - deviator
          /*
          tauxy[0][0] = -sgn(dudx[1]) * abs(rhs[1]);    // dp/dy
          tauxy[0][1] = -sgn(dudx[2]) * abs(rhs[2]);    // dp/dz
          tauxy[1][0] = -sgn(dudx[3]) * abs(rhs[0]);    // dp/dx
          tauxy[1][1] = -sgn(dudx[5]) * abs(rhs[2]);    // dp/dz
          tauxy[2][0] = -sgn(dudx[6]) * abs(rhs[0]);    // dp/dx
          tauxy[2][1] = -sgn(dudx[7]) * abs(rhs[1]);    // dp/dy
          */

          //  x-velocity
          unew[1] = uvec[1] + dt * (rhs[0] + gravity[0]);

          // y-velocity
          unew[2] = uvec[2] + dt * (rhs[1] + gravity[1]);

          // z-velocity
          unew[3] = uvec[3] + dt * (rhs[2] + gravity[2]);

          pi->put_new_state_vars(unew);
        }
      }
    }


  // iterate over hashtable to update state_variables
  HTIterator *it2 = new HTIterator(P_table);
  while ((pi = (Particle *) it2->next()))
    if (pi->is_real())
      pi->update_state_vars();

  // clean up
  delete itr, it2;

  return 0;
}
