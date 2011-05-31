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
mom_update (int myid, HashTable * P_table, HashTable * BG_mesh,
            MatProps * matprops, double dt)
{
  int i, j, k;
  vector < Key > neighs;
  double dx[DIMENSION], xi[DIMENSION], si[DIMENSION], sj[DIMENSION];
  double dx2[DIMENSION];
  double ustar[NO_OF_EQNS], rhs[DIMENSION], unew[NO_OF_EQNS];
  double dwdx[DIMENSION], uvecj[DIMENSION];
  double bndcrd1[DIMENSION], bndcrd2[DIMENSION];
  double divergU, bedfrict[DIMENSION], dudx[DIMSQRD];
  double beta = matprops->bedfrict;
  const double gravity[2] = { 0, -1. };
  double tauxy[DIMENSION];
  Particle *pi; 

  // three point gauss-quadrature points
  double gqp[3] = {-0.7746, 0, 0.7746 };
  double gqw[3] = {0.5556, 0.8889, 0.5556 };

  // friction data
  double intfrict = (1. - matprops->tanintfrict);
  double sinintfrict = matprops->sinintfrict;

  HTIterator *itr = new HTIterator (BG_mesh);
  Bucket *buck;

  // go through bucket table
  while ( (buck = (Bucket *) itr->next ()) )
    if ( buck->is_active() &&  !(buck->is_guest()) )
    {
      Key *neighbors = buck->get_neighbors();
      vector<Key> plist = buck->get_plist();
      vector<Key>::iterator ip;
      for (ip=plist.begin(); ip != plist.end(); ip++)
      {
        pi = (Particle *) P_table->lookup (*ip);
        if ( pi->is_real() )
        {
          for (i = 0; i < DIMENSION; i++)
            xi[i] = *(pi->get_coords() + i);

          // expanded smoothing length for Momentum equation
	  double h = pi->get_smlen();
	  double hi = sqrt2*h;
	  double h3 = 3*h;
	  double supp = 3*hi;

	  const double *uvec = pi->get_state_vars();
	  // density must always be positive
	  assert (uvec[0] > 0);
	  double Vi = 1.0/uvec[0];

	  // reset rhs to zero, for current particle
	  for (i = 0; i < DIMENSION; i++)
          {
	    rhs[i] = 0;
            tauxy[i] = 0;
	    bedfrict[i] = 0;
          }

	  // compute bed-friction
	  if ( buck->get_bucket_type() == MIXED )
          {
            double bslope = -1.0/buck->get_bndnorm(xi);
            double cost   = 1./sqrt(1. + pow(bslope, 2.));
            double sint   = cost*bslope;
            int num_bp = PARTICLE_DENSITY-1;
            for ( j=0; j < num_bp; j++)
            {
              bndcrd1[0] = *(buck->get_bnd_xcrd()+j);
              bndcrd1[1] = buck->get_bndZ(bndcrd1);
              bndcrd2[0] = *(buck->get_bnd_xcrd()+j+1);
              bndcrd2[1] = buck->get_bndZ(bndcrd2);
              for (k=0; k<DIMENSION; k++)
              {
                dx[k] = xi[k] - bndcrd1[k];
		dx2[k]= xi[k] - bndcrd2[k];
              }
              if ( in_support(dx, h3) && in_support(dx2, h3) )
              {
                double x_a = bndcrd1[0];
                double f_a = *(buck->get_f_coef()+j);
                double x_b = bndcrd2[0];
                double f_b = *(buck->get_f_coef()+j);
                // form linear function of f_coef
                double c1 = (f_b - f_a)/(x_b - x_a);
                double c2 = (f_a*x_b - f_b*x_a)/(x_b - x_a);
                  
                // integrate f(z)w(z) over the segment 
  	        double intgval = 0.;
                for (k=0; k<3; k++)
                {
                  double gqx = ((x_b - x_a)/2.)*gqp[k] + (x_a + x_b)/2.;
                  double gqz = buck->get_bndZ(&gqx);
                  si[0] = (xi[0]-gqx)/hi;
                  si[1] = (xi[1]-gqz)/hi;
                  double w = weight(si,hi);
                  double f = c1*gqx + c2;
                  intgval += gqw[k]*f*w*(x_b - x_a)/2.;
                }
                bedfrict[0] += intgval*cost;
                bedfrict[1] += intgval*sint;
              }
            }
          }
	  for (i=0; i < NEIGH_SIZE; i++)
	    if ( *(buck->get_neigh_proc()+i) > -1 )
	    {
	      Bucket *buck_neigh = (Bucket *) BG_mesh->lookup(neighbors[i]);
              if ( !(buck_neigh) && (*(buck->get_neigh_proc()+i)) != myid)
                continue;
              assert(buck_neigh);
	      if ( buck_neigh->get_bucket_type() == MIXED )
              {
		// get dirction cosines
		double bslope = -1.0/buck_neigh->get_bndnorm(xi);
		double cost   = 1./sqrt(1. + pow(bslope, 2.));
		double sint   = cost*bslope;
		// get dimensions of boundary particles

		int num_bp = PARTICLE_DENSITY-1;
		for ( j=0; j < num_bp; j++)
                {
                  bndcrd1[0] = *(buck_neigh->get_bnd_xcrd()+j);
                  bndcrd1[1] = buck_neigh->get_bndZ(bndcrd1);
                  bndcrd2[0] = *(buck_neigh->get_bnd_xcrd()+j+1);
                  bndcrd2[1] = buck_neigh->get_bndZ(bndcrd2);
		  for (k=0; k<DIMENSION; k++)
                  {
		    dx[k] = xi[k] - bndcrd1[k];
		    dx2[k]= xi[k] - bndcrd2[k];
                  }
		  if ( in_support(dx, h3) && in_support(dx2, h3) )
                  {
		    double x_a = bndcrd1[0];
		    double f_a = *(buck_neigh->get_f_coef()+j);
		    double x_b = bndcrd2[0];
		    double f_b = *(buck_neigh->get_f_coef()+j);
		    // form linear function of f_coef
		    double c1 = (f_b - f_a)/(x_b - x_a);
		    double c2 = (f_a*x_b - f_b*x_a)/(x_b - x_a);
                  
		    // integrate f(z)w(z) over the segment 
		    double intgval = 0.;
		    for (k=0; k<3; k++)
                    {
		      double gqx = ((x_b - x_a)/2.)*gqp[k] + (x_a + x_b)/2.;
		      double gqz = buck_neigh->get_bndZ(&gqx);
		      si[0] = (xi[0]-gqx)/hi;
		      si[1] = (xi[1]-gqz)/hi;
		      double w = weight(si,hi);
		      double f = c1*gqx + c2;
		      intgval += gqw[k]*f*w*(x_b - x_a)/2.;
                    }
		    bedfrict[0] += intgval*cost;
		    bedfrict[1] += intgval*sint;
                  }
		}
	      }
	    }
	  if ( isnan(bedfrict[0]) )
	  {
            fprintf(stderr,"ERROR: Getting NaN's for bedfriction values\n");
            return -2;
          }
          pi->put_bedfrict(bedfrict);

	  // list of neighbors
	  neighs = pi->get_neighs ();
	  int no_of_neighs = neighs.size ();
	  for (j = 0; j < no_of_neighs; j++)
	  {
	    Particle *pj = (Particle *) P_table->lookup (neighs[j]);
	    assert (pj);

	    // self contribution is zero as dw(0)=0
	    if (*pi == *pj)
	      continue;
	    double dist = 0;
	    for (i = 0; i < DIMENSION; i++)
	    {
	      dx[i] = xi[i] - *(pj->get_coords () + i);
	      si[i] = dx[i] / hi;
	      dist += dx[i] * dx[i];
            }
	    dist = sqrt (dist);
	    if (in_support (dx, supp))
            {
	      double mj = pj->get_mass ();
	      double Vj = 1.0 / pj->get_density ();

	      // solve the Riemann problem between to particles 
	      riemann_solve (ustar, pi, pj, matprops, dt);
	      if ( ustar[0] <= 0 )
              {
                fprintf(stderr," FATAL ERROR: rho* = %f < 0 \n",ustar[0]);
                return -1;
              }
	      double pstar = matprops->pressure (ustar[0]);

	      // Cij and Dij linear interpolation constants 
	      // variable names used in Inutsuka 2003 paper
	      double Cij = (Vi - Vj) / dist;
	      double Dij = (Vi + Vj) / 2.0;
	      double Vij = pow(pi->get_smlen()*Cij*0.5,2.) 
		+ pow(Dij,2.);

	      // pre-compute weight function derivatives
	      for (k = 0; k < DIMENSION; k++)
		dwdx[k] = d_weight (si, hi, k);

	      // Velocity update 
	      rhs[0] += -2*mj*pstar*Vij*dwdx[0]; 
	      rhs[1] += -2*mj*pstar*Vij*dwdx[1];
 
              // shear tractions
              tauxy[0] += 2*mj*Vij*pstar*dwdx[1]*sinintfrict;
              tauxy[1] += 2*mj*Vij*pstar*dwdx[0]*sinintfrict;
	    }
	  } // end loop over neighs
	  //  calculate divergence 
          for (k=0; k < DIMSQRD; k++)
	    dudx[k] = *(pi->get_d_vel() + k);

	  divergU = 0.;
	  for (k = 0; k < DIMENSION; k++)
	    divergU += dudx[k*DIMENSION +k];

	  // update state variables
	  // density
	  unew[0] = uvec[0] - dt*uvec[0]*divergU;
	  if (unew[0] <= 0 || isnan (unew[0]) )
	  {
	    fprintf(stderr, "FATAL ERROR: Negative/Zero density\n");
            fprintf(stderr, "at (%f, %f) = %f\n",xi[0], xi[1], unew[0]);
	    return (-1);
	  }
          
          tauxy[0] = 0; tauxy[1] = 0;
	  //  x-velocity
	  unew[1] = uvec[1] + dt*(rhs[0] - sgn(dudx[1])*tauxy[0] +
                              gravity[0] + bedfrict[0]);

	  // y-velocity
	  unew[2] = uvec[2] + dt*(rhs[1] - sgn(dudx[2])*tauxy[1] +
                              gravity[1] + bedfrict[1]);

	  pi->put_new_state_vars (unew);
	}
      }
    }

  // iterate over hashtable to update state_variables
  HTIterator *it2 = new HTIterator (P_table);
  while ( (pi = (Particle *) it2->next ()) )
    if (pi->is_real ())
      pi->update_state_vars ();

   // clean up
  delete itr, it2;
  return 0;
}
