
/*
 * =====================================================================================
 *
 *       Filename:  reimann_ivars.cc
 *
 *    Description:  
 *
 *        Created:  04/27/2009 02:52:44 PM EDT
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

#include <cmath>
using namespace std;

#include <properties.h>
#include <particle.h>

#include "constants.h"
#include "sph_header.h"
#include "riemann_solve.h"

#define SECOND_ORDER
double 
riemann_solve (Particle * pi, Particle * pj, 
              MatProps * matprops, double dt)
{

  int i;
  double dU[2];                  // Magnitude of discontinuity
  double U1[2], U2[2];           // density and speed in e_{ij} direction
  double Uleft[2], Uright[2];    // left and right states across discontinuity
  double dUidx[2], dUjdx[2];     // gradients for density and velocity
  double dist, dx[DIMENSION], dir[DIMENSION];
  double Ui[NO_OF_EQNS], Uj[NO_OF_EQNS];   // state variables

  // get state-variables
  for (i = 0; i < NO_OF_EQNS; i++)
  {
    Ui[i] = (*(pi->get_state_vars() + i));
    Uj[i] = (*(pj->get_state_vars() + i));
  }

  // get vector xi - xj, and distance
  dist = 0;
  for (i = 0; i < DIMENSION; i++)
  {
    dx[i] = (*(pi->get_coords() + i) - *(pj->get_coords() + i));
    dist += dx[i] * dx[i];
  }

  // compute n_{ij} unit vector
  dist = sqrt (dist);
  for (i = 0; i < DIMENSION; i++)
    dir[i] = dx[i] / dist;

  /*
  // get gradients of state_vars
  for (i = 0; i < 2; i++)
  {
    dUidx[i] = dot ((double *) pi->get_d_state_vars () + i * DIMENSION, dir);
    dUjdx[i] = dot ((double *) pj->get_d_state_vars () + i * DIMENSION, dir);
  }
  */
  // rotate state-vars along n_{ij}
  U1[0] = Ui[0];
  U1[1] = dot (Ui + 1, dir);

  U2[0] = Uj[0];
  U2[1] = dot (Uj + 1, dir);

  for (i = 0; i < 2; i++)
  {
    dUidx[i] = (U1[i] - U2[i]) / dist;
    dUjdx[i] = (U1[i] - U2[i]) / dist;
  }

  // setup Reimann problem
  double hi = pi->get_smlen();
  double ci = matprops->sound_speed(pi->get_density());
  double hj = pj->get_smlen();
  double cj = matprops->sound_speed(pj->get_density());

  // specific volumes at particle i and j
  double vi = 1.0 / U1[0];
  double vj = 1.0 / U2[0];

  // linear interpolation of density 
  double Cij = (vi - vj) / dist;
  double Dij = (vi + vj) / 2.0;

  // find out specific volume contributions
  // NOTE: Vij's are squared quantities
  double vij = (0.5 * (hi * Cij * hi *Cij)) + (Dij * Dij);

  // sstar is the point where we setup RP
  double sstar = hi * hi * Cij * Dij / (2 * vij);
  double dsi = sstar + ci * dt / 2.0 - dist / 2.0;
  double dsj = sstar - cj * dt / 2.0 + dist / 2.0;

  // polynomial reconstruction to setup RP
  for (i = 0; i < 2; i++)
  {
    Uright[i] = U1[i] + dUidx[i] * dsi;
    Uleft[i] = U2[i] + dUjdx[i] * dsj;
  }

  for (i = 0; i < 2; i++)
    dU[i] = Uright[i] - Uleft[i];

  /* ***********************************
   *    Solve the Riemann problem
   * ***********************************
   */
  double eval[2], temp[2], t1, c;

  // -- using left-side data
  c = matprops->sound_speed (Uleft[0]);
  eval[0] = Uleft[1] - c;
  eval[1] = Uleft[1] + c;
  t1 = Uleft[0] / c;
  temp[0] = 0.5 * (dU[0] - dU[1] * t1);
  temp[1] = 0.5 * (dU[0] + dU[1] * t1);

  double rhol = Uleft [0];
  for (i = 0; i < 2; i++ )
    if ( eval[i] < 0 )
      rhol += temp[i];

  // -- using right-side data
  c = matprops->sound_speed (Uright[0]);
  eval[0] = Uright[1] - c;
  eval[1] = Uright[1] + c;
  t1 = Uright[0] / c;
  temp[0] = 0.5 * (dU[0] - dU[1] * t1);
  temp[1] = 0.5 * (dU[0] + dU[1] * t1);

  double rhor = Uright [0];
  for (i = 0; i < 2; i++ )
    if ( eval[i] > 0 )
      rhor -= temp[i];

  // Although they should be equal, we'll average them anyway
  return  (0.5 * (rhol + rhor));
}
