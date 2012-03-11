
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

/* TODO rotations in 3-D   */

#include <cmath>
using namespace std;

#include <properties.h>
#include <particle.h>

#include "constants.h"
#include "sph_header.h"
#include "riemann_solve.h"

int
riemann_solve(double ustar[], Particle * p1, Particle * p2, MatProps * mp_ptr,
              double dt)
{

  int i, j, N, ierr;
  double eval[NO_OF_EQNS];      // Eigenvalues
  double evec[NO_OF_EQNS][NO_OF_EQNS];  // Eigenvectors
  double dU[NO_OF_EQNS];        // Magnitude of discontinuity
  double alph[NO_OF_EQNS];      // jumps along Eigenvectors
  double Usl[NO_OF_EQNS], Usr[NO_OF_EQNS];      // Rotated State Variables
  double U1[NO_OF_EQNS], U2[NO_OF_EQNS],        // state variables
    sum_left[NO_OF_EQNS], sum_right[NO_OF_EQNS], Uleft[NO_OF_EQNS], Uright[NO_OF_EQNS]; //reimann states
  double dist, dx[DIMENSION], dir[DIMENSION];

  for (i = 0; i < NO_OF_EQNS; i++)
  {
    U1[i] = (*(p1->get_state_vars() + i));
    U2[i] = (*(p2->get_state_vars() + i));
  }
  dist = 0;
  for (i = 0; i < DIMENSION; i++)
  {
    dx[i] = (*(p1->get_coords() + i) - *(p2->get_coords() + i));
    dist += dx[i] * dx[i];
  }
  dist = sqrt(dist);
  for (i = 0; i < DIMENSION; i++)
    dir[i] = dx[i] / dist;

  // rotate in the local coordinate direction
  rotate(U1, dir);
  rotate(U2, dir);

  // setup Reimann problem
  double hi = p1->get_smlen();
  double ci = mp_ptr->sound_speed(p1->get_density());
  double hj = p2->get_smlen();
  double cj = mp_ptr->sound_speed(p2->get_density());

  setup_rp(U1, U2, Uleft, Uright, dist, hi, hj, dt, ci, cj);

  for (i = 0; i < NO_OF_EQNS; i++)
    dU[i] = Uright[i] - Uleft[i];

  /* ***********************************
   *    Solve the Riemann problem
   * ***********************************
   */

  // -- using left-side data
  for (i = 0; i < NO_OF_EQNS; i++)
    sum_left[i] = 0;

  // eigen structure
  double cl = mp_ptr->sound_speed(Uleft[0]);

  ierr = eigen_decomp(Uleft[0], Uleft[1], cl, dU, eval, alph, evec);

  for (i = 0; i < NO_OF_EQNS; i++)
    if (eval[i] < 0)
      for (j = 0; j < NO_OF_EQNS; j++)
        sum_left[j] += alph[i] * evec[j][i];

  // evaluate star values.
  for (i = 0; i < NO_OF_EQNS; i++)
    Usl[i] = Uleft[i] + sum_left[i];

  // -- using right-side data
  for (i = 0; i < NO_OF_EQNS; i++)
    sum_right[i] = 0;

  double cr = mp_ptr->sound_speed(Uright[0]);

  ierr = eigen_decomp(Uright[0], Uright[1], cr, dU, eval, alph, evec);
  for (i = 0; i < NO_OF_EQNS; i++)
    if (eval[i] > 0)
      for (j = 0; j < NO_OF_EQNS; j++)
        sum_right[j] += alph[i] * evec[j][i];

  for (i = 0; i < NO_OF_EQNS; i++)
    Usr[i] = Uright[i] - sum_right[i];

  // Although they should be equal, we'll average them anyway
  for (i = 0; i < DIMENSION; i++)
    ustar[i] = 0.5 * (Usl[i] + Usr[i]);

  return 0;
}
