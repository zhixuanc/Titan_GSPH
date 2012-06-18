/*
 * =====================================================================================
 *
 *       Filename:  eigen.c
 *
 *    Description:  
 *
 *        Created:  03/25/2010 03:21:42 PM EDT
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

#include <stdio.h>
#include <stdlib.h>

#include "constants.h"

void eigen2d (double rho, double u, double c, double dU[],
              double lambda[], double alpha[], double evec[][3])
{
  double t2, t3;

  // C(R, optimize, resultname = "lambda");
  lambda[0] = u - c;
  lambda[1] = u + c;
  lambda[2] = u;

  // C(V, optimize, resultname = "evec");
  t2 = rho/c;
  evec[0][0] = -t2;
  evec[0][1] = t2;
  evec[0][2] = 0;
  evec[1][0] = 1;
  evec[1][1] = 1;
  evec[1][2] = 0;
  evec[2][0] = 0;
  evec[2][1] = 0;
  evec[2][2] = 1;

  //C(del, optimize, resultname = "alpha");
  t3 = c / rho * dU[0];
  alpha[0] = -t3/0.2e1 + dU[1]/0.2e1;
  alpha[1] =  t3/0.2e1 + dU[1]/0.2e1;
  alpha[2] =  dU[2];

  return;
}

/* the code is generated in maple15 */
void eigen3d (double rho, double u, double c, double dU[],
              double lambda[], double alpha[], double evec[][4])
{
  if (NO_OF_EQNS != 4)
  {
    fprintf(stderr,"ERROR: No. of Equations in eigen.c,\
            don't match with No. of eqations in model\n");
    exit(1);
  }

  // eigen-values
  lambda[0] = u;
  lambda[1] = u;
  lambda[2] = u + c;
  lambda[3] = u - c;

  // eigen-vectors
  double t1 = rho / c;
  evec[0][0] = 0;
  evec[0][1] = 0;
  evec[0][2] = t1;
  evec[0][3] = -t1;
  evec[1][0] = 0;
  evec[1][1] = 0;
  evec[1][2] = 1;
  evec[1][3] = 1;
  evec[2][0] = 0;
  evec[2][1] = 1;
  evec[2][2] = 0;
  evec[2][3] = 0;
  evec[3][0] = 1;
  evec[3][1] = 0;
  evec[3][2] = 0;
  evec[3][3] = 0;

  // projection of jump along the eigen-vectors
  t1 = c / rho * dU[0];
  alpha[0] = dU[3];
  alpha[1] = dU[2];
  alpha[2] = t1 / 2.0 + dU[1] / 2.0;
  alpha[3] = -t1 / 2.0 + dU[1] / 2.0;

  return;
}
