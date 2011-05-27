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

void eigen_decomp(double rho, double u, double c, double dU[],
                  double lambda[], double alpha[], double evec[][NO_OF_EQNS])
{
  if (NO_OF_EQNS != 3)
  {
    fprintf(stderr,"ERROR: No. of Equations in eigen.c,\
            don't match with No. of eqations in model\n");
    exit(1);
  }
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
