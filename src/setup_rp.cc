/*
 * =====================================================================================
 *
 *       Filename:  setup_rp.c
 *
 *    Description:  sets up Riemann problem between to particles
 *
 *        Created:  05/04/2008 01:24:19 PM EDT
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

#include <math.h>

#include "constants.h"

void setup_rp (double U1[], double U2[], double Ul[], double Ur[], double dist,
               double h1, double h2, double dt, double ci, double cj)
{
  double dUds;
  int i,ierr;

  // specific volumes at particle i and j
  double vi=1.0/U1[0];
  double vj=1.0/U2[0];
  
  // linear interpolation of density 
  double Cij=(vi-vj)/dist;
  double Dij=(vi+vj)/2.0;

  // find out specific volume contributions
  // NOTE: Vij's are squared quantities
  double vij = pow(h1*Cij*0.5,2)+(Dij*Dij);

  // sstar is the point where we setup RP
  double sstar = h1*h1*Cij*Dij/(2*vij);
  double dsi=sstar + ci*dt/2.0 - dist/2.0;
  double dsj=sstar - cj*dt/2.0 + dist/2.0;

  // Setup Riemann problem at s*
  for ( i=0; i<NO_OF_EQNS; i++)
  {
    dUds=(U1[i]-U2[i])/dist;
    Ur[i]=U1[i] + dUds*dsi;
    Ul[i]=U2[i] + dUds*dsj;
  }
  return;
}
