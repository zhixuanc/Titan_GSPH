
/*
 * =====================================================================================
 *
 *       Filename:  partilce.cc
 *
 *    Description:  
 *
 *        Created:  01/07/2008 01:41:12 PM EST
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
 *        License:  GPL v2
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
#include <iostream>
using namespace std;

#include "particle.h"
#include "constants.h"

// general purpose contstructor
Particle::Particle (unsigned *keyin, double *crd, double m, double h, int ptype)
{
  int i;

  mass = m;
  smlen = h;
  update_delayed = false;
  if (ptype == 1)
    ghost = true;
  else
    ghost = false;
  guest = false;
  reflection = false;

  for (i = 0; i < KEYLENGTH; i++)
    key.key[i] = *(keyin + i);

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = *(crd + i);
    bedfrict[i] = 0;
  }

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0;

  for (i = 0; i < DIMSQRD; i++)
    d_vel[i] = 0;

  state_vars[0] = 1.0;
  return;
}

// Contructor for flux particles from flux source
Particle::Particle (unsigned keyin[], double crd[], double ma, double hl,
                    double vel[])
{
  int i;

  mass = ma;
  smlen = hl;
  ghost = false;
  guest = false;
  reflection = false;
  state_vars[0] = 1.0;

  for (i = 0; i < KEYLENGTH; i++)
    key.key[i] = keyin[i];

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = crd[i];
    state_vars[1 + i] = vel[i];
    bedfrict[i] = 0;
    d_vel[i] = 0;
  }
}

Particle::Particle ()
{
  int i;

  for (i = 0; i < KEYLENGTH; i++)
    key.key[i] = 0;

  mass = 0.;
  smlen = 0.;
  guest = false;
  reflection = false;

  for (i = 0; i < DIMENSION; i++)
  {
    coord[i] = 0.;
    bedfrict[i] = 0.;
  }

  for (i = 0; i < NO_OF_EQNS; i++)
    state_vars[i] = 0.;

  for (i = 0; i < DIMSQRD; i++)
    d_vel[i] = 0.;
}

bool Particle::operator== (const Particle & rhs) const
{
  for (int i = 0; i < KEYLENGTH; i++)
    if (key.key[i] != (rhs.getKey ()).key[i])
      return false;
  return true;
}
