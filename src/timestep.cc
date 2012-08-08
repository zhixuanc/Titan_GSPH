
/*
 * =====================================================================================
 *
 *       Filename:  timestep.cc
 *
 *    Description:  
 *
 *        Created:  06/19/2008 12:40:01 PM EDT
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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cmath>
#include <cassert>
using namespace std;

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <hashtab.h>
#include <particle.h>
#include <properties.h>

#include "constants.h"

double
timestep(HashTable * P_table, MatProps * matprops_ptr)
{
  int i, j, k;
  double xi[DIMENSION];
  double dt, temp;

  dt = 1.0E+10;                 // Initialize to very high value 
  HTIterator *itr = new HTIterator(P_table);
  Particle *p_curr = NULL;

  while ((p_curr = (Particle *) itr->next()))
    if (p_curr->is_real())
    {
      // calc speed of sound through the medium
      double c = matprops_ptr->sound_speed(p_curr->get_density());

      temp = p_curr->get_smlen() / c;
      if (temp < dt)
        dt = temp;
    }

  // delete HT Iterator
  delete itr;

  return (0.2 * dt);
}
