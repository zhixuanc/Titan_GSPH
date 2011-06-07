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
#include <config.h>
#endif

#include <vector>
#include <cmath>
#include <cassert>
using namespace std;

#ifdef HAVE_MPI_H
#include <mpi.h>
#endif

#include <hashtab.h>
#include <particle.h>
#include <properties.h>

#include "constants.h"

double timestep (HashTable *p_table, MatProps *matprops_ptr)
{
  int i,j,k;
  int no_of_neighs;
  double xi[DIMENSION];
  double dt, temp;
  vector<Key> neighs;
  
  dt = 1.0E+10; // Initialize to very high value 
  HTIterator *itr = new HTIterator(p_table);
  Particle *pi;
  while ( (pi = (Particle *) itr->next()) )
  {
    if ( pi->is_real() )
    {
      double c = matprops_ptr->sound_speed(pi->get_density());
      for (i=0; i<DIMENSION; i++)
        xi[i]=*(pi->get_coords()+i);
      neighs=pi->get_neighs();
      no_of_neighs=neighs.size();
      for (j=0; j<no_of_neighs; j++)
      {
        Particle *pj = (Particle *) p_table->lookup(neighs[j]);
        assert(pj);
        if ( *pi == *pj )
          continue;
         
	double d=0;
        for (i=0; i<DIMENSION; i++)
	  d += pow((xi[i]-*(pj->get_coords()+i)),2);

	d = sqrt(d);
	temp = pi->get_smlen()/c;
	if ( temp < dt )
	   dt = temp;
      }
    }
  }

  // delete HT Iterator
  delete itr;
  return (0.2*dt);
}
