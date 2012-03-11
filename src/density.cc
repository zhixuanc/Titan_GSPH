
/*
 * =====================================================================================
 *
 *       Filename:  density.cc
 *
 *    Description:  
 *
 *        Created:  11/14/2008 02:07:27 PM EST
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
#include <cassert>
#include <iostream>
using namespace std;

#include <hashtab.h>

#include "particle.h"
#include "constants.h"
#include "sph_header.h"

void
smooth_density(HashTable * P_table)
{
  int i, j, k, no_of_neighs;
  unsigned jkey[KEYLENGTH];
  double xi[DIMENSION], ds[DIMENSION], s[DIMENSION];
  double tmprho, wnorm, wght, density;
  Key tmpkey;

  vector < Key > neighs;

  // create a Hash Table iterator instance
  HTIterator *itr = new HTIterator(P_table);
  Particle *pi;

  // iterate over the table
  while (pi = (Particle *) itr->next())
  {
    if (pi->is_real())
    {
      for (i = 0; i < DIMENSION; i++)
        xi[i] = (*(pi->get_coords() + i));
      double hi = pi->get_smlen();
      double supp = 3.0 * hi;

      tmprho = 0;
      wnorm = 0;
      neighs = pi->get_neighs();
      no_of_neighs = neighs.size();
      for (j = 0; j < no_of_neighs; j++)
      {
        // get the neighbors
        Particle *pj = (Particle *) P_table->lookup(neighs[j]);

        assert(pj);
        // guests are included ghosts are not
        if (!pj->is_ghost())
        {
          for (i = 0; i < DIMENSION; i++)
            ds[i] = xi[i] - *(pj->get_coords() + i);
          if (in_support(ds, supp))
          {
            for (k = 0; k < DIMENSION; k++)
              s[k] = ds[k] / hi;
            wght = weight(s, hi);
            tmprho += wght * (pj->get_mass());
            wnorm += wght * (pj->get_mass()) / (pj->get_density());
          }
        }
      }
      density = tmprho / wnorm;
      assert(density > 0);
      pi->put_density(density);
    }
  }

  // clean up stuff
  delete itr;

  return;
}
