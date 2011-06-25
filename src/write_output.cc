/*
 * =====================================================================================
 *
 *       Filename:  write_output.cc
 *
 *    Description:  
 *
 *        Created:  03/24/2010 04:20:03 PM EDT
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
# include <config.h>
#endif

#ifdef HAVE_MPI_H
# include <mpi.h>
#endif

#include <cstdio>
#include <cassert>
#include <vector>
using namespace std;

#include <hashtab.h>

#include <properties.h>
#include <particle.h>
#include <outforms.h>
#include "constants.h"
#include "sph_header.h"

void write_output(int myid, int numprocs, 
                  HashTable *P_table, HashTable *BG_mesh, 
                  TimeProps *timeprops, int format)
{

  if ( format&1 )
    write_h5part (myid, numprocs, P_table, timeprops);

  if ( format&2 )
    write_matlab (myid, P_table, BG_mesh, timeprops);

  return;
}



void write_debug_info( int myid, HashTable *P_table, int index)
{

  char fname[20];
  sprintf(fname, "Step%02d%06d.dat", myid, index);
  HTIterator *itr = new HTIterator(P_table);
  Particle *p_curr = NULL;
  while ( (p_curr = (Particle *) itr->next()) )
    if ( (p_curr->getKey().key[0]==1748747027) && 
         (p_curr->is_real()) )
    {
      FILE *fp = fopen(fname, "w");
      vector<Key> neighs = p_curr->get_neighs();
      vector<Key>::iterator itr_p;
      fprintf(fp,"%e, %e, %e, %e, %e\n", 
                 *(p_curr->get_coords()),
                 *(p_curr->get_coords()+1),
                 *(p_curr->get_state_vars()),
                 *(p_curr->get_state_vars()+1),
                 *(p_curr->get_state_vars()+2));
      for (itr_p=neighs.begin(); itr_p!=neighs.end(); itr_p++)
      {
        Particle *p_neigh = (Particle *) P_table->lookup(*itr_p);
        assert(p_neigh);
        if ( *p_curr == *p_neigh ) continue;
        fprintf(fp,"%e, %e, %e, %e, %e\n", 
                 *(p_neigh->get_coords()),
                 *(p_neigh->get_coords()+1),
                 *(p_neigh->get_state_vars()),
                 *(p_neigh->get_state_vars()+1),
                 *(p_neigh->get_state_vars()+2));
      }
      fclose(fp);
    }
}
