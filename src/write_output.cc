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
using namespace std;

#include <hashtab.h>

#include <properties.h>
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
