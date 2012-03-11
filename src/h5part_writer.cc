
/*
 * =====================================================================================
 *
 *       Filename:  h5part_writer.cc
 *
 *    Description:  H5PART is Paraview and Visit readable
 *
 *        Created:  08/18/2010 12:38:59 PM EDT
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

#include <hdf5.h>

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <cstdio>
#include <vector>
using namespace std;

#include <hashtab.h>
#include "constants.h"
#include "particle.h"
#include "sph_header.h"
#include "hdf5calls.h"

#define WRITE_GHOSTS

void
write_h5part(int myid, int numproc, HashTable * P_table, TimeProps * timepros)
{
  int i, j;
  vector < double >x, y, z, Vx, Vy, Vz, rho;
  char filename[18];
  static int step = 0;
  hid_t fp;
  herr_t ierr;

#ifdef PARALLEL_IO
  sprintf(filename, "pvplot_out.h5part");
#else
  sprintf(filename, "pvplot%03d.h5part", myid);
#endif

  if (timepros->ifstart())
    fp = GH5_fopen(filename, 'w');
  else
    fp = GH5_fopen(filename, 'a');

  HTIterator *itr = new HTIterator(P_table);
  Particle *pi = NULL;
  int my_count = 0;

  while ((pi = (Particle *) itr->next()))
  {
#ifdef PARALLEL_IO
    if (!pi->is_guest())
    {
#endif
      rho.push_back(pi->get_density());
      x.push_back(*(pi->get_coords()));
      Vx.push_back(*(pi->get_vel()));
      y.push_back(*(pi->get_coords() + 1));
      Vy.push_back(*(pi->get_vel() + 1));
      my_count++;
      z.push_back(*(pi->get_coords() + 2));
      Vz.push_back(*(pi->get_vel() + 2));
#ifdef PARALLEL_IO
    }
#endif
  }
  char group[10];

  sprintf(group, "Step#%d", step);
  step++;

  hid_t gid = GH5_gopen(fp, group, 'w');

  int start;
  int size = 0;
  int *id_lims = new int[numproc + 1];

  for (i = 0; i < numproc + 1; i++)
    id_lims[i] = 0;

#ifdef PARALLEL_IO
  // make space to recv num of particles on all proc
  int *size_arr = new int[numproc];

  for (i = 0; i < numproc; i++)
    size_arr[i] = 0;

  // get sizes from all around the world 
  MPI_Allgather(&my_count, 1, MPI_INT, size_arr, 1, MPI_INT, MPI_COMM_WORLD);

  // create array index partitons
  id_lims[0] = 0;
  for (i = 1; i < numproc + 1; i++)
    id_lims[i] = id_lims[i - 1] + size_arr[i - 1];
  start = id_lims[myid];

  // get problem size by adding up     
  for (i = 0; i < numproc; i++)
    size += size_arr[i];

  delete[]size_arr;
#else
  size = my_count;
  id_lims[myid] = 0;
  id_lims[myid + 1] = size;
#endif
  int dims[2] = { size, 0 };
  double *buf = new double[my_count];
  int *ibuf = new int[my_count];

  j = 0;
  for (i = id_lims[myid]; i < id_lims[myid + 1]; i++)
  {
    ibuf[j] = i;
    j++;
  }

  // particle ids
  ierr = GH5_Write(gid, "ID", dims, (void *) ibuf, start, my_count, INTTYPE);

  // x-coordinates
  copy(x.begin(), x.end(), buf);
  ierr = GH5_Write(gid, "x", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // y-coordinates
  copy(y.begin(), y.end(), buf);
  ierr = GH5_Write(gid, "y", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // z-coordinates
  copy(z.begin(), z.end(), buf);
  ierr = GH5_Write(gid, "z", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // density of particles
  copy(rho.begin(), rho.end(), buf);
  ierr = GH5_Write(gid, "Rho", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in x-directions
  copy(Vx.begin(), Vx.end(), buf);
  ierr = GH5_Write(gid, "Vx", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in y-directions
  copy(Vy.begin(), Vy.end(), buf);
  ierr = GH5_Write(gid, "Vy", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // Velocity in z-directions
  copy(Vz.begin(), Vz.end(), buf);
  ierr = GH5_Write(gid, "Vz", dims, (void *) buf, start, my_count, DOUBLETYPE);

  // close file and group
  ierr = GH5_gclose(gid);

  // flush output 
  ierr = H5Fflush(fp, H5F_SCOPE_LOCAL);
  ierr = H5Fclose(fp);

  // free memory
  delete[]id_lims;
  delete[]buf;
  delete[]ibuf;
  delete itr;

  return;
}
