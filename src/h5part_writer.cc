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
# include <config.h>
#endif

#ifdef HAVE_HDF5_H
# include <hdf5.h>
# include "hdf5calls.h"
#endif

#ifdef HAVE_MPI_H
# include <mpi.h>
#endif

#include <cstdio>
#include <vector>
using namespace std;

#include <hashtab.h>
#include "constants.h"
#include "particle.h"
#include "sph_header.h"

#define WRITE_GHOSTS

void write_h5part (int myid, int numproc, HashTable *P_table, TimeProps *timepros)
{
  int i, j;
  vector<double> x, y, z, Vx , Vy, Vz, rho;
  char filename[18];
  static int step =0;
  hid_t fp;
#ifdef PARALLEL_IO
  sprintf(filename,"pvplot_out.h5part");
#else
  sprintf(filename,"pvplot%03d.h5part",myid);
#endif

  if ( timepros->ifstart() )
#ifdef PARALLEL_IO
    fp = GH5_fopen_parallel (filename, 'w');
#else
    fp = GH5_fopen( filename, 'w');
#endif
  else
#ifdef PARALLEL_IO
    fp = GH5_fopen_parallel (filename, 'a'); 
#else
    fp = GH5_fopen (filename,'a');
#endif

  HTIterator *itr = new HTIterator(P_table);
  Particle *pi= NULL;
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
      y.push_back(*(pi->get_coords()+1));
      Vy.push_back(*(pi->get_vel()+1));
      my_count++;
#ifdef THREE_D
      z.push_back(*(pi->get_coords()+2));
      Vz.push_back(*(pi->get_vel()+2));
#else
      z.push_back(0);
#endif
#ifdef PARALLEL_IO
    }
#endif
  }
  char group[10];
  sprintf(group,"Step#%d",step);
  step++;

  hid_t gid = GH5_gopen(fp, group, 'w'); 

  int size;
#ifdef PARALLEL_IO
  // make space to recv num of particles on all proc
  int *size_arr = new int [numproc];
  for (i=0; i<numproc; i++)
    size_arr[i] = 0;

  // get sizes from all around the world 
  MPI_Allgather(&my_count, 1, MPI_INT, size_arr, 1, MPI_INT, MPI_COMM_WORLD);

  // create array index partitons
  int *id_lims = new int [numproc+1];
  id_lims[0] = 0;
  for (i=1; i<numproc+1; i++)
    id_lims[i] = id_lims[i-1] + size_arr[i-1];

  // get problem size by adding up     
  for (i=0; i<numproc; i++)
    size += size_arr[i];
#else
  size = my_count;
#endif
  int    dims[2] = {size, 0};
  double *buf = new double[my_count];
  int    *ibuf = new int[my_count];

#ifdef PARALLEL_IO
  j=0;
  int    start = id_lims[myid];
  for (i=id_lims[myid]; i<id_lims[myid+1]; i++)
  {
    ibuf[j] = i;
    j++;
  }
  GH5_par_writedata(gid, "ID",dims, ibuf, start, my_count);
 
  copy(x.begin(), x.end(), buf);
  GH5_par_writedata(gid, "x", dims, buf, start, my_count);

  copy(y.begin(), y.end(), buf);
  GH5_par_writedata(gid, "y", dims, buf, start, my_count);

  copy(rho.begin(), rho.end(), buf);
  GH5_par_writedata(gid, "Rho", dims, buf, start, my_count);

  copy(Vx.begin(), Vx.end(), buf);
  GH5_par_writedata(gid, "Vx", dims, buf, start, my_count);

  copy(Vy.begin(), Vy.end(), buf);
  GH5_par_writedata(gid, "Vy", dims, buf, start, my_count);

# ifdef THREE_D
  copy(z.begin(), z.end(), buf);
  GH5_par_writedata(gid, "z", dims, buf, start, my_count);

  copy(Vz.begin(), Vz.end(), buf);
  GH5_par_writedata(gid, "Vz", dims, buf, start, my_count);
# else
  copy(z.begin(), z.end(), buf);
  GH5_par_writedata(gid, "z", dims, buf, start, my_count);
# endif // THREE_D

  // close file and group
  GH5_gclose(gid);
  GH5_fclose(fp);

// if particler is not built with parallel-io support
#else
  for (i=0; i<size; i++)
    ibuf[i] = i;
  GH5_writedata(gid, "ID",dims, ibuf);
 
  copy(x.begin(), x.end(), buf);
  GH5_writedata(gid, "x", dims, buf);

  copy(y.begin(), y.end(), buf);
  GH5_writedata(gid, "y", dims, buf);

  copy(rho.begin(), rho.end(), buf);
  GH5_writedata(gid, "Rho", dims, buf);

  copy(Vx.begin(), Vx.end(), buf);
  GH5_writedata(gid, "Vx", dims, buf);

  copy(Vy.begin(), Vy.end(), buf);
  GH5_writedata(gid, "Vy", dims, buf);

# ifdef THREE_D
  copy(z.begin(), z.end(), buf);
  GH5_writedata(gid, "z", dims, buf);

  copy(Vz.begin(), Vz.end(), buf);
  GH5_writedata(gid, "Vz", dims, buf);
# else
  copy(z.begin(), z.end(), buf);
  GH5_writedata(gid, "z", dims, buf);
# endif 
  GH5_gclose(gid);
  GH5_fclose(fp);
#endif // PARALLEL_IO

  // free memory
  delete []buf;
  delete []ibuf;
  delete itr;

  return;
}
