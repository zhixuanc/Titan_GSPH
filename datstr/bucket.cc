
/*
 * =====================================================================================
 *
 *       Filename:  bucket.cc
 *
 *    Description:  
 *
 *        Created:  04/11/2010 06:41:40 PM EDT
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

//#define CYLINDERICAL_PILE
#define PARABOLIC_PILE

#include <cstdio>
#include <cassert>
using namespace std;

#include <particle.h>
#include <hilbert.h>
#include <hashtab.h>

#include "bucket.h"

int ineigh[3][3][3] = { {{13, 12, 14}, {10, 9, 11}, {16, 15, 17}},
{{4, 3, 5}, {1, 0, 2}, {7, 6, 8}},
{{22, 21, 23}, {19, 18, 17}, {25, 24, 26}}
};

void surffit (double *, double *, double *);

// constructors
Bucket::Bucket (unsigned *k, double *minx, double *maxx, int buck_type,
                double *elev, int myid, int *nproc, Key * neigh)
{
  int i;
  double coordMatix[9];

  bucket_type = buck_type;
  myprocess = myid;
  guest_flag = 0;
  //active = false;
  active = false;
  particles_type = 0;

  for (i = 0; i < KEYLENGTH; i++)
    key.key[i] = k[i];

  for (i = 0; i < DIMENSION; i++)
  {
    mincrd[i] = minx[i];
    maxcrd[i] = maxx[i];
  }
  for (i = 0; i < NEIGH_SIZE; i++)
  {
    neighbors[i] = neigh[i];
    neigh_proc[i] = nproc[i];
  }

  // boundary points
  double del = (maxcrd[0] - mincrd[0]) / PARTICLE_DENSITY;
  double del2 = del * 0.5;

  for (i = 0; i < PARTICLE_DENSITY; i++)
  {
    if (bucket_type == MIXED)
    {
      bnd_xcrd[i] = mincrd[0] + del2 + i * del;
      bnd_ycrd[i] = mincrd[1] + del2 + i * del;
    }
    else
    {
      bnd_xcrd[i] = 0.;
      bnd_ycrd[i] = 0.;
    }
  }
  // list of particles in the bucket
  particles.clear ();
  new_plist.clear ();

  // fit lower triangle
  coordMatix[0] = mincrd[0];    //x(i)
  coordMatix[1] = mincrd[1];    //y(i)
  coordMatix[2] = 1;
  coordMatix[3] = maxcrd[0];    //x(i+1)
  coordMatix[4] = mincrd[1];    //y(i)
  coordMatix[5] = 1;
  coordMatix[6] = maxcrd[0];    //x(i+1)
  coordMatix[7] = maxcrd[1];    //x(i+1)
  coordMatix[8] = 1;
  surffit (coordMatix, elev, lower_tri);

  // fit upper triangle
  coordMatix[0] = maxcrd[0];    //x(i+1)
  coordMatix[1] = maxcrd[1];    //y(i+1)
  coordMatix[2] = 1;
  coordMatix[3] = mincrd[0];    //x(i)
  coordMatix[4] = maxcrd[1];    //y(i+1)
  coordMatix[5] = 1;
  coordMatix[6] = mincrd[0];    //x(i)
  coordMatix[7] = mincrd[0];    //x(i)
  coordMatix[8] = 1;
  surffit (coordMatix, elev + 1, lower_tri);

} // end constuctor

Bucket::Bucket ()
{
  int i, j;

  myprocess = -1;
  active = false;
  guest_flag = 0;
  bucket_type = 0;
  particles_type = 0;

  for (i = 0; i < DIMENSION; i++)
  {
    mincrd[i] = 0.;
    maxcrd[i] = 0.;
  }

  for (i = 0; i < 4; i++)
  {
    lower_tri[i] = 0.;
    upper_tri[i] = 0.;
  }

  for (i = 0; i < NEIGH_SIZE; i++)
  {
    neigh_proc[i] = -1;
    for (j = 0; j < KEYLENGTH; j++)
      neighbors[i].key[j] = 0;
  }
  particles.clear ();
  new_plist.clear ();
  return;
}

//! boundary normal from point x, if its a boundary bucket
int
Bucket::get_bnd_normal (double pnt[], double normal[]) const
{
  if (bucket_type == MIXED)
  {
    if (abs (get_lower_tri_dist (pnt)) < abs (get_upper_tri_dist (pnt)))
      for (int i = 0; i < 3; i++)
        normal[i] = lower_tri[i];
    else
      for (int i = 0; i < 3; i++)
        normal[i] = upper_tri[i];
    return 0;
  }
  else
  {
    for (int i = 0; i < DIMENSION; i++)
      normal[i] = 0.;
    return -1;
  }
}

Key
Bucket::which_neigh (int dir[]) const
{
  int i, j, k;

  i = dir[0];
  j = dir[1];
  k = dir[2];
  return neighbors[ineigh[i][j][k]];
}

int
Bucket::which_neigh_proc (int dir[]) const
{
  int i, j, k;

  i = dir[0];
  j = dir[1];
  k = dir[2];
  return neigh_proc[ineigh[i][j][k]];
}

void
Bucket::put_neigh_proc (int dir[], int proc)
{
  int i = dir[0], j = dir[1], k = dir[2];

  neigh_proc[ineigh[i][j][k]] = proc;
  return;
}

bool
Bucket::find_neigh_dir (Key neigh, int dir[]) const
{
  int i, j, k;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
      {
        Key comp_key = neighbors[ineigh[i][j][k]];

        if (compare_keys (neigh, comp_key))
        {
          dir[0] = i;
          dir[1] = j;
          dir[2] = k;
          return true;
        }
      }
  return false;
}

/* 
 * Normal Distance between point and boundary:
 * finite value is returned only if the 
 * intersection between normal from point
 * and boundary is within current bucket, 
 * else infinity is returned
 */
double
Bucket::get_bnddist (double pnt[], double intsct[]) const
{
  int i;
  double d1 = get_lower_tri_dist (pnt);
  double d2 = get_upper_tri_dist (pnt);

  if (abs (d1) < abs (d2))
  {
    for (i = 0; i < DIMENSION; i++)
      intsct[i] = pnt[i] - d1 * lower_tri[i];
    return d1;
  }
  else
  {
    for (i = 0; i < DIMENSION; i++)
      intsct[i] = pnt[i] - d2 * upper_tri[i];
    return d2;
  }
}

// put ghost particles in active buckets
int
Bucket::put_ghost_particles (HashTable * P_table,
                             HashTable * BG_mesh, MatProps * matprops)
{
  int i, j, k;
  unsigned key[KEYLENGTH];
  int Down[DIMENSION] = { 0, 0, 1 };
  unsigned keylen = KEYLENGTH;
  double mindom[DIMENSION], maxdom[DIMENSION];
  double crd[DIMENSION], normc[DIMENSION];

  vector < Key > neighparts;
  Key keystr;

  if (have_ghost_particles ())
  {
    fprintf (stderr,
             "Trying to put ghosts in bucket that already have them. \n");
    fprintf (stderr, "%f, %f, %f, %f\n", mincrd[0], maxcrd[0], mincrd[1],
             maxcrd[1]);
    exit (1);
  }

  for (i = 0; i < DIMENSION; i++)
  {
    maxdom[i] = *(P_table->get_maxDom () + i);
  }

  // particle properties
  double ma = matprops->particle_mass;
  double hl = matprops->smoothing_length;
  double dx = hl;
  double dx2 = 0.5 * dx;
  double seed[DIMENSION];

  // get handle to the neigh below
  Key nghkey = which_neigh (Down);
  Bucket *neigh = (Bucket *) BG_mesh->lookup (nghkey);

  assert (neigh);
  int npts = PARTICLE_DENSITY;

  for (i = 0; i < npts; i++)
    for (j = 0; j < npts; j++)
    {
      // get seeds from boundary points
      seed[0] = bnd_xcrd[i];
      seed[1] = bnd_ycrd[j];
      seed[2] = get_bndZ (seed);

      int iend = DIMENSION - 1;

      if (this->contains (seed))
      {
        // put particles
        for (k = 0; k < DIMENSION - 1; k++)
          crd[k] = seed[k];

        for (int ighost = 0; ighost < NUM_GHOST_ROWS; ighost++)
        {
          crd[DIMENSION - 1] = seed[DIMENSION - 1] - dx2 - ighost * dx;

          // create a new particle and add to hash-table
          for (k = 0; k < DIMENSION; k++)
            normc[k] = (crd[k] - mindom[k]) / (maxdom[k] - mindom[k]);

          // generate hashtable keys
          HSFC2d (normc, &keylen, key);
          for (k = 0; k < KEYLENGTH; k++)
            keystr.key[k] = key[k];

          // create new ghost partilce
          Particle *pt = new Particle (key, crd, ma, hl, 1);

          P_table->add (key, pt);

          // now find out which bucket it belongs to
          if (this->contains (crd))
            particles.push_back (keystr);
          else
            neighparts.push_back (keystr);
        }
      }
    }
  if (neighparts.size () > 0)
  {
    neigh->add_ghosts (neighparts);
    if (neigh->get_bucket_type () == UNDERGROUND)
      neigh->set_ghost_particles (true);
    neigh->mark_active ();
  }

  set_ghost_particles (true);
  return 0;
}

/*! Solves Ax = b for 3x3 matrix */
void
surffit (double *f,             //! Matrix of {x(i), y(i), 1}
         double *z,             //! Vector of elevations
         double *poly)          //! polynomial constants
{
  double t1 = (f[4] * f[8]) - (f[5] * f[7]);
  double t2 = (-f[1] * f[8]) + (f[2] * f[7]);
  double t3 = (f[1] * f[5]) - (f[2] * f[4]);
  double t4 = (t2 * f[3]) + (t3 * f[6]) + (f[0] * t1);

  t4 = 0.1e1 / t4;
  poly[0] = (t1 * z[0] + t2 * z[1] + t3 * z[2]) * t4;
  poly[1] = (-(-f[6] * f[5] + f[3] * f[8]) * z[0] +
             (-f[6] * f[2] + f[0] * f[8]) * z[1] -
             (-f[3] * f[2] + f[0] * f[5]) * z[2]) * t4;
  poly[2] = ((-f[6] * f[4] + f[3] * f[7]) * z[0] -
             (-f[6] * f[1] + f[0] * f[7]) * z[1] +
             (-f[3] * f[1] + f[0] * f[4]) * z[2]) * t4;
  return;
}
