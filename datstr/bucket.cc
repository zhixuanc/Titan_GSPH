
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
#include <sph_header.h>
#include "bucket.h"

int ineigh[3][3][3] = { {{13, 12, 14}, {10, 9, 11}, {16, 15, 17}},
{{4, 3, 5}, {1, 0, 2}, {7, 6, 8}},
{{22, 21, 23}, {19, 18, 17}, {25, 24, 26}}
};

// constructors
Bucket::Bucket (unsigned *keyi, double *minx, double *maxx, int buck_type,
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
    key.key[i] = keyi[i];

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

//(i,j+1)_________ (i+1,j+1)
//    |          /|
//    |        /  |
//    |      /    |
//    |    /      |
//    |  /        |
//    |/          |
//(i,j)-----------(i+1,j)
//   
  for (i = 0; i < 4; i++)
    poly[i] = 0.;

  if (bucket_type == MIXED)
  {
    double xcrd[2], ycrd[2];
    xcrd[0] = mincrd[0];    //x(i)
    xcrd[1] = maxcrd[0];    //x(i+1)
    ycrd[0] = mincrd[1];    //y(j)
    ycrd[1] = maxcrd[1];    //y(j+1)
    lsq_surf4 (xcrd, ycrd, elev, poly);
    bool printit = false;
    for (i=0; i<4; i++)
      if (isnan (poly[i]))
        exit (51);
  }
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
    poly[i] = 0.;


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

/*! boundary normal from point \f$\mathbf{x}\f$, provided
 * /f$ \mathbf{x} \in \Gamma^i\f$.
 */
int
Bucket::get_bnd_normal (double pnt[], double normal[]) const
{
  int i;

  // if this isn't a boundary bucket, it shouldn't
  if ( bucket_type != MIXED )
    return 1;
  normal[0] = -(poly[0] + poly[2] * pnt[1]);    // P1 + P3*x
  normal[1] = -(poly[1] + poly[2] * pnt[0]);    // P2 + P3*y
  normal[2] = 1.;

  // normalize
  double d = 0.;
  for (i = 0; i < DIMENSION; i++)
    d += normal[i] * normal[i];

  for (i = 0; i < DIMENSION; i++)
    normal[i] /= sqrt(d);
  return 0;
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
  if (bucket_type == MIXED)
  {
    calc_intersection (pnt , intsct);
    if ( this->contains (intsct) )
    {
      double d = 0;
      for (int i=0; i<DIMENSION; i++)
        d += (pnt[i] - intsct[i]) * (pnt[i] - intsct[i]);
      return sqrt (d);
    }
    else
      return 1.E+10;
  }
  else
    return 1.E+10;
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
    mindom[i] = *(P_table->get_minDom () + i);
  }

  // particle properties
  double ma = matprops->particle_mass;
  double hl = matprops->smoothing_length;
  double del = hl;
  double del2 = 0.5 * del;
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

      if (this->contains (seed))
      {
        // put particles
        crd[0] = seed[0];
        crd[1] = seed[1];

        for (int ighost = 0; ighost < NUM_GHOST_ROWS; ighost++)
        {
          crd[2] = seed[2] - del2 - ighost * del;

          // create a new particle and add to hash-table
          for (k = 0; k < DIMENSION; k++)
            normc[k] = (crd[k] - mindom[k]) / (maxdom[k] - mindom[k]);

          // generate hashtable keys
          HSFC3d (normc, &keylen, key);
          for (k = 0; k < KEYLENGTH; k++)
            keystr.key[k] = key[k];

          // create new ghost partilce
          Particle *pt = new Particle (key, crd, ma, hl, 1);
          // add new particle to hash-table
          P_table->add (key, pt);

          // now find out, which bucket it belongs to
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


/*! iterates for point \f$\mathbf{x}_b\f$ on boundary from point 
 *  \f$\mathbf{x}\f$, outside the boundary, such that 
 *  \f$ \mathbf{x}_b - \mathbf{x} \f$ is normal to the boundary,
 *  given it is in the form : \f$ a x + b y + c x y + d - z = 0 \f$
 *  Newton's method, upto 5 cycles is used.
 */
int 
Bucket::calc_intersection (double * pt, double * xnew) const
{
  register int i, j;
  register double xold[DIMENSION];
  register double pl[4];
  double tmp[DIMENSION];
  double tol = 1.0E-5;
  double err = 0.;

  // hopefully it will staty in registers
  for (i = 0; i < 4; i++)
    pl[i] = poly[i];

  // inital guess
  xnew[0] = pt[0];
  xnew[1] = pt[1];
  xnew[2] = this->get_bndZ (xnew);

  // Newton's method
  for (i=0; i<5; i++)
  {
    // copy x(n+1) to x(n)
    for (j = 0; j < DIMENSION; j++)
      xold[j] = xnew[j];

    // function at xold
    double fn[3] = { (pl[0] * xold[0]) + (pl[1] * xold[1]) +
                     (pl[2] * xold[0] * xold[1]) + pl[3] - xold[2]     ,
                     (xold[0] - pt[0]) + 
                     ((pl[0] + (pl[2] * xold[0])) * (xold[2] - pt[2])) ,
                     (xold[1] - pt[1]) +
                     ((pl[1] + (pl[2] * xold[1])) * (xold[2] - pt[2]))
                   };
 
    // jacobian at xold
    double a1 = pl[0] + pl[2] * xold[1];
    double a2 = pl[1] + pl[2] * xold[0];
    double a3 = pl[2] * (xold[2] - pt[2]);

    /* jacobian = {{ a1, a2, -1},
                   { 1,  a3, a1},
                   { a3, 1,  a2}
                  }
     */
    // --- begin --- maple generated code
    double t1 = (-a3 *  a2 + a1);
    double t2 = ( a2 *  a2);
    double t3 = ( a1 *  a1);
    double t4 = ( a2 *  a1);
    double t5 = t3 + t2 + (-2 * t4 - a3) * a3 + 1;
    t4 = t4 + a3;
    t5 = 1. / t5;
    double t6 = -a2 + a1 * a3;
    tmp[0] = (t1 * fn[0] + (t2 + 1.) * fn[1] - t4 * fn[2]) * t5;
    tmp[1] = (-t6 * fn[0] - t4 * fn[1] + (t3 + 1.0) * fn[2]) * t5;
    tmp[2] = ((-1 + a3 * a3) * fn[0] + t1 * fn[1] - t6 * fn[2]) * t5;
    // --- end ----- maple generated code

    // upate x(n+1)
    for (j = 0; j < DIMENSION; j++)
      xnew[j] = xold[j] - tmp[j];

    // check error
    err = 0;
    for (j = 0; j < DIMENSION; j++)
      err += (xnew[j] - xold[j])*(xnew[j] - xold[j]);
    if (sqrt (err) < tol ) break;
  }
  if ( sqrt (err) > tol ) return 1;
  return 0;
}
