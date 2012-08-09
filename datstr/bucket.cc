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

int ineigh[3][3][3] = {{{13, 12, 14}, {10,  9, 11}, {16, 15, 17}},
                       {{4,   3,  5}, {1,   0,  2}, {7,   6,  8}},
                       {{22, 21, 23}, {19, 18, 20}, {25, 24, 26}}
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

//(i,j+1)_________(i+1,j+1)
//    |           |
//    |           |
//    |           |
//    |           |
//    |           |
//(i,j)-----------(i+1,j)
//   
  for (i = 0; i < 4; i++)
    poly[i] = 0.;

  if (bucket_type == MIXED)
  {

    // transform to local coordinate sys
    double xcrd[2], ycrd[2];
    xcrd[0] = 0;    //x(i)
    xcrd[1] = maxcrd[0] - mincrd[0];    //x(i+1)
    ycrd[0] = 0;    //y(j)
    ycrd[1] = maxcrd[1] - mincrd[1];    //y(j+1)

    // and transform elevs too
    for (i = 0; i < 4; i++)
      elev[i] -= mincrd[2];     

    // fit the surface // 4 pts - 4 constants
    poly_surf (xcrd, ycrd, elev, poly);

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
Bucket::get_bnd_normal (double point[], double normal[]) const
{
  int i;
  // if this isn't a boundary bucket, it shouldn't
  if ( bucket_type != MIXED )
    return 1;

  double pnt[DIMENSION-1];
  for (i = 0; i < DIMENSION-1; i++)
    pnt[i] = point[i] - mincrd[i];
  
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
Bucket::find_neigh_dir (Key keyin, int dir[]) const
{
  int i, j, k;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
      {
        if (neighbors[ineigh[i][j][k]] == keyin)
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
  if (!calc_intersection (pnt , intsct))
  {
    double d = 0;
    for (int i=0; i<DIMENSION; i++)
      d += (pnt[i] - intsct[i]) * (pnt[i] - intsct[i]);
    return sqrt (d);
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
  int Down[3] = { 0, 0, 1 };
  unsigned keylen = KEYLENGTH;
  double mindom[DIMENSION], maxdom[DIMENSION];
  double coord[DIMENSION], normc[DIMENSION];
  double seed[DIMENSION];
  vector<Key> new_ghosts;

  // particle properties
  double ma = matprops->particle_mass;
  double hl = matprops->smoothing_length;
  double del = hl;
  double del2 = 0.5 * del;

#ifdef DEBUG
  if (this->has_ghost_particles ())
  {
    fprintf (stderr,
             "Trying to put ghosts in bucket that already have them. \n");
    fprintf (stderr, "%f, %f, %f, %f\n", mincrd[0], maxcrd[0], mincrd[1],
             maxcrd[1]);
    exit (1);
  }
#endif

  // make sure this is the top-most boundary (MIXED) bucket
  int Up[3] = {0, 0, 2};
  Bucket * upper = (Bucket *) BG_mesh->lookup (which_neigh (Up));
  if ( upper->get_bucket_type () != OVERGROUND )
  {
    fprintf (stderr, "ghost generation is not happening in correct bucket\n");
    fflush (stderr);
    exit (1);
  }

  // get min-max domain for key generation
  for (i = 0; i < DIMENSION; i++)
  {
    maxdom[i] = *(P_table->get_maxDom () + i);
    mindom[i] = *(P_table->get_minDom () + i);
  }

  // generate all the ghosts now, we'll figure out 
  // where to put them later
  int npts = PARTICLE_DENSITY;
  for (i = 0; i < npts; i++)
    for (j = 0; j < npts; j++)
    {
      // get seeds from boundary points
      seed[0] = bnd_xcrd[i];
      seed[1] = bnd_ycrd[j];
      seed[2] = get_bndZ (seed);

      // generate particles
      coord[0] = seed[0];
      coord[1] = seed[1];

      for (int ighost = 0; ighost < NUM_GHOST_ROWS; ighost++)
      {
        coord[2] = seed[2] - del2 - ighost * del;

        // create a new particle and add to hash-table
        for (k = 0; k < DIMENSION; k++)
          normc[k] = (coord[k] - mindom[k]) / (maxdom[k] - mindom[k]);

        // generate hashtable keys
        HSFC3d (normc, &keylen, key);
        Key keystr (key);

        // create new ghost partilce
        Particle *pt = new Particle (key, coord, ma, hl, 1);

#ifdef DEBUG
        // check if particle is a ducplicate
        if (P_table->lookup (key))
          fprintf (stderr,"Particle already exists.\n");
#endif

        // add new particle to hash-table
        P_table->add (key, pt);

        // copy key to vector
        new_ghosts.push_back (keystr);
      }
    }

  // since we started at the top-most boundary
  // bucket, we'll go down one-by-one
  Bucket * buck = this;
  vector<Key> temp;
  vector<Key>::iterator p_itr;
  while ( buck )
  {
    for (p_itr = new_ghosts.begin (); p_itr != new_ghosts.end(); p_itr++)
    {
      Particle *pghost = (Particle *) P_table->lookup (*p_itr);
      for (i = 0; i < DIMENSION; i++)
        coord[i] = *(pghost->get_coords() + i);

      if ( buck->contains (coord) )
        temp.push_back (*p_itr);
    }
    // add ghost particle to the list of particles
    if ( temp.size () > 0 )
    {
      buck->add_ghosts (temp);
      temp.clear ();
    }
    buck->set_ghost_particles (true);
    buck->mark_active ();
    buck = (Bucket *) BG_mesh->lookup (buck->which_neigh (Down));
  }

  return 0;
}


/*! iterates for point \f$\mathbf{x}_b\f$ on boundary from point 
 *  \f$\mathbf{x}\f$, outside the boundary, such that 
 *  \f$ \mathbf{x}_b - \mathbf{x} \f$ is normal to the boundary,
 *  given it is in the form : \f$ a x + b y + c x y + d - z = 0 \f$
 *  Newton's method, upto 5 cycles is used.
 */
int 
Bucket::calc_intersection (double * point, double * xnew) const
{
  register int i, j;
  register double xold[DIMENSION];
  register double pl[4];
  double pt[DIMENSION], tmp[DIMENSION];
  double tol = 1.0E-5;
  double err = 0.;

  // transform point to local coordinates
  for (i = 0; i < DIMENSION; i++)
    pt[i] = point[i] - mincrd[i];

  // hopefully it will stay in registers
  for (i = 0; i < 4; i++)
    pl[i] = poly[i];

  // inital guess
  xnew[0] = 0.5 * (maxcrd[0] - mincrd[0]);
  xnew[1] = 0.5 * (maxcrd[1] - mincrd[1]);
  xnew[2] = pl[0] * xnew[0] + pl[1] * xnew[1] + 
            pl[2] * xnew[0] * xnew[1] + pl[3];

  // Newton's method
  for (i=0; i<8; i++)
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

  // transform back to global coordinates
  for (i = 0; i < DIMENSION; i++)
    xnew[i] += mincrd[i];

  return 0;
}
