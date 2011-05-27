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

int ineigh[3][3]  = {{4, 3, 5},{ 1, 0, 2},{7, 6, 8}};
int Left[2] = { 1, 0 };
int Right[2]= { 2, 0 };
int Down[2] = { 0, 1 };
int Up[2]   = { 0, 2 };

// constructors
Bucket::Bucket(unsigned *k, double *minx, double *maxx, 
           int buck_type, int bnd_type, double* pcoef,
           int myid, int ngprpc[], Key neighc[])
{
  int i;

  bucket_type = buck_type;
  myprocess = myid;
  guest_flag = 0;
  active = false;
  boundary_type = bnd_type;
  real_part_flag = false;
  ghost_part_flag = false;

  for (i=0; i<KEYLENGTH; i++)
    key.key[i]=k[i];

  for (i=0; i<DIMENSION; i++)
  {
    mincrd[i] = minx[i];
    maxcrd[i] = maxx[i];
  }
  for (i=0; i<NEIGH_SIZE; i++)
  {
    neighbors[i]=neighc[i];
    neigh_proc[i]=ngprpc[i];
  }

  // boundary points
  double del = (maxcrd[0]-mincrd[0])/PARTICLE_DENSITY;
  double del2 = del*0.5;
  for (i=0; i<PARTICLE_DENSITY; i++)
  {
    if ( bucket_type == MIXED )
      bnd_xcrd[i] = mincrd[0] + del2 + i*del;
    else
      bnd_xcrd[i] = 0.;
#ifdef THREE_D
    if ( bucket_type == MIXED )
      bnd_ycrd[i] = mincrd[1] + del2 + i*del;
    else 
      bnd_ycrd[i] = 0.;
#endif
  }
  // list of particles in the bucket
  particles.clear();
  new_plist.clear();

  for (i=0; i<3; i++)
    poly[i] = pcoef[i];

} // end constuctor

Bucket::Bucket ()
{
  int i, j;
  myprocess = -1;
  active = false;
  guest_flag = 0;
  bucket_type = 0;
  real_part_flag = false;
  ghost_part_flag = false;

  for (i=0; i < DIMENSION; i++ )
  {
    mincrd[i] = 0.;
    maxcrd[i] = 0.;
  }

#ifdef THREE_D
  for (i=0; i < 4; i++)
#else
  for (i=0; i < 3; i++)
#endif 
    poly[i] = 0.;

  for (i=0; i < NEIGH_SIZE; i++)
  {
    neigh_proc[i] = -1;
    for ( j=0; j < KEYLENGTH; j++ )
      neighbors[i].key[j] = 0;
  }
  particles.clear();
  new_plist.clear();

  return;
}

// Normal to the boundary
double Bucket::get_bndnorm(double pnt[]) const
{
  double dx = pnt[0] - poly[0];
  double dy = pnt[1] - poly[1];
  switch (boundary_type)
  {
    case LINE:
      if ( abs(poly[0]) > 0 )
        return (poly[1]/poly[0]);
      else
        return HUGE_VAL;
    case CIRCLE:
      if ( abs(dx) > 0 )
        return (dy/dx);
      else
        return HUGE_VAL;
  }
  return -1.;
}

Key Bucket::which_neigh(int dir[]) const
{
  int i,j,k;
  i=dir[0];
  j=dir[1];
#ifdef THREE_D
  k=dir[2];
  return neighbors[ineigh[i][j][k]];
#else
  return neighbors[ineigh[i][j]];
#endif 
}

int Bucket::which_neigh_proc(int dir[]) const
{
  int i,j,k;
  i=dir[0];
  j=dir[1];
#ifdef THREE_D
  k=dir[2];
  return neigh_proc[ineigh[i][j][k]];
#endif
  return neigh_proc[ineigh[i][j]];
}

void Bucket::put_neigh_proc(int dir[], int proc)
{
  int i = dir[0], j = dir[1];
#ifdef THREE_D
  int k = dir[2];
  neigh_proc[ineigh[i][j][k]] = proc;
#else
  neigh_proc[ineigh[i][j]] = proc;
#endif
  return;
}

bool Bucket::find_neigh_dir ( Key neigh, int dir[] )
{
  int i, j;
  for (i=0; i<=2; i++ )
    for (j=0; j<=2; j++ )
    {
      Key comp_key = neighbors[ineigh[i][j]];
      if ( compare_keys(neigh, comp_key) )
      {
        dir[0] = i; 
        dir[1] = j;
        return true;
      }
    }
  return false;
}

// Normal Distance between point and boundary:
// finite value is returned only if the 
// intersection between normal from point
// and boundary is within current bucket, 
// else infinity is returned
double Bucket::get_bnddist(double pnt[], double intsct[])
{
  int i;
  double t1, m;
  double dst,c2;
  bool   interior;
  switch (boundary_type)
  {
    case LINE:
      // compute distance from the plane/line with poly coeffs
      dst = poly[DIMENSION];
      for (i=0; i<DIMENSION; i++)
        dst += poly[i]*pnt[i];
      // line constant
      t1 = -dst;
      // get intersection of perpendicular from pnt and plane/line
      for (i=0; i<DIMENSION; i++)
        intsct[i] = pnt[i] + poly[i]*t1;

      // if bucket contains this points
      if ( this->contains (intsct) )
        return (dst);
      else
        return HUGE_VAL;
    case CIRCLE:
      // check if point is interior to the circle
      // if pnt is interior distance is negative
      if ( sqrt(pow(pnt[0]-poly[0],2) 
              + pow(pnt[1]-poly[1],2)) <= poly[2] )
        interior = true;
      else
        interior = false;

      // slope of line between point and center of circle
      m = (pnt[1]-poly[1])/(pnt[0]-poly[0]);
      // solve for intersection of line through center
      // and the circle
      t1 = poly[2]/sqrt(1. + m*m);
      intsct[0] = poly[0] - t1;
      intsct[1] = poly[1] + m*(intsct[0] - poly[0]);
      if ( this->contains(intsct) )
      {
        dst = sqrt(pow(intsct[0]-pnt[0],2) + pow(intsct[1]-pnt[1],2));
        if ( interior )
          return (-dst);
        else
          return dst;
      }
      else
      {
        intsct[0] = poly[0] + t1;
        intsct[1] = poly[1] + m*(intsct[0] - poly[0]);
        if ( this->contains( intsct ) )
        {
          dst = sqrt(pow(intsct[0]-pnt[0],2) + pow(intsct[1]-pnt[1],2));
          if ( interior )
            return (-dst);
          else
            return dst;
        }
        else
          return HUGE_VAL;
      }
    default:
      fprintf(stderr,"ERROR!!! Unknown boundary type in \"Bucket.cc\" \n");
      exit(1);
  }
}

// put ghost particles in active buckets
int Bucket::put_ghost_particles(HashTable *P_table, 
                  HashTable *BG_mesh, MatProps *matprops)
{
  int      i,j,k;
  unsigned key[KEYLENGTH];
  unsigned keylen=KEYLENGTH;
  double   mindom[DIMENSION],maxdom[DIMENSION];
  double   crd[DIMENSION], normc[DIMENSION];
  vector<Key> neighparts;
  Key      keystr;

  if ( ghost_part_flag )
  {
    fprintf(stderr,"Trying to put ghosts in bucket that already have them. \n");
    fprintf(stderr,"%f, %f, %f, %f\n", mincrd[0], maxcrd[0], mincrd[1], maxcrd[1]);
    exit(1);
  }

  for (i=0; i<DIMENSION; i++)
  {
    mindom[i]=*(P_table->get_minDom()+i);
    maxdom[i]=*(P_table->get_maxDom()+i);
  }

  // particle properties
  double ma = matprops->particle_mass;
  double hl = matprops->smoothing_length;
  double dx = hl;
  double dx2 = 0.5*dx;
  double seed[DIMENSION];

  // get handle to the neigh below
  Key nghkey = which_neigh(Down);
  Bucket *neigh = (Bucket *) BG_mesh->lookup(nghkey);
  int npts = PARTICLE_DENSITY;

  for (i=0; i < npts; i++)
#ifdef THREE_D
    for (j=0; j < npts; j++)
#endif
  {
    // get seeds from boundary points
    seed[0] = bnd_xcrd[i];
#ifdef THREE_D
    seed[1] = bnd_ycrd[j];
    seed[2] = get_bndZ(seed);
#else
    seed[1] = get_bndZ(seed);
#endif

    int iend = DIMENSION-1;
    if ( (seed[iend] > mincrd[iend]) &&
         (seed[iend] < maxcrd[iend]) )
    {
      // put particles
      for (k=0; k < DIMENSION-1; k++)
        crd[k] = seed[k];

      for (int ighost=0; ighost<NUM_GHOST_ROWS; ighost++)
      {
        crd[DIMENSION-1] = seed[DIMENSION-1] - dx2 - ighost*dx;

        // create a new particle and add to hash-table
        for (k=0; k<DIMENSION; k++)
          normc[k] = (crd[k]-mindom[k])/(maxdom[k]-mindom[k]);

        // generate hashtable keys
        HSFC2d (normc, &keylen, key);
        for (k=0; k<KEYLENGTH; k++)
          keystr.key[k]=key[k];

        // create new ghost partilce
        Particle *pt = new Particle(key, crd, ma, hl, 1);
        P_table->add(key,pt);

        // now find out which bucket it belongs to
        if ( this->contains(crd) )
          particles.push_back(keystr);
        else
          neighparts.push_back(keystr);
      }
    }
  }
  if ( neighparts.size() > 0 )
  {
    neigh->add_ghosts(neighparts);
    if ( neigh->get_bucket_type() == UNDERGROUND )
      neigh->put_have_ghost_particles(true);
    neigh->mark_active();
  }

  ghost_part_flag = true;
  return 0;
}

