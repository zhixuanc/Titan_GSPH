
/*
 * =====================================================================================
 *
 *       Filename:  init_piles.cc
 *
 *    Description:  
 *
 *        Created:  07/29/2010 11:28:35 AM EDT
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

#include <vector>
#include <iostream>
#include <cassert>
using namespace std;

#include <mpi.h>
#include <hashtab.h>
#include <hilbert.h>
#include <bucket.h>
#include <particle.h>
#include <multiproc.h>
#include <properties.h>

#include "sph_header.h"

//#define CYLINDERICAL_PILE
//#define TRAPIZOIDAL_PILE

double
pile_shape(double crd[], double pcen[], double major, double minor, double hmax)
{
  int i;
  double vec[DIMENSION];
  double d = 0;

  for (i = 0; i < DIMENSION; i++)
    vec[i] = crd[i] - pcen[i];
#if defined CYLINDERICAL_PILE
  d = sqrt((vec[0] * vec[0]) + (vec[1] * vec[1]));
  if (d < major)
    return hmax;
  else
    return 0;
#elif defined TRAPIZOIDAL_PILE
  double slope = 0.79543;

  if ((vec[0] > 0) && (vec[0] < major))
    d = hmax - vec[0] * slope;
  if (d < 0.)
    d = 0;
  return d;
#else
  return (hmax - pow(vec[0] / major, 2) - pow(vec[1] / minor, 2));
#endif
}

void
init_piles(HashTable * P_table, HashTable * BG_mesh,
           PileProps * pileprops, MatProps * matprops, int numproc, int myid)
{

  int i, j, k, ierr;
  unsigned key[KEYLENGTH];
  unsigned keylen = KEYLENGTH;
  double mincrd[DIMENSION], maxcrd[DIMENSION];
  double mindom[DIMENSION], maxdom[DIMENSION];
  double pilecen[DIMENSION], normc[DIMENSION];
  double seed[DIMENSION], pcrd[DIMENSION];
  double bnd1[DIMENSION], bnd2[DIMENSION], poly[DIMENSION + 1];
  Bucket *Curr_buck = NULL;

  // direction indices on upper bucket
  int Up[DIMENSION] = { 0, 0, 2 };
  int num_particle = 0;

  // start putting piles
  double mass = matprops->particle_mass;
  double smlen = matprops->smoothing_length;
  double dx = smlen;
  double dx2 = 0.5 * dx;

  // get min-max domain from hashtable, for key generation
  for (i = 0; i < DIMENSION; i++)
  {
    mindom[i] = *(P_table->get_minDom() + i);
    maxdom[i] = *(P_table->get_maxDom() + i);
  }

  // Pile information
  double major = pileprops->majorrad;
  double minor = pileprops->minorrad;
  double pheight = pileprops->pileheight;
  Key center = pileprops->CenBucket;

  double major_sqrd = major * major;
  double minor_sqrd = minor * minor;

  Bucket * Cen_buck = (Bucket *) BG_mesh->lookup(center);
  double tempout[DIMENSION];
  for ( i = 0; i < DIMENSION; i++)
    tempout[i] = 0.;

#ifdef MULTI_PROC
  if (Cen_buck)
  {
    tempout[0] = pileprops->xCen;
    tempout[1] = pileprops->yCen;
    tempout[2] = Cen_buck->get_bndZ (tempout);
  }
  // sync with all procs
  MPI_Allreduce(tempout, pilecen, DIMENSION, MPI_DOUBLE, 
                MPI_SUM, MPI_COMM_WORLD);
  pileprops->xCen = pilecen[0];
  pileprops->yCen = pilecen[1];
  pileprops->zCen = pilecen[2];
#else
  assert(Cen_buck);
  pilecen[0] = pileprops->xCen;
  pilecen[1] = pileprops->yCen;
  pilecen[2] = Cen_buck->get_bndZ (pilecen);
  pileprops->zCen = pilecen[2];
#endif

  // add particles 
  HTIterator * itr = new HTIterator (BG_mesh);
  Bucket * Bnd_buck;

  while ((Bnd_buck = (Bucket *) itr->next ()))
    if (Bnd_buck->get_bucket_type () == MIXED)
    {
      for (i = 0; i < DIMENSION; i++)
      {
        mincrd[i] = *(Bnd_buck->get_mincrd () + i);
        maxcrd[i] = *(Bnd_buck->get_maxcrd () + i);
      }

      // negative-negative corner
      bnd1[0] = mincrd[0];
      bnd1[1] = mincrd[1];

      // positive-positive corner
      bnd2[0] = maxcrd[0];
      bnd2[1] = maxcrd[1];

      double d1 = pow(bnd1[0] - pilecen[0], 2) + pow(bnd1[1] - pilecen[1], 2);
      double d2 = pow(bnd2[0] - pilecen[0], 2) + pow(bnd1[1] - pilecen[1], 2);
      double d3 = pow(bnd2[0] - pilecen[0], 2) + pow(bnd2[1] - pilecen[1], 2);
      double d4 = pow(bnd1[0] - pilecen[0], 2) + pow(bnd2[1] - pilecen[1], 2);

      double dmax = max(max(d1, d2), max(d3, d4));
      double dmin = min(min(d1, d2), min(d3, d4));
      
      if ((dmax < major_sqrd) || (dmin < major_sqrd))
      {
        int Nx = (int) round((maxcrd[0] - mincrd[0]) / dx);
        int Ny = (int) round((maxcrd[1] - mincrd[1]) / dx);

        for (i = 0; i < Nx; i++)
          for (j = 0; j < Ny; j++)
          {
            seed[0] = mincrd[0] + dx2 + i * dx;
            seed[1] = mincrd[1] + dx2 + j * dx;
            seed[2] = Bnd_buck->get_bndZ(seed);
            double htemp = 0;

            if (Bnd_buck->contains(seed))
              htemp = pile_shape(seed, pilecen, major, minor, pheight);

            if (htemp > 0)
            {
              Curr_buck = Bnd_buck;
              pcrd[0] = seed[0];
              pcrd[1] = seed[1];
              int Nz = (int) round(htemp / dx);

              for (k = 0; k < Nz; k++)
              {
                pcrd[2] = seed[2] + dx2 + k * dx;
                for (int ii = 0; ii < DIMENSION; ii++)
                  normc[ii] = (pcrd[ii] - mindom[ii]) /
                    (maxdom[ii] - mindom[ii]);

                if (!(Curr_buck->contains(pcrd))
                    && (Curr_buck->which_neigh_proc(Up) > -1))
                {
                  Curr_buck =
                    (Bucket *) BG_mesh->lookup(Curr_buck->which_neigh(Up));
                  assert(Curr_buck);
                }
                HSFC3d (normc, &keylen, key);

                // check for duplicates
                if (P_table->lookup(key))
                {
                  fprintf(stderr, "ERROR: Trying to add particle "
                                  "twice on same location.\n");
                  exit(1);
                }
                Particle * pnew = new Particle(key, pcrd, mass, smlen, 0);
               
                // add to hash-table
                P_table->add(key, pnew);
                num_particle++;
                Key tmpkey(key);
                Curr_buck->add_real_particle(tmpkey);
              }
            }
          }
      }
    }

  int ntotal = 0;

#ifdef MULTI_PROC
  MPI_Allreduce (&num_particle, &ntotal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if ( ntotal < 100 )
  {
    if ( myid == 0 )
    {
      fprintf (stderr, "No. of particles = %d\n", ntotal);
      fprintf (stderr, "Not enough particles.\n");
    }
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
#else
  if ( num_particle < 100 )
  {
    fprintf (stderr, "No. of particles = %d\n", ntotal);
    fprintf (stderr, "Not enough particles.\n");
    exit (1);
  }
#endif
 
  //  mark all neighbors with real particles active
  itr->reset();
  while ((Curr_buck = (Bucket *) itr->next()))
  {
    Curr_buck->mark_inactive ();
    Key *nkey = Curr_buck->get_neighbors ();
    for (i = 0; i < NEIGH_SIZE; i++)
      // the guest buckets are not moved yet, 
      // hence only check buckets which belong to current proc
      if (*(Curr_buck->get_neigh_proc () + i) == myid )
      {
        Bucket *neigh = (Bucket *) BG_mesh->lookup (nkey[i]);
        if (neigh->has_real_particles ()) 
        {
          Curr_buck->mark_active ();
          break; 
        }
      }
  }

  delete itr;

  return;
}
