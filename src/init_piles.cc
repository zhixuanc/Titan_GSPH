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
#include <config.h>
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
#include <properties.h>

#include "sph_header.h"

#define CYLINDERICAL_PILE

void check_bnd_pt (double bndpt[], double mincrd[], 
                   double maxcrd[], double bndCoeff[])
{
  if ( bndpt[1] < mincrd[1] )
  {
    bndpt[1] = mincrd[1];
    bndpt[0] = -(bndCoeff[1]*bndpt[1] + bndCoeff[2])/bndCoeff[0];
  }
  else if ( bndpt[1] > maxcrd[1] )
  {
    bndpt[1] = maxcrd[1];
    bndpt[0] = -(bndCoeff[1]*bndpt[1] + bndCoeff[2])/bndCoeff[0];
  }
  return;
}

double  pile_shape (double crd[], double pcen[], 
                    double major, double minor, double hmax)
{
  int i;
  double vec[DIMENSION];
  double d = 0;

  for (i=0; i < DIMENSION; i++ )
    vec[i] = crd[i] - pcen[i];
#if defined CYLINDERICAL_PILE
#  ifdef THREE_D
  d = sqrt((vec[0]*vec[0])+(vec[1]*vec[1]));
#  else
  d = vec[0];
#  endif
  if ( d < major )
    return hmax;
  else 
    return 0;
#elif defined PARABOLOID
#  ifdef THREE_D
  return (hmax - pow(vec[0]/major,2) - pow(vec[1]/minor,2));
#  else
  return (hmax - pow(vec[0]/major,2));
#  endif
#endif
}

void init_piles (HashTable *P_table, HashTable *BG_mesh,
                 PileProps *pileprops, MatProps *matprops,
                 int numproc, int myid)
{

  int      i,j,k,ipile, ierr;
  unsigned key[KEYLENGTH];
  unsigned keylen = KEYLENGTH;
  double   mincrd[DIMENSION], maxcrd[DIMENSION]; 
  double   mindom[DIMENSION], maxdom[DIMENSION];
  double   pilecen[DIMENSION], normc[DIMENSION];
  double   seed[DIMENSION], pcrd[DIMENSION];
  double   bnd1[DIMENSION], bnd2[DIMENSION], poly[DIMENSION+1];
  Bucket   *Curr_bk = NULL;
  Key      tmpkey;
 
#ifdef THREE_D
  int     Up[DIMENSION] = { 0, 0, 2 };
#else
  int     Up[DIMENSION] = { 0, 2 };
#endif

  // start putting piles
  int numpile=pileprops->NumPiles;
  double mass = matprops->particle_mass;
  double smlen = matprops->smoothing_length;
  double dx = smlen;
  double tol = 0.1*dx;
  double dx2 = 0.5*dx;

  // get min-max domain from hashtable, for key generation
  for (i=0; i<DIMENSION; i++)
  {
    mindom[i] = *(P_table->get_minDom()+i);
    maxdom[i] = *(P_table->get_maxDom()+i);
  }
  if ( numpile < 1 )
  {
    cerr<<"ERROR: unable to run as no pile is defined."<<endl;
    exit(0);
  }

  ipile = 0;
  int  root;
  // Pile information
  double major = pileprops->majorrad[ipile];
  double minor = pileprops->minorrad[ipile];
  double pheight = pileprops->pileheight[ipile];
  Key center = pileprops->CenBucket[ipile];

  Bucket *Cen_bk = (Bucket *) BG_mesh->lookup(center);
  pilecen[0] = pileprops->xCen[ipile];
#ifdef THREE_D
  pilecen[1] = pileprops->yCen[ipile];
#endif

  double zCen;
#ifdef MULTI_PROC
  if ( Cen_bk )
  {
    zCen = Cen_bk->get_bndZ(pilecen);
  } 
  else
    zCen = 0;

  // sync with all procs
  MPI_Allreduce( &zCen, &pilecen[DIMENSION-1], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
  assert(Cen_bk);
  pilecen[DIMENSION-1] = Cen_bk->get_bndZ(pilecen);
#endif

  // add particles 
  HTIterator *itr = new HTIterator(BG_mesh);
  Bucket *Bnd_bk;
  while ((Bnd_bk = (Bucket *) itr->next()))
    if ( Bnd_bk->get_bucket_type() == MIXED )
    {
      for (i=0; i < DIMENSION; i++ )
      {
        mincrd[i] = *(Bnd_bk->get_mincrd()+i);
        maxcrd[i] = *(Bnd_bk->get_maxcrd()+i);
      }
      // get bnd polynomial
      for (i=0; i < DIMENSION+1; i++)
        poly[i] = *(Bnd_bk->get_bndcoeff()+i);

      // 1st - boundary point
      bnd1[0]=mincrd[0];
      bnd1[1]=Bnd_bk->get_bndZ(bnd1);
      check_bnd_pt (bnd1, mincrd, maxcrd, poly);

      // 2nd - boundary point
      bnd2[0]=maxcrd[0];
      bnd2[1]=Bnd_bk->get_bndZ(bnd2);
      check_bnd_pt (bnd2, mincrd, maxcrd, poly);

      double d1 = abs(bnd1[0]-pilecen[0]);
      double d2 = abs(bnd2[0]-pilecen[0]);

      if ( (d1 < major) || (d2 < major))
      {
        int Nx = (int) round((maxcrd[0]-mincrd[0])/dx);
        for (i=0; i < Nx; i++ )
        {
          seed[0] = mincrd[0] + dx2 + i*dx;
          seed[1] = Bnd_bk->get_bndZ(seed);
          double htemp = 0;
          if ( Bnd_bk->contains(seed) )
            htemp = pile_shape(seed, pilecen, major, minor, pheight);
          if ( htemp > 0 )
          {
            Curr_bk = Bnd_bk;
            pcrd[0] = seed[0];
            int Nz = (int) round(htemp/dx);
            for (k=0; k < Nz; k++ )
            {
              pcrd[1] = seed[1] + dx2 + k*dx;

              for (int ii=0; ii < DIMENSION; ii++)
                normc[ii] = (pcrd[ii]-mindom[ii])/(maxdom[ii]-mindom[ii]);
              if ( !(Curr_bk->contains(pcrd)) &&
                   (Curr_bk->which_neigh_proc(Up) > -1) )
              {
                  Curr_bk = (Bucket *) BG_mesh->lookup(Curr_bk->which_neigh(Up));
                  assert(Curr_bk);
              }
              HSFC2d(normc, &keylen, key);
              // check for duplicates
              if ( P_table->lookup(key) )
              {
                fprintf(stderr,"Bummer: Trying to add particle twice on same location\n");
                exit(1);
              }
              Particle *pnew = new Particle (key, pcrd, mass, smlen, 0);
              P_table->add(key, pnew);
              tmpkey.copy_key(key);
              Curr_bk->add_real_particle(tmpkey);
            }
          }
        }
      }
    }

  //  mark all neighbors with real particles active
  itr->reset();
  while ( (Curr_bk=(Bucket *) itr->next()) )
    if ( Curr_bk->get_plist().size() > 0 )
    {
      Key *nkey = Curr_bk->get_neighbors();
      for (i=0; i<NEIGH_SIZE; i++)
        if (*(Curr_bk->get_neigh_proc()+i) == myid)
        {
          Bucket *nbucket = (Bucket *) BG_mesh->lookup(nkey[i]);
          nbucket->mark_active();
        }
    }

  delete itr;
  return;
}
