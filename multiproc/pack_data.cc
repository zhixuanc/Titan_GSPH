/*
 * =====================================================================================
 *
 *       Filename:  pack_data.cc
 *
 *    Description:  
 *
 *        Created:  02/03/2011 11:55:04 AM EST
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
#include <iostream>
#include <vector>
using namespace std;

#include <constants.h>
#include <hashtab.h>
#include <bucket.h>
#include <particle.h>
#include "pack_data.h"

// for BSFC repartitioning scheme
void pack_bucket (BucketPack *buckpack, Bucket *sendbuck, int process)
{
  int j, i;
  buckpack->myprocess = process;
  buckpack->bucket_type = sendbuck->bucket_type;
  buckpack->particles_type = sendbuck->particles_type;
  buckpack->activeflag = (int) sendbuck->active;
  for (j=0; j<KEYLENGTH; j++)
    buckpack->key[j] = sendbuck->key.key[j];

  for (i = 0; i < NEIGH_SIZE; i++)
  {
    buckpack->neigh_proc[i] = sendbuck->neigh_proc[i];
    for ( j=0; j<KEYLENGTH; j++)
      buckpack->neighs[i*KEYLENGTH+j] = sendbuck->neighbors[i].key[j];
  }

  vector<Key> particles = sendbuck->particles;
  int psize = (int) particles.size();
  if ( psize > MAX_PARTICLES_PER_BUCKET)
  {
    cerr << "Number of particles exceed Maximum allowable limit." << endl;
    cerr << "Error at line " << __LINE__ <<" in file " << __FILE__ << endl;
    exit (1);
  }

  // pack particle keys within the bucket
  buckpack->NumParticles = psize;
  for ( i=0; i < psize; i++ )
    for ( j=0; j < KEYLENGTH; j++ )
      buckpack->particles[i*KEYLENGTH+j] = particles[i].key[j];

  // bucket upper and lower limits
  for (i=0; i < DIMENSION; i++)
  {
    buckpack->mincrd[i] = sendbuck->mincrd[i];
    buckpack->maxcrd[i] = sendbuck->maxcrd[i];
  }

  // boundary fucntion
  for (i=0; i<4; i++)
  {
    buckpack->poly1[i] = sendbuck->lower_tri[i];
    buckpack->poly2[i] = sendbuck->upper_tri[i];
  }

  // boundary points and value of friction coeficents
  if ( sendbuck->get_bucket_type() == MIXED )
  {
    for (i=0; i < PARTICLE_DENSITY; i++)
    {
      buckpack->bndx[i] = sendbuck->bnd_xcrd[i];
      buckpack->bndy[i] = sendbuck->bnd_ycrd[i];
    }
    for (i = 0;  i < PARTICLE_DENSQRD*DIMENSION; i++)
      buckpack->fric[i] = sendbuck->f_coef[i];
  }
  else
  {
    for (i=0; i < PARTICLE_DENSITY; i++)
    {
      buckpack->bndx[i] = 0.;
      buckpack->bndy[i] = 0.;
    }

    for (i=0; i < PARTICLE_DENSQRD*DIMENSION; i++)
      buckpack->fric[i] = 0.;
  }
  return;
}


void pack_particles (Particle *psend, ParticlePack *pack_array)
{
  int j;
  pack_array->ghost = (int) psend->ghost;
  pack_array->mass  = psend->mass;
  pack_array->smlen = psend->smlen;
   
  for (j=0; j < KEYLENGTH; j++)
    pack_array->key[j] = psend->key.key[j];

  for (j=0; j < DIMENSION; j++)
    pack_array->coords[j] = psend->coord[j];

  for (j=0; j < NO_OF_EQNS; j++)
    pack_array->state_vars[j] = psend->state_vars[j];

  return;
}


void unpack_bucket (BucketPack *recvdBuck, Bucket *buck, int myid)
{

  int i, j;
  vector<Key> plist;

  buck->myprocess = myid;
  buck->bucket_type = recvdBuck->bucket_type;
  buck->particles_type = recvdBuck->particles_type;
  buck->active = (bool) recvdBuck->activeflag;

  for ( i=0; i < KEYLENGTH; i++ )
    buck->key.key[i] = recvdBuck->key[i];

  for ( i=0; i < NEIGH_SIZE; i++ )
  {
    buck->neigh_proc[i] = recvdBuck->neigh_proc[i];
    for ( j=0; j < KEYLENGTH; j++ )
      buck->neighbors[i].key[j]  = recvdBuck->neighs[i*KEYLENGTH+j];
  }

  for ( i=0; i < DIMENSION; i++ )
  {
    buck->mincrd[i] = recvdBuck->mincrd[i];
    buck->maxcrd[i] = recvdBuck->maxcrd[i];
  }

  for ( i=0; i < 4; i++ )
  {
    buck->lower_tri[i] = recvdBuck->poly1[i];
    buck->upper_tri[i] = recvdBuck->poly2[i];
  }

  // if the bucket is a boundary bucket, get unpack boundary points
  // and friction coeficients
  if ( recvdBuck->bucket_type == MIXED )
  {
    for (i=0; i<PARTICLE_DENSITY; i++)
    {
      buck->bnd_xcrd[i] = recvdBuck->bndx[i];
      buck->bnd_ycrd[i] = recvdBuck->bndy[i];
    }
    for (i=0; i<PARTICLE_DENSQRD*DIMENSION; i++)
      buck->f_coef[i]=recvdBuck->fric[i];
  }
  else
  {
    for (i=0; i<PARTICLE_DENSITY; i++)
    {
      buck->bnd_xcrd[i] = 0.;
      buck->bnd_ycrd[i] = 0.;
    }
    for (i=0; i<PARTICLE_DENSQRD*DIMENSION; i++)
      buck->f_coef[i] = 0.;
  }

  // unpack particle keys
  buck->particles.clear();
  int psize = recvdBuck->NumParticles;
  Key tmpkey;
  for ( i=0; i < psize; i++ )
  {
    for ( j=0; j < KEYLENGTH; j++ )
      tmpkey.key[j] = recvdBuck->particles[i*KEYLENGTH+j];
    buck->particles.push_back(tmpkey);
  }
  return;
}

void unpack_particle (ParticlePack *packet, Particle *part)
{
  int i;
  part->ghost = (bool) packet->ghost;

  for ( i=0; i < KEYLENGTH; i++ )
    part->key.key[i] = packet->key[i];
  
  part->mass  = packet->mass;
  part->smlen = packet->smlen;
  
  for ( i=0; i < DIMENSION; i++ )
    part->coord[i] = packet->coords[i];

  for ( i=0; i < NO_OF_EQNS; i++ )
    part->state_vars[i] = packet->state_vars[i];

  return;
}
