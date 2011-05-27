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
  buckpack->bucket_type = sendbuck->get_bucket_type();
  buckpack->boundary_type = sendbuck->get_bndtype();
  buckpack->real_particles = (int) sendbuck->have_real_particles();
  buckpack->ghost_particles = (int) sendbuck->have_ghost_particles();
  buckpack->activeflag = (int) sendbuck->is_active();
  for (j=0; j<KEYLENGTH; j++)
    buckpack->key[j] = (sendbuck->getKey()).key[j];

  for (i = 0; i < NEIGH_SIZE; i++)
  {
    buckpack->neigh_proc[i] = *(sendbuck->get_neigh_proc()+i);
    Key *kptr = sendbuck->get_neighbors();
    for ( j=0; j<KEYLENGTH; j++)
      buckpack->neighs[i*KEYLENGTH+j] = (kptr+i)->key[j];
  }

  vector<Key> my_prtcls = sendbuck->get_plist();
  int psize = (int) my_prtcls.size();
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
      buckpack->particles[i*KEYLENGTH+j] = my_prtcls[i].key[j];

  // bucket upper and lower limits
  for (i=0; i < DIMENSION; i++)
  {
    buckpack->mincrd[i] = *(sendbuck->get_mincrd()+i);
    buckpack->maxcrd[i] = *(sendbuck->get_maxcrd()+i);
  }

  // boundary fucntion
#ifdef THREE_D
  for (i=0; i<4; i++)
#else
  for (i=0; i<3; i++)
#endif
    buckpack->poly[i] = *(sendbuck->get_bndcoeff()+i);   

  // boundary points and value of friction coeficents
  if ( sendbuck->get_bucket_type() == MIXED )
  {
    for (i=0; i < PARTICLE_DENSITY; i++)
    {
      buckpack->bndx[i] = *(sendbuck->get_bnd_xcrd()+i);
#ifdef THREE_D
      buckpack->bndy[i] = *(sendbuck->get_bnd_ycrd()+i);
#endif
    }

#ifdef THREE_D
    for (i=0; i < PARTICLE_DENSQRD; i++)
#else
    for (i=0; i < PARTICLE_DENSITY; i++)
#endif 
      buckpack->fric[i] = *(sendbuck->get_f_coef()+i);
  }
  else
  {
    for (i=0; i < PARTICLE_DENSITY; i++)
    {
      buckpack->bndx[i] = 0.;
#ifdef THREE_D
      buckpack->bndy[i] = 0.;
#endif
    }

#ifdef THREE_D
    for (i=0; i < PARTICLE_DENSQRD; i++)
#else
    for (i=0; i < PARTICLE_DENSITY; i++)
#endif 
      buckpack->fric[i] = 0.;
  } 
  return;
}


void pack_particles (Particle *psend, ParticlePack *pack_array)
{
  int j;
  pack_array->ghost = (int) psend->is_ghost();
  pack_array->mass  = psend->get_mass();
  pack_array->smlen = psend->get_smlen();
   
  for (j=0; j < KEYLENGTH; j++)
    pack_array->key[j] = (psend->getKey()).key[j];

  for (j=0; j < DIMENSION; j++)
    pack_array->coords[j] = *(psend->get_coords()+j);

  for (j=0; j < NO_OF_EQNS; j++)
    pack_array->state_vars[j] = *(psend->get_state_vars()+j);

  return;
}


void unpack_bucket (BucketPack *recvdBuck, Bucket *buck, int myid)
{

  int i, j;
  vector<Key> plist;

  buck->myprocess = myid;
  buck->bucket_type = recvdBuck->bucket_type;
  buck->boundary_type  = recvdBuck->boundary_type;
  buck->real_part_flag = (bool) recvdBuck->real_particles;
  buck->ghost_part_flag = (bool) recvdBuck->ghost_particles;
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

#ifdef THREE_D
  for ( i=0; i < 4; i++ )
#else
  for (i=0; i < 3; i++ )
#endif
    buck->poly[i] = recvdBuck->poly[i];

  // if the bucket is a boundary bucket, get unpack boundary points
  // and friction coeficients
  if ( recvdBuck->bucket_type == MIXED )
  {
    for (i=0; i<PARTICLE_DENSITY; i++)
    {
      buck->bnd_xcrd[i] = recvdBuck->bndx[i];
#ifdef THREE_D
      buck->bnd_ycrd[i] = recvdBuck->bndy[i];
#endif
    }
#ifdef THREE_D
    for (i=0; i<PARTICLE_DENSQRD; i++)
#else
    for (i=0; i<PARTICLE_DENSITY; i++)
#endif
      buck->f_coef[i]=recvdBuck->fric[i];
  }
  else
  {
    for (i=0; i<PARTICLE_DENSITY; i++)
    {
      buck->bnd_xcrd[i] = 0.;
#ifdef THREE_D
      buck->bnd_ycrd[i] = 0.;
#endif
    }
#ifdef THREE_D
    for (i=0; i<PARTICLE_DENSQRD; i++)
#else
    for (i=0; i<PARTICLE_DENSITY; i++)
#endif
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
