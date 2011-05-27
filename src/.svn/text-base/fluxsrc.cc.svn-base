/*
 * =====================================================================================
 *
 *       Filename:  fluxsrc.cc
 *
 *    Description:  
 *
 *        Created:  11/23/2010 02:26:07 PM EST
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
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <hilbert.h>
#include <bucket.h>
#include <particle.h>
#include <sph_header.h>
#include <properties.h>


void init_fluxsrc (HashTable *P_table, HashTable *BG_mesh,
                   MatProps *matprops, FluxProps *fluxprops)
{

  int      i,j,k;
  unsigned key[KEYLENGTH];
  unsigned keylen = KEYLENGTH;
  double   vel[DIMENSION], coord[DIMENSION]; 
  double   mindom[DIMENSION], maxdom[DIMENSION];
  double   seed[DIMENSION], normc[DIMENSION];
  double   fluxcrd[DIMENSION];
  Bucket     *bucketsrc, *bndbucket, *curr_bucket;
  Key      keystr;

  // search directions
  int Left[2] = {1, 0};
  int Down[2] = {0, 1};
  int   Up[2] = {0, 2};

  // particle properties
  double ma = matprops->particle_mass;
  double hl = matprops->smoothing_length;
  double tangvel = fluxprops->tangvel;
  double dx = hl;
  double dx2 = 0.5*dx;
  double plength = 0.1;
  double pheight = 0.012; 
  double cosa  = 0.8746197;
  double sina  = 0.4848096;
  // define parabolic profile
  double c1  = plength;
  double c2  = plength/(pheight*pheight);
  vel[0] = tangvel*cosa;
  vel[1] = tangvel*sina;
  for (i=0; i < DIMENSION; i++)
  {
    mindom[i] = *(P_table->get_minDom()+i);
    maxdom[i] = *(P_table->get_maxDom()+i);
  }

  // get source coodinate and source bucket
  Key  bucketsrckey = fluxprops->bucketsrckey;
  bucketsrc = (Bucket *) BG_mesh->lookup(bucketsrckey);
  assert(bucketsrc);  
  fluxcrd[0] = fluxprops->xSrc;
  fluxcrd[1] = bucketsrc->get_bndZ(fluxcrd);

  int nx = (int) round(plength/dx);
  //int nz = (int) round(pheight/dx);

  bndbucket = bucketsrc;
  for ( i=0; i < nx; i++ )
  {
    seed[0]= fluxcrd[0] - i*dx;
    seed[1] = bndbucket->get_bndZ(seed);
    if ( !bndbucket->contains ( seed ) )
    {
      Key *nkey = bndbucket->get_neighbors();
      for (j=0; j<NEIGH_SIZE; j++)
        if ( *(bndbucket->get_neigh_proc()+j) > -1 )
        {
          Bucket *nbucket = (Bucket *) BG_mesh->lookup(nkey[j]);
          assert(nbucket);
          if ( (nbucket->get_bucket_type () == MIXED) && 
                nbucket->contains(seed) )
          {
            bndbucket = nbucket;
            break;
          }
        }
    }
    if ( !bndbucket->contains(seed) )
    {
      fprintf(stderr,"Can\'t find the correct bucket for seed. Quitting\n");
      exit(1);
    }

    curr_bucket = bndbucket;
    coord[0] = seed[0];
    double  tmpx = coord[0] - fluxcrd[0];
    double  tmpz = sqrt((tmpx + c1)/c2);
    int nz = (int) round(tmpz/dx);
    if ( nz < 2 )
      nz = 2;
    for ( j=0; j < nz; j++ )
    {
      coord[1] = seed[1] + dx2 + j*dx;
      for ( k=0; k < DIMENSION; k++ )
        normc[k] = (coord[k]-mindom[k])/(maxdom[k]-mindom[k]);
      HSFC2d ( normc, &keylen, key );
      Particle *p = new Particle ( key, coord, ma, hl, vel);
      P_table->add(key, p);

      for ( k=0; k < KEYLENGTH; k++ )
        keystr.key[k] = key[k];

      if ( !curr_bucket->contains(coord) )
      {
        Key up = curr_bucket->which_neigh(Up);
        curr_bucket = (Bucket *) BG_mesh->lookup(up);
        if ( !curr_bucket->contains(coord) )
        {
          fprintf(stderr, "Error!!! Can't find bucket for particle\n");
          exit(1);
        }
      }
      curr_bucket->add_real_particle(keystr);
    } 
  }

  return;
}

/* ***********************************************
 *
 *
 *
 * ***********************************************/

void update_fluxsrc (HashTable *P_table, HashTable *BG_mesh, 
                     MatProps *matprops, FluxProps *fluxprops,
                     TimeProps *timeprops)
{
  int      i, j;
  unsigned key[KEYLENGTH];
  unsigned keylen = KEYLENGTH;
  double   vel[DIMENSION], coord[DIMENSION];
  double   mindom[DIMENSION], maxdom[DIMENSION];
  double   virtcrd[DIMENSION], normc[DIMENSION];
  double   fluxcrd[DIMENSION];
  Key      keystr;
  int      Up[2] = { 0, 2 };

  // search directions
  double   ma = matprops->particle_mass;
  double   hl = matprops->smoothing_length;
  double   tangvel = fluxprops->tangvel;
  double   dx = hl;
  double   dx2 = 0.5*dx;
  Bucket     *bucketsrc, *bnd_bucket = NULL, *curr_bucket=NULL;

  vel[0]  = tangvel*0.8746;
  vel[1]  = tangvel*0.4848;

  for (i=0; i < DIMENSION; i++)
  {
    mindom[i] = *(P_table->get_minDom()+i);
    maxdom[i] = *(P_table->get_maxDom()+i);
  }
  Key bucketsrckey = fluxprops->bucketsrckey;
  bucketsrc = (Bucket *) BG_mesh->lookup(bucketsrckey);
  assert(bucketsrc);
  fluxcrd[0] = fluxprops->xSrc;
  fluxcrd[1] = bucketsrc->get_bndZ(fluxcrd);

  int nz = (int) round(0.012/dx);
  if ( nz < 6 )
    nz = 6;

  double time = timeprops->time;
  curr_bucket = bucketsrc;
  for (i=0; i<nz; i++)
  {
    coord[0] = fluxcrd[0];
    coord[1] = fluxcrd[1] + dx2 + i*dx;
    for (j=0; j<DIMENSION; j++)
    {
      virtcrd[j] = coord[j] - time*vel[j];
      normc[j] = (virtcrd[j]-mindom[j])/(maxdom[j]-mindom[j]);
    }
    HSFC2d (normc, &keylen, key);
    Particle *p = new Particle (key, coord, ma, hl, vel);
    P_table->add(key, p);
    
    for (j=0; j < KEYLENGTH; j++)
      keystr.key[j] = key[j];

    if ( !curr_bucket->contains(coord) )
    {
      Key upkey = curr_bucket->which_neigh(Up);
      curr_bucket = (Bucket *) BG_mesh->lookup(upkey);
    }
    curr_bucket->add_real_particle(keystr);
  }
}
