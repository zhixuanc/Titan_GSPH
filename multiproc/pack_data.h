
/*
 * =====================================================================================
 *
 *       Filename:  pack_data.h
 *
 *    Description:  
 *
 *        Created:  01/31/2011 05:28:38 PM EST
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

#ifndef PACK_DATA__H
#  define PACK_DATA__H

#  include <constants.h>

struct PARTICLE_PACK
{
  // Integers
  int ghost;
  unsigned key[DIMENSION];

  //doubles
  double mass;
  double smlen;
  double coords[DIMENSION];
  double state_vars[NO_OF_EQNS];
};
typedef struct PARTICLE_PACK ParticlePack;

struct BUCKET_PACK
{
  // Integers
  int myprocess;
  int activeflag;
  int bucket_type;
  int particles_type;
  int NumParticles;
  int neigh_proc[NEIGH_SIZE];

  // Unsigned integers
  unsigned key[KEYLENGTH];
  unsigned neighs[NEIGH_SIZE * KEYLENGTH];
  unsigned particles[MAX_PARTICLES_PER_BUCKET * KEYLENGTH];

  // Doubles
  double mincrd[DIMENSION];
  double maxcrd[DIMENSION];
  double poly1[4];
  double poly2[4];
  double bndx[PARTICLE_DENSITY];
  double bndy[PARTICLE_DENSITY];
  double fric[PARTICLE_DENSQRD * DIMENSION];
};
typedef struct BUCKET_PACK BucketPack;

#endif // PACK_DATA__H
