
/*
 * =====================================================================================
 *
 *       Filename:  mpi_struct.cc
 *
 *    Description:  
 *
 *        Created:  01/31/2011 05:05:46 PM EST
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

#include <mpi.h>
#include <constants.h>
#include <bnd_image.h>

#include "pack_data.h"

MPI_Datatype BUCKET_TYPE;
MPI_Datatype PARTICLE_TYPE;
MPI_Datatype BND_IMAGE_TYPE;

void
GMFG_new_MPI_Datatype ()
{

  int blockcounts[3];
  MPI_Datatype types[3];
  MPI_Aint displs[3];
  int d;
  BucketPack *buck = new BucketPack;

  blockcounts[0] = 5 + NEIGH_SIZE;
  blockcounts[1] = KEYLENGTH * (1 + NEIGH_SIZE + MAX_PARTICLES_PER_BUCKET);
  blockcounts[2] = (2 * DIMENSION) + 8 + (2 * PARTICLE_DENSITY) +
    (PARTICLE_DENSQRD * DIMENSION);
  MPI_Address (&(buck->myprocess), &displs[0]);
  MPI_Address (&(buck->key[0]), &displs[1]);
  MPI_Address (&(buck->mincrd[0]), &displs[2]);

  types[0] = MPI_INT;
  types[1] = MPI_UNSIGNED;
  types[2] = MPI_DOUBLE;

  for (d = 2; d >= 0; d--)
    displs[d] -= displs[0];

  MPI_Type_struct (3, blockcounts, displs, types, &BUCKET_TYPE);
  MPI_Type_commit (&BUCKET_TYPE);

  //create the 2nd new d_type
  int blockcounts2[3];
  MPI_Datatype types2[3];
  MPI_Aint displs2[3];

  ParticlePack *particlePack = new ParticlePack;

  // 1 int , 2 unsigned, bunch of doubles
  blockcounts2[0] = 1;
  blockcounts2[1] = KEYLENGTH;
  blockcounts2[2] = 2 + DIMENSION + NO_OF_EQNS;

  // get adresses
  MPI_Address (&(particlePack->ghost), &displs2[0]);
  MPI_Address (&(particlePack->key), &displs2[1]);
  MPI_Address (&(particlePack->mass), &displs2[2]);

  types2[0] = MPI_INT;
  types2[1] = MPI_UNSIGNED;
  types2[2] = MPI_DOUBLE;

  for (d = 2; d >= 0; d--)
    displs2[d] -= displs2[0];

  MPI_Type_struct (3, blockcounts2, displs2, types2, &PARTICLE_TYPE);
  MPI_Type_commit (&PARTICLE_TYPE);

  // create 3rd datatype for Boundary Images
  int blockcounts3[3] = { 2, 2 * KEYLENGTH, DIMENSION + NO_OF_EQNS };
  MPI_Datatype type3[3] = { MPI_INT, MPI_UNSIGNED, MPI_DOUBLE };
  MPI_Aint displs3[3];

  // get adresses
  BndImage *bndimage = new BndImage ();

  MPI_Address (&(bndimage->buckproc), &(displs3[0]));
  MPI_Address (&(bndimage->bucket_key), &(displs3[1]));
  MPI_Address (&(bndimage->coord), &(displs3[2]));

  for (d = 2; d >= 0; d--)
    displs3[d] -= displs3[0];

  MPI_Type_struct (3, blockcounts3, displs3, type3, &BND_IMAGE_TYPE);
  MPI_Type_commit (&BND_IMAGE_TYPE);

  //New data types are created at this point
}
