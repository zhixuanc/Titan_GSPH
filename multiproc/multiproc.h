
/*
 * =====================================================================================
 *
 *       Filename:  multiproc.h
 *
 *    Description:  
 *
 *        Created:  03/14/2011 08:55:06 PM EDT
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

#ifndef MULTIPROC__H
#  define MULTIPROC__H

#  include <vector>
using namespace std;

#  include <hashtab.h>
#  include <bnd_image.h>

#  include "buckhead.h"

//! register MPI_structs 
void GMFG_new_MPI_Datatype();

void send_foreign_images(
                          //! My process-id
                          int,
                          //! Total size of MPI pool
                          int,
                          //! Hash-Table of SPH particles
                          HashTable *,
                          //! Hash-Table of Buckets
                          HashTable *,
                          //! STL Vector of ghost reflections
                          vector < BndImage > *);

//! Update data across, inter-proc boundaries
void move_data(
                //! Number of total processes in the mix
                int,
                //! My process rank 
                int,
                //! Hash-Table of SPH particles
                HashTable *,
                //! Hash-Table of Buckets 
                HashTable *
                //! STL vector of ghost reflections
  );

//! Update ghost reflections from other processors
int move_bnd_images(
                     //! My process ID
                     int,
                     //! Total number of process
                     int,
                     //! HashTable of particles
                     HashTable *,
                     //! HashTable of buckets
                     HashTable *,
                     //! STL Vector of Images
                     vector < BndImage >);

//! repartion the domain if load-balance has changed
int repartition(
                 //! STL vector of Partition Table Keys
                 vector < BucketHead > &,
                 //! Hash-Table of SPH particles
                 HashTable *,
                 //! Hash-Table of Buckets
                 HashTable *);

//! delete guest buckets and particles
void delete_guest_buckets(
                           //! Hash-Table of Buckets
                           HashTable *,
                           //! Hash-Table of particles
                           HashTable *);

#endif // MULTIPROC__H
