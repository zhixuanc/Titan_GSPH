/*
 * =====================================================================================
 *
 *       Filename:  bgmesh.h
 *
 *    Description:  
 *
 *        Created:  08/31/2010 02:46:57 PM EDT
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

#ifndef BGMESH__H
#define BGMESH__H

#include <vector>
using namespace std;

#include <hashtab.h>
#include <bnd_image.h>
#include <properties.h>


void update_bgmesh(
    //! Particle Hash-table
    HashTable *, 
    //! Bucket Hash-table
    HashTable *,
    //! Material properties
    MatProps *,
    //! Process ID
    int ,
    //! check if Background Mesh was changed
    int *
   );

void search_bnd_images(
    //! ProcessID
    int ,
    //! Hash-table of particles
    HashTable *,
    //! Hash-table of buckets
    HashTable *,
    //! Vector of Boundary reflections
    vector<BndImage> *,
    //! flag to reset image-table
    int
   );

#endif
