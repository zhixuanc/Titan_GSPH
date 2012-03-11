/*
 * =====================================================================================
 *
 *       Filename:  bucketstr.h
 *
 *    Description:  BucketStruct is defined to get modular data from 
 *                  preprocessor, in HDF5 format.
 *
 *        Created:  07/16/2010 01:42:40 PM EDT
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

#ifndef BUCKETSTR__H
#define BUCKETSTR__H

#include <constants.h>

typedef struct 
{
  unsigned key[KEYLENGTH];
  unsigned neighs[NEIGH_SIZE*KEYLENGTH];
  int      buckettype;
  int      boundary;
  int      myproc;
  int      neigh_proc[NEIGH_SIZE];
  double   xcoord[2];
  double   ycoord[2];
  double   zcoord[2];
  double   elev[4];

} BucketStruct;

#endif // BUCKETSTR__H
