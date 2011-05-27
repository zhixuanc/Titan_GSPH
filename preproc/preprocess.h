/*
 * =====================================================================================
 *
 *       Filename:  preprocess.h
 *
 *    Description:  
 *
 *        Created:  03/27/2010 05:40:45 PM EDT
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

#ifndef PREPROCESS__H
#define PREPROCESS__H

#include <vector>
using namespace std;

#include <buckstr.h>

struct buckhead {
  int index;
  int proc;
  double xcrd;
  unsigned segment[KEYLENGTH];

  // constructor
  buckhead (int i, unsigned key[], double crd)
  {
    index = i;
    xcrd = crd;
    for (i=0; i < KEYLENGTH; i++ )
      segment[i] = key[i];
  }
};
typedef struct buckhead BuckHead;

// 2-D inclined plane
void GIS_get_elevation ( int, double , double *);

//! Write Background mesh and particle data to HDF5 file
void createfunky(
                 //! Numboer of processes in a multiproc run
                 int ,
                 //! Size of Hash Table constants array
                 int,
                 //! Hash Table Constants
                 double *,
                 //! Size of array containing Background grid data
                 int ,
                 //! Background grid data
                 BucketStruct *
                );

//! Generate key based on location
void determine_the_key (
                //! Normalized coordinates
                double norm_coord[],
                //! Keylength
                unsigned keylen,
                //! key to be generated
                unsigned key[],
                //! Max key
                unsigned maxkey[],
                //! Min key
                unsigned minkey[]
                );

#ifdef __cplusplus
extern "C" {
#endif

void dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);

#ifdef __cplusplus
}
#endif

#endif // PREPROCESS__H
