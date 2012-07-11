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
#include <iostream>
using namespace std;

#include <buckstr.h>

struct ColumnHead 
{
  int xind, yind, proc;
  unsigned key[KEYLENGTH];

  // constructor
  ColumnHead (int i, int j, unsigned keyi[])
  {
    xind = i;
    yind = j;
    proc = 0;
    for (i = 0; i < KEYLENGTH; i++)
      key[i] = keyi[i];
  }

  bool operator < (const ColumnHead & rhs) const
  {
    if ( key[0] < rhs.key[0] )
      return true;
    else if ( key[0] > rhs.key[0] )
      return false;
    else if ( key[1] < rhs.key[1] )
      return true;
    else 
      return false;
  }
};

//! Write Background mesh and particle data to HDF5 file
void createfunky(
                 //! Numboer of processes in a multiproc run
                 int ,
                 //! Size of Hash Table constants array
                 int,
                 //! Hash Table Constants
                 double *,
                 //! Background grid data
                 vector <BucketStruct> &, 
                 //! partition table keys
                 vector <unsigned> &
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

#endif // PREPROCESS__H
