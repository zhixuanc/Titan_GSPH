
/*
 * =====================================================================================
 *
 *       Filename:  funky2txt.cc
 *
 *    Description:  
 *
 *        Created:  03/10/2011 04:03:33 PM EST
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
#include <cstdlib>
using namespace std;

#include <buckstr.h>
#include "hdf5calls.h"

#define DATA "/Buckets"

int
main(int argc, char **argv)
{

  int i;
  char filename[20];

  if (argc != 2)
  {
    printf("USAGE: funky2ascii <proc number>\n");
    exit(1);
  }

  int myid = atoi(argv[1]);

  sprintf(filename, "funky%04d.h5", myid);

  hsize_t dims[2] = { 0, 0 };
  hid_t fp = GH5_fopen(filename, 'r');
  int num = GH5_getsize(fp, DATA, dims);

  num = (int) dims[0];

  BucketStruct *buckets = new BucketStruct[num];

  GH5_read_grid_data(fp, DATA, buckets);
  GH5_fclose(fp);

  sprintf(filename, "grid%04d.dat", myid);
  FILE *f2 = fopen(filename, "w");

  for (i = 0; i < num; i++)
    fprintf(f2, "%e, %e, %e, %e\n",
            buckets[i].xcoord[0], buckets[i].zcoord[0],
            buckets[i].xcoord[1], buckets[i].zcoord[1]);

  fclose(f2);
  delete[]buckets;
  return 0;
}
