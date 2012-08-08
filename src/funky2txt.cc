
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <vector>
#include <map>
using namespace std;

#ifdef MULTI_PROC
#include <mpi.h>
#endif

#include <buckstr.h>
#include "hdf5calls.h"

#define DATA "/Buckets"

// define vertex
struct Vertex 
{
  int index;
  int key;
  double coord[3];

  Vertex (int i, double *xyz)
  {
    index = i;
    for (i = 0; i < 3; i++)
      coord[i] = xyz[i];
  }
};

int genkey (double * coord, double * min, double * max)
{
  int i, key;
  int constant[3] = {1000000, 10000, 100};
  key = 0;
  for (i = 0; i < 3; i++)
    key += (int) ((coord[i] - min[i])/(max[i] - min[i]) * constant[i]);
  return key;
}

int main(int argc, char **argv)
{

  int i, key;
  char filename[20];
  double coord[3], mincrd[3], maxcrd[3];

#ifdef MULTI_PROC
  i = MPI_Init (&argc, &argv);
#endif

  if (argc != 2)
  {
    printf("USAGE: funky2ascii <proc number>\n");
    exit(1);
  }

  int myid = atoi(argv[1]);

  sprintf(filename, "funky%04d.h5", myid);

  hsize_t dims[2] = { 0, 0 };
  hid_t fh5 = GH5_fopen(filename, 'r');
  int num = GH5_getsize(fh5, DATA, dims);

  num = (int) dims[0];

  BucketStruct *buckets = new BucketStruct[num];

  GH5_read_grid_data(fh5, DATA, buckets);
  GH5_fclose(fh5);

  sprintf (filename, "grid%04d.txt", myid);
  FILE * fp = fopen (filename, "w");
  
  int num_vert = 0;
  vector <Vertex> vert_table;
  for (i = 0; i < num; i++)
  {
    // vertex 0 
    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    vert_table.push_back (Vertex (num_vert++, coord));        

    // vertex 1
    if ( buckets[i].neigh_proc[22] != myid )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[0];
      vert_table.push_back ( Vertex (num_vert++, coord));
    }

    // vertex 2
    if ( buckets[i].neigh_proc[25] != myid )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[0];
      vert_table.push_back (Vertex (num_vert++, coord));
    }

    // vertex 3
    if ( buckets[i].neigh_proc[16] != myid )
    {
      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[0];
      vert_table.push_back ( Vertex (num_vert++, coord));
    }

    // vertex 4
    if ( buckets[i].neigh_proc[14] != myid )
    {
      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[1];
      vert_table.push_back (Vertex (num_vert++, coord));
    }

    // vertex 5
    if ( buckets[i].neigh_proc[23] != myid )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[1];
      vert_table.push_back ( Vertex (num_vert++, coord));
    }

    // vertex 6
    if ( buckets[i].neigh_proc[26] != myid )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      vert_table.push_back ( Vertex(num_vert++, coord));
    }

    // vertex 7
    if ( buckets[i].neigh_proc[17] != myid )
    {
      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      vert_table.push_back (Vertex (num_vert++, coord));
    }

    if ( mincrd[0] > buckets[i].xcoord[0] )
      mincrd[0] = buckets[i].xcoord[0];

    if ( mincrd[1] > buckets[i].ycoord[0] )
      mincrd[1] = buckets[i].ycoord[0];

    if ( mincrd[2] > buckets[i].zcoord[0] )
      mincrd[2] = buckets[i].zcoord[0];

    if ( maxcrd[0] < buckets[i].xcoord[1] )
      maxcrd[0] = buckets[i].xcoord[1];

    if ( maxcrd[1] < buckets[i].ycoord[1] )
      maxcrd[1] = buckets[i].ycoord[1];

    if ( maxcrd[2] < buckets[i].zcoord[1] )
      maxcrd[2] = buckets[i].zcoord[1];

  }

  // write number of vertices to file
  fprintf (fp, "%d\n", (int) vert_table.size());

  // store vertex-indices into hash-map
  map <int, int> search_table;
  vector <Vertex>::iterator itr;
  for (itr = vert_table.begin(); itr != vert_table.end(); itr++)
  {
    key = genkey (itr->coord, mincrd, maxcrd);
    itr->key = key;
    search_table[key] = itr->index;
    fprintf (fp, "%f, %f, %f\n", itr->coord[0], itr->coord[1], itr->coord[2]);
  } 

  /* create face connectivity data */
  int face[4];
  double cube[8][3];
  for ( i = 0; i < num; i++)
  {
    // everyone writes negative side faces
    // begin -- face 0 [3, 2, 1, 0]
    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[1];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[0] = search_table[key];

    coord[0] = buckets[i].xcoord[1];
    coord[1] = buckets[i].ycoord[1];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[1] = search_table[key];

    coord[0] = buckets[i].xcoord[1];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[2] = search_table[key];

    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[3] = search_table[key];
    fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    // end face -- 0

    // begin face -- 1 [0 1 5 4]
    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[0] = search_table[key];

    coord[0] = buckets[i].xcoord[1];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[1] = search_table[key];

    coord[0] = buckets[i].xcoord[1];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[1];
    key = genkey (coord, mincrd, maxcrd);
    face[2] = search_table[key];

    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[1];
    key = genkey (coord, mincrd, maxcrd);
    face[3] = search_table[key];
    fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    // end face 1

    // begin face 2 [1 2 6 5]
    if ( buckets[i].neigh_proc[22] != myid )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[0];
      key = genkey (coord, mincrd, maxcrd);
      face[0] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[0];
      key = genkey (coord, mincrd, maxcrd);
      face[1] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[2] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[3] = search_table[key];
      fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    }
    // end face 2
 
    // begin face 3 [2, 3, 7, 6]
    if ( buckets[i].neigh_proc[16] )
    {
      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[0];
      key = genkey (coord, mincrd, maxcrd);
      face[0] = search_table[key];

      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[0];
      key = genkey (coord, mincrd, maxcrd);
      face[1] = search_table[key];

      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[2] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[3] = search_table[key];
      fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    }
    // end face 3

    // begin face 4 [3, 0, 4, 7]
    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[1];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[0] = search_table[key];

    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[0];
    key = genkey (coord, mincrd, maxcrd);
    face[1] = search_table[key];

    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[0];
    coord[2] = buckets[i].zcoord[1];
    key = genkey (coord, mincrd, maxcrd);
    face[2] = search_table[key];

    coord[0] = buckets[i].xcoord[0];
    coord[1] = buckets[i].ycoord[1];
    coord[2] = buckets[i].zcoord[1];
    key = genkey (coord, mincrd, maxcrd);
    face[3] = search_table[key];
    fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    // end face 4

    // begin face 5 [4 5 6 7]
    if ( buckets[i].neigh_proc[14] != myid )
    {
      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[0] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[0];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[1] = search_table[key];

      coord[0] = buckets[i].xcoord[1];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[2] = search_table[key];

      coord[0] = buckets[i].xcoord[0];
      coord[1] = buckets[i].ycoord[1];
      coord[2] = buckets[i].zcoord[1];
      key = genkey (coord, mincrd, maxcrd);
      face[3] = search_table[key];
      fprintf (fp, "%d, %d, %d, %d\n", face[0], face[1], face[2], face[3]);
    }
    // end face 4
  }

#ifdef MULTI_PROC
  MPI_Finalize ();
#endif

  fclose (fp);
  delete [] buckets;
  return 0;
}
