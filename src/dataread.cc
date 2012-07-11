
/*
 * =====================================================================================
 *
 *       Filename:  dataread.cc
 *
 *    Description:  
 *
 *        Created:  03/15/2010 01:49:34 PM EDT
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

#define deg2rad(A)  ((A)*(0.01745329252))

#include <iostream>
#include <string>
#include <limits>
using namespace std;

#include <hdf5.h>
#include "hdf5calls.h"

#ifdef MULTI_PROC
#  include <mpi.h>
#endif

#include <hashtab.h>
#include <bucket.h>
#include <buckstr.h>
#include <buckhead.h>
#include <particle.h>
#include <properties.h>

#include "sph_header.h"

int
Read_Data(MatProps * matprops, TimeProps * timeprops,
          PileProps * pileprops, FluxProps * fluxprops, int *format)
{

  int i;
  double rotang, intfrict, bedfrict;

  ifstream inD1 ("scale.data", ios::in);
  if (inD1.good())
  {
    inD1 >> matprops->LENGTH_SCALE;
    inD1 >> matprops->GRAVITY_SCALE;
  }
  else
  {
    matprops->LENGTH_SCALE = 1.;
    matprops->GRAVITY_SCALE = 1.;
  }

  ifstream inD2 ("simulation.data", ios::in);

  if (inD2.fail())
  {
    cerr << "ERROR: Can't find \"simulation.data\" input file." << endl;
    exit(1);
  }
  int numpiles;

  // Pile properties
  // read data for all the piles
  double temp;
  double len_scale = matprops->LENGTH_SCALE;
  inD2 >> temp;
  pileprops->pileheight = temp / len_scale;
  inD2 >> temp;
  pileprops->xCen = temp / len_scale;
  inD2 >> temp;
  pileprops->yCen = temp / len_scale;
  inD2 >> temp;
  pileprops->majorrad = temp / len_scale;
  inD2 >> temp;
  pileprops->minorrad = temp / len_scale;
  inD2 >> rotang;
  rotang = deg2rad (rotang);
  pileprops->cosrot = cos(rotang);
  pileprops->sinrot = sin(rotang);

  double time_scale = 1;
  // simulation time properties
  inD2 >> timeprops->max_time;
  inD2 >> timeprops->max_steps;
  inD2 >> timeprops->timeoutput;
  inD2 >> *(format);
  timeprops->TIME_SCALE = time_scale;
  timeprops->ndtimeoutput = timeprops->timeoutput / time_scale;
  timeprops->ndmax_time = timeprops->max_time / time_scale;

  // material properties
  inD2 >> matprops->P_CONSTANT;
  inD2 >> matprops->GAMMA;
  inD2 >> matprops->smoothing_length;
  inD2 >> intfrict;
  inD2 >> bedfrict;
  matprops->particle_mass = pow(matprops->smoothing_length, DIMENSION);
  matprops->intfrict = deg2rad (intfrict);
  matprops->bedfrict = bedfrict;
  matprops->tanintfrict = tan(matprops->intfrict);
  matprops->sinintfrict = sin(matprops->intfrict);
  inD2.close();

  // Flux source data
  ifstream influx("fluxsource.data", ios::in);

  if (influx.good())
  {
    fluxprops->have_src = true;
    influx >> fluxprops->xSrc;
    influx >> fluxprops->starttime;
    influx >> fluxprops->stoptime;
    influx >> fluxprops->tangvel;
  }
  influx.close();
  return 0;
}

int
Read_Grid(HashTable ** P_table, HashTable ** BG_mesh,
          vector < BucketHead > & partition_tab, MatProps * matprops,
          PileProps * pileprops, FluxProps * fluxprops)
{
  int No_of_Buckets;
  int BG_TABLE_SIZE = 100000;
  int P_TABLE_SIZE = 400000;
  double keyrange[2], mindom[DIMENSION], maxdom[DIMENSION];
  unsigned key[KEYLENGTH];
  Key tempkey;
  double hvars[6], min_crd[DIMENSION], max_crd[DIMENSION];
  char filename[14];
  int Down[DIMENSION] = { 0, 0, 1 };

  BucketStruct *bucket;
  Key neigh_keys[NEIGH_SIZE];
  int neigh_proc[NEIGH_SIZE];
  double elev[4];
  int i, j, k, myid;

  // infinity
  double infty;
  if (numeric_limits < double >::has_infinity)
    infty = numeric_limits < double >::infinity();
  else
    infty = HUGE_VAL;

  // Read Hash-table related constants
  myid = 0;
#ifdef MULTI_PROC
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif
  sprintf(filename, "funky%04d.h5", myid);
  hid_t fp = GH5_fopen_serial(filename, 'r');

  // Read Hash table constants
  GH5_readdata(fp, "/hashtable_constants", hvars);
  mindom[0] = hvars[0];         // min x
  maxdom[0] = hvars[1];         // max x
  mindom[1] = hvars[2];         // min y
  maxdom[1] = hvars[3];         // max y
  mindom[2] = hvars[4];         // min z
  maxdom[2] = hvars[5];         // max z

  // Create two new Hash-tables
  for (i = 0; i < DIMENSION; i++)
  {
    mindom[i] /= matprops->LENGTH_SCALE;
    maxdom[i] /= matprops->LENGTH_SCALE;
  }

  // create hash-table for back-ground mesh
  *BG_mesh = new HashTable (BG_TABLE_SIZE, 2017, mindom, maxdom);

  // get the size of BG Mesh
  hsize_t dims[2];
  int num = GH5_getsize(fp, "/Buckets", dims);

  No_of_Buckets = (int) dims[0];

  // allocate memory for buckets
  bucket = new BucketStruct[No_of_Buckets];

  // read BG Mesh data
  GH5_read_grid_data(fp, "/Buckets", bucket);

  for (i = 0; i < No_of_Buckets; i++)
  {
    // hash-table keys
    for (j = 0; j < KEYLENGTH; j++)
      key[j] = bucket[i].key[j];

    // min coordinates
    min_crd[0] = bucket[i].xcoord[0];
    min_crd[1] = bucket[i].ycoord[0];
    min_crd[2] = bucket[i].zcoord[0];

    // max coordinates
    max_crd[0] = bucket[i].xcoord[1];
    max_crd[1] = bucket[i].ycoord[1];
    max_crd[2] = bucket[i].zcoord[1];

    if (bucket[i].myproc != myid)
    {
      fprintf(stderr, "ERROR: Input data is not correct. Aborting.\n");
      fprintf(stderr, "myid = %d, data_proc = %d\n", myid, bucket[i].myproc);
#ifdef MULTI_PROC
      MPI_Abort(MPI_COMM_WORLD, myid);
#else
      exit(1);
#endif
    }

    for (j = 0; j < NEIGH_SIZE; j++)
    {
      for (k = 0; k < KEYLENGTH; k++)
        tempkey.key[k] = bucket[i].neighs[j * KEYLENGTH + k];
      neigh_keys[j] = tempkey;
      neigh_proc[j] = bucket[i].neigh_proc[j];
    }

    // Check if current bucket has flux source location
    if (fluxprops->have_src)
    {
      double xsrc = fluxprops->xSrc;

      if ((xsrc > min_crd[0]) && (xsrc < max_crd[0])
          && (bucket[i].buckettype == 2))
      {
        for (k = 0; k < KEYLENGTH; k++)
          tempkey.key[k] = bucket[i].key[k];
        fluxprops->bucketsrckey = tempkey;
        // shift fluxsrc a little bit
        double dx = matprops->smoothing_length;
        int nn = (int) round((xsrc - min_crd[0]) / dx);

        fluxprops->xSrc = min_crd[0] + nn * dx + dx / 2.;
      }
    }
    // Check if current bucket has center of any of the piles
    double xcen = pileprops->xCen;
    double ycen = pileprops->yCen;
    if ((bucket[i].buckettype == 2) &&
        (xcen > min_crd[0]) && (xcen <= max_crd[0]) &&
        (ycen > min_crd[1]) && (ycen <= max_crd[1]))
    {
      for (k = 0; k < KEYLENGTH; k++)
        tempkey.key[k] = key[k];
      pileprops->CenBucket = tempkey;
      pileprops->xCen = 0.5 * (min_crd[0] + max_crd[0]);
      pileprops->yCen = 0.5 * (min_crd[1] + max_crd[1]);
    }

    // buckettype flag
    int btflag;
    for (j = 0; j < 4; j++)
      elev[j] = 0.;

    switch (bucket[i].buckettype)
    {
    case 1:
      btflag = UNDERGROUND;
      break;
    case 2:
      btflag = MIXED;
      for (j = 0; j < 4; j++)
        elev[j] = bucket[i].elev[j];
      break;
    case 3:
      btflag = OVERGROUND;
      break;
    default:
      fprintf(stderr,"ERROR: Unknown buckettype flag.\nCheck the preoprocessor\n");
      exit(1);
    }

    // create a new bucket
    Bucket * buck = new Bucket(key, min_crd, max_crd, btflag, elev,
                               myid, neigh_proc, neigh_keys);
    (*BG_mesh)->add(key, buck);
  }

  // free memory
  delete [] bucket;

  // clear out the partition table
  partition_tab.clear();
 
  // get the size of BG Mesh
  GH5_getsize (fp, "/partition_table", dims);
  if ((dims[0] <= 0))
  {
    fprintf (stderr,"Error! unable to read partition table.\n");
    exit (1);
  }

  double center[2];
  unsigned * part_keys = new unsigned [dims[0]];
  int incr = 2 + KEYLENGTH;
  GH5_readdata (fp, "/partition_table", part_keys);
  for (i = 0; i < dims[0]; i += incr)
  {
    Bucket * buck = (Bucket *) (*BG_mesh)->lookup (part_keys + i + 2);
    if (buck->which_neigh_proc (Down) != -1)
    {
      fprintf (stderr, "ERROR: Partition table do not have correct keys.\n");
      exit (1);
    }
    partition_tab.push_back (BucketHead (part_keys + i, part_keys + i + 2));
  }

  // Create hash-table for particles 
  *P_table = new HashTable(P_TABLE_SIZE, 2017, mindom, maxdom);
  return 0;
}
