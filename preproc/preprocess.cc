/*
 * =====================================================================================
 *
 *       Filename:  preprocess.cc
 *
 *    Description:  
 *
 *        Created:  03/26/2010 01:43:08 PM EDT
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
# include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>
using namespace std;

#include <unistd.h>
#include <hilbert.h>
#include <GisApi.h>

#include "buckstr.h"
#include "preprocess.h"

#define MAXKEY ~(0u)

void usage()
{
  cerr <<"Usage :"<< "./preprocess <nproc> <smooth-length> <gis database> <gis location>"
                  << " <gis mapset> <gis map>" << endl;
  exit(0);
}

void get_plane_coef(double x[], double y[], double z[], double poly[])
{
  double temp1 = z[1] - z[0];
  double temp2 = x[1] - x[0];
  double temp3 = sqrt(pow(temp1,2)+pow(temp2,2));
  poly[0] = -temp1/temp3;
  poly[1] = temp2/temp3;
  poly[2] = (temp1*x[0]-temp2*z[0])/temp3;
  return;
}


bool compare_keys ( const ColumnHead  & buck1, const ColumnHead & buck2)
{

  if      ( buck1.key[0] < buck2.key[0] )
    return true;
  else if ( buck1.key[0] > buck2.key[0] )
    return false;
  else if ( buck1.key[1] < buck2.key[1] )
    return true;
  else if ( buck1.key[1] > buck2.key[1] )
    return false;
  
   cerr << "ERROR: Two keys match exactly, expand keylength" << endl;
   exit(1);
}

int main(int argc, char *argv[])
{

  int i, j, k, l;
  int np, kk;
  double xcrd[2],ycrd[2], zcrd[2],cntr[DIMENSION],smlen;
  unsigned keylen=KEYLENGTH;
  unsigned key[KEYLENGTH], max_kc[KEYLENGTH], min_kc[KEYLENGTH];
  double mindom[DIMENSION], maxdom[DIMENSION], polyconst[DIMENSION+1];
  string gis_db, gis_location, gis_mapname, gis_mapset;


  // hashtable constants
  double htvars[8];
  double keyrange[2];

  if ( argc != 7 )
    usage();

  np = atoi(argv[1]);
  smlen = atof(argv[2]);
  gis_db = argv[3];
  gis_location = argv[4];
  gis_mapset = argv[5];
  gis_mapname = argv[6];
   
  i = Initialize_GIS_data (argv[3], argv[4], argv[5], argv[6]);
  if ( i!=0 )
  {
    cerr << "Problem with GIS data." << endl;
    exit(1);
  }

  double del=PARTICLE_DENSITY*smlen;
  double resolution = del;

  // get domain limits from GIS
  // min
  Get_xmin(resolution, &mindom[0]);
  Get_ymin(resolution, &mindom[1]);
  Get_elev_min(resolution, &mindom[2]);
  // max
  Get_xmax(resolution, &maxdom[0]);
  Get_ymax(resolution, &maxdom[1]);
  Get_elev_max(resolution, &maxdom[2]);

  // max number of buckets along each directions
  int nx = (int) ((maxdom[0]-mindom[0])/(del)) + 1;
  int ny = (int) ((maxdom[1]-mindom[1])/(del)) + 1;
  int nz = (int) ((maxdom[2]-mindom[2])/(del)) + 1;

  maxdom[0] = mindom[0] + nx*del;
  maxdom[1] = mindom[1] + ny*del;
  maxdom[2] = mindom[2] + nz*del;

  // fill hashtable variables
  htvars[0] = keyrange[0];
  htvars[1] = keyrange[1];
  htvars[2] = mindom[0];
  htvars[3] = maxdom[0];
  htvars[4] = mindom[1];
  htvars[5] = maxdom[1];
  htvars[6] = mindom[2];
  htvars[7] = maxdom[2];


  int Nbucket = nx*ny*nz;
  for (i=0; i<KEYLENGTH; i++)
  {
    min_kc[i] =MAXKEY;
    max_kc[i] =0;
  }

  vector<ColumnHead> partition_table;

  // data-structure for back-ground mesh
  // is a 2-d array of linkded lists
  // using <STL List> for the purpose
  BucketStruct ***bgmesh = new BucketStruct**[nx];
  for (i=0; i<nx; i++)
  {
    bgmesh[i] = new BucketStruct* [ny];

    for (j=0; j<ny; j++)
      bgmesh[i][j] = new BucketStruct [nz];
  }

  int n[4];
  double elev[4];
  for (i=0; i<nx; i++)
  {
    xcrd[0] = mindom[0] + i*del;
    xcrd[1] = mindom[0] + (i+1)*del;
    for (j=0; j<ny; j++)
    {
      ycrd[0] = mindom[1] + j*del;
      ycrd[1] = mindom[1] + (j+1)*del;

      Get_elevation(resolution, xcrd[0], ycrd[0], elev);
      Get_elevation(resolution, xcrd[1], ycrd[0], elev+1);
      Get_elevation(resolution, xcrd[1], ycrd[1], elev+2);
      Get_elevation(resolution, xcrd[0], ycrd[1], elev+3);

      for (k=0; k<4; k++)
        n[k] = (int) ((elev[k]-mindom[2])/del);

      // sort n1 , n2
      int nmin = min(min(n[0], n[1]), min(n[2], n[3]));
      int nmax = max(max(n[0], n[1]), max(n[2], n[3]));

      for (k=0; k<nz; k++)
      {
        zcrd[0] = mindom[2] + k*del;
        zcrd[1] = mindom[2] + (k+1)*del;

        for (l=0; l<2; l++)
        {
          bgmesh[i][j][k].xcoord[l] = xcrd[l];
          bgmesh[i][j][k].ycoord[l] = ycrd[l];
          bgmesh[i][j][k].zcoord[l] = zcrd[l];
        }

        if ( k < nmin )
        {
          bgmesh[i][j][k].buckettype = 1;
          for (l=0; l<DIMENSION+1; l++)
            bgmesh[i][j][k].elev[l] = 0;
        }
        else if ( k > nmax )
        {
          bgmesh[i][j][k].buckettype = 3;
          for (l=0; k<DIMENSION+1; l++)
            bgmesh[i][j][k].elev[l] = 0;
        }
        else
        {
          bgmesh[i][j][k].buckettype = 2;
          for (l=0; l<DIMENSION+1; l++)
            bgmesh[i][j][k].elev[l]=elev[l];
        }
        
        // generate hash-key for bucket
        double normc[DIMENSION];
        cntr[0] = (xcrd[0]+xcrd[1])*0.5;
        cntr[1] = (ycrd[0]+ycrd[1])*0.5;
        cntr[2] = (zcrd[0]+zcrd[1])*0.5;
        for ( l=0; l<DIMENSION; l++)
          normc[l]=(cntr[l]-mindom[l])/(maxdom[l]-mindom[l]);

        // determine key and update min-max keys
        determine_the_key (normc, keylen, key, max_kc, min_kc);
        for (l=0; l<KEYLENGTH; l++)
          bgmesh[i][j][k].key[l] = key[l];
        if ( k==0 )
          partition_table.push_back(ColumnHead(i, j, key));
      }
    }
  }

  // determin neighbors
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++)
      {
        int ncount=0;
        for (int kk=k-1; kk<k+2; kk++)
          for (int jj=j-1; jj<j+2; jj++)
            for (int ii=i-1; ii<i+2; ii++)
            {
              if ((ii>-1)&&(ii<nx) &&
                  (jj>-1)&&(jj<ny) &&
                  (kk>-1)&&(kk<nz))
              {
                bgmesh[i][j][k].neighs[ncount++]=bgmesh[ii][jj][kk].key[0];
                bgmesh[i][j][k].neighs[ncount++]=bgmesh[ii][jj][kk].key[1];
              }
              else
              {
                bgmesh[i][j][k].neighs[ncount++]=0;
                bgmesh[i][j][k].neighs[ncount++]=0;
              }
            } 
      }
  // order buckets according to keys
  sort(partition_table.begin(), partition_table.end(), compare_keys);

  int bucks_per_proc = (nx*ny/np);
  for (i=0; i < nx*ny; i++ )
  {
    int myid = i/bucks_per_proc;
    if ( myid >= np )  myid = np-1;
    j = partition_table[i].xind;
    k = partition_table[i].yind;
    partition_table[i].proc = myid;
    for (l=0; l<nz; l++)
      bgmesh[j][k][l].myproc=myid;
  }

  for (int iproc=0; iproc<np; iproc++)
  {
    vector<BucketStruct> proc_bucks;
    for (i=0; i<nx; i++)
      for (j=0; j<ny; j++)
        for (k=0; k<nz; k++)
          if ( bgmesh[i][j][k].myproc == iproc )
            proc_bucks.push_back(bgmesh[i][j][k]);
    createfunky (iproc, 8, htvars, &(proc_bucks));
    proc_bucks.clear();
  }

  // Create write initial data to HDF5 file
  cout << "Total "<< Nbucket <<", "<< bucks_per_proc <<" buckets per proc" << endl;
  cout << "it's ready to run"<<endl;


  // clean up 
  for (i=0; i<nx; i++)
  {
    for (j=0; j<ny; j++)
      delete [] bgmesh[i][j];

    delete [] bgmesh[i];
  }
  delete [] bgmesh;
  return 0;
}



void determine_the_key(double norm_coord[], unsigned nkey, unsigned key[],
                             unsigned ma[], unsigned mi[])
{
  // call Hilbert's Space Filling Curve
  HSFC2d (norm_coord, &nkey, key);
  
  // min-max keys
  if(key[0]>ma[0] || (key[0]==ma[0] && key[1]>ma[1])) {ma[0]=key[0]; ma[1]=key[1];}
  if(key[0]<mi[0] || (key[0]==mi[0] && key[1]<mi[1])) {mi[0]=key[0]; mi[1]=key[1];}

  return;
}
