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
#include <fstream>
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
  else //if ( buck1.key[1] > buck2.key[1] )
    return false;
  
  cerr << "ERROR: Two keys match exactly, expand keylength" << endl;
  cerr << buck1.key[0] <<", " << buck1.key[1] << " == "
          << buck2.key[0] <<", " << buck2.key[1] << endl;
  exit(1);
}

int main(int argc, char *argv[])
{

  int i, j, k, l;
  int np, kk;
  double xcrd[2],ycrd[2], zcrd[2],cntr[DIMENSION],smlen;
  unsigned keylen=KEYLENGTH;
  unsigned key[KEYLENGTH], key2[2];
  double mindom[DIMENSION], maxdom[DIMENSION], polyconst[DIMENSION+1];
  string gis_db, gis_location, gis_mapname, gis_mapset;


  // hashtable constants
  double htvars[6];

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
  Get_xmax (resolution, &maxdom[0]);
  Get_ymax (resolution, &maxdom[1]);
  Get_elev_max (resolution, &maxdom[2]);

  /* expand z-span a little bit */
  // shift min-z by 1.5 times mesh-resolution
  mindom[2] -= 1.5 * del;

  // read "simulation.data" file
  ifstream inSim ("simulation.data", ios::out);
  if ( inSim.fail () )
  {
    cerr << "Failed to read \"simulation.data\" file." << endl;
    exit (1);
  }

  // increase maxdom by pileheight + 1.5 time mesh-resolution
  double pheight, xcen, ycen, zcen;
  inSim >> pheight;
  inSim >> xcen;
  inSim >> ycen;
  inSim.close ();

  // get elevation at the pile-center
  Get_elevation (resolution, xcen, ycen, &zcen);

  if ( maxdom[2] - mindom[2] < 1.5 * pheight )
    maxdom[2] = mindom[2] + pheight + 1.5 * del;
  else
    maxdom[2] += 1.5 * del;
 
  // max number of buckets along each directions
  int nx = (int) ((maxdom[0]-mindom[0])/(del)) + 1;
  int ny = (int) ((maxdom[1]-mindom[1])/(del)) + 1;
  int ncol = (int) (pheight / del + 0.5) + 2;
  int nz = (int) ((maxdom[2]-mindom[2])/(del)) + 1 + ncol;

  maxdom[0] = mindom[0] + nx*del;
  maxdom[1] = mindom[1] + ny*del;
  maxdom[2] = mindom[2] + nz*del;

  // fill hashtable variables
  htvars[0] = mindom[0];
  htvars[1] = maxdom[0];
  htvars[2] = mindom[1];
  htvars[3] = maxdom[1];
  htvars[4] = mindom[2];
  htvars[5] = maxdom[2];


  int Nbucket = nx*ny*nz;
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
  for (i = 0; i < nx; i++)
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

      for (k = 0; k < nz; k++)
      {
        zcrd[0] = mindom[2] + k*del;
        zcrd[1] = mindom[2] + (k+1)*del;

        for (l=0; l<2; l++)
        {
          bgmesh[i][j][k].xcoord[l] = xcrd[l];
          bgmesh[i][j][k].ycoord[l] = ycrd[l];
          bgmesh[i][j][k].zcoord[l] = zcrd[l];
        }

        if ((k < nmin - 1) || (k > nmax + ncol))
        {
          bgmesh[i][j][k].buckettype = 0;
          for (l=0; l<DIMENSION+1; l++)
            bgmesh[i][j][k].elev[l] = 0;
        }
        else if ( k == nmin - 1)
        {
          bgmesh[i][j][k].buckettype = 1;
          for (l=0; l<DIMENSION+1; l++)
            bgmesh[i][j][k].elev[l] = 0;
        }
        else if ( k > nmax )
        {
          bgmesh[i][j][k].buckettype = 3;
          for (l=0; l<DIMENSION+1; l++)
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

        // determine key 
        HSFC3d (normc, & keylen, key);
  
        // put key value to BG_mesh structure
        for (l=0; l<KEYLENGTH; l++)
          bgmesh[i][j][k].key[l] = key[l];

        // find key value for partition bucket
        if ( k==0 )
        {
          HSFC2d (normc, & keylen, key2);
          ColumnHead temp_head (i, j, key2);
          partition_table.push_back(temp_head);
        }
      }
    }
  }

  // order buckets according to keys
  sort(partition_table.begin(), partition_table.end());

  int bucks_per_proc = (nx*ny/np);
  int nxny = nx * ny;
  for (i=0; i < nxny; i++ )
  {
    int myid = i/bucks_per_proc;
    if ( myid >= np )  myid = np-1;
    partition_table[i].proc = myid;
    j = partition_table[i].xind;
    k = partition_table[i].yind;
    for (l = 0; l < nz; l++)
      if ( bgmesh[j][k][l].buckettype )
        bgmesh[j][k][l].myproc = myid;
      else
        bgmesh[j][k][l].myproc = -1;
  }
  
  // determine neighbors
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
      {
        int ncount=0;
        int npcount=0;
        for (int ii=i-1; ii<i+2; ii++)
          for (int jj=j-1; jj<j+2; jj++)
            for (int kk=k-1; kk<k+2; kk++)
            {
              if ( ! bgmesh[i][j][k].buckettype )
                continue;

              if ((ii>-1)&&(ii<nx) &&
                  (jj>-1)&&(jj<ny) &&
                  (kk>-1)&&(kk<nz))
              {
                bgmesh[i][j][k].neighs[ncount++] = bgmesh[ii][jj][kk].key[0];
                bgmesh[i][j][k].neighs[ncount++] = bgmesh[ii][jj][kk].key[1];
                if ((ii==i) && (jj==j) && (kk==k))
                  bgmesh[i][j][k].neigh_proc[npcount++] = -2;
                else
                  bgmesh[i][j][k].neigh_proc[npcount++] = bgmesh[ii][jj][kk].myproc;
              }
              else
              {
                bgmesh[i][j][k].neighs[ncount++]=0;
                bgmesh[i][j][k].neighs[ncount++]=0;
                bgmesh[i][j][k].neigh_proc[npcount++] = -1;
              }
            }
      }

  int kcount;
  vector <ColumnHead> :: iterator c_itr = partition_table.begin ();
  for (int iproc = 0; iproc < np; iproc++)
  {
    vector<BucketStruct> proc_bucks;
    vector <unsigned> partition_keys;
    while ( c_itr->proc == iproc )
    {
      i = c_itr->xind;
      j = c_itr->yind;
      kcount = 0;
      for (k = 0; k < nz; k++)
        if ( bgmesh[i][j][k].buckettype > 0)
        {
          if (bgmesh[i][j][k].myproc == iproc)
          {
            proc_bucks.push_back(bgmesh[i][j][k]);
            kcount++;
          }
          else
          {
            cerr << "Error: proc-ids don't match" << endl;
            exit (1);
          }
        }
        // copy partition table keys
        for (k = 0; k < KEYLENGTH; k++)
          partition_keys.push_back (c_itr->key[k]);

        // search for the first-bucket
        for (l = 0; l < nz; l++)
          if (bgmesh[i][j][l].buckettype == 1)
            break;

        // copy key of the first bucket in the column
        for (k = 0; k < KEYLENGTH; k++)
          partition_keys.push_back (bgmesh[i][j][l].key[k]);

        // advance the iterator
        c_itr++;
        if ( c_itr == partition_table.end () )
          break;
    }
    createfunky (iproc, 6, htvars, proc_bucks, partition_keys);
    proc_bucks.clear();
  }

  // Create write initial data to HDF5 file
  cout << "Total "<< Nbucket <<", "<< bucks_per_proc * nz <<" buckets per proc" << endl;
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
  HSFC3d (norm_coord, &nkey, key);
  
  // min-max keys
  if(key[0]>ma[0] || (key[0]==ma[0] && key[1]>ma[1])) {ma[0]=key[0]; ma[1]=key[1];}
  if(key[0]<mi[0] || (key[0]==mi[0] && key[1]<mi[1])) {mi[0]=key[0]; mi[1]=key[1];}

  return;
}
