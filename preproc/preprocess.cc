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
#include <algorithm>
#include <cassert>
using namespace std;

#include <unistd.h>
#include <getopt.h>
#include <hilbert.h>
#include <buckstr.h>
#include "preprocess.h"
#include "gisapi.h"

#ifdef THREE_D
#  define create_ind(I,J,NY,K,NZ) \
                     ((I)*((NY)*(NZ))) \
                     +((J)*(NZ)) + (K)
#else
#  define create_ind(I,J,NY) ((I)*(NY))+(J)
#endif

#define MAXKEY ~(0u)

void usage()
{
  cerr <<"Usage :"<< " preprocess [options] "<< endl;
  cerr <<" Options:" <<endl;
  cerr <<"  -p, --nproc= <no of processes>" <<endl;
  cerr <<"  -s, --smlen= <smoothing length>" <<endl; 
  cerr <<"  -m, --gismap= <Gis Raster Map>" <<endl; 
  exit(0);
}

void get_line_coef(double x[], double z[], double poly[])
{
  double temp1 = z[1] - z[0];
  double temp2 = x[1] - x[0];
  double temp3 = sqrt(pow(temp1,2)+pow(temp2,2));
  poly[0] = -temp1/temp3;
  poly[1] = temp2/temp3;
  poly[2] = (temp1*x[0]-temp2*z[0])/temp3;
  return;
}


bool compare_keys ( const BuckHead  & buck1, const BuckHead & buck2)
{

  if      ( buck1.xcrd < buck2.xcrd )
    return true;
  else 
    return false;
}

int main(int argc, char *argv[])
{

  int i, j, k, l;
  int np, kk;
  double xcrd[2],zcrd[2],cntr[2],smlen;
  unsigned keylen=2;
  unsigned key[2], max_kc[2], min_kc[2], max_kp[2], min_kp[2];
  double mindom[2], maxdom[2], polyconst[3];
  char *gispath;
  int  polytype;


  // hashtable constants
  double htvars[6];
  double keyrange[2];

  if ( argc < 2 )
    usage();

  // Invokes ctor `GetOpt (int argc, char **argv, 
  //                       char *optstring);'
  int c;
  int long_opt_index = 0;
  int longval;
  char *myarg;
  struct option long_options[] =
  { // long options array CAUTION: case-senstive //
    { "gismap", 1, &longval, 'g' },
    { "nproc", 1, &longval, 'p' },
    { "smlen", 1, &longval, 's' },
    { 0, 0, 0, 0 },      /* terminating zoros */	
  };

    while ((c = getopt_long(argc, argv, "g::p:s:h", long_options, &long_opt_index)) != -1) 
    {
      switch (c) 
      {
        case 'g':   /* long_opt_index does not make sense for these */
          gispath = optarg;
          break;
        case 'p':
          np = atoi(optarg);
          break;
        case 's':
          smlen = atof(optarg);
          break;
        case 0:     /* this is returned for long options with option[i].flag set (not NULL). */
                    /* long_opt_index is now relevant */
          switch (longval) 
          {
            case 'g':
              gispath = optarg;
              break;
            case 'p':
              np = atoi(optarg);
              break;
            case 's':
              myarg = optarg;
              smlen = atof(optarg);
              break;
          }
          break;
        case 'h': 
          usage ();
          break;
        default:
          printf("You, lamah!\n");
       }
    }


  // keyrange for hashtable
  keyrange[1]= pow(2.0,32) + 1;
  keyrange[0]= keyrange[1]/np;

  double del=PARTICLE_DENSITY*smlen;
  Get_xmin(&mindom[0]);
  Get_ymin(&mindom[1]);
  Get_xmax(&maxdom[0]);
  Get_ymax(&maxdom[1]);

  // find out ghost particles
  // x_edge = xCen + Radus * cos(a)
  double xedge = xCen + Radius/sqrt(1 + pow(Slope,2));
  int    nedge = (int) ((xedge - mindom[0])/del);
  int nx = (int) ((maxdom[0]-mindom[0])/(del)) + 1;
  int ny = (int) ((maxdom[1]-mindom[1])/(del)) + 1;

  mindom[0] = xedge - (nedge*del);
  maxdom[0] = mindom[0] + nx*del;
  maxdom[1] = mindom[1] + ny*del;

  // fill hashtable variables
  htvars[0] = keyrange[0];
  htvars[1] = keyrange[1];
  htvars[2] = mindom[0];
  htvars[3] = maxdom[0];
  htvars[4] = mindom[1];
  htvars[5] = maxdom[1];


  int Nbucket = nx*ny;
  // Initialize min-max keys
  for (i=0; i<2; i++)
  {
    min_kc[i] = min_kp[i] = MAXKEY;
    max_kc[i] = max_kp[i] = 0;
  }

  vector<BuckHead> repart_keys;
  BucketStruct *bgmesh = new BucketStruct [Nbucket];
  int **neigh_index    = new int * [Nbucket];
  int **repart_inds    = new int * [nx];
  for ( i=0; i < Nbucket; i++ )
    neigh_index[i] = new int [NEIGH_SIZE];
  for ( i = 0; i < nx; i++ )
    repart_inds[i] = new int [ny];


  int ind=0;
  for (i=0; i<nx; i++)
  {
    int n1, n2;
    int pt1=0, pt2=0;
    double zcen;
    xcrd[0] = mindom[0] + i*del;
    xcrd[1] = mindom[0] + (i+1)*del;
    zcrd[0] = GIS_get_elevation(xcrd[0], &pt1);
    zcrd[1] = GIS_get_elevation(xcrd[1], &pt2);
    n1 = (int) ((zcrd[0] - mindom[1])/del);
    n2 = (int) ((zcrd[1] - mindom[1])/del);

    polytype = (pt1 | pt2);
    switch (polytype)
    {
      case LINE:
        get_line_coef (xcrd, zcrd, polyconst);
        break;
      case CIRCLE:
        zcen  = Slope*xCen + Intcpt;
        polyconst[0] = xCen;
        polyconst[1] = zcen;
        polyconst[2] = Radius;
        break;
      case HYBRID:
        if ( pt1 == LINE )
        {
          zcrd[1] = Slope*xcrd[1] + Intcpt;
          n2 = (int) ((zcrd[1] - mindom[1])/del);
        }
        else
        {
          zcrd[0] = Slope*xcrd[0] + Intcpt;
          n1 = (int) ((zcrd[0] - mindom[1])/del);
        }
        get_line_coef (xcrd, zcrd, polyconst);
        polytype = LINE;
        break;
      default:
        fprintf(stderr,"Something gone wrong in Get_elevation\n");
        exit(1);
    }
    for (j=0; j<ny; j++)
    {
      ind = create_ind(i,j,ny);
      if ( (ind>-1) && (ind<Nbucket) )
      {
        zcrd[0] = mindom[1] + j*del;
        zcrd[1] = mindom[1] + (j+1)*del;

        for (l=0; l<2; l++)
        {
          bgmesh[ind].xcoord[l] = xcrd[l];
          bgmesh[ind].zcoord[l] = zcrd[l];
        }

        repart_inds[i][j]=ind;
        // find neighbors
        kk=0;
        int neigh_ind;
        for (int ii=i-1; ii<=i+1; ii++)
          for (int jj=j-1; jj<=j+1; jj++)
          {
            neigh_ind = create_ind(ii, jj, ny);
            if ((ii < 0)||(ii >= nx))
              neigh_ind = -1;
            if ((jj < 0)||(jj >= ny))
              neigh_ind = -1;
            if ((ii == i)&&(jj == j))
              neigh_ind = -2;
            neigh_index[ind][kk++]=neigh_ind;
          }

        // sort n1 , n2
        int nmin = min(n1, n2);
        int nmax = max(n1, n2);

        if ( j < nmin )
        {
          bgmesh[ind].boundary = 0;
          bgmesh[ind].buckettype = 1;
          for (k=0; k<3; k++)
            bgmesh[ind].poly[k] = 0;
        }
        else if ( j > nmax )
        {
          bgmesh[ind].boundary = 0;
          bgmesh[ind].buckettype = 3;
          for (k=0; k<3; k++)
            bgmesh[ind].poly[k] = 0;
        }
        else
        {
          bgmesh[ind].boundary = 1;
          bgmesh[ind].bndtype  = polytype;
          bgmesh[ind].buckettype = 2;
          for (k=0; k<3; k++)
            bgmesh[ind].poly[k] = polyconst[k];
        }
        
        // generate hash-key for bucket
        double normc[2];
        cntr[0] = (xcrd[0]+xcrd[1])/2;
        cntr[1] = (zcrd[0]+zcrd[1])/2;
        for ( k=0; k<DIMENSION; k++)
          normc[k]=(cntr[k]-mindom[k])/(maxdom[k]-mindom[k]);

        // determine key and update min-max keys
        determine_the_key (normc, keylen, key, max_kc, min_kc);
        for (k=0; k<KEYLENGTH; k++)
          bgmesh[ind].key[k] = key[k];
      
        if ( j == 0 )
        {
          BuckHead tmp(i, key, 0.5*(xcrd[0]+xcrd[1]));
          repart_keys.push_back(tmp);
        }
      }
    }
  }

  // order buckets according to keys
  sort(repart_keys.begin(), repart_keys.end(), compare_keys);

  int bucks_per_proc = (nx/np);
  for (i=0; i < nx; i++ )
  {
    int myid = i/bucks_per_proc;
    if ( myid >= np )  myid = np-1;
    k = repart_keys[i].index;
    repart_keys[i].proc = myid;
    for ( j=0; j < ny; j++ )
    {
      bgmesh[repart_inds[k][j]].myproc = myid;
    }
  }

  for ( i=0; i < Nbucket; i++ )
    for (j =0; j < NEIGH_SIZE; j++ )
    {
      kk = neigh_index[i][j];
      if ( kk > -1 )
      {
        bgmesh[i].neigh_proc[j] = bgmesh[kk].myproc;
        for ( k=0; k < KEYLENGTH; k++ )
          bgmesh[i].neighs[j*KEYLENGTH+k] = bgmesh[kk].key[k];
      }
      else
        bgmesh[i].neigh_proc[j] = kk;
    }

  // Create write initial data to HDF5 file
  createfunky (np, 6, htvars, Nbucket, bgmesh);
  cout << "Total "<< Nbucket <<", "<< bucks_per_proc <<" buckets per proc" << endl;
  cout << "it's ready to run"<<endl;


  // clean up 
  delete [] bgmesh;
  for (i=0; i < Nbucket; i++)
    delete [] neigh_index[i];
  for (i=0; i < nx; i++)
    delete [] repart_inds[i];
  delete [] neigh_index;
  delete [] repart_inds;

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
