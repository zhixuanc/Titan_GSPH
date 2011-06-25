/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description:  
 *   14th of March 2011:
 *      Patition table is added as an after-thought. It is a STL
 *      vector of Hash-Table keys of buckets with z = zmin. Hence
 *      instead of partitioning the whole domain, we partition its
 *      projection along x-y plane.
 *
 *******************************************************************
 * $Id: repartition_BSFC.C,v 1.2 2003/11/25 22:13:04 kdalbey Exp $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include <cassert>
#include <vector>
using namespace std;

#include <hashtab.h>
#include <bucket.h>

#include "exvar.h"
#include "buckhead.h"
#include "repartition_BSFC.h"

// computes standard deviation 
double calc_std_deviation (double mean, double *load_arr, int size)
{
  int i;
  double variance = 0;
  for (i=0; i<size; i++)
  {
    double temp = load_arr[i] - mean;
    variance += temp*temp;
  }
  return sqrt(variance/(double) size);
}

int repartition(vector<BucketHead> *PartitionTable , HashTable* P_table, 
                HashTable* BG_mesh)
{
  int i, j, k;                           /* local variables */
  int num_local_objects;                 /* the number of objects this processor owns */
  double total_weight = 0;                /* the amount of work for all objects */
  double dev0, dev1;

#ifdef THREE_D
  int Up[DIMENSION] = {0, 0, 2};
  int Down[DIMENSION] = {0, 0, 1};
#else
  int Up[DIMENSION] = {0, 2};
  int Down[DIMENSION] = {0, 1};
#endif

  int myid, numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if ( numprocs < 2 )
    return 0;

  //printf("proc %d has entered the repartitioning scheme...\n",myid);
  
  /* get application data (number of objects, ids, weights, and coords */
  vector<double> weights;
  vector<BucketHead>::iterator ibuck;
  for ( ibuck=PartitionTable->begin(); ibuck != PartitionTable->end(); ibuck++)
  {
    Bucket *Curr_buck = (Bucket *) BG_mesh->lookup(ibuck->get_bucket());
    assert(Curr_buck);
    double iwght = 0.;
    do 
    {
      vector<Key> particles = Curr_buck->get_plist();
      vector<Key>::iterator p_itr;
      for (p_itr=particles.begin(); p_itr!=particles.end(); p_itr++)
      {
        Particle *p_curr = (Particle *) P_table->lookup(*p_itr);
        if ( !p_curr->is_ghost() )
          iwght += 1.;
      }
      Curr_buck = (Bucket*) BG_mesh->lookup(Curr_buck->which_neigh(Up));
    } while ((Curr_buck->which_neigh_proc(Up)) != -1 );

    weights.push_back(iwght);
    total_weight += iwght;
  }
  double *work_per_proc = new double [numprocs];
  MPI_Allgather(&total_weight, 1, MPI_DOUBLE, 
                work_per_proc, 1, MPI_DOUBLE, MPI_COMM_WORLD);

  double global_weight =0;
  for (i=0; i<numprocs; i++)
    global_weight += work_per_proc[i];

  // Allocate memeory for SFC_VERTICES array
  num_local_objects = (int ) PartitionTable->size();
  BSFC_VERTEX_PTR sfc_vertices = new BSFC_VERTEX [num_local_objects];

  // fill up the sfc_vertices array which stores 
  // all the necessary info about the load-balancing objects
  for (j=0; j < num_local_objects; j++)
  {
    sfc_vertices[j].lb_weight =(float) weights[j];
    Key buck_key = (*PartitionTable)[j].get_bucket();
    for(k=0;k<KEYLENGTH;k++)
        sfc_vertices[j].sfc_key[k] = buck_key.key[k];
  }    

  // get total number of objects for all procs
  int   num_global_objects;
  MPI_Allreduce(&num_local_objects, &num_global_objects, 1, 
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // local-size at each proc
  int  *proc_size_arr = new int [numprocs];
  for ( i=0; i < numprocs; i++ )
    proc_size_arr[i] =0;
  MPI_Allgather(&num_local_objects, 1, MPI_INT, 
                proc_size_arr, 1, MPI_INT, MPI_COMM_WORLD);

  // create array to store indices of objects in global array
  int *ind_arr = new int [numprocs+1];
  ind_arr[0] = 0;
  for ( i=1; i < numprocs+1; i++ )
  {
    ind_arr[i] =0;
    for ( j=0; j < i; j++ )
      ind_arr[i] += proc_size_arr[j];
  }
  
  // create local array for lb_weights and corresponding procs
  double *loc_wght_arr = new double [num_global_objects];
  for (i=0; i < num_global_objects; i++)
    loc_wght_arr[i] = 0.;

  for (i=ind_arr[myid], j=0; i < ind_arr[myid+1]; i++, j++ )
    loc_wght_arr[i] = sfc_vertices[j].lb_weight;

  double *gl_wght_arr = new double [num_global_objects];
  double *gl_work_arr = new double [numprocs];
  int    *gl_proc_arr = new int    [num_global_objects];

  for ( i=0; i < num_global_objects; i++ )
    gl_proc_arr[i] = 0;

  for ( i=0; i < numprocs; i++ )
    gl_work_arr[i] = 0.;

  // create global weight array
  MPI_Allreduce(loc_wght_arr, gl_wght_arr, num_global_objects,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double avg_work = global_weight/((double) numprocs);
  double my_ideal_share = avg_work*1.025;

  // calculcate initial imbalance
  dev0 = calc_std_deviation (avg_work, work_per_proc, numprocs);

  double my_work = 0;
  double residual_work = global_weight;
  double temp;

  // do successive backword sweeps
  int iter = numprocs -1;
  int start_proc = numprocs-1;
  int start_pos  = num_global_objects-1;
  while ( iter > 0 )
  {
    int curr_proc = start_proc;
    i = start_pos;
    while ( (curr_proc >= 0) && (i >= 0) )
    {
      temp = my_work + gl_wght_arr[i];
      if ( temp < my_ideal_share )
      {
        my_work = temp;
        residual_work -= gl_wght_arr[i];
        gl_proc_arr[i] = curr_proc;
        gl_work_arr[curr_proc] = my_work;
        i--;
      }
      else
      {
        my_work  = 0;
        curr_proc--;
      }
    }
    if ( residual_work > 0 )
    {
      start_proc--;
      // find new start position
      for (i=start_pos; i > 0; i--)
        if ( gl_proc_arr[i] == start_proc )
          break;

      // re-calculate the ideal_share 
      start_pos = i;
      temp = 0;
      for (i=start_pos; i >= 0; i--)
        temp += gl_wght_arr[i];
      residual_work = temp;
      my_ideal_share = (temp/(double) (start_proc+1))*1.025;
      my_work = 0;
      iter--;
    }
    else
      break;
  }

  // get imbalance 
  dev1 = calc_std_deviation (avg_work, gl_work_arr, numprocs);

  // if can't do better than current distribution return
  if (dev0 <= dev1) 
  {
    // clean up
    delete [] proc_size_arr;
    delete [] ind_arr;
    delete [] loc_wght_arr;
    delete [] gl_wght_arr;
    delete [] gl_proc_arr;
    delete [] gl_work_arr;
    delete [] work_per_proc;
    delete [] sfc_vertices;
    return 0;
  }

  for (i=ind_arr[myid], j=0; i < ind_arr[myid+1]; i++, j++ )
    sfc_vertices[j].destination_proc = gl_proc_arr[i];

  // update current proc info
  for (i=0; i < num_local_objects; i++)
  {
    Key bkey = (*PartitionTable)[i].get_bucket();
    Bucket *Curr_buck = (Bucket *) BG_mesh->lookup(bkey);
    assert(Curr_buck);
    Curr_buck->put_myprocess(sfc_vertices[i].destination_proc);
    do
    {
      Curr_buck = (Bucket *) BG_mesh->lookup(Curr_buck->which_neigh(Up));
      Curr_buck->put_myprocess(sfc_vertices[i].destination_proc);
    } while ( Curr_buck->which_neigh_proc(Up) != -1 );
  }
 
  // a function name can't be more self-explanatory
  BSFC_update_and_send_elements(myid, numprocs, P_table, BG_mesh);

  PartitionTable->clear();
  // update the repartion info
  Bucket *buck = NULL;
  HTIterator *itr = new HTIterator(BG_mesh);
  while ((buck = (Bucket *) itr->next()))
    if ( buck->which_neigh_proc(Down) == -1 )
    {
      Key bkey = buck->getKey();
      double xx = (*buck->get_mincrd()+*buck->get_maxcrd())*0.5;
      BucketHead bhead(bkey, xx);
      PartitionTable->push_back(bhead);
    }

  // sort bucket-heads
  sort (PartitionTable->begin(), PartitionTable->end());

  // clean up
  delete [] proc_size_arr;
  delete [] ind_arr;
  delete [] loc_wght_arr;
  delete [] gl_wght_arr;
  delete [] gl_proc_arr;
  delete [] gl_work_arr;
  delete [] work_per_proc;
  delete [] sfc_vertices;
  delete itr;
  return 2;
}
