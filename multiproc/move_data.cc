
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
 * Description: most of this function is taken from move_data of
 *              TITAN, with functionality for migrating particles
 *              along with elements (called Buckets in SPH code)
 *
 *******************************************************************
 * $Id: move_data.C,v 1.1.1.1 2003/08/13 19:26:11 sorokine Exp $ 
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <cassert>
#include <vector>
using namespace std;

#include <constants.h>
#include <hashtab.h>
#include <bucket.h>
#include <particle.h>
#include <bnd_image.h>

#include "exvar.h"
#include "multiproc.h"
#include "pack_data.h"
#include "repartition_BSFC.h"

void
move_data (int nump, int myid, HashTable * P_table, HashTable * BG_mesh)
{
  int i, j, ierr;
  const int NEW = 1;
  const int OLD = -1;

  // if number of procs == 1 , don't waste time here
  if (nump < 2)
    return;

  // __DEBUG_LOG_FILE__
  // char logfile[15];
  // sprintf(logfile, "debug%03d.log", myid);
  // FILE *fp = fopen(logfile, "w");
  // __DEBUG_LOG_FILE__

  // send_info array
  int *check_proc = new int[nump];
  int *send_info = new int[2 * nump];

  for (i = 0; i < 2 * nump; i++)
    send_info[i] = 0;

  /* count how many buckets we should send and receive from other procs */
  HTIterator *itr = new HTIterator (BG_mesh);
  Bucket *buck;

  while ((buck = (Bucket *) itr->next ()))
    if (buck->is_active () && !buck->is_guest ())
    {
      const int *neigh_proc = buck->get_neigh_proc ();

      vector < Key > plist = buck->get_plist ();

      for (i = 0; i < nump; i++)
        check_proc[i] = 0;
      // find out number of buckets and particles to send-recv
      for (i = 0; i < NEIGH_SIZE; i++)
        if ((neigh_proc[i] > -1) &&
            (neigh_proc[i] != myid) && (check_proc[neigh_proc[i]] == 0))
        {
          check_proc[neigh_proc[i]] = 1;
          send_info[2 * neigh_proc[i]]++;
          send_info[2 * neigh_proc[i] + 1] += plist.size ();
        }
    }
  send_info[2 * myid] = 0;      // don't need to send info to myself
  send_info[2 * myid + 1] = 0;  // ditto for particles

  /* 
   * A proc only receives from the proc it sends to,
   * but the recv_size may (most likely will) differ from
   * send_size. Hence we have to make extra MPI calls to get
   * recv_sizes.
   */
  int *recv_info = new int[2 * nump];

  for (i = 0; i < 2 * nump; i++)
    recv_info[i] = 0;
  int tag1 = 3210;
  MPI_Request *req1 = new MPI_Request[2 * nump];

  for (i = 0; i < nump; i++)
    if (send_info[2 * i] > 0)
    {
      // post non-blocking receives
      MPI_Irecv ((recv_info + 2 * i), 2, MPI_INT, i, tag1, MPI_COMM_WORLD,
                 (req1 + i));
      // post all the sends
      MPI_Isend ((send_info + 2 * i), 2, MPI_INT, i, tag1, MPI_COMM_WORLD,
                 (req1 + nump + i));
    }

  MPI_Status status;

  // wait for all the sends to finish
  for (i = 0; i < nump; i++)
    if (send_info[2 * i] > 0)
      ierr = MPI_Wait ((req1 + nump + i), &status);

  // Wait for all the recvs to finish 
  for (i = 0; i < nump; i++)
    if (send_info[2 * i] > 0)
      ierr = MPI_Wait ((req1 + i), &status);

  // clean up allocated memory
  delete[]req1;

  /* Mark particles in ghost-buckets "old" (if present) */
  itr->reset ();
  while ((buck = (Bucket *) itr->next ()))
    if (buck->is_guest ())
    {
      vector < Key > plist = buck->get_plist ();
      vector < Key >::iterator p_itr;
      for (p_itr = plist.begin (); p_itr != plist.end (); p_itr++)
      {
        Particle *p_old = (Particle *) P_table->lookup (*p_itr);

        p_old->put_new_old (OLD);
      }
    }

  /* post receives */
  // size of data to be received
  int recv_count[2] = { 0, 0 };
  for (i = 0; i < nump; i++)
  {
    recv_count[0] += recv_info[2 * i];
    recv_count[1] += recv_info[2 * i + 1];
  }

  // allocate space for data to be received
  MPI_Request *recv_req = new MPI_Request[2 * nump];
  BucketPack *recv_buckets = new BucketPack[recv_count[0]];
  ParticlePack *recv_particles = new ParticlePack[recv_count[1]];
  int counter_recv[2] = { 0, 0 };
  int buck_tag = 22674, part_tag = 32491;       // random tags

  for (i = 0; i < nump; i++)
  {
    if (recv_info[2 * i] != 0)
    {
      j =
        MPI_Irecv ((recv_buckets + counter_recv[0]), recv_info[2 * i],
                   BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                   (recv_req + 2 * i));
      counter_recv[0] += recv_info[2 * i];
    }
    if (recv_info[2 * i + 1] != 0)
    {
      j =
        MPI_Irecv ((recv_particles + counter_recv[1]), recv_info[2 * i + 1],
                   PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                   (recv_req + 2 * i + 1));
      counter_recv[1] += recv_info[2 * i + 1];
    }
  }                             /* done with receives */

  /* put (GHOST) elements to be moved in the proper arrays */
  // size of data to be sent
  int send_count[2] = { 0, 0 };
  for (i = 0; i < nump; i++)
  {
    send_count[0] += send_info[2 * i];
    send_count[1] += send_info[2 * i + 1];
  }
  // Allocate space for send data
  BucketPack *send_buckets = new BucketPack[send_count[0]];
  ParticlePack *send_particles = new ParticlePack[send_count[1]];
  int *counter_send_proc = new int[2 * nump];

  counter_send_proc[0] = 0;
  counter_send_proc[1] = 0;
  for (i = 1; i < nump; i++)
  {
    counter_send_proc[2 * i] =
      counter_send_proc[2 * (i - 1)] + send_info[2 * (i - 1)];
    counter_send_proc[2 * i + 1] =
      counter_send_proc[2 * (i - 1) + 1] + send_info[2 * (i - 1) + 1];
  }

  /* Pack buckets and particles that are being sent over to other procs */
  itr->reset ();
  while ((buck = (Bucket *) itr->next ()))
    if (buck->is_active () && !buck->is_guest ())
    {
      const int *neigh_proc = buck->get_neigh_proc ();

      // set check proc to zero
      for (i = 0; i < nump; i++)
        check_proc[i] = 0;

      // pack data to be sent over to neighs
      for (i = 0; i < NEIGH_SIZE; i++)
        if ((neigh_proc[i] > -1) &&
            (neigh_proc[i] != myid) && (check_proc[neigh_proc[i]] == 0))
        {
          check_proc[neigh_proc[i]] = 1;
          pack_bucket ((send_buckets + counter_send_proc[2 * neigh_proc[i]]),
                       buck, myid);
          counter_send_proc[2 * neigh_proc[i]]++;

          // pack particles
          vector < Key > plist = buck->get_plist ();
          if (plist.size () > 0)
          {
            vector < Key >::iterator ip;
            for (ip = plist.begin (); ip != plist.end (); ip++)
            {
              Particle *psend = (Particle *) P_table->lookup (*ip);

              assert (psend);
              pack_particles (psend,
                              (send_particles +
                               counter_send_proc[2 * neigh_proc[i] + 1]));
              counter_send_proc[2 * neigh_proc[i] + 1]++;
            }
          }
        }
    }

  //first nump is for buckets, 2nd nump is for particles
  MPI_Request *send_req = new MPI_Request[2 * nump];
  int counter[2] = { 0, 0 };
  for (i = 0; i < nump; i++)
  {
    if (send_info[2 * i] != 0)
    {
      j = MPI_Isend ((send_buckets + counter[0]), send_info[2 * i],
                     BUCKET_TYPE, i, buck_tag, MPI_COMM_WORLD,
                     (send_req + 2 * i));
      counter[0] += send_info[2 * i];
    }
    if (send_info[2 * i + 1] != 0)
    {
      j = MPI_Isend ((send_particles + counter[1]), send_info[2 * i + 1],
                     PARTICLE_TYPE, i, part_tag, MPI_COMM_WORLD,
                     (send_req + 2 * i + 1));
      counter[1] += send_info[2 * i + 1];
    }
  }

  int add_counter = 0, update_counter = 0;
  int count = 0;

  // unpack particles after receiving them
  for (i = 0; i < nump; i++)
    if (recv_info[2 * i + 1] != 0)
    {
      j = MPI_Wait ((recv_req + 2 * i + 1), &status);
      for (j = 0; j < recv_info[2 * i + 1]; j++)
      {
        Particle *pcurr =
          (Particle *) P_table->lookup ((recv_particles + count)->key);
        // if the particle doesn't exist on this proc, create a new one
        if (!pcurr)
        {
          Particle *newp = new Particle ();

          unpack_particle ((recv_particles + count), newp);
          newp->put_new_old (NEW);
          newp->put_guest_flag (true);
          P_table->add (newp->getKey (), newp);
        }
        // if it alread exists, copy the variables
        else
        {
          unpack_particle ((recv_particles + count), pcurr);
          pcurr->put_new_old (NEW);
          pcurr->put_guest_flag (true);
        }
        count++;
      }
    }
  // unpack buckets 
  // All the particle-packets from proc (i) should arrive before we start
  // unpacking buckets-packets from proc (i), otherwise there will an extra 
  // overhead of deleting and re-creating particles that existed already
  count = 0;
  for (i = 0; i < nump; i++)
    if (recv_info[2 * i] != 0)
    {
      ierr = MPI_Wait ((recv_req + 2 * i), &status);
      for (j = 0; j < recv_info[2 * i]; j++)
      {
        Bucket *bcurr =
          (Bucket *) (BG_mesh->lookup ((recv_buckets + count)->key));
        if (!bcurr)
        {
          // this bucket doesn't exist on this proc
          Bucket *new_buck = new Bucket ();

          unpack_bucket ((recv_buckets + count), new_buck,
                         (recv_buckets + count)->myprocess);
          BG_mesh->add (new_buck->getKey (), new_buck);
          new_buck->put_guest_flag (true);
          add_counter++;
        }
        else
        {
          // bucket is present on this proc 
          // delete old particles 
          vector < Key > plist = bcurr->get_plist ();
          vector < Key >::iterator ip;
          for (ip = plist.begin (); ip != plist.end (); ip++)
          {
            Particle *p_old = (Particle *) P_table->lookup (*ip);

            if (p_old->get_new_old () == OLD)
            {
              P_table->remove (*ip);
              delete p_old;
            }
          }
          unpack_bucket ((recv_buckets + count), bcurr,
                         (recv_buckets + count)->myprocess);
          update_counter++;
        }
        count++;
      }
    }

  // __DEBUG_LOG_FILE__
  // fclose(fp);
  // __DEBUG_LOG_FILE__

  // make sure all the sends are completed
  for (i = 0; i < 2 * nump; i++)
    if (send_info[i] != 0)
      ierr = MPI_Wait ((send_req + i), &status);

  // clean up
  delete[]check_proc;
  delete[]recv_buckets;
  delete[]recv_particles;
  delete[]recv_info;
  delete[]recv_req;

  delete[]send_req;
  delete[]send_buckets;
  delete[]send_particles;
  delete[]send_info;
  delete[]counter_send_proc;
  delete itr;

  return;
}

/* delete the ghost elements that were put in the element hashtable */
void
delete_guest_buckets (HashTable * BG_mesh, HashTable * P_table)
{
  int delete_counter = 0;
  HTIterator *itr = new HTIterator (BG_mesh);
  Bucket *buck;

  while ((buck = (Bucket *) itr->next ()))
    if (buck->is_guest ())
    {
      vector < Key > plist = buck->get_plist ();
      vector < Key >::iterator ip;
      // delete particles in the bucket
      for (ip = plist.begin (); ip != plist.end (); ip++)
      {
        Particle *pdel = (Particle *) P_table->lookup (*ip);

        P_table->remove (*ip);
        delete pdel;
      }

      // now remove the bucket
      BG_mesh->remove (buck->getKey ());
      delete buck;

      delete_counter++;
    }

  delete itr;

  return;
}
