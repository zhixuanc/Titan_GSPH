
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
#include <vector>
#include <cassert>
using namespace std;

#include <mpi.h>

#include <hashtab.h>
#include <particle.h>
#include <bucket.h>
#include <bnd_image.h>
#include <exvar.h>

int
move_bnd_images (int myid, int nump, HashTable * P_table,
                 HashTable * BG_mesh, vector < BndImage > & Image_table)
{

  if (nump < 2)
    return 0;

  int i, j, k;
  double uvec[NO_OF_EQNS];

  int img_tag = 3152;
  int *send_info = new int[nump];
  int *recv_info = new int[nump];
  int *recv_info2 = new int[nump];

  for (i = 0; i < nump; i++)
  {
    send_info[i] = 0;
    recv_info[i] = 0;
    recv_info2[i] = 0;
  }

  /* if the bucket that contains reflected image belongs to 
   * foreing process, current proc needs to receive data from
   * that proc
   */
  
  vector < BndImage >::iterator i_img;
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)
      recv_info2[i_img->buckproc]++;

  /* calulate the amount of data that needs to be sent */
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->partproc != myid)
      send_info[i_img->partproc]++;

  // communicate send_info to update recv_info
  MPI_Request *req = new MPI_Request[2 * nump];
  int tag0 = 5132;

  /* start receiving data from foreign procs */
  for (i = 0; i < nump; i++)
    if (recv_info2[i] > 0)
      MPI_Irecv ((recv_info + i), 1, MPI_INT, i, tag0, MPI_COMM_WORLD,
                 (req + i));

  /* start sending data to foreign procs */
  for (i = 0; i < nump; i++)
    if (send_info[i] > 0)
      MPI_Isend ((send_info + i), 1, MPI_INT, i, tag0, MPI_COMM_WORLD,
                 (req + nump + i));

  MPI_Status status;

  // wait for sends to finish
  for (i = 0; i < nump; i++)
    if (recv_info2[i] > 0)
      k = MPI_Wait ((req + nump + i), &status);

  // wait for recvs to finish
  for (i = 0; i < nump; i++)
    if (send_info[i] > 0)
      k = MPI_Wait ((req + i), &status);

  delete[]req;

  // check data integrity
  for (i = 0;  i < nump; i++)
    if ( recv_info[i] != recv_info2[i] )
    {
      fprintf(stderr,"data sizes don't match %s : %d\n", __FILE__, __LINE__);
      exit (1);
    }

  //  allocate recv_buffer
  BndImage **recv_buf = new BndImage *[nump];

  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
      recv_buf[j] = new BndImage[recv_info[j]];

  // post reveive calls
  MPI_Request *img_recv_req = new MPI_Request[nump];
  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
      MPI_Irecv (recv_buf[j], recv_info[j], BND_IMAGE_TYPE, j,
                 img_tag, MPI_COMM_WORLD, (img_recv_req + j));

  // allocate memory
  int *img_counter = new int[nump];
  BndImage **send_buf = new BndImage *[nump];

  for (j = 0; j < nump; j++)
  {
    if (send_info[j] > 0)
      send_buf[j] = new BndImage[send_info[j]];
    img_counter[j] = 0;
  }

  // pack Images
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->partproc != myid)
    {
      j = i_img->partproc;
      send_buf[j][img_counter[j]] = *i_img;
      img_counter[j]++;
    }

  // send buffer
  MPI_Request *img_send_req = new MPI_Request[nump];

  for (j = 0; j < nump; j++)
    if (send_info[j] > 0)
      MPI_Isend (send_buf[j], send_info[j], BND_IMAGE_TYPE, j,
                 img_tag, MPI_COMM_WORLD, (img_send_req + j));

  // Wait and update ghost particles
  for (j = 0; j < nump; j++)
    if (recv_info[j] > 0)
    {
      MPI_Wait ((img_recv_req + j), &status);
      for (k = 0; k < recv_info[j]; k++)
      {
        Particle *pghost =
          (Particle *) P_table->lookup (recv_buf[j][k].ghost_key);
        assert (pghost);
        for (int i2 = 0; i2 < NO_OF_EQNS; i2++)
          uvec[i2] = recv_buf[j][k].state_vars[i2];
        pghost->put_state_vars (uvec);
      }
    }

  // Wait for sends to finish
  for (j = 0; j < nump; j++)
    if (send_info[j] > 0)
      MPI_Wait ((img_send_req + j), &status);

  delete[]img_send_req;
  delete[]img_recv_req;
  delete[]img_counter;
  for (j = 0; j < nump; j++)
  {
    if (send_info[j] > 0)
      delete[]send_buf[j];
    if (recv_info[j] > 0)
      delete[]recv_buf[j];
  }
  delete[]send_buf;
  delete[]recv_buf;
  delete[]send_info;
  delete[]recv_info;

  return 0;
}

// send reflctions that belong to neighboring partitions
void
send_foreign_images (int myid, int numprocs, HashTable * P_table,
                     HashTable * BG_mesh, vector < BndImage > & Image_table)
{

  if (numprocs < 2)
    return;

  int i, j;
  int dir[DIMENSION];
  double coord[DIMENSION], tmpc[DIMENSION];
  double refc[DIMENSION], intsct[DIMENSION];
  unsigned ghst_key[KEYLENGTH], buck_key[KEYLENGTH];

  int *send_info = new int[numprocs];

  for (i = 0; i < numprocs; i++)
    send_info[i] = 0;

  // count number of BoundaryImages to send to other procs
  vector < BndImage >::iterator i_img;
  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)
      send_info[i_img->buckproc]++;

  // allocate send buffer
  BndImage **sendbuf = new BndImage *[numprocs];

  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      sendbuf[i] = new BndImage[send_info[i]];
  int tag = 929;

  // get receive information
  int *recv_info = new int[numprocs];

  for (i = 0; i < numprocs; i++)
    recv_info[i] = 0;

  // gather recv_info
  MPI_Alltoall (send_info, 1, MPI_INT, recv_info, 1, MPI_INT, MPI_COMM_WORLD);

  // allocate receive buffer
  BndImage **recvbuf = new BndImage *[numprocs];

  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      recvbuf[i] = new BndImage[recv_info[i]];

  // post receive 
  MPI_Request *recv_req = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      MPI_Irecv (recvbuf[i], recv_info[i], BND_IMAGE_TYPE,
                 i, tag, MPI_COMM_WORLD, &recv_req[i]);

  // pack and send images
  int *count = new int[numprocs];

  for (i = 0; i < numprocs; i++)
    count[i] = 0;

  for (i_img = Image_table.begin (); i_img != Image_table.end (); i_img++)
    if (i_img->buckproc != myid)
    {
      j = i_img->buckproc;
      sendbuf[j][count[j]] = *i_img;
      count[j]++;
    }
  // post sends
  MPI_Request *send_req = new MPI_Request[numprocs];

  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      MPI_Isend (sendbuf[i], send_info[i], BND_IMAGE_TYPE,
                 i, tag, MPI_COMM_WORLD, &send_req[i]);

  // add to local Image_table after receives
  MPI_Status status;

  for (i = 0; i < numprocs; i++)
    if (recv_info[i] != 0)
    {
      j = MPI_Wait (&recv_req[i], &status);
      for (j = 0; j < recv_info[i]; j++)
        Image_table.push_back (recvbuf[i][j]);
    }

  // clean up
  for (i = 0; i < numprocs; i++)
    if (recv_info[i] > 0)
      delete[]recvbuf[i];
  delete[]recvbuf;
  delete[]recv_req;
  delete[]recv_info;

  // Wait for sends to finish
  for (i = 0; i < numprocs; i++)
    if (send_info[i] != 0)
      j = MPI_Wait (&send_req[i], &status);

  for (i = 0; i < numprocs; i++)
    if (send_info[i] > 0)
      delete[]sendbuf[i];
  delete[]sendbuf;
  delete[]send_info;
  delete[]send_req;
  delete[]count;
}
