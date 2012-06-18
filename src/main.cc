
/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  
 *
 *        Created:  01/07/2008 01:57:30 PM EST
 *         Author:  Dinesh Kumar (dkumar), dkumar@buffalo.edu
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
#include <vector>
using namespace std;

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef MULTI_PROC
#  include <mpi.h>
#  include <multiproc.h>
#endif

#include <hashtab.h>
#include <bgmesh.h>
#include <bnd_image.h>
#include <properties.h>
#include <buckhead.h>

#include "sph_header.h"
#include "particler.h"

int
main(int argc, char **argv)
{
  int i, j, ierr = 0;
  double dt;
  int added_ghosts;
  MatProps *matprops = new MatProps();
  TimeProps *timeprops = new TimeProps();
  PileProps *pileprops = new PileProps();
  FluxProps *fluxprops = new FluxProps();
  HashTable *P_table, *BG_mesh;

  vector < BucketHead > partition_table;
  vector < BndImage > Image_table;
  int format = 0;
  int adapt = 0;
  int lost = 0;
  int lostsum = 0;
  double loc_data[2], glob_data[2];
  int myid, numprocs;

#ifdef MULTI_PROC
  double start, finish;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  GMFG_new_MPI_Datatype();
  start = MPI_Wtime();
#else
  numprocs = 1;
  myid = 0;
#endif

  //Read input data
  if (Read_Data(matprops, timeprops, pileprops, fluxprops, &format) != 0)
  {
    cerr << "ERROR: Can't read input data\n";
    exit(1);
  }

  //read initial particle distribution
  if (Read_Grid (&P_table, &BG_mesh, &partition_table, matprops, pileprops,
       fluxprops) != 0)
  {
    cerr << "ERROR: Can't read Initial grid\n";
    exit(1);
  }

  // add piles
  init_piles(P_table, BG_mesh, pileprops, matprops, numprocs, myid);

  // update mesh
  update_bgmesh(P_table, BG_mesh, matprops, myid, &added_ghosts);

#ifdef MULTI_PROC
  // wait till initialization has finished
  MPI_Barrier (MPI_COMM_WORLD);

  // Initial repartition
  repartition (partition_table, P_table, BG_mesh);

  // move data
  move_data (numprocs, myid, P_table, BG_mesh);
#endif

  // search mirror imgaes of ghost particles into boundary
  search_bnd_images(myid, P_table, BG_mesh, Image_table, 1);

#ifdef MULTI_PROC
  // send reflections that belong to other procs
  send_foreign_images(myid, numprocs, P_table, BG_mesh, Image_table);
#endif

  // apply boundary conditions
  apply_bcond(myid, P_table, BG_mesh, matprops, Image_table);

  // Write inital configuration
  write_output (myid, numprocs, P_table, BG_mesh, timeprops, format);

  while (!timeprops->ifend())
  {
    // search and update neighbors
    search_neighs(myid, P_table, BG_mesh);

    // calculate time-step
    dt = timestep(P_table, matprops);

#ifdef MULTI_PROC
    // get-minimum time-step for multiproc run
    loc_data[0] = dt;
    loc_data[1] = (double) -ierr;
    MPI_Allreduce(loc_data, glob_data, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = glob_data[0];
    if (glob_data[1] != 0)
      MPI_Abort(MPI_COMM_WORLD, ierr);
#endif

    timeprops->incrtime(&dt);
    if (myid == 0)
      cout << "Time-step: " << timeprops->step << " dt=" << dt
        << " time=" << timeprops->timesec() << endl;

#ifdef MULTI_PROC
    // ghost need to be updated only before updating momentum
    move_bnd_images(myid, numprocs, P_table, BG_mesh, Image_table);
#endif

    // velocity gradients for density update
    if (calc_gradients(P_table) != 0)
      ierr = 1;

    // update momentum
    if (mom_update(myid, P_table, BG_mesh, matprops, timeprops) != 0)
    {
      cerr << "Momentum update failed on proc" << myid <<
        " at time-step : " << timeprops->step << endl;
      cerr << "Check outfile for proc " << myid << " for errors." << endl;
      ierr = 2;
    } 

#ifdef MULTI_PROC
    // update guests on all procs
    move_data(numprocs, myid, P_table, BG_mesh);
#endif

    // smooth out density oscillations (if any)
    smooth_density(P_table);

#ifdef MULTI_PROC
    // update guests on all procs
    move_data(numprocs, myid, P_table, BG_mesh);
#endif

    // calc friction coef at the solid boundary
    calc_f_coef(myid, P_table, BG_mesh, matprops);

    // update particle positions
    adapt = update_pos(myid, P_table, BG_mesh, fluxprops, timeprops, &lost);

    if ( fluxprops->have_src && timeprops->addmaterial() )
      update_fluxsrc (P_table, BG_mesh, matprops, fluxprops, timeprops);

    // update mesh
    if (adapt)
    {
      update_bgmesh(P_table, BG_mesh, matprops, myid, &added_ghosts);
      if (added_ghosts)
        search_bnd_images(myid, P_table, BG_mesh, Image_table, 0);
    }

#ifdef MULTI_PROC
    // update guests as density has changed since last update
    move_data(numprocs, myid, P_table, BG_mesh);

    /* DYNAMIC LOAD BALANCING */

    if ((numprocs > 1) && (timeprops->step % 1000 == 0))
    {
      // remove guest buckets and particles
      delete_guest_buckets (BG_mesh, P_table);

      // repartition the domain
      i = repartition (partition_table, P_table, BG_mesh);

      // send new guests
      move_data (numprocs, myid, P_table, BG_mesh);

      // search ghost refections
      search_bnd_images (myid, P_table, BG_mesh, Image_table, 1);

      // send reflections that belong to other procs
      send_foreign_images (myid, numprocs, P_table, BG_mesh, Image_table);
    }
#endif

    // apply boundary conditions
    ierr = apply_bcond (myid, P_table, BG_mesh, matprops, Image_table);

    // write output if needed
    //write_output(myid, numprocs, P_table, BG_mesh, timeprops, format);
    if (timeprops->ifoutput())
      write_output (myid, numprocs, P_table, BG_mesh, timeprops, format);

  }
#ifdef MULTI_PROC
  move_data (numprocs, myid, P_table, BG_mesh);
#endif

#ifdef MULTI_PROC
  finish = MPI_Wtime();
  MPI_Reduce (&lost, &lostsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Finalize ();
#endif
  if (myid == 0)
  {
    printf("A total of %d SPH Particles were lost.\n", lostsum);
    double walltime = finish - start;
    int hours = (int) (walltime / 3600.);
    int mins = (int) ((walltime - hours * 3600) / 60);
    double secs = walltime - (hours * 3600) - (mins * 60);

    printf ("Computation time for a %d proc run was %d:%02d:%f\n",
           numprocs, hours, mins, secs);
  }

  return 0;
}
