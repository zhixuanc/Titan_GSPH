/*
 * =====================================================================================
 *
 *       Filename:  particler.h
 *
 *    Description:  
 *
 *        Created:  03/22/2010 11:51:36 PM EDT
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

#ifndef PARTICLER__H
#define PARTICLER__H

#include <vector>
#include <list>
using namespace std;

#include <hashtab.h>
#include <bnd_image.h>
#include <multiproc.h>

//! calc f_coef
int  calc_f_coef(
                 int ,
                 HashTable *,
                 HashTable *,
                 MatProps *
                );

//! Calculate velocity gradients
int calc_gradients(
                    HashTable * //! Hashtable of SPH partilces
                   );

//! Apply boundary conditions
int apply_bcond (
                 int        ,  //! Process ID
                 HashTable *,  //! Hashtable of SPH partilces
                 HashTable *,  //! Background mesh
                 MatProps *,    //! Structure of material properties 
                 list <BndImage> & //! table of ghost images
                );

//! Put initial piles. All piles are Elliptical-parablpoids
void init_piles(
                HashTable *,   //! HashTable of SPH partilces
                HashTable *,   //! Hashtable containing cells of backgound mesh
                PileProps *,   //! Structure of pile properties data
                MatProps *,    //! Structure of material properties data
                int ,          //! Number of procs
                int            //! my process id
               );
//! Momentum update function. This where most of the work is done
int mom_update(
               int ,           //! my Proc ID
               HashTable *,    //! HashTable of SPH partilces
               HashTable *,    //! HashTable of Mesh elements
               MatProps *,     //! Material properties
               TimeProps *     //! struct for simulation props
              );

//! Read simulation data \f${\it i.e.}\f$ simulation time, output format, pile pros etc
int Read_Data (
               MatProps *,     //! Structure containg material mroperties
               TimeProps *,    //! Structure containg time mroperties
               PileProps *,    //! Structure containg pile information
               FluxProps *,
               int *           //! output file format
              );

//! Read background grid created by preprocessor
int Read_Grid (
               HashTable **,   //! Pointer to Hashtable for background mesh
               HashTable **,   //! Pointer to HashTable for partilces
               vector<BucketHead> &,  //! Vector of sorted partition table keys
               MatProps  *,    //! Structure containing material properties
               PileProps *,    //! Structure containing initial piles data
               FluxProps *,    //! Structure cotaining pros of flux source
               int ,           //! my process id
               int ,           //! number of total processes
               int *           //! array of flags, for cummnication with other procs
              );

//! Search neighbors for their proximity
int  search_neighs(
                   int myid   , //! My processor ID
                   HashTable *, //! HashTable of SPH partilces
                   HashTable *  //! HashTable of cells of background mesh
                  );

//! Smooth density (low-pass filter)
void smooth_density (
                     HashTable * //! HashTable of SPH partilces
                    );

//! Calculate time increment, depending upon CFL condition
double timestep(
                HashTable *,    //! HashTable of SPH partilces
                MatProps *      //! Structure of material properties
               );

//! Update particle positions and their relationship with background Mesh
int update_pos(
               int ,         //! my proc id
               HashTable *,  //! HashTable of SPH particles
               HashTable *,  //! HashTable of cells of background mesh
               FluxProps *,  //! Flux properties struct
               TimeProps *,  //! Time properties struct
               int *         //! Pointer to number of particles removed
              );

//! Write output, in specified format
void write_output(
                  int ,        //! my process id
                  int ,        //! total number of procs
                  HashTable *, //! HashTable of SPH particles
                  HashTable *, //! HashTable of Background Cells
                  vector<BucketHead> &, //! table of partitioned keys
                  TimeProps *, //! Structure time properties
                  int     //! file format, {\it i.e.} hdf5, tecplot etc
                 );

// flux
void init_fluxsrc (
                   HashTable *,
                   HashTable *,
                   MatProps  *,
                   FluxProps *
                  );

void update_fluxsrc (
                    HashTable *,
                    HashTable *,
                    MatProps  *,
                    FluxProps *,
                    TimeProps *
                   );

#endif // PARTICLER__H
