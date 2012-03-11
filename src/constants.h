/*
 * =====================================================================================
 *
 *       Filename:  constants.h
 *
 *    Description:  
 *
 *        Created:  03/10/2009 05:03:40 PM EDT
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

#ifndef CONSTANT__H
#define CONSTANT__H

const int NUM_NODES=4;
const int DIMENSION=3;
const int KEYLENGTH=2;
const int NO_OF_EQNS=4;
const int DIMSQRD=9;
const int NOEQxDIM=12;
const int NEIGH_SIZE=27;
const int MAX_PARTICLES_PER_BUCKET=1000;

// Directions
const int XDIR=0;
const int YDIR=1;
const int ZDIR=2;

// Boundary Conditions
const int NONE    =0;
const int DIRCHLET=1;
const int NEUMANN =2;

// Ghost particles
const int NUM_GHOST_ROWS=5;
// Number of particles per cell per dimension
const int PARTICLE_DENSITY=5;
const int PARTICLE_DENSQRD=25;


// Bucket TYPES
const int UNDERGROUND = 0xA;
const int MIXED       = 0xB;
const int OVERGROUND  = 0xC;
const float  LOAD_BALANCE_TOLERANCE = 1.001;
const double PI = 3.14159265358979;
const double TINY = 1.0E-08;

#endif // CONSTANT__H
