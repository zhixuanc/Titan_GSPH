/*
 * =====================================================================================
 *
 *       Filename:  riemann_solve.h
 *
 *    Description:  
 *
 *        Created:  03/24/2010 01:35:04 PM EDT
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

#ifndef RIEMANN_SOLVE__H
#define RIEMANN_SOLVE__H

#include "properties.h"
#include "particle.h"

/*!
 *  reimann_solve() computes Riemann invariants at a 
 *  particle-particle interface.
 */
int riemann_solve(
    //! Riemann invariant density at intersection of two particles
    double *ustar,
    //! particle P_{i}
    Particle *pi,
    //! particle P_{j}
    Particle *pj,
    //! Material Properties
    MatProps *matprops,
    //! time-step dt
    double dt
    );

/*!
 *  setup_rp() calculates the density weighted interface for 
 *  a pair of particles and Interpolates state variables to
 *  setup a Reimann problem at the interface.
 */
void setup_rp (
    //! state variables at point i
    double *Ui,      
    //! state variables at point j
    double *Uj,
    //! states projected at the left inter-particle interface
    double *U_left,
    //! states projected at the right inter-particle interface
    double *U_right,
    //! inter-particle distance
    double distance,
    //! smoothing length at point i
    double hi,                  
    //! smoothing length at point j
    double hj, 
    //! time step
    double dt,        
    //! speed of sound at i
    double ci,
    //! speed of sound at j
    double cj
    );

extern "C" 
{

/*!
 * eigen_decomp() solves Riemann problem, setup by setup_rp(). 
 * Reimann invariants are stored on Ustar. The Reimann Soltion
 * is very much problem specific. Any change in constitutive model
 * requires this function to be modified. 
 */

int eigen_decomp (
    //! density 
    double rho,      
    //! velocity in tangetial direction
    double ux,
    //! Velocity of Sound in the medium
    double c,
    //! jumps in left and right states
    double du[],
    //! eigen-values returned by the function
    double eval[],
    //! jumps resolved along eigen-vectors
    double alph[],
    //! eigen-vectors returned by the function
    double evec[][NO_OF_EQNS]
    );
}

#endif
