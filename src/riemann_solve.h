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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#define eigen_decomp eigen1d

/*!
 *  reimann_solve() computes Riemann invariants at a 
 *  particle-particle interface.
 */
double riemann_solve(
       //! particle P_{i}
       Particle *pi,
       //! particle P_{j}
       Particle *pj,
       //! Material Properties
       MatProps *matprops,
       //! time-step dt
       double dt);

extern "C" 
{

/*!
 * eigen () solves Riemann problem, setup by setup_rp(). 
 * Reimann invariants are stored on Ustar. The Reimann Soltion
 * is very much problem specific. Any change in constitutive model
 * requires this function to be modified. 
 */
int eigen2d (
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
    double evec[][3]);

int eigen3d (
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
    double evec[][4]);

}// extern "C"

#endif // RIEMANN_SOLVE__H
