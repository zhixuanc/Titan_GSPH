
/*
 * =====================================================================================
 *
 *       Filename:  sph_header.h
 *
 *    Description:  
 *
 *        Created:  03/10/2009 05:00:42 PM EDT
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
#ifndef SPH_HEADER__H
#  define SPH_HEADER__H

#  include <cmath>
#  include <ctime>
#  include <cstdlib>
#  include <iostream>
using namespace std;

#  include <hashtab.h>
#  include <properties.h>

#  include "constants.h"

/*
 *  Walltime
 */
inline double
walltime()
{
  return (double) clock() / (double) CLOCKS_PER_SEC;
}

/*
 *  sgn
 */
inline double
sgn(double x)
{
  if (x < 0)
    return -1.0;
  else
    return 1.0;
}

/*! 
 * in_support() checks if the neigboring particles in within support or not.
 */
inline bool
in_support(
            //! cartesian projections of inter particle distance
            const double *dx,
            //! support
            double h)
{
  for (int i = 0; i < DIMENSION; i++)
    if (abs(*(dx + i)) > h)
      return false;
  return true;
};

/*! 
 * weight() computes and returns the weight of a point, based on its
 * location. Gussian weight function is used. 
 */
double weight(
               //!  s = (x_i-x_j) / h
               double *s,
               //! smoothing length 
               double h);

/*!
 * d_weight computes and returns the derivative of Gaussian weight
 * function and returns its value. Similar to weight(), function is 
 * overloaded.
 */
double d_weight(
                 //! sx = (x-xj)/h
                 double *s,
                 //! smoothing length 
                 double h,
                 //! { 0, 1, 2 }: direction of differentiation 
                 int dir);

/*!
 * rotate() rotates the vector of state variables anti-clockwise. 
 * It involves rotation of velocity vector and rotation of stress-tensor.
 * It is used to rotate velocity etc. to inter-particle local coordinate system
 */
void rotate(
             //! the vector of state variables
             double *u,
             //! consines of rotation
             double *cosines);

/*! 
 * Reflects vector \f$\overline{u}\f$ in plane with \f$\hat{n}\f$
 * normal
 */
void reflect(
              //! incident vector
              double *,
              //! reflection
              double *,
              //! normal to plane of reflection
              double *);

/*! 
 * dot product of vector A and vector B
 */
double dot(
            //! Vector A
            double *,
            //! Vecotr B
            double *);

/*!
 *  Rotation of Vector.
 *  For any higher order multiplication BLAS etc is preferred
 */
void matrix_vec_mult(
                      //! Matix A ( assumed square matrix )
                      double *A,
                      // matirx leading dim, i.e. N
                      int N,
                      //! Vector b
                      double *b,
                      //! c = A*b
                      double *c);

/*!
 *  For very small dimensional matericies 
 *  For any higher order multiplication BLAS etc is preferred
 */
void matrix_matrix_mult(
                         //! Matix A ( N x P )
                         double *A,
                         //! Matix A leading dim, i.e. N
                         int N,
                         //! Matrix A lagging dim i.e. P
                         int P,
                         //! Matrix B ( P x M )
                         double *B,
                         //! Matrix B lagging dim, i.e. M
                         int M,
                         //! C = A*B
                         double *C);

/*!
 *  Solve Liner equations
 *  Used to compute velocity gradients
 */
void linsolve(
               //! Matix A
               double *A,
               //! Matix Dimension, i.e. (N x N)
               int N,
               //! RHS ( set of vectors, in columns )
               double *b,
               //! No of vectors in the set, (N x M)
               int M,
               //! solution (N x M
               double *d);

void sph_exit(
               //! Error code
               int);

#endif //SPH_HEADER__H
