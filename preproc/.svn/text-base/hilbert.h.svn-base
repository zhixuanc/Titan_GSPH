/*
 * =====================================================================================
 *
 *       Filename:  hilbert.h
 *
 *    Description:  
 *
 *        Created:  03/15/2010 01:25:11 PM EDT
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

#ifndef HILBERT_H
#define HILBERT_H

/* 2D Hilbert Space-filling curve */

#ifdef __cplusplus
extern "C" {
#endif

void hsfc2d(
  unsigned   coord[] , /* IN: Normalized integer coordinates */
  unsigned * nkey ,    /* IN: Word length of key */
  unsigned   key[]    /* OUT: space-filling curve key */
  );


/* 3D Hilbert Space-filling curve */

void hsfc3d(
  unsigned   coord[] , /* IN: Normalized integer coordinates */
  unsigned * nkey ,    /* IN: Word length of 'key' */
  unsigned   key[]     /* OUT: space-filling curve key */
  );

/* API for 2-D Hilbert Space Filling Curve */
void HSFC2d (
  double     coord[] , /* IN: Normalized floating point coordinates */
  unsigned * nkey ,    /* IN: Word length of key */
  unsigned   key[]     /* OUT: space-filling curve key */
  );

/* API for 3-D Hilbert Space Filling Curve */
void HSFC3d (
  double     coord[] , /* IN: Normalized floating point coordinates */
  unsigned * nkey ,    /* IN: Word length of key */
  unsigned   key[]    /* OUT: space-filling curve key */
  );

#ifdef __cplusplus
}
#endif

#endif // HILBERT_H
