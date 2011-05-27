/*
 * =====================================================================================
 *
 *       Filename:  inc_plane.cc
 *
 *    Description:  
 *
 *        Created:  07/26/2010 03:15:45 PM EDT
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

#include <cmath>
#include <constants.h>
#include "gisapi.h"

double abs ( double a)
{
  return (a > 0 ? a : -a);
}

double GIS_get_elevation (double x, int *polytype)
{
  double z = Slope * x + Intcpt;
  *polytype = LINE;
  return 0.2;
}

int Get_xmin(double *crd)
{
  *crd = -0.02;
  return 0;
}

int Get_ymin(double *crd)
{
  *crd = -0.01;
  return 0;
}

int Get_xmax(double *crd)
{
  *crd = 1.5;
  return 0;
}

int Get_ymax (double *crd)
{
  *crd = 1.4;
  return 0;
}
