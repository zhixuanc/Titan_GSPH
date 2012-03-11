
/*
 * =====================================================================================
 *
 *       Filename:  properties.h
 *
 *    Description:  
 *
 *        Created:  03/23/2010 02:01:20 PM EDT
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

#ifndef PRPERTIES__H
#  define PRPERTIES__H

#  include <cmath>
using namespace std;

#  include <hashtab.h>

struct TimeProps
{
  //! current non-dimensional time
  double time;

  //! current non-dimensional time-step
  double dtime;

  //! maximum simulation time in sec
  double max_time;

  //! non-dimensional max simulation time
  double ndmax_time;

  //! time ineterval after which output is to be written
  double timeoutput;

  //! non-dimensional timeoutput
  double ndtimeoutput;

  //! time-period beteen matieral injection
  double fluxrate_time;

  //! material injection flag
  bool mat_add_time;

  //! material add counter
  int imatin;

  //! current time-step
  int step;

  //! count of previous outputs
  int ioutput;

  //! maximum time-steps allowed
  int max_steps;

  //! Scale to Normalize the time
  double TIME_SCALE;

  //! constructor to allocate default values
    TimeProps ()
  {
    max_time = 1.0;
    timeoutput = 0.02;
    TIME_SCALE = 1.0;
    ndmax_time = max_time = 1.;
    ndtimeoutput = 0.02;
    ioutput = 0;
    imatin = 0;
    step = 0;
    max_steps = 10;
  }

  void incrtime (double *dt)
  {
    // first reduce dt to hit output or end time "exactly"
    if ((time + (*dt)) > ndtimeoutput)
      (*dt) = ndtimeoutput - time;
    if ((time + (*dt)) > ndmax_time)
      (*dt) = ndmax_time - time;
    // then increment time
    time += (*dt);
    dtime = (*dt);
    step++;
  }

  bool ifstart ()
  {
    return (step == 0);
  }                             //before first time step
  bool iffirst ()
  {
    return (step == 1);
  }                             //at first time step
  bool ifend ()
  {
    return ((time >= ndmax_time) || (step > max_steps));
  }

  bool ifoutput ()
  {
    if (time >= ndtimeoutput)
    {
      ioutput++;                //using ioutput eliminates roundoff
      ndtimeoutput = ((ioutput + 1) * timeoutput) / TIME_SCALE;
      return true;
    }
    else
      return false;
  }

  bool addmaterial ()
  {
    if (mat_add_time)
    {
      mat_add_time = false;
      return true;
    }
    return false;
  }

  void chunktime (int *hours, int *minutes, double *seconds)
  {
    double dimtime = time * TIME_SCALE;

    *hours = ((int) dimtime) / 3600;
    *minutes = (((int) dimtime) % 3600) / 60;
    *seconds = dimtime - (double) (*hours * 3600 + *minutes * 60);
  }

  double timesec ()
  {
    return (time * TIME_SCALE);
  }

};

struct MatProps
{
  //! phi_{int}, the internal friction angle (must be GREATER than the bedfriction angle)
  double intfrict;

  //! tan(phi_{int}), tangent of the internal friction angle
  double tanintfrict;

  //! sin(phi_{int}), sine of the internal friction angle
  double sinintfrict;

  //! phi_{bed}, the bed friction angle, must be LESS than the internal friction
  //! angle and should be greater than about 8 degrees 
  double bedfrict;

  //! tan(phi_{bed}), tangent of the bed friction angle
  double tanbedfrict;

  //! constant in equation of state
  double P_CONSTANT;

  //! slope limiting stuff
  double GAMMA;

  //! SPH smoothing length
  double smoothing_length;

  //! mass of individual particle
  double particle_mass;

  //! length scaling factor
  double LENGTH_SCALE;

  //! gravity scaling factor
  double GRAVITY_SCALE;

  //! normaliszed initial density
  double Rho0;

  //! ZERO Approximation
  double TINY;

  //! constructor allocates default properties 
    MatProps ()
  {
    intfrict = 0.3491;          // 20 deg
    tanintfrict = tan (intfrict);
    bedfrict = 0.1745;          // 10 deg
    tanbedfrict = tan (bedfrict);
    Rho0 = 1;
    P_CONSTANT = 1.;
    GAMMA = 7.;
    LENGTH_SCALE = 1.;
    GRAVITY_SCALE = 9.81;
    TINY = 1.0E-06;
  }

  double pressure (double rho)
  {
    double n = GAMMA;
    double k = P_CONSTANT;
    double tmp = rho / Rho0;

    return k * (pow (tmp, n) - 1);
  }
  double sound_speed (double rho)
  {
    double n = GAMMA;
    double p = pressure (rho);
    double k = P_CONSTANT;
    double c = rho <= TINY ? 0. : sqrt (n * (p + k) / rho);

    return c;
  }
};

struct PileProps
{
  //! Number of piles
  int NumPiles;

  //! Maximum pile height
  double *pileheight;

  //!array holding x coordinate of pile center
  double *xCen;

  //!array holding y coordinate of pile center
  double *yCen;

  //!array holding y coordinate of pile center
  double *zCen;

  //! Array of HT Keys of Buckets containing center
  Key *CenBucket;

  //!array holding the major (x before rotation) radius
  double *majorrad;

  //!array holding the minor (y before rotation) radius
  double *minorrad;

  //!array holding the cosine of the rotation angle
  double *cosrot;

  //!array holding the sine of the rotation angle;
  double *sinrot;

  //! this constuctor initializes the number of piles to zero.
    PileProps ()
  {
    NumPiles = 0;
  }

  //! function allocates space for the pile data
  void allocpiles (int nump)
  {
    NumPiles = nump;
    pileheight = new double[nump];
    xCen = new double[nump];
    yCen = new double[nump];
    zCen = new double[nump];

    CenBucket = new Key[nump];
    majorrad = new double[nump];
    minorrad = new double[nump];
    cosrot = new double[nump];
    sinrot = new double[nump];
  }

  //! this function deallocates the dynamically out array members of the PileProps structure
  ~PileProps ()
  {
    if (NumPiles > 0)
    {
      delete[]pileheight;
      delete[]xCen;
      delete[]yCen;
      delete[]zCen;
      delete[]majorrad;
      delete[]minorrad;
      delete[]cosrot;
      delete[]sinrot;
    }
  }
};

struct FluxProps
{

  bool have_src;

  //! X coord of flux source
  double xSrc;

  //! source bucket
  Key bucketsrckey;

  //! start time
  double starttime;

  //! stop time
  double stoptime;

  //! tangenital velocity function
  double tangvel;

    FluxProps ()
  {
    have_src = false;
  }
};

#endif // PRPERTIES__H
