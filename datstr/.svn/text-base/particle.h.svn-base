/*
 * =====================================================================================
 *
 *       Filename:  partilce.h
 *
 *    Description:  
 *
 *        Created:  01/07/2008 12:41:11 PM EST
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

#ifndef PARTICLE__H
#define PARTICLE__H

#include <vector>
using namespace std;

#include <constants.h>
#include <hashtab.h>
#include <pack_data.h>


class Particle
{
  // friend functions in repartition
  friend void unpack_particle (ParticlePack *, Particle *);

  private:

    //! Data Structure properties
    Key key;

    //! mass of the each particle 
    double mass;

    //! support of each particle
    double smlen;

    //! current (normalized) postion of each particle
    double coord[DIMENSION];

    //! Bedfriction
    double bedfrict[DIMENSION];

    //! state_vars are normalized velocities, stresses etc
    double state_vars[NO_OF_EQNS];

    //! new_state_vars stores updated state variables
    double new_state_vars[NO_OF_EQNS];

    //! derivatives of 
    double d_vel[DIMSQRD];

    //! gravity vector
    double gravity[DIMENSION];

    //!  neighbors
    vector<Key> neighs;

    //! if the particle is real, guest or a ghost
    bool ghost, guest;
    
    //! update delayed
    bool update_delayed;

    //! new_old for ghost updates across process boundaries
    int new_old;

    //! aleardy searched reflection, if ghost, if flag is up
    bool reflection;
  public:
    //! Constructors
    //! Base constructor
    Particle ();

    //! Constructor for piles and boundary ghosts
    Particle (
              //! hash-table key
              unsigned *,
              //! particle coodinates
              double *, 
              //! particle mass
              double , 
              //! smoothing length
              double , 
              //! ghost/real flag
              int
              //! 
             );
    //! Constructor 2 for flux source
    Particle (
              //! hash-table key
              unsigned *,
              //! particle coordinates
              double *,
              //! particle mass
              double ,
              //! particle smoothing length
              double ,
              //! initial velocties
              double *
             );

    //! get hash-table key 
    Key getKey() const 
    { 
      return key; 
    }

    //! get particle density
    double get_density() const 
    { 
      return state_vars[0]; 
    }

    //! get particle mass
    double get_mass() const 
    { 
      return mass; 
    }

    //! get smoothing length of current paricle
    double get_smlen() const 
    { 
      return smlen; 
    }

    //!
    const double * get_coords() const 
    { 
      return coord; 
    }

    //!
    const double * get_vel() const 
    { 
      return (state_vars+1); 
    }

    //!
    const double * get_bedfrict () const 
    { 
      return bedfrict; 
    }

    //!
    const double * get_state_vars() const 
    { 
      return state_vars; 
    }

    //!
    const double * get_new_state_vars() const 
    { 
      return new_state_vars; 
    }

    //!
    const double * get_d_vel () const 
    { 
      return d_vel; 
    }

    //!
    const double * get_gravity() const 
    { 
      return gravity; 
    }

    //!
    vector<Key> get_neighs () const 
    { 
      return neighs; 
    }

    //! get new_old info
    int get_new_old ()
    {
      return new_old;
    }

    //! check if particle is real
    bool is_real() const 
    { 
      return ((!ghost) && (!guest));
    }

    //! check if particle is ghost
    bool is_ghost() const 
    { 
      return ghost; 
    }

    //!
    bool is_guest () const  
    { 
      return guest;
    }

    //!
    bool is_not_updated () const  
    { 
      return update_delayed; 
    }

    //! check if already searched the ghost reflection 
    bool have_reflection()
    {
      return reflection;
    }
   
    //! update density value
    void put_density(double den)
    { 
      state_vars[0]=den; 
    }

    //! update smoothing length
    void put_smlen (double h)
    {
      smlen=h;
    }

    //! update neighbor information
    void put_neighs(vector<Key> n)
    {
      neighs = n;
    }

   //! get delayed update flag
    void put_update_delayed (bool val)
    {
      update_delayed = val;
    }

    //! update bed-friction coefficents
    void put_bedfrict (double *bfric)
    {
      for (int i=0; i<DIMENSION; i++)
        bedfrict[i] = *(bfric + i);
    }

    //! update new state variables 
    void put_new_state_vars(double u[])
    {
      for (int i=0; i<NO_OF_EQNS; i++)
        new_state_vars[i]=u[i];
    }

    //! put state_vars. This should only be uesed for ghost particles
    void put_state_vars(double u[])
    {
      for (int i=0; i<NO_OF_EQNS; i++)
        state_vars[i]= u[i];
    }

    //! update state_vars
    void update_state_vars()
    {
      for ( int i=0; i<NO_OF_EQNS; i++)
        state_vars[i]=new_state_vars[i];
    }

    // update slopes
    void put_d_vel(double *du)
    {
      for (int i=0; i<DIMSQRD; i++)
        d_vel[i]= *(du+i);
    }

    // update partilce positions
    void put_coords(double *x)
    {
      for (int i=0; i<DIMENSION; i++)
        coord[i]= *(x+i);
    }

    //! update guest info
    void put_guest_flag (bool val)
    {
      guest = val;
    }

    //! update new_old info
    void put_new_old (int info)
    {
      new_old = info;
    }

    //! set reflection flag
    void set_have_reflection (bool val)
    {
      reflection = val;
    }

    // Operator overloads
    bool operator==(const Particle &rhs) const;
};


#endif // PARTICLE__H
