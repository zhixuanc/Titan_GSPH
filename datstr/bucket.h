
/*
 * =====================================================================================
 *
 *       Filename:  bucket.h
 *
 *    Description:  
 *
 *        Created:  04/11/2010 05:46:27 PM EDT
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

#ifndef BUCKET__H
#  define BUCKET__H

#  include <vector>
using namespace std;

#  include <hashtab.h>
#  include <constants.h>
#  include <properties.h>
#  include <pack_data.h>

const int FIRST_BIT_UP = 0x1;
const int SECND_BIT_UP = 0x2;

// Bucket is a unit of background mesh
class Bucket
{
  // friends in repartition
  friend void pack_bucket (BucketPack *, Bucket *, int);
  friend void unpack_bucket (BucketPack *, Bucket *, int);

private:
    Key key;
  double lb_weight;
  double mincrd[DIMENSION];
  double maxcrd[DIMENSION];
  double bnd_xcrd[PARTICLE_DENSITY];
  double bnd_ycrd[PARTICLE_DENSITY];
  double f_coef[PARTICLE_DENSQRD * DIMENSION];
  double poly[4];
  bool active;
  bool guest_flag;
  int newold;
  int particles_type;
  int myprocess;
  int bucket_type;
  int neigh_proc[NEIGH_SIZE];
  Key neighbors[NEIGH_SIZE];
  vector < Key > particles;
  vector < Key > new_plist;

public:
  //! Contructors
  //! base constructor
    Bucket ();

  //! constructor called from Read_data();
    Bucket (
             //! Hashtable key of the Bucket
             unsigned *,
             //! Minimum coordinates
             double *,
             //! Maximum coordinates
             double *,
             //! Bucket type ( over-ground, underground or boundary)
             int,
             //! boundary elevations in order [i,j; i+1,j; i+1,j+1; i,j+1]
             double *,
             //! my procees id
             int,
             //! Neighbor process info
             int *,
             //! Array of neighboring buckets HT keys
             Key *);

  //! change process id (only called from repartition)
  void put_myprocess (int myid)
  {
    myprocess = myid;
  }

  //! put repartition weights
  void put_lb_weight (double wght)
  {
    lb_weight = wght;
  }

  //! mark bucket as guest, i.e. belongs to different proc
  void put_guest_flag (int fl)
  {
    guest_flag = fl;
  }

  //! add a particle to the bucket
  void add_particle (Key pk)
  {
    new_plist.push_back (pk);
  }

  //! put a list of particles
  void put_plist (vector < Key > pl)
  {
    particles = pl;
  }

  //! put a list of new particles
  void put_new_plist (vector < Key > pl)
  {
    new_plist.insert (new_plist.end (), pl.begin (), pl.end ());
  }

  //! update particles to new particles
  void update_particles ()
  {
    particles = new_plist;
    new_plist.clear ();
  }

  //! empty particle vector
  void empty_plist ()
  {
    particles.clear ();
    new_plist.clear ();
  }

  //! puts ghosts in underground buckets
  void add_ghosts (vector < Key > ghosts)
  {
    particles.insert (particles.end (), ghosts.begin (), ghosts.end ());
  }

  //! mark the bucket active
  void mark_active ()
  {
    active = true;
  }

  //! mark bucket inactive
  void mark_inactive ()
  {
    active = false;
  }

  //! check if point is contained in bucket or not
  bool contains (double pnt[]) const
  {
    for (int i = 0; i < DIMENSION; i++)
      if ((pnt[i] < mincrd[i]) || (pnt[i] >= maxcrd[i]))
        return false;
    return true;
  }

  //! put have real particles flags
  void set_real_particles (bool flag)
  {
    if (flag)
      particles_type |= FIRST_BIT_UP;
    else
      particles_type &= SECND_BIT_UP;
  }

  //! put have ghosts
  void set_ghost_particles (bool flag)
  {
    if (flag)
      particles_type |= SECND_BIT_UP;
    else
      particles_type &= FIRST_BIT_UP;
  }

  //! put newold information
  void put_new_old (int info)
  {
    newold = info;
  }

  // update friction_coefs
  void put_f_coef (double *fcoef)
  {
    int npts = PARTICLE_DENSQRD * DIMENSION;

    for (int i = 0; i < npts; i++)
      f_coef[i] = *(fcoef + i);
  }

  // access methods
  //! Access HT key of current bucket
  Key getKey () const
  {
    return key;
  }

  //! get my process id
  int get_myprocess ()
  {
    return myprocess;
  }

  //! Access minimum coordinates of current bucket
  const double *get_mincrd () const
  {
    return mincrd;
  }

  //! Access maximum coordinates of current bucket
  const double *get_maxcrd () const
  {
    return maxcrd;
  }

  //! Access bucket-type
  int get_bucket_type ()
  {
    return bucket_type;
  }

  //! check if bucket belongs to different proc
  int is_guest ()
  {
    return guest_flag;
  }


  //! get newold info
  int get_new_old ()
  {
    return newold;
  }

  //! get repartition weights
  double get_lb_weight () const
  {
    return lb_weight;
  }

  //! Compare hash-table keys for equality
  bool compare_keys (Key k1, Key k2) const
  {
    for (int i = 0; i < KEYLENGTH; i++)
      if (k1.key[i] != k2.key[i])
        return false;
    return true;
  }

  //! Value of elevation z(x,y) using linear interpolation
  double get_bndZ (double x[]) const
  {
    return (poly[0] * x[0] + poly[1] * x[1] + 
            poly[2] * x[0] * x[1] + poly[3]);
  }

  // get array of x-boundary points
  const double *get_bnd_xcrd () const
  {
    return bnd_xcrd;
  }

  // get array of y-boundary points
  const double *get_bnd_ycrd () const
  {
    return bnd_ycrd;
  }

  // get array of friction coefficients
  double get_f_coef (int i, int j, int k) const
  {
    return f_coef[i * (PARTICLE_DENSITY + DIMENSION) + j * DIMENSION + k];
  }

  // update neigh_proc
  void put_neigh_proc (int *, int);

  /*! 
   * distance of point from the boundary, 
   * also find point of intersection with perpendicular
   */

  //! Check is the bucket is active/inactive
  bool is_active () const
  {
    return active;
  }

  //! check if bucket has any real particles
  bool has_real_particles () const
  {
    return (particles_type & FIRST_BIT_UP);
  }

  //! check if bucket has ghost particles
  bool has_ghost_particles () const
  {
    return (particles_type & SECND_BIT_UP);
  }

  //! get list of particles in the current bucket
  vector < Key > get_plist ()
  {
    return particles;
  }

  //!  array of neighboring buckets
  Key *get_neighbors ()
  {
    return neighbors;
  }

  //! get neigh_proc info. neigh_proc also tells if neigh is boundary
  const int *get_neigh_proc () const
  {
    return neigh_proc;
  }

  //! get HT key of the neighbor is  up,down etc direction
  Key which_neigh (int dir[]) const;

  //! get neigh_proc info in up, down etc direction
  int which_neigh_proc (int dir[]) const;

  //! find out direction of the neighbor
  bool find_neigh_dir (Key k, int dir[]) const;

  /*!
   *  get boundary normal and return 0 if the bucket is a boundary bucket
   *  else return 1
   *  @param pnt : coordinates of any point, ususally a SPH particle
   *  @param normal : cosines of normal to the boundary
   *  @return 
   *  0 if current bucket is a boundary bucket, 1 otherwise
   */
  int get_bnd_normal (double * pnt, double * normal) const;

  //!  get distance of point x fom the boundary, 
  double get_bnddist (double * pnt, double * intsct) const;

  /*! calculate intersection of line and the boundary, such
   * that the line is normal to boundary at pt. of intersection
   */
  int calc_intersection (double * pnt, double * intsct) const;

  //*! Add particles for initial piles. 
  void add_real_particle (Key k)
  {
    active = true;
    set_real_particles (true);
    particles.push_back (k);
  }

  //! Add ghost particles to the bucket
  int put_ghost_particles (
    //! Particle HashTable
    HashTable *,
    //! Background Mesh
    HashTable *,
    //! Material properties
    MatProps *);

  void Copy_data (
    //! void pointer to datastream
    void *,
    //! HashTable of particles
    HashTable *,
    //! current process id
    int);
};

#endif // BUCKET__H
