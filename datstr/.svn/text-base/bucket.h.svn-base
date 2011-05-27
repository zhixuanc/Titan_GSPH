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
#define BUCKET__H

#include <vector>
using namespace std;

#include <hashtab.h>
#include <constants.h>
#include <properties.h>
#include <pack_data.h>

// boundary point structure
//typedef struct
//{
//  double coord[DIMENSION];   
//  double f_coef;
//} BndPnt;



// Bucket is a unit of background mesh
class Bucket
{
  // friends in repartition
  friend void unpack_bucket (BucketPack *, Bucket *, int);

  private:
    Key         key;
    double      mincrd[DIMENSION];
    double      maxcrd[DIMENSION];
    double      bnd_xcrd[PARTICLE_DENSITY];
#ifdef THREE_D
    double      poly[4];
    double      bnd_ycrd[PARTICLE_DENSITY];
    double      f_coef[PARTICLE_DENSQRD];
#else
    double      poly[3];
    double      f_coef[PARTICLE_DENSITY];
#endif
    double      lb_weight;
    bool        real_part_flag;
    bool        ghost_part_flag;
    bool        active;
    bool        guest_flag;
    int         boundary_type;
    int         bucket_type;
    int         neigh_proc[NEIGH_SIZE];
    int         myprocess;
    Key         neighbors[NEIGH_SIZE];
    vector<Key> particles;
    vector<Key> new_plist;

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
         //! Bucket type ( over-ground, underground or mixed)
         int ,
         //! Boundary type ( LINE or other shape )
         int ,
         //! polynomial coefs that define boundary (if its is boundary bucket)
         double *,
         //! my procees id
         int ,
         //! Neighbor process info
         int *,
         //! Array of neighboring buckets HT keys
         Key * 
         );


    //! change process id (only called from repartition)
    void put_myprocess (int mypid)
    {
      myprocess = mypid;
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

    //! get repartition weights
    double get_lb_weight ()
    {
      return lb_weight;
    }

    //! Compare hash-table keys for equality
    bool compare_keys (Key k1, Key k2)
    {
      for (int i=0; i < KEYLENGTH; i++ )
        if ( k1.key[i] != k2.key[i] )
          return false;
      return true;
    }

    //! add a particle to the bucket
    void add_particle(Key pk)
    { 
      new_plist.push_back(pk); 
    }

    //! put a list of particles
    void put_plist (vector<Key> pl)
    { 
      particles = pl;
    }

    //! put a list of new particles
    void put_new_plist(vector<Key> pl)    
    { 
      new_plist.insert(new_plist.end(), pl.begin(), pl.end());
    }

    //! update particles to new particles
    void update_particles()               
    { 
      particles = new_plist;
      new_plist.clear();
    }

    //! empty particle vector
    void empty_plist()
    {
      particles.clear();
      new_plist.clear();
    }

    //! puts ghosts in underground buckets
    void add_ghosts(vector<Key> ghosts)
    { 
      particles.insert(particles.end(), ghosts.begin(), ghosts.end());
    }   

    //! mark the bucket active
    void mark_active()
    { 
      active = true;
    }

    //! mark bucket inactive
    void mark_inactive()
    { 
      active = false;
    }

    //! check if point is contained in bucket or not
    bool contains(double pnt[])
    {
      for (int i=0; i<DIMENSION; i++)
        if ((pnt[i] < mincrd[i]) || (pnt[i] > maxcrd[i]))
          return false;
      return true;
    }

    //! put have real particles flags
    void put_have_real_particles (bool flag)
    { 
      real_part_flag = flag;
    }

    //! put have ghosts
    void put_have_ghost_particles (bool flag)
    {
      ghost_part_flag = flag;
    }

    // access methods 
    //! Access HT key of current bucket
    Key getKey() const
    {  
      return key;
    }

    //! get my process id
    int get_myprocess()
    {
      return myprocess;
    }

    //! boundary type ... line or circle?
    int get_bndtype () const
    {
      return boundary_type;
    }

    //! Access minimum coordinates of current bucket
    const double * get_mincrd() const
    { 
      return mincrd;
    }

    //! Access maximum coordinates of current bucket
    const double * get_maxcrd() const
    { 
      return maxcrd;
    }

    //! Access boundary normal (if it is a boundary bucket
    double get_bndnorm(double []) const;

    //! Access bucket-type
    int get_bucket_type()
    { 
      return bucket_type; 
    }

    //! check if bucket belongs to different proc
    int is_guest()
    {
      return guest_flag;
    }
    //! Get cofficients of boundary equatons
    const double * get_bndcoeff() const 
    { 
      return poly; 
    }

    // Value of elevation z(x,y) ... 
    // using linear interpolation
    double get_bndZ (double *x)
    {

#ifdef THREE_D
      return (-(poly[0]*(*x) + poly[1]*(*(x+1)) + poly[3])/poly[2]);
#else
      switch ( boundary_type )
      {
        case LINE:
          return ((-poly[0]*(*x) - poly[2])/poly[1]);
        case CIRCLE:
          return (poly[1] + sqrt(pow(poly[2],2) - pow(*x-poly[0],2)));
      }
#endif
      return -1.;
    }

    // get array of x-boundary points
    const double * get_bnd_xcrd() const
    {
      return bnd_xcrd;
    }

    // get array of y-boundary points
#ifdef THREE_D
    const double * get_bnd_ycrd() const
    {
      return bnd_ycrd;
    }
#endif
      
    // get array of friction coefficients
    const double * get_f_coef() const
    {
      return f_coef;
    }

    // update friction_coefs
    void put_f_coef ( double *fcoef )
    {
#ifdef THREE_D
       int npts = PARTICLE_DENSQRD;
#else
       int npts = PARTICLE_DENSITY;
#endif
       for (int i=0; i<npts; i++)
         f_coef[i] = *(fcoef+i);
    }

    // update neigh_proc
    void put_neigh_proc (int *, int );

    /*! 
     * distance of point from the boundary, 
     * also find point of intersection with perpendicular
     */
    double get_bnddist(double *, double *);

    //! Check is the bucket is active/inactive
    bool is_active() const
    { 
      return active;
    }

    //! check if bucket has any real particles
    bool have_real_particles () const
    { 
      return real_part_flag;
    }


    //! check if bucket has ghost particles
    bool have_ghost_particles ()
    {
      return ghost_part_flag;
    }

    //! get list of particles in the current bucket
    vector<Key> get_plist()
    { 
      return particles;
    }

    //!  array of neighboring buckets
    Key * get_neighbors()
    { 
      return neighbors;
    }

    //! get neigh_proc info. neigh_proc also tells if neigh is boundary
    const int * get_neigh_proc() const
    { 
      return neigh_proc;
    }

    //! get HT key of the neighbor is  up,down etc direction
    Key which_neigh( int dir[] ) const;

    //! get neigh_proc info in up, down etc direction
    int which_neigh_proc( int dir[] ) const;

    //! find out direction of the neighbor
    bool find_neigh_dir (Key k, int d[]);

    //*! Add particles for initial piles. 
    void add_real_particle(Key k)
    {
       active = true;
       real_part_flag = true;
       particles.push_back(k);
    }
 
    //! Add ghost particles to the bucket
    int put_ghost_particles (
                //! Particle HashTable
                HashTable *,
                //! Background Mesh
                HashTable *,
                //! Material properties
                MatProps *
              );

    //! Mark buckets that contains real particles
    void add_extra_ghosts (
                //! Particle HashTable
                HashTable *,
                //! Background Mesh
                HashTable *,
                //! Material properties
                MatProps *
              );
    void Copy_data (
                //! void pointer to datastream
                void *,
                //! HashTable of particles
                HashTable *,
                //! current process id
                int 
              );
};


#endif // BUCKET__H
