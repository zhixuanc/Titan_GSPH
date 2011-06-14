/*
 * =====================================================================================
 *
 *       Filename:  calc_fcoef.cc
 *
 *    Description:  
 *
 *        Created:  10/21/2010 01:34:08 PM EDT
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

#include <vector>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <bucket.h>
#include <particle.h>
#include <properties.h>

#include "constants.h"
#include "sph_header.h"

int calc_f_coef(int myid, HashTable *P_table, 
                HashTable *BG_mesh, MatProps *matprops)
{
  int i, j, k, l;
  int count;
  double dx[DIMENSION], sj[DIMENSION], vel[DIMENSION];
  double coord[DIMENSION];
  double bnd_xcrd[PARTICLE_DENSITY];
  double bnd_ycrd[PARTICLE_DENSITY];
  double f_coef[PARTICLE_DENSITY];
  double hj, mj, rhoj;
  vector<Key> plist;
  vector<Key>::iterator p_itr;

  double beta = matprops->bedfrict;
  HTIterator *itr = new HTIterator(BG_mesh);
  Bucket *curr_bucket;
  while ( curr_bucket = (Bucket *) itr->next() )
    if ( curr_bucket->is_active() && !(curr_bucket->is_guest()) &&
         (curr_bucket->get_bucket_type() == MIXED) )
    {
      // initalize 
      for (i=0; i<PARTICLE_DENSITY; i++)
      {
        bnd_xcrd[i] = *(curr_bucket->get_bnd_xcrd()+i);
        f_coef[i] = 0;
      }

      count=0;
      for (i=0; i<PARTICLE_DENSITY; i++)
      {
        double wnorm  = 0.;
        // coords of bnd_pnt
        coord[0] = bnd_xcrd[0];
        coord[1] = curr_bucket->get_bndZ(coord);

        double slope = -1./curr_bucket->get_bndnorm(coord);
        double cost = 1.0/sqrt(1 + pow(slope,2.));
        double sint = slope*cost;

        // look into thyself
        plist = curr_bucket->get_plist();
        for (p_itr=plist.begin(); p_itr!=plist.end(); p_itr++)
        {
          Particle *pj = (Particle *) P_table->lookup(*p_itr);
          hj = pj->get_smlen();
          mj = pj->get_mass();
          rhoj = pj->get_density();
          for (l=0; l<DIMENSION; l++)
          {
            dx[l] = coord[l] - *(pj->get_coords()+l);
            sj[l] = dx[l]/hj;
            vel[l] = *(pj->get_vel()+l);
          }
          double w=weight(sj,hj);
          // velocity parallel to the boundary
          double vel_dot_n = vel[0]*cost + vel[1]*sint;
          f_coef[count] -= mj*beta*vel_dot_n*w/rhoj;
            wnorm  += mj*w/rhoj;
        }
        // look into neighbors
        Key *neighbors = curr_bucket->get_neighbors();
        for ( j=0; j<NEIGH_SIZE; j++ )
          if ( *(curr_bucket->get_neigh_proc()+j) > -1 )
          {
            Bucket *neigh = (Bucket *) BG_mesh->lookup(neighbors[j]);
            if ( !(neigh) && 
                 (*(curr_bucket->get_neigh_proc()+j) != myid) )
              continue;

            // skip bucket totally below the boundary
            if ( neigh->get_bucket_type() == UNDERGROUND )
              continue;

            plist = neigh->get_plist();
            // get contribution of each particle in the support
            for ( k=0; k < plist.size(); k++)
            {
              Particle *pj = (Particle *) P_table->lookup(plist[k]);
              if ( pj->is_real() )
              {
                hj = pj->get_smlen();
                mj = pj->get_mass();
                rhoj = pj->get_density();
                for (l=0; l<DIMENSION; l++)
                { 
                  dx[l] = coord[l] - *(pj->get_coords()+l);
                  vel[l] = *(pj->get_vel()+l);
                  sj[l] = dx[l]/hj;
                }
                double w=weight(sj,hj);
                // velocity parallel to the boundary
                double vel_dot_n = vel[0]*cost + vel[1]*sint;
                f_coef[count] -= mj*beta*vel_dot_n*w/rhoj;
                wnorm  += mj*w/rhoj;
              }
            }
          }
        if ( wnorm > TINY )
          f_coef[count] /= wnorm;
        else
          f_coef[count] = 0.;
        count++;
      }
      curr_bucket->put_f_coef(f_coef);
    }
  return 0;
}
