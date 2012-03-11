
/*
 * =====================================================================================
 *
 *       Filename:  update_pos.cc
 *
 *    Description:  
 *
 *        Created:  04/02/2010 11:54:04 AM EDT
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

int
update_pos(int myid, HashTable * P_table, HashTable * BG_mesh,
           FluxProps * fluxprops, TimeProps * timeprops, int *lost)
{

  int dir[DIMENSION];
  double pos[DIMENSION], coord[DIMENSION], vel[DIMENSION];
  double mincrd[DIMENSION], maxcrd[DIMENSION];
  double smooth_vel[DIMENSION], vj[DIMENSION], wnorm;
  double dx[DIMENSION], s[DIMENSION], junk[DIMENSION];
  unsigned pkey[KEYLENGTH];
  int i, j;
  int adapt = 0;
  double bndnorm;
  const double v_coef = 0.5;
  double dt = timeprops->dtime;
  double fluxsrc[2];
  bool add_material_check, add_material;

  // coordinate of flux source
  if (fluxprops->have_src)
  {
    fluxsrc[0] = fluxprops->xSrc;
    add_material = true;
  }
  else
    add_material = false;

  // update particle positions
  HTIterator *itr = new HTIterator(P_table);
  Particle *p;

  while (p = (Particle *) itr->next())
    // guest particles are moved as well as some of them
    // will move to different partitions
    if (!(p->is_ghost()))
    {
      // velocity and coodinates
      for (i = 0; i < DIMENSION; i++)
      {
        vel[i] = *(p->get_vel() + i);
        coord[i] = *(p->get_coords() + i);
      }

      // update particle positions
      for (i = 0; i < DIMENSION; i++)
        pos[i] = coord[i] + dt * vel[i];
      p->put_coords(pos);
    }

  // move-in and move out particles from buckets
  vector < Key > plist;
  vector < Key >::iterator p_itr;
  vector < Key > my_realp, my_particles;

  HTIterator *it2 = new HTIterator(BG_mesh);
  Bucket *curr_bucket, *buck_save;

  while (curr_bucket = (Bucket *) it2->next())
    if (curr_bucket->is_active())
    {
      my_realp.clear();
      my_particles.clear();

      // get bucket limits
      for (i = 0; i < DIMENSION; i++)
      {
        mincrd[i] = *(curr_bucket->get_mincrd() + i);
        maxcrd[i] = *(curr_bucket->get_maxcrd() + i);
      }

      // on/off  add_material_check flag
      add_material_check = false;
      if ((fluxprops->have_src) && (curr_bucket->get_bucket_type() == MIXED) &&
          (fluxsrc[0] > mincrd[0]) && (fluxsrc[0] < maxcrd[0]))
        add_material_check = true;

      // get list of particles in the bucket
      plist = curr_bucket->get_plist();
      for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
      {
        Particle *p_curr = (Particle *) P_table->lookup(*p_itr);

        assert(p_curr);
        if (!(p_curr->is_ghost()))
        {
          for (i = 0; i < DIMENSION; i++)
            pos[i] = *(p_curr->get_coords() + i);

          // compute distance from boundary, 
          // ... if current particle is real
          double bnddist = 1 / TINY;

          if (p_curr->is_real())
          {
            if (curr_bucket->get_bucket_type() == MIXED)
              bnddist = curr_bucket->get_bnddist(pos, junk);
            for (i = 0; i < NEIGH_SIZE; i++)
              if (*(curr_bucket->get_neigh_proc() + i) > -1)
              {
                Key tmpkey = *(curr_bucket->get_neighbors() + i);
                Bucket *neigh = (Bucket *) BG_mesh->lookup(tmpkey);

                if (!(neigh) && (*(curr_bucket->get_neigh_proc() + i) != myid))
                  continue;
                assert(neigh);
                if (neigh->get_bucket_type() == MIXED)
                {
                  double temp = neigh->get_bnddist(pos, junk);

                  if (temp < bnddist)
                    bnddist = temp;
                }
              }
          }
          // check if need to add more particles
          if (add_material_check &&
              abs(fluxsrc[0] - pos[0]) < p_curr->get_smlen())
            add_material = false;

          // if particle crosses boundary -- remove it
          if (bnddist < -TINY)
          {
            P_table->remove(*p_itr);
            delete p_curr;

            lost++;
            continue;
          }
          // check where particle is going, 
          // ... if going anywhere at all
          bool left_curr_bucket = false;

          for (i = 0; i < DIMENSION; i++)
          {
            if (pos[i] < mincrd[i])
            {
              dir[i] = 1;
              left_curr_bucket = true;
            }
            else if (pos[i] >= maxcrd[i])
            {
              dir[i] = 2;
              left_curr_bucket = true;
            }
            else
              dir[i] = 0;
          }

          // determine, to which bucket particle migrating
          // there are 3 special cases, that are worth considering, i.e. 
          //    1) a native particle enters a guest bucket
          //    2) a guest particle  enters native bucket
          //    3) leaves a guest bucket and leaves domain
          if (left_curr_bucket)
          {
            if (!p_curr->is_guest())    // should always have neighbors
            {
              Key neigh_key = curr_bucket->which_neigh(dir);
              Bucket *neigh = (Bucket *) BG_mesh->lookup(neigh_key);

              assert(neigh->contains(pos));

              // if real particle moves into guest bucket
              // delete it
              if (neigh->is_guest())
              {
                P_table->remove(*p_itr);
                delete p_curr;

                continue;
              }
              else              // if neighbor is a native ... 
              {
                // if an empty bucket got particle ...
                //   1 -> turn of have_realp flag
                //   2 -> turn on adapt flag
                if (!neigh->have_real_particles())
                {
                  neigh->set_real_particles(true);
                  adapt = 1;
                }
                neigh->add_particle(*p_itr);
              }
            }
            else if (curr_bucket->which_neigh_proc(dir) == myid)
            {
              //case 2 : guest particle enter native bucket
              Key neigh_key = curr_bucket->which_neigh(dir);
              Bucket *neigh = (Bucket *) BG_mesh->lookup(neigh_key);

              assert(neigh->contains(pos));

              // the particle is coming from partition not 
              // on current proc
              p_curr->put_guest_flag(false);
              neigh->add_particle(*p_itr);
              if (!neigh->have_real_particles())
              {
                neigh->set_real_particles(true);
                adapt = 1;
              }
            }
            else                // a guest particles leaves curent domain
            {
              // case 3: guest particle leave the current domain
              P_table->remove(*p_itr);
              delete p_curr;

              continue;
            }
          }
          else
            my_realp.push_back(*p_itr);
        }
        else
          my_particles.push_back(*p_itr);
      }
      if (!my_realp.empty())
        my_particles.insert(my_particles.end(), my_realp.begin(),
                            my_realp.end());
      else
        curr_bucket->set_real_particles(false);

      // update new list 
      curr_bucket->put_new_plist(my_particles);
    }

  if (add_material)
    timeprops->mat_add_time = true;

  // Update particle lists in the buckets
  it2->reset();
  while (curr_bucket = (Bucket *) it2->next())
  {
    if ((curr_bucket->is_active()) &&
        (curr_bucket->get_bucket_type() != UNDERGROUND))
    {
      curr_bucket->update_particles();
      int numofp = curr_bucket->get_plist().size();

      if (numofp > 75)
      {
        int nfl = 0, ngh = 0;

        vector < Key > pl = curr_bucket->get_plist();
        for (i = 0; i < pl.size(); i++)
        {
          Particle *p = (Particle *) P_table->lookup(pl[i]);

          fprintf(stderr, "%f, %f, %d\n", *(p->get_coords()),
                  *(p->get_coords() + 1), p->is_real());
          if (p->is_real())
            nfl++;
          else
            ngh++;
        }
        Key ck = curr_bucket->getKey();

        fprintf(stderr, "Problem in (%u, %u) bucket.\n \
                        No. of particles exceeded threshold: %d\n", ck.key[0], ck.key[1], numofp);
        fprintf(stderr, "real particles: %d, ghost particles: %d\n", nfl, ngh);
        fprintf(stderr, " Min coords: %e, %e\n", *(curr_bucket->get_mincrd()),
                *(curr_bucket->get_mincrd() + 1));
        fprintf(stderr, " Max coords: %e, %e\n", *(curr_bucket->get_maxcrd()),
                *(curr_bucket->get_maxcrd() + 1));
        fprintf(stderr, " ... at %s: %d\n\n", __FILE__, __LINE__);
        return 13;
      }
    }
  }

  // clean up
  delete itr, it2;

  return adapt;
}
