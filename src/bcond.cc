
/*
 * =====================================================================================
 *
 *       Filename:  bcond.cc
 *
 *    Description:  
 *
 *        Created:  08/04/2010 02:25:33 PM EDT
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
#include <iostream>
using namespace std;

#include <hashtab.h>
#include <bucket.h>
#include <particle.h>
#include <bnd_image.h>

#include "constants.h"
#include "particler.h"
#include "sph_header.h"

int
apply_bcond(int myid, HashTable * P_table, HashTable * BG_mesh,
            MatProps * matprops, vector < BndImage > *Image_table)
{

  int i, j, k;
  unsigned tmpkey[KEYLENGTH];
  double coord[DIMENSION];
  double refc[DIMENSION], supp;
  double intsct[DIMENSION], bnddist;
  double normal[DIMENSION + 1];
  double uvec2[NO_OF_EQNS];
  double uvec[NO_OF_EQNS], state_vars[NO_OF_EQNS], wnorm, smlen;
  double dx[DIMENSION], s[DIMENSION];
  int Up[DIMENSION] = { 0, 0, 2 };
  Key *neighbors;

  vector < Key > plist;
  vector < Key >::iterator p_itr;
  vector < BndImage >::iterator i_img;

  for (i_img = Image_table->begin(); i_img != Image_table->end(); i_img++)
    if (i_img->buckproc == myid)
    {
      // reflection coordinates
      for (i = 0; i < DIMENSION; i++)
        refc[i] = i_img->coord[i];

      // reset variables
      bnddist = 1.0E10;
      wnorm = 0.;
      for (i = 0; i < NO_OF_EQNS; i++)
        uvec[i] = 0.;

      // get hold of bucket containing the image
      Bucket *buck = (Bucket *) BG_mesh->lookup(i_img->bucket_key);

      assert(buck);
      if (buck->get_bucket_type() == MIXED)
      {
        bnddist = buck->get_bnddist(refc, intsct);
        for (i = 0; i <= DIMENSION; i++)
          buck->get_bnd_normal(refc, normal);
      }
      // search neighbors for particles within 3-h neighborhood
      // go through bucket that has the particle
      plist = buck->get_plist();
      for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
      {
        Particle *pj = (Particle *) P_table->lookup(*p_itr);

        assert(pj);
        if (!pj->is_ghost())
        {
          double h = pj->get_smlen();

          supp = 3 * h;
          for (j = 0; j < DIMENSION; j++)
            dx[j] = refc[j] - *(pj->get_coords() + j);

          // if Particle(j) is in 3-h of reflection ...
          if (in_support(dx, supp))
          {
            for (j = 0; j < DIMENSION; j++)
              s[j] = dx[j] / h;
            for (j = 0; j < NO_OF_EQNS; j++)
              state_vars[j] = *(pj->get_state_vars() + j);
            double w = weight(s, h);
            double mj = pj->get_mass();

            uvec[0] += mj * w;
            uvec[1] += mj * w * state_vars[1] / state_vars[0];
            uvec[2] += mj * w * state_vars[2] / state_vars[0];
            uvec[3] += mj * w * state_vars[3] / state_vars[0];
            wnorm += mj * w / state_vars[0];
          }
        }
      }
      // now search neighboring buckets ...
      neighbors = buck->get_neighbors();
      for (i = 0; i < NEIGH_SIZE; i++)
        if (*(buck->get_neigh_proc() + i) > -1)
        {
          Bucket *buck_neigh = (Bucket *) BG_mesh->lookup(neighbors[i]);

          if (!(buck_neigh) && (*(buck->get_neigh_proc() + i) != myid))
            continue;
          assert(buck_neigh);
          // distance from the boundary
          if (buck_neigh->get_bucket_type() == MIXED)
          {
            double temp = buck_neigh->get_bnddist(refc, intsct);

            if (temp < bnddist)
              buck_neigh->get_bnd_normal(refc, normal);
          }

          // search buckets for real particles in 3h neighborhood
          if (buck_neigh->have_real_particles())
          {
            plist = buck_neigh->get_plist();
            for (p_itr = plist.begin(); p_itr != plist.end(); p_itr++)
            {
              Particle *pj = (Particle *) P_table->lookup(*p_itr);

              assert(pj);
              if (!pj->is_ghost())
              {
                double h = pj->get_smlen();

                supp = 3 * h;
                for (j = 0; j < DIMENSION; j++)
                  dx[j] = refc[j] - *(pj->get_coords() + j);

                // if Particle(j) is in 3-h of reflection ...
                if (in_support(dx, supp))
                {
                  for (j = 0; j < DIMENSION; j++)
                    s[j] = dx[j] / h;
                  for (j = 0; j < NO_OF_EQNS; j++)
                    state_vars[j] = *(pj->get_state_vars() + j);
                  double w = weight(s, h);
                  double mj = pj->get_mass();

                  uvec[0] += mj * w;
                  uvec[1] += mj * w * state_vars[1] / state_vars[0];
                  uvec[2] += mj * w * state_vars[2] / state_vars[0];
                  uvec[3] += mj * w * state_vars[3] / state_vars[0];
                  wnorm += mj * w / state_vars[0];
                }
              }
            }
          }
        }

      if (wnorm > 0)
      {
        for (i = 0; i < NO_OF_EQNS; i++)
          uvec[i] /= wnorm;
        reflect(uvec, uvec2, normal);
      }
      else
      {
        uvec[0] = 1.;
        for (i = 0; i < DIMENSION; i++)
          uvec[i] = 0.;
      }
      if (i_img->partproc == myid)
      {
        Particle *p_ghost = (Particle *) P_table->lookup(i_img->ghost_key);

        assert(p_ghost);
        p_ghost->put_state_vars(uvec2);
        p_ghost->put_update_delayed(false);
      }
      else
      {
        for (i = 0; i < NO_OF_EQNS; i++)
          i_img->state_vars_interp[i] = uvec[i];
      }
    }

  // now update particles which have no reflections
  /*
     HTIterator *p_itr = new HTIterator (P_table);
     Particle *gp;
     while ( gp = (Particle *) p_itr->next() )
     if ( gp->is_ghost () && gp->is_not_updated () )
     {
     for (i=0; i<DIMENSION; i++)
     coord[i] = *(gp->get_coords()+i);

     double hi = gp->get_smlen();
     supp = 3*hi;
     wnorm = 0;
     for (i=0; i<NO_OF_EQNS; i++)
     uvec[i] = 0;

     vector<Key> neighs = gp->get_neighs();
     for (i=0; i<neighs.size(); i++)
     {
     Particle *pj = (Particle *) P_table->lookup(neighs[i]);
     if ( !pj )
     {
     fprintf(stderr,"Problem at particle: %e, %e \n",
     coord[0], coord[1]);
     return 1;
     }

     if ( pj->is_ghost() && (!pj->is_not_updated()) )
     {
     for (j=0; j<DIMENSION; j++)
     {
     dx[j] = coord[j] - *(pj->get_coords()+j);
     s[j]  = dx[j]/hi;
     }
     if ( in_support ( dx, supp ) )
     {
     double mj = pj->get_mass();
     double w = weight(s, hi);
     for (j=0; j<NO_OF_EQNS; j++)
     state_vars[j] = *(pj->get_state_vars()+j);
     uvec[0] += w*mj;
     uvec[1] += mj*state_vars[1]*w/state_vars[0];
     uvec[2] += mj*state_vars[2]*w/state_vars[0];
     wnorm   += mj*w/state_vars[0];
     }
     }
     }
     if ( abs(wnorm) > TINY )
     {
     for (i=0; i<NO_OF_EQNS; i++)
     uvec[i] /= wnorm;
     }
     else
     {
     uvec[0] = 1.0;
     uvec[1] = 0;
     uvec[2] = 0;
     }
     gp->put_state_vars(uvec);
     }
     // clean up
     delete p_itr;
   */

  return 0;
}
