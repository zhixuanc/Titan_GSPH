
/*
 * =====================================================================================
 *
 *       Filename:  mesh_update.cc
 *
 *    Description:  
 *
 *        Created:  08/31/2010 02:26:37 PM EDT
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

#include <cassert>
#include <vector>
using namespace std;

#include <hashtab.h>
#include <properties.h>
#include <bucket.h>
#include <particle.h>

void
update_bgmesh (HashTable * P_table, HashTable * BG_mesh,
               MatProps * matprops, int myid, int *added_ghosts)
{
  int i, j;
  double dx = matprops->smoothing_length;
  Bucket *curr_bucket, *neigh;

  *added_ghosts = 0;
  HTIterator *itr = new HTIterator (BG_mesh);

  itr = new HTIterator (BG_mesh);
  while ((curr_bucket = static_cast < Bucket * >(itr->next ())))
    if ((curr_bucket->get_bucket_type () == MIXED) &&
        (curr_bucket->is_active ()) &&
        (!curr_bucket->is_guest ()) && (!curr_bucket->have_ghost_particles ()))
    {
      Key *neighbors = curr_bucket->get_neighbors ();
      bool put_ghosts = false;

      for (i = 0; i < NEIGH_SIZE; i++)
        if (*(curr_bucket->get_neigh_proc () + i) == myid)
        {
          Bucket *buck_neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);

          if (buck_neigh->have_real_particles ())
            put_ghosts = true;
        }
        else if (*(curr_bucket->get_neigh_proc () + i) > -1)
          put_ghosts = true;

      // put ghost particles if any neighbor has real particles
      if (put_ghosts)
      {
        int status =
          curr_bucket->put_ghost_particles (P_table, BG_mesh, matprops);
        *added_ghosts = 1;

        // make sure not sure to mark current bucket visited even if 
        // the area between bucket and bounary is too small to hold 
        // even a single ghost particle. This will avoid duplicating ghosts
        // during next adaptation
        curr_bucket->set_ghost_particles (true);
      }
    }

  // visit every bucket 
  itr->reset ();
  while ((curr_bucket = static_cast < Bucket * >(itr->next ())))
  {
    // mark it inactive to start with
    curr_bucket->mark_inactive ();

    // if any neighbor as any real particle, mark current bucket active
    Key *neigh_buckets = curr_bucket->get_neighbors ();

    for (i = 0; i < NEIGH_SIZE; i++)
      if (*(curr_bucket->get_neigh_proc () + i) == myid)
      {
        neigh = (Bucket *) BG_mesh->lookup (neigh_buckets[i]);
        if (neigh->get_plist ().size () > 0)
          curr_bucket->mark_active ();
      }
    // if bucket is still inactive, delete any ghost particles it has
    if (!curr_bucket->is_active ())
    {
      vector < Key > plist = curr_bucket->get_plist ();
      vector < Key >::iterator ip;
      if (plist.size () > 0)
      {
        for (ip = plist.begin (); ip != plist.end (); ip++)
        {
          Particle *p_del = (Particle *) P_table->lookup (*ip);

          P_table->remove (p_del->getKey ());
          delete p_del;
        }
        curr_bucket->empty_plist ();
      }
    }
  }

  delete itr;

  return;
}
