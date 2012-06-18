
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#ifdef HAVE_MPI_H
#  include <mpi.h>
#endif

#include <vector>
#include <iostream>
using namespace std;

#include <constants.h>
#include <hashtab.h>

#include "bucket.h"

#include <particle.h>
#include <bnd_image.h>
#include <exvar.h>


void
search_bnd_images (int myid, HashTable * P_table, HashTable * BG_mesh,
                   vector < BndImage > & Image_table, int reset)
{
  bool found_bucket;
  int i, j;
  int dir[DIMENSION];
  double coord[DIMENSION], tmpc[DIMENSION];
  double refc[DIMENSION], intsct[DIMENSION];
  unsigned ghst_key[KEYLENGTH], buck_key[KEYLENGTH];

  // always reset after repatition
  if (reset)
  {
    Image_table.clear ();
    HTIterator *p_itr = new HTIterator (P_table);
    Particle *p_temp = NULL;

    while ((p_temp = (Particle *) p_itr->next ()))
      p_temp->set_reflection (false);
    delete p_itr;
  }

  Bucket *buck = NULL, *buck2 = NULL;
  HTIterator *itr = new HTIterator (BG_mesh);

  while ((buck = (Bucket *) itr->next ()))
    if (buck->have_ghost_particles () && (!buck->is_guest ()))
    {
      vector < Key > plist = buck->get_plist ();
      vector < Key >::iterator p_itr;
      const int * neigh_proc = buck->get_neigh_proc ();
      const Key * neighbors = buck->get_neighbors ();

      if (plist.size () > 0)
        for (p_itr = plist.begin (); p_itr != plist.end (); p_itr++)
        {
          Particle *p_ghost = (Particle *) P_table->lookup (*p_itr);
          if (p_ghost->is_ghost ())
          {
            // coordinate of ghost particle
            for (i = 0; i < DIMENSION; i++)
              coord[i] = *(p_ghost->get_coords () + i);

            // initalize distance and state-variables
            double bnddist = 1.0E+10;     // very large number

            // if particle is ghost
            if (! p_ghost->have_reflection ())
            {
              // first look in current bucket, 
              // ... if current bucket has boundary information
              double bnddist = 1.0E+10;   // some large value
              if (buck->get_bucket_type () == MIXED)
              {
                buck2 = buck;
                bnddist = buck->get_bnddist (coord, intsct);
              }

              // then check neighboring buckets, and compare
              // ... perpendicular distance
              for (i = 0; i < NEIGH_SIZE; i++)
                if (neigh_proc[i] > -1)
                {
                  Bucket * buck_neigh = (Bucket *) BG_mesh->lookup (neighbors[i]);

                  // skip neighbor if it belongs to foreign process
                  if (!(buck_neigh) && (neigh_proc[i] != myid))
                    continue;

                  // involve only if neighbor bucket is a boundary bucket
                  if (buck_neigh->get_bucket_type () == MIXED)
                  {
                    double tmpdist = buck_neigh->get_bnddist (coord, tmpc);
                    if (tmpdist < bnddist)
                    {
                      bnddist = tmpdist;
                      buck2 = buck_neigh;
                      for (j = 0; j < DIMENSION; j++)
                        intsct[j] = tmpc[j];
                    }
                  }
                }

              /* 
               *  Now, we have found the samllest perpendicular distance
               *  and bucket which contains the nearest boundary,
               *  we can find the mirror image, i.e. a point 
               *  along normal at 2*bnddist from current particle coords
               */
              for (j = 0; j < DIMENSION; j++)
                refc[j] = 2 * intsct[j] - coord[j];

              // find the bucket containing reflection
              do
              {
                found_bucket = true;
                for (j = 0; j < DIMENSION; j++)
                  if (refc[j] < *(buck2->get_mincrd () + j))
                  {
                    dir[j] = 1;
                    found_bucket = false;
                  }
                  else if (refc[j] > *(buck2->get_maxcrd () + j))
                  {
                    dir[j] = 2;
                    found_bucket = false;
                  }
                  else
                    dir[j] = 0;

                if (! found_bucket)
                buck2 = (Bucket *) BG_mesh->lookup (buck2->which_neigh (dir));
              }
              while (! found_bucket);

              for (j = 0; j < KEYLENGTH; j++)
              {
                ghst_key[j] = p_ghost->getKey ().key[j];
                buck_key[j] = buck2->getKey ().key[j];
              }
              BndImage bndimg (buck2->get_myprocess (), myid, buck_key, 
                               ghst_key, refc);
              Image_table.push_back (bndimg);
              p_ghost->set_reflection (true);
            }
          }
        }
    }

  delete itr;
}
