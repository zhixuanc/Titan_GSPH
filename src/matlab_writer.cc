
/*
 * =====================================================================================
 *
 *       Filename:  matlab_writer.cc
 *
 *    Description:  
 *
 *        Created:  08/18/2010 12:48:11 PM EDT
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

#if (HAVE_CONFIG_H)
#  include <config.h>
#endif

#include <cstdio>
using namespace std;

#include <hashtab.h>
#include <constants.h>
#include <particle.h>
#include <bucket.h>
#include <outforms.h>

void
write_matlab(int myid, HashTable * P_table, HashTable * BG_mesh,
             TimeProps * timeprops)
{
  char file1[20], file2[20];
  static int icount = 0;
  double mincrd[DIMENSION], maxcrd[DIMENSION];

  sprintf(file1, "active%03d%04d.dat", myid, icount);
  sprintf(file2, "prtcls%03d%04d.dat", myid, icount);
  icount++;

  FILE *fp = fopen(file1, "w");
  FILE *f2 = fopen(file2, "w");

  HTIterator *itr = new HTIterator(BG_mesh);
  Bucket *buck;

  while ((buck = (Bucket *) itr->next()))
    if (buck->is_active())
    {
      for (int i = 0; i < DIMENSION; i++)
      {
        mincrd[i] = *(buck->get_mincrd() + i);
        maxcrd[i] = *(buck->get_maxcrd() + i);
      }
      int buckettype = buck->get_bucket_type();

      fprintf(fp, "%e, %e, %e, %e, %d\n", mincrd[0], mincrd[1],
              maxcrd[0], maxcrd[1], buckettype);
      vector < Key > plist = buck->get_plist();
      vector < Key >::iterator ip;
      for (ip = plist.begin(); ip != plist.end(); ip++)
      {
        Particle *p = (Particle *) P_table->lookup(*ip);
        int ghost = (int) p->is_real();
        const double *coord = p->get_coords();
        const double *state_vars = p->get_state_vars();

        fprintf(f2, "%e, %e, ,%d, %e, %u\n", *(coord), *(coord + 1), ghost,
                *(state_vars), p->get_neighs().size());
      }
    }
  fclose(fp);
  fclose(f2);
}
