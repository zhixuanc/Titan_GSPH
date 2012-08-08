
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
#include <vector>
#include <cassert>
using namespace std;

#include <hashtab.h>
#include <constants.h>
#include <particle.h>
#include <buckhead.h>
#include <bucket.h>
#include <outforms.h>
#include <hdf5calls.h>

void
write_matlab(int myid, HashTable * P_table, HashTable * BG_mesh,
             TimeProps * timeprops, vector <BucketHead>  & partition_table)
{
    int i, j;
    int Up[3] = {0, 0, 2};
    char file1[25];
    static int icount = 0;
    double mincrd[2], maxcrd[2];
    Bucket * buck2 = NULL;

    sprintf(file1, "buckets%03d%06d.h5", myid, timeprops->step);
    icount++;

    hid_t fp = GH5_fopen_serial (file1, 'w');

    int size = (int) partition_table.size ();
    double * xcoord = new double [size * 4];
    double * ycoord = new double [size * 4];
    int * active = new int [size];
   
    vector <BucketHead>::iterator itr;
    j = 0;
    for (itr = partition_table.begin(); itr != partition_table.end(); itr++)
    {
        Bucket * buck = (Bucket *) BG_mesh->lookup ( (unsigned *)itr->get_buck_head ());
        assert (buck);

        active[j] = 0;
        buck2 = buck;
        while (buck2->get_bucket_type () != OVERGROUND)
        {
            
            buck2 = (Bucket *) BG_mesh->lookup (buck2->which_neigh (Up));
            if ( buck2->has_real_particles ())
            {
                active[j] += 4;
                break;
            }
        }

        for (i = 0; i < 2; i++)
        {
            mincrd[i] = *(buck->get_mincrd()+i);
            maxcrd[i] = *(buck->get_maxcrd()+i);
        }

        // corner 0
        xcoord[4*j] = mincrd[0];
        ycoord[4*j] = mincrd[1];

        // corner 1
        xcoord[4*j +1] = maxcrd[0];
        ycoord[4*j +1] = mincrd[1];

        // corner 2
        xcoord[4*j + 2] = maxcrd[0];
        ycoord[4*j + 2] = maxcrd[1];

        // corner 3
        xcoord[4*j + 3] = mincrd[0];
        ycoord[4*j + 3] = maxcrd[1];

        if ( buck->is_active () )
            active[j] += 1;

        // add 1 for paritcles
        if ( buck->has_ghost_particles ())
            active[j]++;

        // increment the counter
        j++;
    }

    int dims1[2] = {size, 4};
    int dims2[2] = {size, 1};
    GH5_WriteS (fp, "/xcoord", dims1, (void *) xcoord, 0, 0, DOUBLETYPE);
    GH5_WriteS (fp, "/ycoord", dims1, (void *) ycoord, 0, 0, DOUBLETYPE);
    GH5_WriteS (fp, "/active", dims2, (void *) active, 0, 0, INTTYPE);
    GH5_fclose (fp);

    delete [] xcoord;
    delete [] ycoord;
    delete [] active;
    return;
}
