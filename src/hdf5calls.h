/*
 * =====================================================================================
 *
 *       Filename:  hdf5calls.h
 *
 *    Description:  
 *
 *        Created:  03/20/2010 06:34:29 PM EDT
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

#ifndef HDF5CALL_H
#define HDF5CALL_H

#include <hdf5.h>
#include <buckstr.h>

hid_t GH5_fopen (const char *filename, char mode);
hid_t GH5_fopen_parallel (const char *filename, char mode);
inline herr_t GH5_fclose(hid_t fp) { return(H5Fclose(fp)); }
inline herr_t GH5_gclose(hid_t id) { return(H5Gclose(id)); }

int GH5_getsize (hid_t fp, const char *fullpath, hsize_t *dims);
hid_t GH5_gopen (hid_t fp, const char *name, char mode);
hid_t GH5_createdataset (hid_t gid, hid_t dataspace, const char *dsetname, unsigned type);
hid_t GH5_cellstruct ();

// Overload data read functions
herr_t GH5_readdata(hid_t fp, const char *fullpath, double *);
herr_t GH5_readdata(hid_t fp, const char *fullpath, int *);
herr_t GH5_readdata(hid_t fp, const char *fullpath, unsigned *);
herr_t GH5_read_grid_data(hid_t fp, const char *path, BucketStruct *c);


// Overloaded data write functions
herr_t GH5_writedata (hid_t fp, const char *path, int *dims, double *buf);
herr_t GH5_writedata (hid_t fp, const char *path, int *dims, int *buf);
herr_t GH5_writedata (hid_t fp, const char *path, int *dims, unsigned *buf);
herr_t GH5_write_grid_data (hid_t fp, const char *path, int size, BucketStruct *c);

// write parallel data
herr_t GH5_par_writedata (hid_t fp, const char *path, int *dims, double *buf, int, int);
herr_t GH5_par_writedata (hid_t fp, const char *path, int *dims, int *buf, int, int);
herr_t GH5_par_writedata (hid_t fp, const char *path, int *dims, unsigned *buf, int, int);

#endif // HDF5CALL_H
