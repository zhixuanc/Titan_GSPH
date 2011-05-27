#! /bin/bash
. /util/Modules/3.1.6/init/bash
module load hdf/5-1.8.6-impi
module load intel-mpi
export CXX=mpiicpc
export CC=mpiicc
export CPPFLAGS=-DMPICH_IGNORE_CXX_SEEK
make clean
make distclean
./configure --disable-3d --enable-parallel --enable-debug --with-hdf5=$HDF5 
make && make install
