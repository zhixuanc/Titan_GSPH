#! /bin/bash
module purge
module load intel-mpi
module load hdf-intel-mpi
export CXX=mpiicpc
export CC=mpiicc
export CPPFLAGS=-DMPICH_IGNORE_CXX_SEEK
make clean
make distclean
./configure --enable-parallel --enable-parallel-IO --with-hdf5=$HDF5 
make && make install
