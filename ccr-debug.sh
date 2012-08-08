#! /bin/bash
module purge
module load use.own
module load hdf18-intel-mpi
export CXX=mpicxx
export CC=mpicc
export CPPFLAGS=-DMPICH_IGNORE_CXX_SEEK
make clean
make distclean
./configure --enable-paralel --disable-parallel-IO --enable-debug --with-hdf5=$HDF5 
make && make install
