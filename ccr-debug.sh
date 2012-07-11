#! /bin/bash
module purge
module load hdf/5-1.8.6-impi
export CXX=mpiicpc
export CC=mpiicc
export CPPFLAGS=-DMPICH_IGNORE_CXX_SEEK
make clean
make distclean
./configure --enable-paralel --disable-parallel-IO --enable-debug --with-hdf5=$HDF5 
make && make install
