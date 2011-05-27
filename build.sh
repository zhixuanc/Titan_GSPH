#! /bin/bash

make clean
make distclean
./configure
make && make install
