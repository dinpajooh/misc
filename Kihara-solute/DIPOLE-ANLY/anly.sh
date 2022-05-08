#!/bin/bash

BD=$PWD
cd $BD

module load intel
make
rm -rf *.o
rm -rf *.mod

./Dipole-anly.e < input



exit
