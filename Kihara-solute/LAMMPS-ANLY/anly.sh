#!/bin/bash

BD=$PWD
cd $BD

lammps=$HOME/mylammps/build/lmp

$lammps < in.anly > OUT-anly

exit
