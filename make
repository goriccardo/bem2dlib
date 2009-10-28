#!/usr/bin/env sh
gfortran -O2 -Wall -llapack -c bem2d.f90
gfortran -O2 -Wall -llapack -o circle2d circle2d.f90 bem2d.o
