#!/usr/bin/env sh
gfortran -O2 -Wall -c geomcircle.f90
gfortran -O2 -Wall -c integrals.f90
gfortran -O2 -Wall -llapack -c bem2d.f90
gfortran -O2 -Wall -llapack -o circle2d circle2d.f90 geomcircle.o integrals.o bem2d.o
