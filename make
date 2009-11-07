#!/usr/bin/env sh
gfortran -O2 -Wall -c fieldgrid.f90
gfortran -O2 -Wall -c geomcircle.f90
gfortran -O2 -Wall -c geomwing1.f90
gfortran -O2 -Wall -c integrals.f90
gfortran -O2 -Wall -c bem2d.f90
gfortran -O2 -Wall -c calcpres.f90
gfortran -O2 -Wall -llapack -o circle2d circle2d.f90 fieldgrid.o geomcircle.o integrals.o bem2d.o
gfortran -O2 -Wall -llapack -o wing2d wing2d.f90 fieldgrid.o geomwing1.o integrals.o bem2d.o calcpres.o
