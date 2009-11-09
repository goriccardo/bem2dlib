#!/usr/bin/env sh
gfortran -O2 -Wall -llapack -o circle2d circle2d.f90 fieldgrid.f90 geomcircle.f90 integrals.f90 bem2d.f90 geom.f90
gfortran -O2 -Wall -llapack -o wing2d wing2d.f90 pressure.f90 geom.f90 fieldgrid.f90 velocity.f90 geomwing.f90 integrals.f90 bem2d.f90
