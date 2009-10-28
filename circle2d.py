#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bem2d import circlegeom, srfmatbc, bcondvel, solvephi, circlefld, \
                  calcfield, fldmatbc, fldmatbcv, calcfieldv
from numpy import array, zeros, savetxt, eye
#from enthought.tvtk.api import tvtk

nelem = 20
radius = 1.
u = array([-1.,0.])
nhx, nhy = 100, 75
xnode = circlegeom(nelem,radius)
B, C = srfmatbc(xnode)
chisrf = bcondvel(xnode, u)
phisrf = solvephi(B,C,chisrf,nelem)
A = 0.5*eye(nelem)+C
xfield = circlefld(radius, nhx, radius*4., nhy, radius*3.)
Bx, By, Cx, Cy = fldmatbcv(xfield, xnode)
fldvel = calcfieldv(phisrf, chisrf, Bx, By, Cx, Cy)


"""
xyz = zeros((fldphi.size,3))
xyz[:,:2] = xfield
grid = tvtk.StructuredGrid(dimensions=(nhx*2,nhy*2,1))
grid.points = xyz
grid.point_data.scalars = fldphi
grid.point_data.scalars.name = 'Potential'
w = tvtk.XMLStructuredGridWriter(input=grid, file_name='grid.vts')
w.write()
"""

