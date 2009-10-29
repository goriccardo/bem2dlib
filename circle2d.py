#!/usr/bin/env python
# -*- coding: utf-8 -*-

def main():

  from bempy2d import circlegeom, srfmatbc, bcondvel, solvephi, circlefld, \
                      calcphifld, fldmatbc, fldmatbcv, calcfieldv, fieldgrid
  from numpy import array, zeros, savetxt, eye, linspace

  nelem = 100
  radius = 1.
  u = array([-1.,0.])
  xmin, xmax = -4,4
  ymin, ymax = -3,3
  nx, ny = 500, 500
  xnode = circlegeom(nelem,radius)
  B, C = srfmatbc(xnode)
  chisrf = bcondvel(xnode, u)
  phisrf = solvephi(B,C,chisrf,nelem)
  xfield = fieldgrid(xmin,xmax,nx,ymin,ymax,ny)
  Bf, Cf = fldmatbc(xfield, xnode)
  phifld = calcphifld(phisrf,chisrf,Bf,Cf)

  xcont = linspace(xmin,xmax,nx)
  ycont = linspace(ymin,ymax,ny)
  Z = phifld.reshape((nx,ny)).T

  plotphicont(xcont,ycont,Z)


def plotphicont(X,Y,Z):
  from pylab import subplot, show, contourf, colorbar, title
  spl = subplot(111)
  cpl = contourf(X,Y,Z,255)
  colorbar(cpl)
  spl.set_aspect('equal','box')
  title(r'Field potential $\varphi$')
  show()

if( __name__ == '__main__'):
  main()

#Bx, By, Cx, Cy = fldmatbcv(xfield, xnode)
#fldvel = calcfieldv(phisrf, chisrf, Bx, By, Cx, Cy)


"""
xyz = zeros((fldphi.size,3))
xyz[:,:2] = xfield
grid = tvtk.StructuredGrid(dimensions=(nx*2,ny*2,1))
grid.points = xyz
grid.point_data.scalars = fldphi
grid.point_data.scalars.name = 'Potential'
w = tvtk.XMLStructuredGridWriter(input=grid, file_name='grid.vts')
w.write()
"""

