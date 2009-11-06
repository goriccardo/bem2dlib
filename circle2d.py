#!/usr/bin/env python
# -*- coding: utf-8 -*-

def main():
  from bempy2d import circlegeom, srfmatbc, bcondvel, solvephi, circlefld, \
                      calcphifld, fldmatbc, fldmatbcv, calcvelfld, fieldgrid
  from numpy import array, zeros, savetxt, eye, linspace
  from pylab import plot, show, grid
  nelem = 100
  radius = 1.
  u = array([-1.,0.])
  xmin, xmax = -4.,  4.
  ymin, ymax = -3.,  3.
  nx, ny = 100, 100
  xysrf = zeros((nelem+1,2))
  xnode = circlegeom(nelem,radius)
  xysrf[:nelem,:] = xnode[:,:2]
  xysrf[nelem,:] = xnode[nelem-1,:2]
  plot(xysrf[:,0],xysrf[:,1])
  B, C = srfmatbc(xnode)
  chisrf = bcondvel(xnode, u)
  phisrf = solvephi(B,C,chisrf,nelem)
  print phisrf
  xfield = fieldgrid(xmin,xmax,nx,ymin,ymax,ny)
  print xfield
#  plot(xfield[:,0],xfield[:,1])
#  show()
#  return 0
  Bf, Cf = fldmatbc(xfield, xnode)
  phifld = calcphifld(phisrf,chisrf,Bf,Cf)

  #Bxf, Byf, Cxf, Cyf = fldmatbcv(xfield,xnode)
  #velfld = calcvelfld(phisrf,chisrf,Bxf,Byf,Cxf,Cyf)

  xcont = linspace(xmin,xmax,nx)
  ycont = linspace(ymin,ymax,ny)
  Z = phifld.reshape((ny,nx))

  plotphicont(xcont,ycont,Z)
#  plotvelfld(xfield,velfld)


def plotphicont(X,Y,Z):
  from pylab import subplot, show, contourf, contour, colorbar, title
  spl = subplot(111)
  cpl = contour(X,Y,Z,100)
  colorbar(cpl)
  spl.set_aspect('equal','box')
  title(r'Field potential $\varphi$')
  show()


def plotvelfld(XY,VF):
  from pylab import quiver, show
  quiver(XY[:,0],XY[:,1],VF[:,0],VF[:,1])
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

