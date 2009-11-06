#!/usr/bin/env python
# -*- coding: utf-8 -*-

def main():
  from bempy2d import geomwing1, srfmatbc, bcondvel, solvephi, circlefld, \
                      calcphifld, fldmatbc, fldmatbcv, calcvelfld, fieldgrid
  from numpy import array, zeros, savetxt, eye, linspace
  from pylab import plot, show, grid
  nelem = 100
  chord = 1.
  thickness = 0.2
  u = array([-1.,0.])
  xysrf = zeros((nelem+1,2))
  xnode = geomwing1(nelem,chord,thickness)
  xysrf[:nelem,:] = xnode[:,:2]
  xysrf[nelem,:] = xnode[0,:2]
  #print xysrf
  plot(xysrf[:,0],xysrf[:,1])
  """
  B, C = srfmatbc(xnode)
  chisrf = bcondvel(xnode, u)
  phisrf = solvephi(B,C,chisrf,nelem)
  print phisrf
  """
  grid()
  show()

if( __name__ == '__main__'):
  main()
