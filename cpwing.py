#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from bempy2d import srfmatbc, bcondvel, solvephi, collocation, \
                    calcpres, geomwing, normals, calchalfcl, calccp
from numpy import array, zeros, savetxt, eye, linspace, sqrt, \
                  arctan2, cos, pi, sin
from numpy.linalg import norm
from pylab import plot, show, grid, figure, subplot, quiver, xlim, \
                  draw, savefig, ylim, xticks, yticks, arange, title, \
                  xlabel, ylabel, text

def main():
    nelem = 151
    chord = 1.
    thick = 0.2
    u = array([-15.,0.])
    plotcp(nelem, chord, thick, u)
    grid()
    show()

def plotcp(nelem, chord, thick, u):
    #Geometry
    xnode = geomwing(nelem, chord, thick); TEat1 = 1
    dt = 1.
    #Boundary conditions
    chisrf = bcondvel(xnode, u)
    #Integral equation solution
    B, C = srfmatbc(xnode)
    phisrf = solvephi(B,C,chisrf)
    cpoint = collocation(xnode)
    #Pressure figure
    spl = subplot(111)
    spl.set_aspect('equal','box')
    plotgeom(xnode)
    #Pressure calculation and plotting
    U = u.reshape((1,2))
    cp = calccp(xnode,TEat1,dt,U,phisrf,chisrf)
    cl = calchalfcl(xnode,TEat1,dt,U,phisrf,chisrf)
    print cl
    plot(cpoint[:nelem/2+1,0],cp[:nelem/2+1,0])
    title(r'Stationary pressure coefficient, $c_L = %9.3f$' % cl[0])
    ylabel(r'$c_p$', size=18)
    xlabel(r'$x/c$', size=18)
    xticks(arange(-0.2,1.3,0.2))
    yticks(arange(-0.8,1.3,0.2))


def plotgeom(xnode):
    from pylab import zeros, plot
    nelem = len(xnode)
    xysrf = zeros((nelem+1,2))
    xysrf[:nelem,:] = xnode[:,:]
    xysrf[nelem,:] = xnode[0,:]
    plot(xysrf[:,0],xysrf[:,1],'k')

if( __name__ == '__main__'):
    main()
