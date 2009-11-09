#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

def main():
    from bempy2d import circlegeom, srfmatbc, bcondvel, solvephi, \
                        circlefld, calcphifld, fldmatbc, fldmatbcv, \
                        calcvelfld, fieldgrid, calcdphidx, calcsrfvel, \
                        collocation, calcpres, geomwing, normals
    from numpy import array, zeros, savetxt, eye, linspace, sqrt, \
                      arctan2, cos, pi, sin
    from numpy.linalg import norm
    from pylab import plot, show, grid, figure, subplot, quiver, xlim
    #Geometry
    figure(1)
    nelem = 100
    R = 1.
    chord = 1.
    thick = 0.2
    ntstep = 10
    Ttot = 1.
    dt = Ttot/(ntstep-1.)
    #xnode = circlegeom(nelem,R); TEat1 = 0
    xnode = geomwing(nelem, chord, thick); TEat1 = 1
    plotgeom(xnode)
    #Boundary conditions
    u = array([-1.5,0.])
    chisrf = bcondvel(xnode, u)
    #Integral equation solution
    B, C = srfmatbc(xnode)
    phisrf = solvephi(B,C,chisrf)
    #Surface velocity calculation
    v = calcsrfvel(xnode,TEat1,phisrf, chisrf)
    #for i in xrange(nelem):
    #    v[i,0] += norm(u)
    cpoint = collocation(xnode)
    n = normals(xnode)
    #Pressure figure
    figure(2)
    spl = subplot(111)
    spl.set_aspect('equal','box')
    plotgeom(xnode)
    #Pressure calculation and plotting
    pres = calcpres(xnode,TEat1,dt,u,phisrf,chisrf)
    #plot(cpoint[:,0], cpoint[:,1], 'o')
    presplot = cpoint+pres*n
    #plot(presplot[:,0],presplot[:,1])
    plot(linspace(chord,0,nelem/2),pres[:nelem/2],'g')
    presn = abs(pres)*n
    col = zeros(nelem)
    for i in xrange(nelem):
        if pres[i] > 0.:
            col[i] = 1.
    quiver(cpoint[:,0], cpoint[:,1],presn[:,0],presn[:,1], col, units='dots', width=0.8)
    xlim(-1,chord+1)
    grid()
    #Analitical for circle
    #TH = linspace(pi/nelem,2*pi-pi/nelem,nelem)
    #presan = 0.5*norm(u)**2*(1-4*sin(TH)**2)
    #presplotan = cpoint+presan.reshape((nelem,1))*n
    #plot(presplotan[:,0],presplotan[:,1])
    figure(1)
    #Field definition
    xmin, xmax = -2.,  2.
    ymin, ymax = -1.5,  1.5
    nx, ny = 100, 100
    xfield = fieldgrid(xmin,xmax,nx,ymin,ymax,ny)
    #r = sqrt(xfield[:,0]**2+xfield[:,1]**2)
    #th = arctan2(xfield[:,1], -xfield[:,0])
    Bf, Cf = fldmatbc(xfield, xnode)
    phifld = calcphifld(phisrf,chisrf,Bf,Cf)
    #phifld += r*u[0]*cos(th)
    #for i in xrange(nx*ny):
    #    if r[i] < R:
    #        phifld[i] = 0.
    xcont = linspace(xmin,xmax,nx)
    ycont = linspace(ymin,ymax,ny)
    Z = phifld.reshape((ny,nx))
    plotphicont(xcont,ycont,Z)
    plotvelfld(cpoint, v)
    show()

def plotphicont(X,Y,Z):
    from pylab import subplot, show, contourf, contour, colorbar, title
    spl = subplot(111)
    cpl = contour(X,Y,Z,100)
    colorbar(cpl)
    spl.set_aspect('equal','box')
    title(r'Field potential $\varphi$')

def plotgeom(xnode):
    from pylab import zeros, plot
    nelem = len(xnode)
    xysrf = zeros((nelem+1,2))
    xysrf[:nelem,:] = xnode[:,:]
    xysrf[nelem,:] = xnode[0,:]
    plot(xysrf[:,0],xysrf[:,1],'k')

def plotvelfld(XY,VF):
    from pylab import quiver, show
    quiver(XY[:,0],XY[:,1],VF[:,0],VF[:,1], units='dots', width=2.)

if( __name__ == '__main__'):
  main()
