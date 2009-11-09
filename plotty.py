#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from bempy2d import circlegeom, srfmatbc, bcondvel, solvephi, \
                    circlefld, calcphifld, fldmatbc, fldmatbcv, \
                    calcvelfld, fieldgrid, calcdphidx, calcsrfvel, \
                    collocation, calcpres, geomwing, normals, \
                    calchalflift

def main():
    from numpy import array, zeros, savetxt, eye, linspace, sqrt, \
                      arctan2, cos, pi, sin
    from numpy.linalg import norm
    from pylab import plot, show, grid, figure, subplot, quiver, xlim, \
                      ion, ioff, draw, savefig, ylim
    #Geometry
    figure(1)
    nelem = 151
    R = 1.
    chord = 1.
    thick = 0.2
    ntstep = 20
    Ttot = 1.
    dt = Ttot/(ntstep-1.)
    #xnode = circlegeom(nelem,R); TEat1 = 0
    xnode = geomwing(nelem, chord, thick); TEat1 = 1
    plotgeom(xnode)
    #Boundary conditions
    u = array([-1.,0.])
    uarr = norm(u)*sinspeedprof(ntstep,Ttot)
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
    #ion()
    spl = subplot(111)
    spl.set_aspect('equal','box')
    plotgeom(xnode)
    #Pressure calculation and plotting
    pres, lift = solvebem(xnode, uarr, TEat1, dt)
    for t in xrange(ntstep):
        hru2 = 0.5*norm(uarr[t,:])**2
        if hru2 > 0.001:
            pres[:,t] /= hru2
            lift[t] /= hru2
    #pres = calcpres(xnode,TEat1,dt,u,phisrf,chisrf)
    #plot(cpoint[:,0], cpoint[:,1], 'o')
    #presplot = cpoint+pres*n
    #plot(presplot[:,0],presplot[:,1])
    xlim(-1,chord+1)
    grid()
    ylim(-1,1)
    for t in xrange(ntstep):
        plot(linspace(chord,0,nelem/2+1),pres[:nelem/2+1,t],'g')
        #if t == 0:
        #    line, = plot(linspace(chord,0,nelem/2),pres[:nelem/2,t],'g')
        #else:
        #    line.set_ydata(pres[:nelem/2,t])
        ylim(-1,1)
        presn = abs(pres[:,t].reshape((nelem,1)))*n
        col = zeros(nelem)
        #print pres[:,t]
        for i in xrange(nelem):
            if pres[i,t] > 0.:
                col[i] = 1.
        quiver(cpoint[:,0], cpoint[:,1],presn[:,0],presn[:,1], col, units='dots', width=0.8)
        #if t == 0:
        #    qul = quiver(cpoint[:,0], cpoint[:,1],presn[:,0],presn[:,1], col, units='dots', width=0.8)
        #else:
        #    qul.set_UVC(presn[:,0],presn[:,1], col)
        #draw()
        #savefig('pres%d' % t)
    #ioff()
    #Analitical for circle
    #TH = linspace(pi/nelem,2*pi-pi/nelem,nelem)
    #presan = 0.5*norm(u)**2*(1-4*sin(TH)**2)
    #presplotan = cpoint+presan.reshape((nelem,1))*n
    #plot(presplotan[:,0],presplotan[:,1])
    #-----------LIFT------------
    figure(3)
    print lift
    plot(linspace(0,Ttot,ntstep),lift)
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

def sinspeedprof(NTstep,Ttot):
    from numpy import cos, pi, linspace, zeros
    tarr = linspace(0,2*pi,NTstep)
    uarr = zeros((NTstep,2))
    uarr[:,0] = (1-cos(tarr))/2.
    return uarr

def solvebem(xnode, uarr, TEat1, dt):
    from numpy import zeros
    ne = xnode.shape[0]
    nt = uarr.shape[0]
    B,C = srfmatbc(xnode)
    phi = zeros((ne,nt))
    chi = phi.copy()
    for i in xrange(nt):
        chi[:,i] = bcondvel(xnode, uarr[i,:])
        phi[:,i] = solvephi(B,C,chi[:,i])
    pres = calcpres(xnode,TEat1,dt,uarr,phi,chi)
    lift = calchalflift(xnode,TEat1,dt,uarr,phi,chi)
    return pres, lift

if( __name__ == '__main__'):
    main()

