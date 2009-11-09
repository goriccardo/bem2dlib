#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from bempy2d import circlegeom, srfmatbc, bcondvel, solvephi, \
                    circlefld, calcphifld, fldmatbc, fldmatbcv, \
                    calcvelfld, fieldgrid, calcdphidx, calcsrfvel, \
                    collocation, calcpres, geomwing, normals, \
                    calchalflift, calccp, calchalfcl

from numpy import array, zeros, savetxt, eye, linspace, sqrt, \
                  arctan2, cos, pi, sin, sqrt
from numpy.linalg import norm
from pylab import plot, show, grid, figure, subplot, quiver, xlim, \
                  draw, savefig, ylim, xticks, yticks, arange, title, \
                  xlabel, ylabel, text, legend

def main():
    nelem = 151
    chord = 1.
    thick = 0.2
    xnode = geomwing(nelem, chord, thick); TEat1 = 1
    cpoint = collocation(xnode)
    ntstep = 41
    Ttot = 1
    oar = [10,20,30,40]
    T = linspace(0.,Ttot,ntstep)
    uarr = sqrt(2)*sinspeedprof(ntstep,Ttot)
    dt = Ttot/(ntstep-1.)
    pres ,lift = solvebem(xnode, uarr, TEat1, dt)
    #PLOT SPEED
    figure(1)
    plotvel(uarr, Ttot, oar)
    #PLOT GEOMETRY and PRESSURE
    figure(2)
    spl = subplot(111)
    spl.set_aspect('equal','box')
    for t in oar:
        plot(cpoint[:nelem/2+1,0],pres[:nelem/2+1,t],label=r'$t=%9.2f$' % T[t])
    legend(loc=0)
    plotgeom(xnode)
    title(r'Pressure, $u=\sqrt{2}$')
    grid()
    ylabel(r'$\frac{p-p_\infty}{\rho}$', size=18)
    xlabel(r'$x/c$', size=18)
    xticks(arange(-0.2,1.3,0.2))
    yticks(arange(-0.8,1.3,0.2))
    #Pressure
    figure(3)
    print lift
    plotcl(lift, Ttot)
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
    tarr = linspace(0,pi,NTstep)
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

def plotvel(uarr,Ttot,oar):
    nt = uarr.shape[0]
    T = linspace(0.,Ttot,nt)
    plot(T,uarr[:,0])
    for i in oar:
        print array([uarr[i,0]])
        plot(array([T[i]]),array([uarr[i,0]]),'bo')
        text(array([T[i]-0.1]),array([uarr[i,0]]),r'$t=%9.2f$' % T[i])
    title(r'Speed diagram, $u_\infty=\frac{1}{2}\left(1-\cos\frac{\pi t}{T}\right)$')
    xlabel(r'$t$',size=18)
    ylabel(r'$u_\infty$',size=18)
    grid()

def plotcl(cl, Ttot):
    nt = len(cl)
    plot(linspace(0.,Ttot,nt),cl)
    title('Lift')
    ylabel(r'$L$', size=18)
    xlabel(r'$t$', size=18)
    grid()

if( __name__ == '__main__'):
    main()

