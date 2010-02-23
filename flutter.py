#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bempy2d import ematrixwing

from math import sqrt
from numpy.linalg import inv, eigvals
from pylab import plot, show, grid, figure, xlabel, ylabel
from numpy import zeros, reshape, dot, linspace, array, pi, real, imag

def main():
    # DATI
    a = 3
    Omega = 0.8
    O2 = Omega**2
    xie = 0.1
    xig = 0.5
    ra = sqrt(0.25)
    Nelem = 119
    NWake = (Nelem+1)//2*10
    rt = 0.1

    M = array([[1, xig],[xig, ra**2+2*xig*xie-xie**2]])
    K = array([[O2, O2*xie],[O2*xie, ra**2+O2*xie**2]])

    I = complex(0,1)
    karr = linspace(0.2,1,100)

    z1 = []
    z2 = []
    kg = []
    for k in karr:
        E = ematrixwing(Nelem, NWake, rt, I*k)/(2.*pi)
        #print E
        A = -dot(inv(K),(M + 1./(a*k**2)*E))
        zz = eigvals(A)
        z1.append(zz[0])
        z2.append(zz[1])
    z1 = array(z1)
    z2 = array(z2)
    figure(1)
    plot(karr,real(z1))
    plot(karr,real(z2))
    xlabel('k')
    ylabel('Re(z)')
    grid()
    figure(2)
    plot(karr,imag(z1))
    plot(karr,imag(z2))
    xlabel('k')
    ylabel('Im(z)')
    #figure(3)
    #plot(real(z1),imag(z1))
    #plot(real(z2),imag(z2))
    #xlabel('Re(z)')
    #ylabel('Im(z)')
    grid()
    show()

if __name__ == '__main__':
    main()
