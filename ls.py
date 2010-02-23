#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bempy2d import ematrixwing

from math import sqrt
from numpy import zeros, reshape, dot, linspace, array, pi
from numpy.linalg import inv
from scipy.optimize import leastsq

def func(emr,parr,E):
    # em = coefficienti [ E0 E1 E2 EB EA EC ]
    assert len(emr) == 44
    assert len(E) == len(parr)*4

    emr = reshape(emr,(22,2))
    em = map(lambda x: complex(*x), emr)
    assert len(em) == 22
    E_B = reshape(em[12:16],(2,2))
    E_C = reshape(em[18:],(2,2))
    E_BAC = zeros(4*len(parr))
    for i in xrange(len(parr)):
        p = parr[i]
        E_A = zeros((2,2))
        E_A[0,0] = p+em[16]
        E_A[1,1] = p+em[17]
        E_A = inv(E_A)
        E_BAC[i*4:(i+1)*4] = dot(E_B,dot(E_A,E_C)).ravel()
        for j in xrange(4):
            E_BAC[i*4+j] += em[j] + p*em[j+4] + p**2*em[j+8]
    ret = E/(2.*pi) - E_BAC
    ret = array(map(toarr,ret)).ravel()

    assert len(ret) == len(parr)*8
    return ret

def toarr(c):
    return array([c.real,c.imag])

def Elong(Nelem,NWake,rt,parr):
    E = []
    for p in parr:
        E.extend(ematrixwing(Nelem, NWake, rt, p).ravel().tolist())
    return E

def singles(E):
    E0 = reshape(E[0:4],(2,2))
    E1 = reshape(E[4:8],(2,2))
    E2 = reshape(E[8:12],(2,2))
    EB = reshape(E[12:16],(2,2))
    EA = zeros((2,2),dtype=complex)
    EA[0,0] = E[16]
    EA[1,1] = E[17]
    EC = reshape(E[18:22],(2,2))
    return E0, E1, E2, EB, EA, EC

def main():
    I = complex(0,1)
    Nelem = 119
    NWake = (Nelem+1)//2*10
    rt = 0.1
    parr = I*linspace(0.1,50,100)
    E = Elong(Nelem,NWake,rt,parr)
    em0 = zeros(44)+0.5
    #print func(em0,parr,E)
    Eout, success = leastsq(func, em0, args=(parr,E), maxfev=10000)
    print "SUCCESS = ", success
    Eout.resize((len(Eout)/2,2))
    Eout = map(lambda x: complex(*x), Eout)
   
    E0, E1, E2, EB, EA, EC = singles(Eout)
    print singles(Eout)
    
    # DATI
    a = 3
    Omega = 0.8
    xie = 0.1
    xig = 0.5
    ra = sqrt(0.25)
    
    O2 = Omega**2
    M = array([[1, xig],[xig, ra**2+2*xig*xie-xie**2]])
    K = array([[O2, O2*xie],[O2*xie, ra**2+O2*xie**2]])

if __name__ == '__main__':
    main()
