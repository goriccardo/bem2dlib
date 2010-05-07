!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2dlap
      IMPLICIT NONE
      integer, parameter :: Nup = 100
      integer, parameter :: Nelem = Nup*2+1
      integer, parameter :: Nfreq = 2
      integer, parameter :: NWake = Nup*10
      integer :: i, j
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(Nelem,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1.D0, Thick = 0.12D0, uscalar = 1.D0
      real(kind=8), parameter :: Freq = 0.05D0, VelAmpl = 0.1D0, alphaAmpl = 10.D0
      real(kind=8) :: UHoriz
      real(kind=8) :: alpha = 0.D0
      real(kind=8), dimension(2), parameter :: Xo = (/0.25,0./)
      complex(kind=8), dimension(Nfreq) :: s
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(Nelem,2) :: CPoint
      real(kind=8), dimension(NWake,2) :: XWnode
      complex(kind=8), dimension(Nelem,Nfreq) :: ChiLap, PhiLap
      complex(kind=8), dimension(Nfreq,2) :: Us
!     real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      complex(kind=8), dimension(Nelem,Nfreq,2) :: srfvelLap
      complex(kind=8), dimension(Nelem,2*Nfreq-1) :: presLap, cpLap
      complex(kind=8), dimension(2*Nfreq-1) :: clLap
      real(kind=8), dimension(Nelem, NWake) :: D
!     Time step
      real(kind=8) :: DT, CalcDT
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      logical :: TEat1 = .true.
!     The program starts here!
      !call GeomWing(Nelem, Xnode, Chord, Thick)
      call GeomNACA00xx(NUp, Chord, Thick, Xnode)
      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
      DT = CalcDT(Nelem, Xnode, UHoriz)
      write(*,*) "Timestep = ",DT
!     Zero Frequency
!      call BCondStraightLap(Nelem, Xnode, uscalar, alpha, Nfreq, s, Us, Chilap)
!     Oscillation
      call BCondOscilLap(Nelem, Xnode, uscalar, alpha, VelAmpl, freq, NFreq, s, Us, ChiLap)
!     Rotation
!      call BCondRotLap(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, freq, NFreq, s, Us, ChiLap)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call SolvePhiLap(Nelem, B, C, Nwake, D, Nfreq, s, DT, ChiLap, PhiLap)
      call calcSrfVelLap(Nelem, Xnode, TEat1, Nfreq, phiLap, chiLap, srfvelLap)
      call CalcPresLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, presLap)
      call CalcCpLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, cpLap)
      call CalcClLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, clLap)
      call Collocation(Nelem, Xnode, Cpoint)
      open(unit=30, file='phisurflap')
      open(unit=31, file='chisurflap')
      open(unit=32, file='velsurflap')
      open(unit=33, file='presurflap')
      open(unit=34, file='cllap')
      open(unit=35, file='cplap')
      write(34,1001) clLap
      do i = 1,Nelem
       write(30,1001) PhiLap(i,:)
       write(31,1001) ChiLap(i,:)
       write(33,1001) Cpoint(i,:), presLap(i,:)
       write(35,1001) cpLap(i,:)
      end do
      do i = 1,NFreq
       do j = 1, Nelem
        write(32,1001) Cpoint(j,:), srfVelLap(j,i,:)
       end do
      end do
      close(unit=30)
      close(unit=31)
      close(unit=32)
      close(unit=34)
      close(unit=35)
 1001 FORMAT('',100(F15.8))
END PROGRAM
