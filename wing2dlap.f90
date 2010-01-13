!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2dlap
      IMPLICIT NONE
      integer, parameter :: Nelem = 119
      integer, parameter :: Nfreq = 2
      integer, parameter :: NWake = (Nelem+1)/2*10
      integer :: i
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(Nelem,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1.D0, Thick = 0.1D0, uscalar = 1.D0
      real(kind=8), parameter :: Freq = 0.05D0, VelAmpl = 0.1D0, alphaAmpl = 10.D0
      real(kind=8) :: UHoriz
      real(kind=8) :: alpha = 5.D0
      real(kind=8), dimension(2) :: Xo = (/0.25,0./)
      complex(kind=8), dimension(Nfreq) :: s
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(NWake,2) :: XWnode
      complex(kind=8), dimension(Nelem,Nfreq) :: ChiLap, PhiLap
      complex(kind=8), dimension(Nfreq,2) :: Us
!     real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      complex(kind=8), dimension(Nelem,Nfreq,2) :: srfvelLap
      complex(kind=8), dimension(Nelem,2*Nfreq-1) :: presLap
      complex(kind=8), dimension(2*Nfreq-1) :: clLap
      real(kind=8), dimension(Nelem, NWake) :: D
!     Time step
      real(kind=8) :: DT, CalcDT
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      logical :: TEat1 = .true.
!     The program starts here!
      call GeomWing(Nelem, Xnode, Chord, Thick)
      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
      DT = CalcDT(Nelem, Xnode, UHoriz)
      write(*,*) "Timestep = ",DT
!     Oscillation
      call BCondOscilLap(Nelem, Xnode, uscalar, alpha, VelAmpl, freq, NFreq, s, Us, ChiLap)
!     Rotation
!      call BCondRotLap(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, freq, NFreq, s, Us, ChiLap)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call SolvePhiLap(Nelem, B, C, Nwake, D, Nfreq, s, DT, ChiLap, PhiLap)
      call calcSrfVelLap(Nelem, Xnode, TEat1, Nfreq, phiLap, chiLap, srfvelLap)
      call CalcPresLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, presLap)
      call CalcClLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, clLap)
      open(unit=30, file='phisurflap')
      open(unit=31, file='chisurflap')
      open(unit=32, file='velsurflap')
      open(unit=33, file='presurflap')
      open(unit=34, file='cllap')
      write(34,1001) clLap
      do i = 1,Nelem
       write(30,1001) PhiLap(i,:)
       write(31,1001) ChiLap(i,:)
       write(33,1001) presLap(i,:)
      end do
      do i = 1,NFreq
       write(32,1001) srfVelLap(1,i,:)
      end do
      close(unit=30)
      close(unit=31)
      close(unit=32)
      close(unit=34)
 1001 FORMAT('',100(F15.8))
END PROGRAM
