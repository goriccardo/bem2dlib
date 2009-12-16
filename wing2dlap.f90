!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2dlap
      IMPLICIT NONE
      integer, parameter :: Nelem = 9
      integer, parameter :: Nfreq = 3
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
!     real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      real(kind=8), dimension(Nelem, NWake) :: D
!     Time step
      real(kind=8) :: DT, CalcDT
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
!     The program starts here!
      call GeomWing(Nelem, Xnode, Chord, Thick)
      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
      DT = CalcDT(Nelem, Xnode, UHoriz)
      write(*,*) "Timestep = ",DT
!     Oscillation
!      call BCondOscilLap(Nelem, Xnode, uscalar, alpha, VelAmpl, freq, NFreq, s, ChiLap)
!     Rotation
      call BCondRotLap(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, freq, NFreq, s, ChiLap)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call SolvePhiLap(Nelem, B, C, Nwake, D, Nfreq, s, DT, ChiLap, PhiLap)
      open(unit=30, file='phisurflap')
      open(unit=31, file='chisurflap')
      do i = 1,Nelem
       write(30,1001) PhiLap(i,:)
       write(31,1001) ChiLap(i,:)
      end do
      close(unit=30)
      close(unit=31)
 1001 FORMAT('',100(F15.8))
END PROGRAM
