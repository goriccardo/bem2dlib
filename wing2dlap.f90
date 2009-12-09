!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2dlab
      IMPLICIT NONE
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      integer, parameter :: Nelem = 9
      integer, parameter :: NWake = (Nelem+1)/2*10
      integer :: i
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(Nelem,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1., Thick = 0.1D0, uscalar = 1., Freq = 0.1D0, Ampl = 0.15D0
      real(kind=8) :: UHoriz
      real(kind=8) :: alpha = 0.
      real(kind=8), parameter :: w = 2.*PI*Freq
      complex(kind=8), parameter :: s = dcmplx(0.,Freq)
      real(kind=8), dimension(2) :: U
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(NWake,2) :: XWnode
      complex(kind=8), dimension(Nelem) :: ChiLap, PhiLap
      complex(kind=8), dimension(Nelem,Nelem) :: DRS
!     real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      real(kind=8), dimension(Nelem, NWake) :: D
!     Time step
      real(kind=8) :: DT, CalcDT
!     The program starts here!
      call GeomWing(Nelem, Xnode, Chord, Thick)
      UHoriz = Uscalar*dcos(alpha/360.D0*PI)
      DT = CalcDT(Nelem, Xnode, UHoriz)
      call BodyRotation(Uscalar, alpha, U)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call MatDRS(Nelem, NWake, D, DT, s, DRS)
      call BCondLap(Nelem, Xnode, Ampl, ChiLap)
      call SolvePhiLap(Nelem, B, C, DRS, ChiLap, PhiLap)
      do i = 1,Nelem
       write(*,1001) cdabs(PhiLap(i))
      end do
 1001 FORMAT('',F15.6)
END PROGRAM
