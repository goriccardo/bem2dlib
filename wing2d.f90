!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      integer, parameter :: Nelem = 9
      integer, parameter :: NTime = 900
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
      real(kind=8), dimension(NTime,2) :: Ut
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(NWake,2) :: XWnode
      real(kind=8), dimension(Nelem,2) :: SrfVelEnd, CPoint
      real(kind=8), dimension(Nelem, NTime) :: Chit, Phit
      real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      real(kind=8), dimension(Nelem, NWake) :: D
      real(kind=8), dimension(NWake, NTime) :: DPhiW
      real(kind=8), dimension(NTime) :: cl
      real(kind=8), dimension(Nelem,NTime) :: cp
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      logical :: TEat1 = .true.
!     Time step
      real(kind=8) :: DT, CalcDT
!     The program starts here!
      call GeomWing(Nelem, Xnode, Chord, Thick)
      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
      DT = CalcDT(Nelem, Xnode, UHoriz)
      write(*,*) "Timestep = ",DT
      write(*,*) "Total time = ", DT*NTime
!     Going straight
!      call BCondStraight(Nelem, Xnode, Uscalar, alpha, Ntime, Ut, Chit)
!     Moving up and down...
!      call BCondOscil(Nelem, Xnode, Uscalar, alpha, VelAmpl, Freq, DT, Ntime, Ut, Chit)
!     Rotating
      call BCondRot(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, Freq, DT, Ntime, Ut, ChiT)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call SolvePhiTime(Nelem, B, C, NWake, D, NTime, Chit, Phit, DPhiW)
      call CalcCl(Nelem, Xnode, TEat1, NTime, DT, Ut, Phit, Chit, cl)
      call CalcCp(Nelem, Xnode, TEat1, NTime, dt, Ut, phit, Chit, cp)
      call CalcSrfVel(Nelem, Xnode, TEat1, Phit(:,1), Chit(:,1), srfvelend)
      call CalcDPhiBody(Nelem, NTime, PhiT, DPhiBody)
      do i = 1,NTime
       open(unit=20, file='cpfile')
       open(unit=21, file='utime')
       open(unit=22, file='chisurf')
       write(20,*) cp(:,i)
       write(21,*) Ut(i,:)
       write(22,*) Chit(:,i)
      end do
      close(unit=22)
      close(unit=21)
      close(unit=20)
      open(unit=16, file='dphi')
      do i = 2,NTime
       write(16,*) DPhiBody(:,i), DPhiW(:,i-1)
      end do
      close(unit=16)
      call Collocation(Nelem, Xnode, Cpoint)
      open(unit=15, file='endvel')
      do i = 1,Nelem
       write(15,*) Cpoint(i,:), srfvelend(i,:)
      end do
      close(unit=15)
!      open(unit=19, file='cpoints')
      do i = 1,Nelem/2
       write(19,*) Cpoint(Nelem/2-i+1,:)
      end do
      do i = 1,NWake
       write(19,*) XWnode(i,:)
      end do
      close(unit=19)
!     Save results
      call SaveCL(NTime, cl)
      call SavePhi(Nelem, NTime, PhiT)
      call SaveXWnode(NWake, XWnode)
END PROGRAM

SUBROUTINE SAVECL(NWake, CL)
      integer, intent(IN) :: NWake
      real(kind=8), dimension(NWake), intent(IN) :: CL
      OPEN(UNIT=14, FILE="cltime")
      WRITE(14,1001) CL
      CLOSE(UNIT=14)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

SUBROUTINE SAVEXWNODE(NWake, XWNODE)
      integer, intent(IN) :: NWake
      integer :: i
      real(kind=8), dimension(NWake,2), intent(IN) :: XWnode
      OPEN(UNIT=12, FILE="xwnode")
      do i = 1,NWake 
       WRITE(12,1001) XWnode(i,:)
      end do
      CLOSE(UNIT=12)
 1001 FORMAT('',2(F15.6),2X)
END SUBROUTINE

!Save the phi on the surface
SUBROUTINE SAVEPHI(N, NTime, PHIT)
      IMPLICIT NONE
      integer, intent(IN) :: N, NTime
      integer :: i
      real(kind=8), dimension(N,NTime), intent(IN) :: PhiT
      OPEN(UNIT=13, FILE="phisurf")
      do i = 1,NTime
       WRITE(13,1001) PhiT(:,i)
      end do
      CLOSE(UNIT=13)
 1001 FORMAT('',1000(F15.6),2X)
END SUBROUTINE


