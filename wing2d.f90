!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      integer, parameter :: Nelem = 9
      integer, parameter :: NTime = 300
      integer, parameter :: NWake = (Nelem+1)/2*10
      integer :: i
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(Nelem,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1., Thick = 0.1
      real(kind=8) :: UHoriz, uscalar = 1., Freq = 0.1, Ampl = 0.15
      real(kind=8) :: alpha = 4.
      real(kind=8), dimension(2) :: U
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
      UHoriz = Uscalar*dcos(alpha/360.D0*PI)
      DT = CalcDT(Nelem, Xnode, UHoriz)
      call BodyMoveSin(NTime, Ampl, Freq, Uscalar, alpha, Ut)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      do i = 1,NTime
       call BCondVel(Nelem, Xnode, U, Chit(:,i))
      end do
      call SolvePhiTime(Nelem, B, C, NWake, D, NTime, Chit, Phit, DPhiW)
      call CalcCl(Nelem, Xnode, TEat1, NTime, DT, Ut, Phit, Chit, cl)
      call CalcCp(Nelem, Xnode, TEat1, NTime, dt, Ut, phit, chit, cp)
      do i = 1,NTime
       open(unit=20, file='cpfile')
       open(unit=21, file='utime')
       write(20,*) cp(:,i)
       write(21,*) Ut(i,:)
      end do
      close(unit=21)
      close(unit=20)
      call CalcSrfVel(Nelem, Xnode, TEat1, Phit(:,1), Chit(:,1), srfvelend)
      call CalcDPhiBody(Nelem, NTime, PhiT, DPhiBody)
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
      open(unit=19, file='cpoints')
      do i = 1,Nelem/2
       write(19,*) Cpoint(Nelem/2-i+1,:)
      end do
      do i = 1,NWake
       write(19,*) XWnode(i,:)
      end do
      close(unit=19)
!     Save results
      call SaveCL(NTime, cl)
      call SavePhi(Nelem, NWake, PhiT)
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


