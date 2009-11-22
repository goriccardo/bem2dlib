!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      integer, parameter :: Nelem = 29
      integer, parameter :: NTime = 50
      integer :: i
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(NELEM,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1.
      real(kind=8) :: uscalar = 1.
      real(kind=8) :: alpha = 6.
      real(kind=8), dimension(2) :: U
      real(kind=8), dimension(NTime,2) :: Ut
      real(kind=8) :: Thick = 0.3
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(NTime,2) :: XWnode
      real(kind=8), dimension(Nelem, Ntime) :: Chit, Phit, D
      real(kind=8), dimension(NTime) :: DPHIW
      real(kind=8), dimension(NTime) :: cl
      logical :: TEat1 = .true.
!     Time step
      real(kind=8) :: DT, CalcDT
!     The program starts here!
      call BodyRotation(Uscalar, alpha, U)
      Ut(:,1) = U(1)
      Ut(:,2) = U(2)
      call GeomWing(Nelem, Xnode, Chord, Thick)
      DT = CalcDT(Nelem, Xnode, Uscalar)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NTime, XWnode)
      call SrfMatBCD(Nelem, Xnode, NTime, XWnode, B, C, D)
      do i = 1,NTime
       call BCondVel(Nelem, XNODE, U, Chit(:,i))
      end do
      call SolvePhiTime(Nelem, B, C, NTime, D, Chit, Phit, DPhiW)
      call CalcCl(Nelem, Xnode, TEat1, NTime, DT, Ut, Phit, Chit, cl)
!     Save results
      call SaveCL(NTime, cl)
      call SavePhi(Nelem, NTime, PhiT)
      call SaveXWnode(NTime, XWNODE)
      call SaveDPhiW(NTime, DPHIW)
END PROGRAM


SUBROUTINE SAVEDPHIW(NTime, DPHIW)
      integer, intent(IN) :: NTime
      real(kind=8), dimension(NTime), intent(IN) :: DPHIW
      OPEN(UNIT=11, FILE="dphiw")
      WRITE(11,1001) DPHIW
      CLOSE(UNIT=11)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

SUBROUTINE SAVECL(NTime, CL)
      integer, intent(IN) :: NTime
      real(kind=8), dimension(NTime), intent(IN) :: CL
      OPEN(UNIT=14, FILE="cltime")
      WRITE(14,1001) CL
      CLOSE(UNIT=14)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

SUBROUTINE SAVEXWNODE(NTime, XWNODE)
      integer, intent(IN) :: NTime
      integer :: i
      real(kind=8), dimension(NTime,2), intent(IN) :: XWnode
      OPEN(UNIT=12, FILE="xwnode")
      do i = 1,NTime 
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
      real(kind=8), dimension(N,NTime), intent(IN) :: PHIT
      OPEN(UNIT=13, FILE="phisurf")
      do i = 1,NTime
       WRITE(13,1001) PhiT(:,i)
      end do
      CLOSE(UNIT=13)
 1001 FORMAT('',1000(F15.6),2X)
END SUBROUTINE


