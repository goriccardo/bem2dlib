!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      integer, parameter :: Nelem = 9
      integer :: i
!     Vector of nodes global coordinates (x,y)
      real(kind=8), dimension(NELEM,2) :: Xnode
!     Circle radius
      real(kind=8), parameter :: Chord = 1.
      real(kind=8) :: uscalar = 1.
      real(kind=8) :: alpha = 6.
      real(kind=8), dimension(2) :: U
      real(kind=8) :: Thick = 0.3
      integer, parameter :: NTime = 10
!     Potential and normal wash on the surface
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(NTime,2) :: XWnode
      real(kind=8), dimension(Nelem, Ntime) :: Chit, Phit, D
      real(kind=8), dimension(NTime) :: DPHIW
!     Time step
      real(kind=8) :: DT, CalcDT
!     The program starts here!
      CALL BodyRotation(Uscalar, alpha, U)
      CALL GeomWing(Nelem, Xnode, Chord, Thick)
      DT = CalcDT(Nelem, Xnode, Uscalar)
      CALL WakeGrid(Nelem, Xnode, Uscalar, DT, NTime, XWnode)
      CALL SrfMatBCD(Nelem, Xnode, NTime, XWnode, B, C, D)
      do i = 1,NTime
       CALL BCondVel(Nelem, XNODE, U, Chit(:,i))
      end do
      CALL SolvePhiTime(Nelem, B, C, NTime, D, Chit, Phit, DPhiW)
!     Save results
      CALL SavePhi(Nelem, NTime, PhiT)
      CALL SaveXWnode(NTime, XWNODE)
      CALL SaveDPhiW(NTime, DPHIW)
END PROGRAM


SUBROUTINE SAVEDPHIW(NTime, DPHIW)
      integer, intent(IN) :: NTime
      real(kind=8), dimension(NTime), intent(IN) :: DPHIW
      OPEN(UNIT=11, FILE="dphiw")
      WRITE(11,1001) DPHIW
      CLOSE(UNIT=11)
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
