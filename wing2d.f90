!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      integer, parameter :: Nelem = 9
      integer :: i
!     Vector of nodes global coordinates (x,y)
      REAL(KIND=8), DIMENSION(NELEM,2) :: Xnode
!     Circle radius
      REAL(KIND=8), PARAMETER :: Chord = 1.
      REAL(KIND=8) :: uscalar = 1.
      REAL(KIND=8) :: alpha = 0.
      REAL(KIND=8), DIMENSION(2) :: U
      REAL(KIND=8) :: Thick = 0.3
      INTEGER, PARAMETER :: NTime = 10
!     Potential and normal wash on the surface
      REAL(KIND=8), DIMENSION(Nelem,Nelem) :: B, C
      REAL(KIND=8), dimension(NTime,2) :: XWnode
      real(kind=8), dimension(nelem, ntime) :: Chit, Phit, D
      real(kind=8), dimension(NTime) :: DPHIW
!     Time step
      REAL(KIND=8) :: DT, CalcDT
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
      INTEGER, INTENT(IN) :: NTime
      REAL(KIND=8), DIMENSION(NTime), INTENT(IN) :: DPHIW
      OPEN(UNIT=11, FILE="dphiw")
      WRITE(11,1001) DPHIW
      CLOSE(UNIT=11)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

SUBROUTINE SAVEXWNODE(NTime, XWNODE)
      INTEGER, INTENT(IN) :: NTime
      integer :: i
      REAL(KIND=8), DIMENSION(NTime,2), INTENT(IN) :: XWNODE
      OPEN(UNIT=12, FILE="xwnode")
      do i = 1,NTime 
       WRITE(12,1001) XWNODE(i,:)
      end do
      CLOSE(UNIT=12)
 1001 FORMAT('',2(F15.6),2X)
END SUBROUTINE

!Save the phi on the surface
SUBROUTINE SAVEPHI(N, NTime, PHIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NTime
      integer :: i
      REAL(KIND=8), DIMENSION(N,NTime), INTENT(IN) :: PHIT
      OPEN(UNIT=13, FILE="phisurf")
      do i = 1,NTime
       WRITE(13,1001) PhiT(:,i)
      end do
      CLOSE(UNIT=13)
 1001 FORMAT('',1000(F15.6),2X)
END SUBROUTINE
