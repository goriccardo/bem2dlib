!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!A circle in a _potential_ flow
PROGRAM wing2d
      IMPLICIT NONE
      INTEGER, PARAMETER :: NELEM = 51
!     Vector of nodes global coordinates (x,y)
      REAL(KIND=8), DIMENSION(NELEM,2) :: XNODE
!     Circle radius
      REAL(KIND=8), PARAMETER :: CHORD = 1.
      REAL(KIND=8) :: uscalar = 1.
      REAL(KIND=8) :: alpha = 6.
      REAL(KIND=8), DIMENSION(2) :: U
      REAL(KIND=8) :: T = 0.3
      INTEGER, PARAMETER :: TIMESTEP = 10.
!     Potential and normal wash on the surface
      REAL(KIND=8), DIMENSION(NELEM) :: phi, chi
      REAL(KIND=8), DIMENSION(NELEM,NELEM) :: B, C
      REAL(KIND=8), DIMENSION(TIMESTEP,2) :: XWNODE
      real(kind=8), dimension(nelem, timestep) :: phit
      real(kind=8), dimension(TIMESTEP,TIMESTEP) :: DPHIW
!     Time length of the acceleration and end-speed
      REAL(KIND=8), PARAMETER :: TIME = 1.
      REAL(KIND=8), PARAMETER :: VB = 5.
!     The program starts here!
      CALL BODYROTATION(uscalar, alpha, u)
      CALL GEOMWING(NELEM, XNODE, CHORD, T)
      CALL WAKEGRID(NELEM, TIMESTEP, XNODE, USCALAR, XWNODE)
      CALL WAKE(NELEM, TIMESTEP, PHIT, DPHIW)
      CALL SRFMATBC(NELEM, XNODE, B, C)
      CALL BCONDVEL(NELEM, XNODE, U, CHI)
      CALL SOLVEPHI(NELEM, B, C, PHI, CHI)
!     Save results
      CALL SAVEPHI(NELEM, PHI)
      CALL SAVEXWNODE(TIMESTEP, XWNODE)
      CALL SAVEDPHIW(TIMESTEP, DPHIW)
END PROGRAM


SUBROUTINE SAVEDPHIW(TIMESTEP, DPHIW)
      INTEGER, INTENT(IN) :: TIMESTEP
      REAL(KIND=8), DIMENSION(TIMESTEP,TIMESTEP), INTENT(IN) :: DPHIW
      OPEN(UNIT=11, FILE="dphiw")
      WRITE(11,1001) DPHIW
      CLOSE(UNIT=11)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

SUBROUTINE SAVEXWNODE(TIMESTEP, XWNODE)
      INTEGER, INTENT(IN) :: TIMESTEP
      REAL(KIND=8), DIMENSION(TIMESTEP,2), INTENT(IN) :: XWNODE
      OPEN(UNIT=12, FILE="xwnode")
      WRITE(12,1001) XWNODE
      CLOSE(UNIT=12)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE

!Save the phi on the surface
SUBROUTINE SAVEPHI(N, PHI)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: PHI
      OPEN(UNIT=13, FILE="phisurf")
      WRITE(13,1001) PHI
      CLOSE(UNIT=13)
 1001 FORMAT('',F15.6,2X)
END SUBROUTINE
