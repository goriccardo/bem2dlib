!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!A circle in a _potential_ flow
PROGRAM wing2d
       IMPLICIT NONE
       INTEGER, PARAMETER :: NELEM = 51
!      Vector of nodes global coordinates (x,y)
       REAL(KIND=8), DIMENSION(NELEM,2) :: XNODE
!      Circle radius
       REAL(KIND=8), PARAMETER :: CHORD = 1., T = 0.3
       REAL(KIND=8), DIMENSION(2) :: U = (/-1.,0./)
!      Potential and normal wash on the surface
       REAL(KIND=8), DIMENSION(NELEM) :: phi, chi
       REAL(KIND=8), DIMENSION(NELEM,NELEM) :: B, C
!      The program starts here!
       CALL GEOMWING(NELEM, XNODE, CHORD, T)
       CALL SRFMATBC(NELEM, XNODE, B, C)
       CALL BCONDVEL(NELEM, XNODE, U, CHI)
       CALL SOLVEPHI(NELEM, B, C, PHI, CHI)
       CALL SAVEPHI(NELEM, PHI, "phisurf")
END PROGRAM

!Save the phi on the surface
SUBROUTINE SAVEPHI(N, PHI, FNAME)
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: PHI
      CHARACTER :: FNAME*128
      OPEN(UNIT=11, FILE=FNAME)
      WRITE(11,1) PHI
      CLOSE(UNIT=11)
   1  FORMAT('',F15.6,2X)
END SUBROUTINE

