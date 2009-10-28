!A circle in a _potential_ flow
PROGRAM circle2d
       IMPLICIT NONE
       INTEGER, PARAMETER :: NELEM = 50
!      Vector of nodes global coordinates (x,y)
       REAL, DIMENSION(NELEM,2) :: XNODE
!      Circle radius
       REAL, PARAMETER :: R = 2.
       REAL, DIMENSION(2) :: U = (/-1.,0./)
!      Potential and normal wash on the surface
       REAL, DIMENSION(NELEM) :: phi, chi
       REAL, DIMENSION(NELEM,NELEM) :: B, C
       CALL CIRCLEGEOM(NELEM, XNODE, R)
       CALL SRFMATBC(NELEM, XNODE, B, C)
       CALL BCONDVEL(NELEM, XNODE, U, CHI)
       CALL SOLVEPHI(NELEM, B, C, PHI, CHI)
       CALL SAVEPHI(NELEM, PHI)
END PROGRAM

!Save the phi on the surface
SUBROUTINE SAVEPHI(NELEM, PHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL, DIMENSION(NELEM), INTENT(IN) :: PHI
      OPEN(UNIT=11, FILE="phisurf")
      WRITE(11,1) PHI
      CLOSE(UNIT=11)
   1  FORMAT('',F15.6,2X)
END SUBROUTINE
