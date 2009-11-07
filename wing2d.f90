!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!A circle in a _potential_ flow
PROGRAM wing2d
       IMPLICIT NONE
       INTEGER, PARAMETER :: NELEM = 50
!      Vector of nodes global coordinates (x,y)
       REAL(KIND=8), DIMENSION(NELEM,2) :: XNODE
!      Circle radius
       REAL(KIND=8), PARAMETER :: CHORD = 1.
       REAL(KIND=8), DIMENSION(2) :: U = (/-1.,0./)
       REAL(KIND=8) :: T = 1.
!      Potential and normal wash on the surface
       REAL(KIND=8), DIMENSION(NELEM) :: phi, chi
       REAL(KIND=8), DIMENSION(NELEM,NELEM) :: B, C
!      Time length of the acceleration and end-speed
       REAL(KIND=8), PARAMETER :: TIME = 1.
       REAL(KIND=8), PARAMETER :: VB = 5.
       INTEGER :: TIMESTEP = 10.
!      Pressure
       REAL(KIND=8), DIMENSION(NELEM) :: PRES
!      Field grid size limits and vector
       INTEGER, PARAMETER :: NX = 100, NY = 100
       REAL(KIND=8), PARAMETER :: XMIN = -4., XMAX = 4., YMIN = -3., YMAX = 3.
       REAL(KIND=8), DIMENSION(NY*NY,2) :: XFIELD
!      Field matrices Bf and Cf
       REAL(KIND=8), DIMENSION(NX*NY,NELEM) :: BF, CF
!      Field potential vector
       REAL(KIND=8), DIMENSION(NX*NY) :: PHIF
       INTEGER :: NFIELD
       NFIELD = NX*NY
!      The program starts here!
       CALL GEOMWING(NELEM, XNODE, CHORD, T)
       CALL SRFMATBC(NELEM, XNODE, B, C)
       CALL BCONDVEL(NELEM, XNODE, U, CHI)
       CALL SOLVEPHI(NELEM, B, C, PHI, CHI)
       CALL SAVEPHI(NELEM, PHI, "phisurf")
       CALL CALCPRES(NELEM, XNODE, TIME, TIMESTEP, VB, PHI, CHI, PRES)
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

