!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!The formula is y = sqrt(x)*(1-x)
!  NELEM    Figure of elements (odd) (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  T        Thickness                (IN)
!  C        Chord                    (IN)
!It use an iterative method to make the elements of the same lenght
SUBROUTINE GEOMWING(NELEM, XNODE, C, T)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL(KIND=8), INTENT(IN) :: T, C
      REAL(KIND=8), DIMENSION(NELEM,2), INTENT(OUT) :: XNODE
      REAL(KIND=8), DIMENSION(NELEM/2) :: DS, DX, PHI
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I, J, NITER = 15
      REAL(KIND=8) :: X, Y, SSQ, SIGMA, MEAN !, THETA
      !First and mid nodes
      XNODE(1,1) = C
      XNODE(1,2) = 0.
      XNODE(NELEM/2+1,1) = 0.
      XNODE(NELEM/2+1,2) = 0.
      XNODE(NELEM/2+2,1) = 0.
      XNODE(NELEM/2+2,2) = 0.
      !Symmetric profile
      DX(:) = C/REAL(NELEM/2,8)
      DO I = 1,NITER
       DO J = 1,NELEM/2
        IF (J .LT. NELEM/2) THEN
         X = C-SUM(DX(:J))
         XNODE(J+1,1) = X
         Y = T / DSQRT(16.D0/27.D0) * DSQRT(X/C)*(1.-X/C)
         XNODE(J+1,2) = Y
         XNODE(NELEM-J+1,1) = X
         XNODE(NELEM-J+1,2) = -Y
        END IF
        DS(J) = DSQRT( (XNODE(J+1,1)- XNODE(J,1))**2 + (XNODE(J+1,2)-XNODE(J,2))**2 )
        PHI(J) = DATAN2( XNODE(J,2) - XNODE(J+1,2) , XNODE(J,1) - XNODE(J+1,1) )
       END DO
       !Mean and stdev
       MEAN = SUM(DS)/DFLOAT(NELEM/2)
       SSQ = 0.
       DO J = 1,NELEM/2
        SSQ = SSQ + (DS(J)-MEAN)**2
       END DO
       SIGMA = DSQRT( 1./DFLOAT(NELEM/2-1) * SSQ )
       !Iteration
       DO J = 0,NELEM/2-2
        DX(NELEM/2-J) = MEAN*DCOS(PHI(NELEM/2-J))
       END DO
       DX(1) = C - SUM(DX(2:NELEM/2))
      END DO
      XNODE(NELEM/2+1,2) = 0.01*MEAN/2.
      XNODE(NELEM/2+2,2) = -0.01*MEAN/2.
END SUBROUTINE


subroutine geomNACA00xx(Nelem, c, t, Xnode)
      implicit none
      integer, intent(IN) :: Nelem
      real(kind=8), intent(IN) :: c, t
      real(kind=8), dimension(Nelem*2+1,2), intent(OUT) :: Xnode
      integer, parameter :: Ncoef = 5
      real(kind=8) :: Ka = 5.D0
      real(kind=8), dimension(Ncoef) :: Coef = (/0.2969D0,-0.1260D0,-0.3516D0,0.2843D0,-0.1015D0/)
      real(kind=8), dimension(Ncoef) :: ECoef = (/0.5D0,1.D0,2.D0,3D0,4.D0/)
      call geomPolySymWing(Ncoef, Ka, Coef, Ecoef, Nelem, c, t, Xnode) 
end subroutine


!Four digit NACA profile:
! c is the chord lenght
! m is the maximum camber (100 m is the first digit)
! p is the location of maximum camber (10 p is the second digit)
! t is the thickness (100 t are the last two digits)
!All values are relative to the chord
!Example: NACA 2312 with chord = 1 is
!  call geomNACAxxx(Nelem, 1, 0.02, 0.3, 0.12, Xnode)
! subroutine geomNACAxxxx(Nelem, c, m, p, t, Xnode)
!       implicit none
!       integer, intent(IN) :: Nelem
!       real(kind=8), intent(IN) :: c, m, p, t
!       real(kind=8), dimension(Nelem*2+1,2), intent(OUT) :: Xnode
!       integer, parameter :: Ncoef = 5
!       real(kind=8) :: Ka = 5.D0
!       real(kind=8), dimension(Ncoef) :: Coef = (/0.2969D0,-0.1260D0,-0.3516D0,0.2843D0,-0.1015D0/)
!       real(kind=8), dimension(Ncoef) :: ECoef = (/0.5D0,1.D0,2.D0,3D0,4.D0/)
!       
! end subroutine


!Create the mesh for a symmetrical wing described by a polynomial
!       __
!      \
! T*Ka*/__ coef(i)*(x/C)^ecoef(i)
!    i=1..Ncoef
!
! Ncoef: number of coefficient for the poly
! Ka:    Coefficient
! Coef:  coefficients
! Ecoef: exponents of coefficients
! Nelem: Number of segments for a side
! c:     Chord
! t:     Thickness
! Xnode: 2*Nelem+1 mesh
subroutine geomPolySymWing(Ncoef, Ka, Coef, Ecoef, Nelem, c, t, Xnode)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM, Ncoef
      REAL(KIND=8), INTENT(IN) :: Ka, T, C
      REAL(KIND=8), DIMENSION(NELEM*2+1,2), INTENT(OUT) :: XNODE
      REAL(KIND=8), DIMENSION(NELEM+1) :: DS, DX, PHI
      REAL(KIND=8), PARAMETER :: PI = 4.D0*DATAN(1.D0)
      real(kind=8), dimension(Ncoef), intent(in) :: Coef
      real(kind=8), dimension(Ncoef), intent(in) :: ECoef
      INTEGER :: I, J, NITER = 15
      REAL(KIND=8) :: X, Y, SSQ, SIGMA, MEAN, XOC
      !First and mid nodes
      XNODE(1,1) = C
      XNODE(1,2) = 0.
      XNODE(NELEM+1,1) = 0.
      XNODE(NELEM+1,2) = 0.
      XNODE(NELEM+2,1) = 0.
      XNODE(NELEM+2,2) = 0.
      !Symmetric profile
      DX(:) = C/DBLE(NELEM)
      DO I = 1,NITER
       DO J = 1,NELEM
        IF (J .LT. NELEM) THEN
         X = C-SUM(DX(:J))
         XNODE(J+1,1) = X
         xoc = X/C
         Y = t*ka*sum(coef*xoc**ecoef)
         XNODE(J+1,2) = Y
         XNODE(2*NELEM-J+2,1) = X
         XNODE(2*NELEM-J+2,2) = -Y
        END IF
        DS(J) = DSQRT( (XNODE(J+1,1)- XNODE(J,1))**2 + (XNODE(J+1,2)-XNODE(J,2))**2 )
        PHI(J) = DATAN2( XNODE(J,2) - XNODE(J+1,2) , XNODE(J,1) - XNODE(J+1,1) )
       END DO
       !Mean and stdev
       MEAN = SUM(DS)/DFLOAT(NELEM)
       SSQ = 0.
       DO J = 1,NELEM
        SSQ = SSQ + (DS(J)-MEAN)**2
       END DO
       SIGMA = DSQRT( 1./DFLOAT(NELEM-1) * SSQ )
       !Iteration
       DO J = 0,NELEM-2
        DX(NELEM-J) = MEAN*DCOS(PHI(NELEM-J))
       END DO
       DX(1) = C - SUM(DX(2:NELEM))
      END DO
      XNODE(NELEM+1,2) = 0.01*MEAN/2.
      XNODE(NELEM+2,2) = -0.01*MEAN/2.
end subroutine

