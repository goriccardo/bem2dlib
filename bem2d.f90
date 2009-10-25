!Calculate the position of the nodes
!  NELEM    Figure of elements  (IN)
!  XNODE    Nodes vector        (OUT)
!  R        Radius              (IN)
SUBROUTINE CIRCLEGEOM(NELEM, XNODE, R)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL, INTENT(IN) :: R
      REAL, DIMENSION(NELEM,2), INTENT(OUT) :: XNODE
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL :: DTH
      DTH = 2.*PI/REAL(NELEM)
      DO I = 1,NELEM
       XNODE(I,1) = R*COS(DTH*REAL(I-1))
       XNODE(I,2) = R*SIN(DTH*REAL(I-1))
      END DO
END SUBROUTINE


!Gives the inode for the ielem
!  IDXNODE must be 1 or 2 (first or second node of the element)
INTEGER FUNCTION NODOFEL(NELEM, IELEM, IDXNODE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM, IELEM, IDXNODE
      IF (IDXNODE .EQ. 1) THEN
        NODOFEL = IELEM
        RETURN
      END IF
      IF (IELEM .EQ. NELEM) THEN
        NODOFEL = 1
        RETURN
      END IF
      NODOFEL = IELEM + 1
END FUNCTION


!The integral of G(XI,XJ) from XJM to XJP, nodes of the J element
REAL FUNCTION BIJ(XI, XJM, XJP)
      IMPLICIT NONE
      REAL, DIMENSION(2) :: XI, XJM, XJP, EX, VXP, VXM
      REAL :: XM, YM, XP, YP, RP, RM
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      EX = XJP - XJM
      EX = EX/SQRT(EX(1)**2 + EX(2)**2)
      VXP = (XJP - XI)
      VXM = (XJM - XI)
      XP = VXP(1)*EX(1) + VXP(2)*EX(2)
      XM = VXM(1)*EX(1) + VXM(2)*EX(2)
!     YP and YM should be equals
      YP =-VXP(1)*EX(2) + VXP(2)*EX(1)
      YM =-VXM(1)*EX(2) + VXM(2)*EX(1)
      RP = SQRT(XP**2 + YP**2)
      RM = SQRT(XM**2 + YM**2)
      BIJ = (XP*LOG(RP) - XP + YP*ATAN2(XP,YP) - XM*LOG(RM) + XM - YM*ATAN2(XM,YM))/(2.*PI)
      RETURN
END FUNCTION


!The integral of dG/dn(XS,XJ) from XJM to XJP, nodes of the J element
REAL FUNCTION CIJ(XS, XJM, XJP)
      IMPLICIT NONE
      REAL, DIMENSION(2) :: XS, XJM, XJP, EX, VXP, VXM
      REAL :: XM, XP, YS
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      EX = XJP - XJM
      EX = EX/SQRT(EX(1)**2 + EX(2)**2)
!     Vector X_+
      VXP = (XS - XJP)
!     Vector X_-
      VXM = (XS - XJM)
      XP = -VXP(1)*EX(1) - VXP(2)*EX(2)
      XM = -VXM(1)*EX(1) - VXM(2)*EX(2)
      YS = -VXP(1)*EX(2) + VXP(2)*EX(1)
      CIJ = (ATAN2(XP,YS) - ATAN2(XM,YS))/(2.*PI)
      RETURN
END FUNCTION


!Calculate the matrix B and C
!  N        Figure of elements  (IN)
!  B        Matrix B            (OUT)
!  C        Matrix C            (OUT)
!  XNODE    Nodes vector        (IN)
SUBROUTINE SRFMATBC(N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL, DIMENSION(N,N), INTENT(OUT) :: B, C
      REAL, DIMENSION(2) :: XI
      REAL :: CIJ, BIJ
      INTEGER :: I, J, NODOFEL
      DO I = 1, N
       !Midpoint of I element
       XI = (XNODE(NODOFEL(N,I,1),:) + XNODE(NODOFEL(N,I,2),:))/2.
       C(I,I) = 0.
       DO J = 1, N
        B(I,J) = BIJ(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:))
        IF (I .NE. J) THEN
          C(I,J) = CIJ(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:))
        END IF
       END DO
      END DO
END SUBROUTINE


!Impose boundary conditions for a moving circle with u velocity
!in the body frame of reference
SUBROUTINE BCONDVEL(N, XNODE, U, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL, DIMENSION(2), INTENT(IN) :: U
      REAL, DIMENSION(N), INTENT(OUT) :: CHI
      REAL, DIMENSION(2) :: T
      INTEGER :: I, NODOFEL
      DO I = 1, N
       T = XNODE(NODOFEL(N,I,2),:) - XNODE(NODOFEL(N,I,1),:)
       T = T/SQRT(T(1)**2 + T(2)**2)
       CHI(I) = U(1)*T(2) - U(2)*T(1)
      END DO
END SUBROUTINE


!Solve the (NxN) linear system (0.5*I-B)*phi = C*chi
SUBROUTINE SOLVEPHI(N, B, C, PHI, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N,N), INTENT(IN) :: B, C
      REAL, DIMENSION(N), INTENT(OUT) :: PHI
      REAL, DIMENSION(N), INTENT(IN) :: CHI
      INTEGER :: I, INFO
      REAL, DIMENSION(N,N) :: A
      REAL, DIMENSION(N) :: RHS, IPIV
      A = C
      DO I = 1, N
       A(I,I) = A(I,I) + 0.5
      END DO
      RHS = MATMUL(B,CHI)
      CALL SGESV(N, 1, A, N, IPIV, RHS, N, INFO)
      IF (INFO .NE. 0) THEN
        WRITE(*,*) "ERROR IN LINEAR SYSTEM"
      END IF
      PHI = RHS
END SUBROUTINE


!Calculate B and C matrices on the field points X
SUBROUTINE FLDMATBC(NX, X, N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL, DIMENSION(NX,2), INTENT(IN) :: X
      REAL, DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL, DIMENSION(NX,N), INTENT(OUT) :: B, C
      REAL, DIMENSION(2) :: XI
      REAL :: BIJ, CIJ
      INTEGER :: I, J, NODOFEL
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        B(I,J) = BIJ(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:))
        C(I,J) = CIJ(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:))
       END DO
      END DO
END SUBROUTINE


!Calculate the potential in the field
!  PHISRF   On surface                          (in)
!  CHISRF   On surface                          (in)
!  B, C     In the field                        (in)
!  PHIFLD   In the field                        (out)
SUBROUTINE CALCFIELD(N, PHISRF, CHISRF, NX, B, C, PHIFLD)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NX
      REAL, DIMENSION(N), INTENT(IN) :: PHISRF, CHISRF
      REAL, DIMENSION(NX,N) :: B, C
      REAL, DIMENSION(NX), INTENT(OUT) :: PHIFLD
      PHIFLD = MATMUL(C, PHISRF) + MATMUL(B, CHISRF)
END SUBROUTINE


!  R        Circle radius
!  NHX      Half of points on X figure
!  DELTAX   Size on X
!  NHY      Half of points on Y figure
!  DELTAY   Size on Y
SUBROUTINE CIRCLEFLD(R, NHX, DELTAX, NHY, DELTAY, XFIELD)
      IMPLICIT NONE
      REAL, INTENT(IN) :: R, DELTAX, DELTAY
      INTEGER, INTENT(IN) :: NHX, NHY
      REAL, DIMENSION(NHX*NHY*4,2), INTENT(OUT) :: XFIELD
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      REAL :: X, Y, DX, DY
      INTEGER :: IY, IX, NX, NY
      NX = NHX*2
      NY = NHY*2
      X = -DELTAX/2.
      DX = DELTAX/REAL(NX-1)
      DO IY = 1, NHY
       X = -DELTAX/2.
       DO IX = 1, NHX
        IF (ABS(X) .GE. R) THEN
          DY = (DELTAY)/REAL(NY-1)
        ELSE
          DY = (DELTAY - 2.*SIN(ACOS(ABS(X)/R)))/REAL(NY-1)
        END IF
        Y = -DELTAY/2. + (IY-1)*DY
!       Bottom left
        XFIELD((IY-1)*NX + IX,:) = (/X, Y/)
!       Bottom right
        XFIELD((IY)*NX - IX + 1,:) = (/-X, Y/)
!       Top left
        XFIELD(NHX*NHY*4 - (IY)*NX + IX,:) = (/X, -Y/)
!       Top right
        XFIELD(NHX*NHY*4 - (IY-1)*NX - IX + 1,:) = (/-X, -Y/)
        X = X + DX
       END DO
      END DO
END SUBROUTINE
