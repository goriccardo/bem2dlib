!Calculate the position of the nodes
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  R        Radius                   (IN)
SUBROUTINE CIRCLEGEOM(NELEM, XNODE, R)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL, INTENT(IN) :: R
      REAL, DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL :: DTH, THETA
      DTH = 2.*PI/REAL(NELEM)
      DO I = 1,NELEM
       THETA = DTH*REAL(I-1)
       XNODE(I,1) = R*COS(THETA)
       XNODE(I,2) = R*SIN(THETA)
       XNODE(I,3) = THETA+DTH*0.5
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


!Local frame of reference
SUBROUTINE LOCFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      IMPLICIT NONE
      REAL, DIMENSION(2), INTENT(IN) :: XI, XJM, XJP
      REAL, INTENT(OUT) :: XM, XP, YS, RM, RP
      REAL, DIMENSION(2) :: EX, VXP, VXM
      EX = XJP - XJM
      EX = EX/SQRT(EX(1)**2 + EX(2)**2)
!     Vector X_+
      VXP = (XJP - XI)
!     Vector X_-
      VXM = (XJM - XI)
!     XP is always greater than XM
      XP = VXP(1)*EX(1) + VXP(2)*EX(2)
      XM = VXM(1)*EX(1) + VXM(2)*EX(2)
!     YS is positive on the inside and negative on the outside
      YS = VXP(1)*EX(2) - VXP(2)*EX(1)
      RP = SQRT(XP**2 + YS**2)
      RM = SQRT(XM**2 + YS**2)
END SUBROUTINE


!BIJ is the integral of G(XI,XJ) from XJM to XJP, nodes of the J element
!CIJ is the integral of dG/dn(XI,XJ) from XJM to XJP, nodes of the J element
SUBROUTINE BCIJ(XI, XJM, XJP, BIJ, CIJ)
      IMPLICIT NONE
      REAL, DIMENSION(2), INTENT(IN) :: XI, XJM, XJP
      REAL, INTENT(OUT) :: BIJ, CIJ
      REAL :: XM, XP, YS, RM, RP
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      CALL LOCFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      BIJ = (XP*LOG(RP) - XP + YS*ATAN2(XP,YS) - XM*LOG(RM) + XM - YS*ATAN2(XM,YS))/(2.*PI)
      CIJ = (ATAN2(XP,YS) - ATAN2(XM,YS))/(2.*PI)
END SUBROUTINE


SUBROUTINE BCIJV(XI, XJM, XJP, TH, BIJX, BIJY, CIJX, CIJY)
      IMPLICIT NONE
      REAL, DIMENSION(2), INTENT(IN) :: XI, XJM, XJP
      REAL, INTENT(IN) :: TH
      REAL, INTENT(OUT) :: BIJX, BIJY, CIJX, CIJY
      REAL :: XM, XP, YS, RM, RP
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      CALL LOCFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      BIJX = (ATAN2(YS,XM)*COS(TH) + LOG(RM)*SIN(TH) - ATAN2(YS,XP)*COS(TH) - LOG(RP)*SIN(TH))/(2.*PI)
      BIJY = (ATAN2(YS,XM)*SIN(TH) - LOG(RM)*COS(TH) - ATAN2(YS,XP)*SIN(TH) + LOG(RP)*COS(TH))/(2.*PI)
      CIJX = ((XM*COS(TH)+YS*SIN(TH))/(RM**2) - (XP*COS(TH)+YS*SIN(TH))/(RP**2))/(2.*PI)
      CIJY = ((XM*SIN(TH)-YS*COS(TH))/(RM**2) - (XP*SIN(TH)-YS*COS(TH))/(RP**2))/(2.*PI)
END SUBROUTINE


!Calculate the matrix B and C
!  N        Figure of elements  (IN)
!  B        Matrix B            (OUT)
!  C        Matrix C            (OUT)
!  XNODE    Nodes vector        (IN)
SUBROUTINE SRFMATBC(N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N,3), INTENT(IN) :: XNODE
      REAL, DIMENSION(N,N), INTENT(OUT) :: B, C
      REAL, DIMENSION(2) :: XI
      REAL :: CIJ, BIJ
      INTEGER :: I, J, NODOFEL
      DO I = 1, N
       !Midpoint of I element
       XI = (XNODE(NODOFEL(N,I,1),:2) + XNODE(NODOFEL(N,I,2),:2))/2.
       C(I,I) = 0.
       DO J = 1, N
        CALL BCIJ(XI, XNODE(NODOFEL(N,J,1),:2), XNODE(NODOFEL(N,J,2),:2),BIJ,CIJ)
        B(I,J) = BIJ
        IF (I .NE. J) THEN
          C(I,J) = CIJ
        END IF
       END DO
      END DO
END SUBROUTINE


!Impose boundary conditions for a moving circle with u velocity
!in the body frame of reference
SUBROUTINE BCONDVEL(N, XNODE, U, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL, DIMENSION(N,3), INTENT(IN) :: XNODE
      REAL, DIMENSION(2), INTENT(IN) :: U
      REAL, DIMENSION(N), INTENT(OUT) :: CHI
      REAL, DIMENSION(2) :: T
      INTEGER :: I, NODOFEL
      DO I = 1, N
       T = XNODE(NODOFEL(N,I,2),:2) - XNODE(NODOFEL(N,I,1),:2)
       T = T/SQRT(T(1)**2 + T(2)**2)
       CHI(I) = U(1)*T(2) - U(2)*T(1)
      END DO
END SUBROUTINE


!Solve the (NxN) linear system (0.5*I - C)*phi = C*chi
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
      REAL, DIMENSION(N,3), INTENT(IN) :: XNODE
      REAL, DIMENSION(NX,N), INTENT(OUT) :: B, C
      REAL, DIMENSION(2) :: XI
      REAL :: BIJ, CIJ
      INTEGER :: I, J, NODOFEL
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJ(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:), BIJ, CIJ)
        B(I,J) = BIJ
        C(I,J) = CIJ
       END DO
      END DO
END SUBROUTINE


!Calculate B and C matrices on the field points X for velocity calculation
SUBROUTINE FLDMATBCV(NX, X, N, XNODE, BX, BY, CX, CY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL, DIMENSION(NX,2), INTENT(IN) :: X
      REAL, DIMENSION(N,3), INTENT(IN) :: XNODE
      REAL, DIMENSION(NX,N), INTENT(OUT) :: BX, BY, CX, CY
      REAL, DIMENSION(2) :: XI
      REAL :: BIJX, BIJY, CIJX, CIJY
      INTEGER :: I, J, NODOFEL
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJV(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:), XNODE(J,3), BIJX, BIJY, CIJX, CIJY)
        BX(I,J) = BIJX
        BY(I,J) = BIJY
        CX(I,J) = CIJX
        CY(I,J) = CIJY
       END DO
      END DO
END SUBROUTINE


!Calculate the velocity in the field
!  PHISRF   On surface                          (in)
!  CHISRF   On surface                          (in)
!  B, C     In the field                        (in)
!  PHIFLD   In the field                        (out)
SUBROUTINE CALCFIELDV(N, PHISRF, CHISRF, NX, BX, BY, CX, CY, VELFLD)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NX
      REAL, DIMENSION(N), INTENT(IN) :: PHISRF, CHISRF
      REAL, DIMENSION(NX,N) :: BX, BY, CX, CY
      REAL, DIMENSION(NX,2), INTENT(OUT) :: VELFLD
      VELFLD(:,1) = MATMUL(BX, CHISRF) - MATMUL(CX, PHISRF)
      VELFLD(:,2) = MATMUL(BY, CHISRF) - MATMUL(CY, PHISRF)
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
      PHIFLD = MATMUL(B, CHISRF) - MATMUL(C, PHISRF)
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
