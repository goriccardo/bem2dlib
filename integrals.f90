!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE


!Local frame of reference (LFOR)
! ------ INPUT ------
!   XS      point X* in the global fram of reference (GFOR)
!   XJM     first point (X-) of the element in the GFOR
!   XJP     second point (X+) of the element in the GFOR
! ------ OUTPUT ------
!   XM      abscissa of the first point (X-) of the element in the LFOR
!   XP      abscissa of the second point (X+) of the element in the LFOR
!   YS      ordinate of the point XS (X*) in the LFOR
!   RM      distance between XS (X*) and XJM (X-)
!   RP      distance between XS (X*) and XJP (X+)
SUBROUTINE LOCFR(XS, XJM, XJP, XM, XP, YS, RM, RP)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: XS, XJM, XJP
      REAL(KIND=8), INTENT(OUT) :: XM, XP, YS, RM, RP
      REAL(KIND=8), DIMENSION(2) :: EX, VXP, VXM
      EX = XJP - XJM
      EX = EX/DSQRT(EX(1)**2 + EX(2)**2)
!     Vector X_+
      VXP = (XJP - XS)
!     Vector X_-
      VXM = (XJM - XS)
!     XP is always greater than XM
      XP = VXP(1)*EX(1) + VXP(2)*EX(2)
      XM = VXM(1)*EX(1) + VXM(2)*EX(2)
!     YS is positive on the inside and negative on the outside
      YS = VXP(1)*EX(2) - VXP(2)*EX(1)
      RP = DSQRT(XP**2 + YS**2)
      RM = DSQRT(XM**2 + YS**2)
END SUBROUTINE


!BIJ is the integral of G(XI,XJ) from XJM to XJP, nodes of the J element
!CIJ is the integral of dG/dn(XI,XJ) from XJM to XJP, nodes of the J element
SUBROUTINE BCIJ(XI, XJM, XJP, BIJ, CIJ)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: XI, XJM, XJP
      REAL(KIND=8), INTENT(OUT) :: BIJ, CIJ
      REAL(KIND=8) :: XM, XP, YS, RM, RP
      REAL(KIND=8) :: AYS, S
      REAL(KIND=8), PARAMETER :: PI = REAL(4.,8)*DATAN(REAL(1.,8))
      CALL LOCFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      AYS = DABS(YS)
      S = DSIGN(1.D0,YS)
      BIJ = (XP*DLOG(RP) - XP - XM*DLOG(RM) + XM + YS*DATAN2(S*XP,AYS) - YS*DATAN2(S*XM,AYS))/(REAL(2.,8)*PI)
      CIJ = -(DATAN2(S*XP,AYS) - DATAN2(S*XM,AYS))/(REAL(2.,8)*PI)
END SUBROUTINE


SUBROUTINE BCIJV(XI, XJM, XJP, BIJX, BIJY, CIJX, CIJY)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: XI, XJM, XJP
      REAL(KIND=8) :: TH
      REAL(KIND=8), INTENT(OUT) :: BIJX, BIJY, CIJX, CIJY
      REAL(KIND=8) :: XM, XP, YS, RM, RP
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      CALL LOCFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      TH = DATAN2(XJP(2)-XJM(2),XJP(1)-XJP(1)) + PI/2.
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
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(N,N), INTENT(OUT) :: B, C
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: CIJ, BIJ
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


!Calculate B and C matrices on the field points X
SUBROUTINE FLDMATBC(NX, X, N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL(KIND=8), DIMENSION(NX,2), INTENT(IN) :: X
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(NX,N), INTENT(OUT) :: B, C
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: BIJ, CIJ
      INTEGER :: I, J, NODOFEL
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJ(XI, XNODE(NODOFEL(N,J,1),:2), XNODE(NODOFEL(N,J,2),:2), BIJ, CIJ)
        B(I,J) = BIJ
        C(I,J) = CIJ
       END DO
      END DO
END SUBROUTINE


!Calculate B and C matrices on the field points X for velocity calculation
SUBROUTINE FLDMATBCV(NX, X, N, XNODE, BX, BY, CX, CY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL(KIND=8), DIMENSION(NX,2), INTENT(IN) :: X
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(NX,N), INTENT(OUT) :: BX, BY, CX, CY
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: BIJX, BIJY, CIJX, CIJY
      INTEGER :: I, J, NODOFEL
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJV(XI, XNODE(NODOFEL(N,J,1),:), XNODE(NODOFEL(N,J,2),:), BIJX, BIJY, CIJX, CIJY)
        BX(I,J) = BIJX
        BY(I,J) = BIJY
        CX(I,J) = CIJX
        CY(I,J) = CIJY
       END DO
      END DO
END SUBROUTINE

