!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

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
