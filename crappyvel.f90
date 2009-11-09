
!Calculate B and C matrices on the field points X for velocity calculation
SUBROUTINE FLDMATBCV(NX, X, N, XNODE, BX, BY, CX, CY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL(KIND=8), DIMENSION(NX,2), INTENT(IN) :: X
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(NX,N), INTENT(OUT) :: BX, BY, CX, CY
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: BIJX, BIJY, CIJX, CIJY
      INTEGER :: I, J, NOE
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJV(XI, XNODE(NOE(N,J,0),:), XNODE(NOE(N,J,1),:), BIJX, BIJY, CIJX, CIJY)
        BX(I,J) = BIJX
        BY(I,J) = BIJY
        CX(I,J) = CIJX
        CY(I,J) = CIJY
       END DO
      END DO
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
