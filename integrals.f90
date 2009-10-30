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
      REAL, DIMENSION(2), INTENT(IN) :: XS, XJM, XJP
      REAL, INTENT(OUT) :: XM, XP, YS, RM, RP
      REAL, DIMENSION(2) :: EX, VXP, VXM
      EX = XJP - XJM
      EX = EX/SQRT(EX(1)**2 + EX(2)**2)
!     Vector X_+
      VXP = (XJP - XS)
!     Vector X_-
      VXM = (XJM - XS)
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
      BIJ = (XP*LOG(RP) - XP + YS*ATAN2(XP,ABS(YS)) - XM*LOG(RM) + XM - YS*ATAN2(XM,ABS(YS)))/(2.*PI)
      CIJ = -(ATAN2(XP,ABS(YS)) - ATAN2(XM,ABS(YS)))/(2.*PI)
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
        CALL BCIJ(XI, XNODE(NODOFEL(N,J,1),:2), XNODE(NODOFEL(N,J,2),:2), BIJ, CIJ)
        B(I,J) = BIJ
        C(I,J) = CIJ
       END DO
      END DO
END SUBROUTINE

