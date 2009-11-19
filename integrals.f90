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


!DIJ is the integral of dG/dn(XI,XJ) on the wake
!  XI is on the body
!  XJ is on the wake
double precision function DIJ(XI, XJM, XJP)
      IMPLICIT NONE
      real(kind=8), dimension(2), intent(IN) :: XI, XJM, XJP
      real(kind=8) :: XM, XP, YS, RM, RP
      real(kind=8) :: AYS, S
      real(kind=8), parameter :: PI = dble(4)*datan(1.D0)
      call LocFR(XI, XJM, XJP, XM, XP, YS, RM, RP)
      AYS = dabs(YS)
      S = dsign(1.D0,YS)
      DIJ = (datan2(S*XP,AYS) - datan2(S*XM,AYS))/(dble(2)*PI)
end function


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


!Calculate the matrix B and C
!  N        Figure of elements  (IN)
!  XNODE    Nodes vector        (IN)
!  B        Matrix B            (OUT)
!  C        Matrix C            (OUT)
SUBROUTINE SRFMATBC(N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(N,N), INTENT(OUT) :: B, C
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: CIJ, BIJ
      INTEGER :: I, J, NOE
      DO I = 1, N
       !Midpoint of I element
       XI = (XNODE(NOE(N,I,0),:) + XNODE(NOE(N,I,1),:))/2.
       C(I,I) = 0.
       DO J = 1, N
        CALL BCIJ(XI, XNODE(NOE(N,J,0),:), XNODE(NOE(N,J,1),:),BIJ,CIJ)
        B(I,J) = BIJ
        IF (I .NE. J) THEN
          C(I,J) = CIJ
        END IF
       END DO
      END DO
END SUBROUTINE


!Calculate the matrices B, C and D
!  N        Figure of elements            (IN)
!  XNODE    Nodes vector                  (IN)
!  NW       Number of wake elements       (IN)
!  XWNODE   Wake's nodes vector           (IN)
!  B        Matrix B                      (OUT)
!  C        Matrix C                      (OUT)
!  D        Matrix D (wake doublet int.)  (OUT)
subroutine SrfMatBCD(N, Xnode, NW, XWnode, B, C, D)
      IMPLICIT NONE
      integer, intent(IN) :: N, NW
      real(kind=8), dimension(N,2), intent(IN) :: Xnode
      real(kind=8), dimension(NW,2), intent(IN) :: XWnode
      real(kind=8), dimension(N,N), intent(OUT) :: B, C
      real(kind=8), dimension(N,NW), intent(OUT) :: D
      real(kind=8), dimension(2) :: XI
      real(kind=8) :: CIJ, BIJ, DIJ
      integer :: i, j, NOE
      do i = 1, N
       !Midpoint of the i-th element
       XI = (Xnode(NOE(N,i,0),:) + Xnode(NOE(N,i,1),:))/2.
       C(i,i) = 0.
       do j = 1, N
        call BCIJ(XI, Xnode(NOE(N,j,0),:), Xnode(NOE(N,j,1),:),BIJ,CIJ)
        B(i,j) = BIJ
        if (i .ne. j) then
          C(i,j) = CIJ
        end if
       end do
       do j = 1, NW-1
        D(i,j) = DIJ(XI, XWnode(j,:), XWnode(j+1,:))
       end do
      end do
end subroutine


!Calculate B and C matrices on the field points X
SUBROUTINE FLDMATBC(NX, X, N, XNODE, B, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX, N
      REAL(KIND=8), DIMENSION(NX,2), INTENT(IN) :: X
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(NX,N), INTENT(OUT) :: B, C
      REAL(KIND=8), DIMENSION(2) :: XI
      REAL(KIND=8) :: BIJ, CIJ
      INTEGER :: I, J, NOE
      DO I = 1, NX
       XI = X(I,:)
       DO J = 1, N
        CALL BCIJ(XI, XNODE(NOE(N,J,0),:2), XNODE(NOE(N,J,1),:2), BIJ, CIJ)
        B(I,J) = BIJ
        C(I,J) = CIJ
       END DO
      END DO
END SUBROUTINE
