!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Impose boundary conditions for a moving circle with u velocity
!in the air frame of reference
SUBROUTINE BCONDVEL(N, XNODE, U, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N,2), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: U
      REAL(KIND=8), DIMENSION(N), INTENT(OUT) :: CHI
      REAL(KIND=8), DIMENSION(2) :: T
      INTEGER :: I, NOE
      DO I = 1, N
       T = XNODE(NOE(N,I,1),:2) - XNODE(NOE(N,I,0),:2)
       T = T/SQRT(T(1)**2 + T(2)**2)
       CHI(I) = U(1)*T(2) - U(2)*T(1)
      END DO
END SUBROUTINE


!Solve the (NxN) linear system (0.5*I - C)*phi = C*chi
SUBROUTINE SOLVEPHI(N, B, C, PHI, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N,N), INTENT(IN) :: B, C
      REAL(KIND=8), DIMENSION(N), INTENT(OUT) :: PHI
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: CHI
      INTEGER :: I, INFO
      REAL(KIND=8), DIMENSION(N,N) :: A
      REAL(KIND=8), DIMENSION(N) :: RHS, IPIV
      A = -C
      DO I = 1, N
       A(I,I) = A(I,I) + 0.5
      END DO
      RHS = MATMUL(B,CHI)
      CALL DGESV(N, 1, A, N, IPIV, RHS, N, INFO)
      IF (INFO .NE. 0) THEN
        WRITE(*,*) "ERROR IN LINEAR SYSTEM"
      END IF
      PHI = RHS
END SUBROUTINE


SUBROUTINE SolvePhiTime(N, B, C, NWake, D, NTime, ChiTime, PhiTime, DPhiW)
      IMPLICIT NONE
      integer, intent(IN) :: N, NWake, NTime
      real(kind=8), dimension(N,N), intent(IN) :: B, C
      real(kind=8), dimension(N,NWake), intent(IN) :: D
      real(kind=8), dimension(N,NTime), intent(IN) :: ChiTime
      real(kind=8), dimension(NWake,NTime), intent(OUT) :: DPhiW
      real(kind=8), dimension(N,NTime), intent(OUT) :: PhiTime
      integer :: I, INFO
      real(kind=8), dimension(N,N) :: A, ATEMP
      real(kind=8), dimension(N) :: RHS, IPIV
      DPhiW(:,:) = 0.
      A = -C
      do i = 1, N
       A(i,i) = A(i,i) + DBLE(0.5)
      end do
      do i = 1,NTime
       ATemp = A
       RHS = MATMUL(B,ChiTime(:,i)) + MATMUL(D, DPhiW(:,i))
       CALL DGESV(N, 1, ATEMP, N, IPIV, RHS, N, INFO)
       IF (INFO .NE. 0) THEN
         WRITE(*,*) "ERROR IN LINEAR SYSTEM"
       END IF
       PhiTime(:,i) = RHS
       CALL Wake(N, NWake, NTime, i, PhiTime, DPhiW)
      end do
END SUBROUTINE


!Calculate the potential in the field
!  PHISRF   On surface                          (in)
!  CHISRF   On surface                          (in)
!  B, C     In the field                        (in)
!  PHIFLD   In the field                        (out)
SUBROUTINE CALCPHIFLD(N, PHISRF, CHISRF, NX, B, C, PHIFLD)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NX
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: PHISRF, CHISRF
      REAL(KIND=8), DIMENSION(NX,N) :: B, C
      REAL(KIND=8), DIMENSION(NX), INTENT(OUT) :: PHIFLD
      PHIFLD = MATMUL(B, CHISRF) + MATMUL(C, PHISRF)
END SUBROUTINE
