!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE


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


!Impose boundary conditions for a moving circle with u velocity
!in the body frame of reference
SUBROUTINE BCONDVEL(N, XNODE, U, CHI)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL(KIND=8), DIMENSION(N,3), INTENT(IN) :: XNODE
      REAL(KIND=8), DIMENSION(2), INTENT(IN) :: U
      REAL(KIND=8), DIMENSION(N), INTENT(OUT) :: CHI
      REAL(KIND=8), DIMENSION(2) :: T
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


!Calculate the velocity in the field
!  PHISRF   On surface                          (in)
!  CHISRF   On surface                          (in)
!  B, C     In the field                        (in)
!  PHIFLD   In the field                        (out)
SUBROUTINE CALCVELFLD(N, PHISRF, CHISRF, NX, BX, BY, CX, CY, VELFLD)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, NX
      REAL(KIND=8), DIMENSION(N), INTENT(IN) :: PHISRF, CHISRF
      REAL(KIND=8), DIMENSION(NX,N) :: BX, BY, CX, CY
      REAL(KIND=8), DIMENSION(NX,2), INTENT(OUT) :: VELFLD
      VELFLD(:,1) = MATMUL(BX, CHISRF) + MATMUL(CX, PHISRF)
      VELFLD(:,2) = MATMUL(BY, CHISRF) + MATMUL(CY, PHISRF)
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
