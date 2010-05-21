subroutine clBemNACA00xx(Nup, alpha, chord, t, cl)
      IMPLICIT NONE
      integer, intent(IN) :: Nup
      integer :: Nelem
      integer, parameter :: Nfreq = 1
      integer :: NWake
      real(kind=8), intent(IN) :: chord, t
      real(kind=8), intent(OUT) :: cl
!     Vector of nodes global coordinates (x,y)
      real(kind=8), allocatable :: Xnode(:,:)
!     Circle radius
      real(kind=8), parameter :: uscalar = 1.D0
      real(kind=8) :: UHoriz
      real(kind=8), intent(IN) :: alpha
      complex(kind=8), dimension(Nfreq) :: s
!     Potential and normal wash on the surface
      real(kind=8), allocatable :: B(:,:), C(:,:), D(:,:)
      real(kind=8), allocatable :: XWnode(:,:)
      complex(kind=8), allocatable :: ChiLap(:,:), PhiLap(:,:)
      complex(kind=8), allocatable :: Us(:,:)
!     real(kind=8), dimension(Nelem/2, NTime) :: DPhiBody
      complex(kind=8) :: clLap
!     Time step
      real(kind=8) :: DT, CalcDT
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      logical :: TEat1 = .true.
      Nelem = 2*Nup+1
      Nwake = Nup*10
      allocate(Xnode(Nelem, 2))
      allocate(B(Nelem,Nelem),C(Nelem,Nelem),D(Nelem,NWake))
      allocate(XWnode(NWake,2))
      allocate(ChiLap(Nelem,Nfreq), PhiLap(Nelem,Nfreq))
      allocate(Us(Nfreq,2))
!     The program starts here!
      call GeomNACA00xx(NUp, chord, t, Xnode)
      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
      DT = CalcDT(Nelem, Xnode, UHoriz)
      call BCondStraightLap(Nelem, Xnode, uscalar, alpha, Nfreq, s, Us, Chilap)
      call WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call SolvePhiLap(Nelem, B, C, Nwake, D, Nfreq, s, DT, ChiLap, PhiLap)
      call CalcClLap(Nelem, Xnode, TEat1, Nfreq, s, Us, phiLap, chiLap, clLap)
      cl = cdabs(clLap)
      deallocate(Xnode, B, C, D, XWnode, ChiLap, PhiLap, Us)
end subroutine

