      module constants
      implicit none

      integer :: nq1, nq2, nq3
      integer :: nqtot
      integer :: nq_out
      integer :: nqBZ
      integer :: nnu
      integer, allocatable :: nstar(:)
      integer :: nr     ! number of input r points (el)
      integer :: nsig     ! number of smearing patterns
      real(8) :: weight

      real(8) :: avec1(3), avec2(3), avec3(3)      
      real(8) :: bvec1(3), bvec2(3), bvec3(3)
      real(8), parameter :: pi=3.14159265358979d0
      real(8), parameter :: Ry2eV=13.6058d0
      real(8), parameter :: Kayser2eV=1d0/8065d0
      end module
