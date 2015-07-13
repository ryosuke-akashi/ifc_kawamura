      module matrices
      implicit none

      real(8), allocatable :: qcrys(:,:), qcart(:,:)
      real(8), allocatable :: rcrys(:,:), rcart(:,:)
      real(8), allocatable :: hr_re(:,:,:), hr_im(:,:,:)
      real(8), allocatable :: wr(:)
      real(8), allocatable :: mass(:)
      real(8), allocatable :: q_out(:,:)
      real(8), allocatable :: w_out_sq(:,:), w_out(:,:)

      complex(8), allocatable :: dynmat(:,:,:)  ! dynmat D(q)
      complex(8), allocatable :: ifc(:,:,:)  ! dynmat D(R)
      complex(8), allocatable :: dynmat_out(:,:,:)  ! dynmat D(q)

      complex(8), allocatable :: modevec_in(:,:,:)  ! input mode vector generated from dynmat
      real(8), allocatable :: gamma_nu(:,:,:)  ! input linewidth gamma(q,nu)
      complex(8), allocatable :: gamma_ab(:,:,:,:)  ! linewidth matrix gamma(q,a,b)
      complex(8), allocatable :: gamma_r(:,:,:,:)  ! linewidth matrix gamma(r,a,b)
      complex(8), allocatable :: gamma_ab_out(:,:,:,:)  ! linewidth matrix gamma(q,a,b)
      real(8), allocatable :: gamma_nu_out(:,:,:)  ! input linewidth gamma(q,nu)
      end module 
