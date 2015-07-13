      subroutine diag_herm_sub(ndim,mat,eig,vec)
      !! return the eigenvalues and eigenvectors of an hermitian matrix
      !! leaving the matrix unchanged
      !! note: vec(:,i): i-th eigenvector
      implicit none
      integer, intent(in)  :: ndim   !! dimension
      complex(8), intent(in)  :: mat(ndim, ndim)   !! matrix
      real(8), intent(out) :: eig(ndim)   !! eigenval
      complex(8), intent(out) :: vec(ndim,ndim)   !! eigenvec

      integer lwork, info   
      real(8), allocatable :: rwork(:)      
      complex(8), allocatable :: work(:) 

      vec(:,:) = mat(:,:)

      lwork = 3*ndim-2
      allocate(work(lwork))
      allocate(rwork(lwork))
      call zheev('V', 'U', ndim, vec(:,:), ndim, eig(:),    &
     &           work,lwork, rwork,info)

      deallocate(work,rwork)
      return
      end subroutine
