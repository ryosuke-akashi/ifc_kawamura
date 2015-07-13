      program test_diag
      implicit none

      integer nmat, i
      complex(8) mat(3,3)
      real(8) eigen(3)
      complex(8) vec(3,3)
      nmat=3

      mat(1,1)=(0d0,0d0)
      mat(1,2)=(-1d0,0d0)
      mat(1,3)=(0d0,0d0)
      mat(2,1)=(-1d0,0d0)
      mat(2,2)=(0d0,0d0)
      mat(2,3)=(0d0,0d0)
      mat(3,1)=(0d0,0d0)
      mat(3,2)=(0d0,0d0)
      mat(3,3)=(-1d0,0d0)

      CALL diag_herm_sub(nmat,mat(:,:),eigen(:),                        &
     &                   vec(:,:) )
      write(*,*)"eigen", eigen(:)
      DO i=1,nmat
       write(*,*)"vec", vec(:,i)
      ENDDO

      end program
