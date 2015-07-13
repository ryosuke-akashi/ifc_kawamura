      subroutine cryst_to_cart(nvec, vec_in, vec_out, switch)
      USE constants, ONLY : avec1, avec2, avec3,                       &
     &                       bvec1, bvec2, bvec3,pi
      implicit none

      integer, intent(in) :: nvec
      real(8), intent(in) :: vec_in(3,nvec)
      real(8), intent(out) :: vec_out(3,nvec)
      integer, intent(in) :: switch
      !! if switch =  1 : vector in reciprocal space 
      !! if switch = -1 : vector in real space 

      !! LOCAL 
      integer iv
      !!/LOCAL 
      write(*,*)"nvec",nvec
      IF(switch.eq.1) THEN
       DO iv = 1, nvec
        write(*,*)"iv=",iv
        vec_out(:,iv) = bvec1(:)*vec_in(1,iv)                           &
     &                + bvec2(:)*vec_in(2,iv)                           &
     &                + bvec3(:)*vec_in(3,iv)
       ENDDO
       vec_out(:,:)=2d0*pi*vec_out(:,:)
      ELSE IF(switch.eq.-1) THEN
       DO iv = 1, nvec
        vec_out(:,iv) = avec1(:)*vec_in(1,iv)                           &
     &                + avec2(:)*vec_in(2,iv)                           &
     &                + avec3(:)*vec_in(3,iv)
       ENDDO
      ELSE
       write(*,*)"ERROR: cryst_to_cart switch"
       stop
      ENDIF

      RETURN
      end subroutine
