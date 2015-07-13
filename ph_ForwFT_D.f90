      subroutine ph_ForwFT_D()
      USE constants, ONLY : nq_out, nr, nnu,weight
      USE flags, ONLY : la2F
      USE matrices, ONLY : rcart, wr, ifc,dynmat_out,q_out,             &
     &                     gamma_r,gamma_ab_out
      implicit none


      !! LOCAL
      integer, parameter :: nsig=10
      real(8) Dr_tmp_re(nnu,nnu,nr)
      real(8) Dr_tmp_im(nnu,nnu,nr)
      real(8) Dq_tmp_re(nnu,nnu)
      real(8) Dq_tmp_im(nnu,nnu)
      integer ir, iq, inu, jnu,isig
      real(8) prod
      real(8) cos_prod, sin_prod


      real(8) gr_tmp_re(nnu,nnu,nsig,nr)
      real(8) gr_tmp_im(nnu,nnu,nsig,nr)
      real(8) gq_tmp_re(nnu,nnu,nsig)
      real(8) gq_tmp_im(nnu,nnu,nsig)
      !!/LOCAL

      !! FT D(nu,nu,nr) --> D(nu,nu,nq)
      Dr_tmp_re(:,:,:)=0d0
      Dr_tmp_im(:,:,:)=0d0
      DO ir=1, nr
      DO inu=1, nnu
      DO jnu=1, nnu
       Dr_tmp_re(jnu,inu,ir)=dble(ifc(jnu,inu,ir))
       Dr_tmp_im(jnu,inu,ir)=dimag(ifc(jnu,inu,ir))
      ENDDO
      ENDDO
      ENDDO

      IF(la2F.eq.1)THEN
       allocate(gamma_ab_out(nnu,nnu,nsig,nq_out))
       gr_tmp_re(:,:,:,:)=0d0
       gr_tmp_im(:,:,:,:)=0d0
       DO ir=1, nr
       DO isig=1,nsig
       DO inu=1, nnu
       DO jnu=1, nnu
        gr_tmp_re(jnu,inu,isig,ir)=dble(gamma_r(jnu,inu,isig,ir))
        gr_tmp_im(jnu,inu,isig,ir)=dimag(gamma_r(jnu,inu,isig,ir))
       ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF

      DO iq=1,nq_out
       Dq_tmp_re(:,:)=0d0
       Dq_tmp_im(:,:)=0d0
       IF(la2F.eq.1)THEN
        gq_tmp_re(:,:,:)=0d0
        gq_tmp_im(:,:,:)=0d0
       ENDIF
       DO ir=1,nr
        prod    =  rcart(1,ir)*q_out(1,iq)                             &
     &           + rcart(2,ir)*q_out(2,iq)                             &
     &           + rcart(3,ir)*q_out(3,iq)
        cos_prod=dcos(prod)
        sin_prod=dsin(prod)
        DO inu=1,nnu
         DO jnu=1,nnu
          Dq_tmp_re(jnu,inu)= Dq_tmp_re(jnu,inu)                        &
     &                     +(cos_prod*Dr_tmp_re(jnu,inu,ir)          &
     &                       -sin_prod*Dr_tmp_im(jnu,inu,ir))        &
     &                       /dble(wr(ir))
          Dq_tmp_im(jnu,inu)= Dq_tmp_im(jnu,inu)                        &
     &                     +(sin_prod*Dr_tmp_re(jnu,inu,ir)          &
     &                       +cos_prod*Dr_tmp_im(jnu,inu,ir))        &
     &                       /dble(wr(ir))
          IF(la2F.eq.1)THEN
           DO isig=1,nsig
            gq_tmp_re(jnu,inu,isig)= gq_tmp_re(jnu,inu,isig)            &
     &                     +(cos_prod*gr_tmp_re(jnu,inu,isig,ir)        &
     &                       -sin_prod*gr_tmp_im(jnu,inu,isig,ir))      &
     &                       /dble(wr(ir))
            gq_tmp_im(jnu,inu,isig)= gq_tmp_im(jnu,inu,isig)            &
     &                     +(sin_prod*gr_tmp_re(jnu,inu,isig,ir)        &
     &                       +cos_prod*gr_tmp_im(jnu,inu,isig,ir))      &
     &                       /dble(wr(ir))
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDDO
!! AKASHI old -150430
!       Dq_tmp_re(:,:)=Dq_tmp_re(:,:)/dble(weight)
!       Dq_tmp_im(:,:)=Dq_tmp_im(:,:)/dble(weight)
!!/ AKASHI old -150430
       DO inu=1,nnu
        DO jnu=1,nnu
         dynmat_out(jnu,inu,iq) =                                       &
     &           dcmplx(Dq_tmp_re(jnu,inu),Dq_tmp_im(jnu,inu))
         IF(la2F.eq.1)THEN
          DO isig=1, nsig
           gamma_ab_out(jnu,inu,isig,iq) =                              &
     &     dcmplx(gq_tmp_re(jnu,inu,isig),                           &
     &            gq_tmp_im(jnu,inu,isig))
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      write(*,*)"gamma_ab_out"
      DO inu=1,nnu
       write(*,*)
       write(*,'(6e20.10)')(gamma_ab_out(jnu,inu,1,26),jnu=1,nnu)
      ENDDO
      !!/FT D(nu,nu,nrp) --> D(nu,nu,nq)
      end subroutine
