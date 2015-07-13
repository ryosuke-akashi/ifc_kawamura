      subroutine ph_InvFT_D()
      USE constants, ONLY : nqtot,nr,nnu,nsig
      USE flags, ONLY : la2F
      USE matrices, ONLY: dynmat, qcart, rcart, ifc,                    &
     &                    gamma_ab,gamma_r
      implicit none

      !! LOCAL
      integer isig
      real(8) Dq_tmp_re(nnu,nnu,nqtot)
      real(8) Dq_tmp_im(nnu,nnu,nqtot)
      real(8) Dr_tmp_re(nnu,nnu)
      real(8) Dr_tmp_im(nnu,nnu)
      real(8) prod
      real(8) cos_prod, sin_prod
      integer ir, iq, inu, jnu

      real(8) gq_tmp_re(nnu,nnu,nsig,nqtot)
      real(8) gq_tmp_im(nnu,nnu,nsig,nqtot)
      real(8) gr_tmp_re(nnu,nnu,nsig)
      real(8) gr_tmp_im(nnu,nnu,nsig)
      !!/LOCAL


      !! IFT D(nu,nu,nq) --> D(nu,nu,nr)
      Dq_tmp_re(:,:,:)=0d0
      Dq_tmp_im(:,:,:)=0d0

      DO iq=1, nqtot
      DO inu=1, nnu
      DO jnu=1, nnu
       Dq_tmp_re(jnu,inu,iq)=dble(dynmat(jnu,inu,iq))
       Dq_tmp_im(jnu,inu,iq)=dimag(dynmat(jnu,inu,iq))
      ENDDO
      ENDDO
      ENDDO

      DO iq=1, nqtot
       write(*,*)
       write(*,*)"Dq_tmp_re"
      DO inu=1, nnu
       write(*,'(256f18.8)')Dq_tmp_re(inu,:,iq)
      ENDDO
      ENDDO
      DO iq=1, nqtot
       write(*,*)
       write(*,*)"Dq_tmp_im"
      DO inu=1, nnu
       write(*,'(256f18.8)')Dq_tmp_im(inu,:,iq)
      ENDDO
      ENDDO

      IF(la2F.eq.1)THEN
       allocate(gamma_r(nnu,nnu,nsig,nr))
       DO iq=1, nqtot
       DO isig=1, nsig
       DO inu=1, nnu
       DO jnu=1, nnu
        gq_tmp_re(jnu,inu,isig,iq)=dble(gamma_ab(jnu,inu,isig,iq))
        gq_tmp_im(jnu,inu,isig,iq)=dimag(gamma_ab(jnu,inu,isig,iq))
       ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF

      DO ir=1,nr
       Dr_tmp_re(:,:)=0d0
       Dr_tmp_im(:,:)=0d0
       IF(la2F.eq.1)THEN
        gr_tmp_re(:,:,:)=0d0
        gr_tmp_im(:,:,:)=0d0
       ENDIF
       DO iq=1,nqtot
        prod    =  rcart(1,ir)*qcart(1,iq)                            &
     &           + rcart(2,ir)*qcart(2,iq)                            &
     &           + rcart(3,ir)*qcart(3,iq) 
        cos_prod=dcos(prod)
        sin_prod=dsin(prod)
        DO inu=1,nnu
         DO jnu=1,nnu
          IF(inu.eq.1 .and. jnu.eq.1)THEN
           write(*,*)"Dr_tmp_re(jnu,inu)",Dr_tmp_re(jnu,inu)
           write(*,*)"Dr_tmp_im(jnu,inu)",Dr_tmp_im(jnu,inu)
           write(*,'(a15,3f18.8)')"prod,con,sin",prod,                  &
     &                           dcos(prod), dsin(prod)
          ENDIF
          Dr_tmp_re(jnu,inu)= Dr_tmp_re(jnu,inu)                        &
     &                     +cos_prod*Dq_tmp_re(jnu,inu,iq)            &
     &                     +sin_prod*Dq_tmp_im(jnu,inu,iq)
          Dr_tmp_im(jnu,inu)= Dr_tmp_im(jnu,inu)                        &
     &                     +cos_prod*Dq_tmp_im(jnu,inu,iq)            &
     &                     -sin_prod*Dq_tmp_re(jnu,inu,iq)
          IF(la2F.eq.1)THEN
           DO isig=1,nsig
            gr_tmp_re(jnu,inu,isig)= gr_tmp_re(jnu,inu,isig)            &
     &                       +cos_prod*gq_tmp_re(jnu,inu,isig,iq)       &
     &                       +sin_prod*gq_tmp_im(jnu,inu,isig,iq)
            gr_tmp_im(jnu,inu,isig)= gr_tmp_im(jnu,inu,isig)            &
     &                       +cos_prod*gq_tmp_im(jnu,inu,isig,iq)       &
     &                       -sin_prod*gq_tmp_re(jnu,inu,isig,iq)
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDDO
!! AKASHI test 150430
       Dr_tmp_re(:,:)=Dr_tmp_re(:,:)/dble(nqtot)
       Dr_tmp_im(:,:)=Dr_tmp_im(:,:)/dble(nqtot)
       IF(la2F.eq.1)THEN
        gr_tmp_re(:,:,:)=gr_tmp_re(:,:,:)/dble(nqtot)
        gr_tmp_im(:,:,:)=gr_tmp_im(:,:,:)/dble(nqtot)
       ENDIF
!!/AKASHI test 150430
       DO inu=1,nnu
        DO jnu=1,nnu
         ifc(jnu,inu,ir) =                                              &
     &            dcmplx(Dr_tmp_re(jnu,inu),Dr_tmp_im(jnu,inu))
         IF(la2F.eq.1)THEN
          DO isig=1, nsig
           gamma_r(jnu,inu,isig,ir)=                                    &
     &          dcmplx(gr_tmp_re(jnu,inu,isig),gr_tmp_im(jnu,inu,isig))
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO


      end subroutine
