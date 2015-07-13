      subroutine ph_ifc_out()
      USE constants, ONLY : nr, nnu, nsig
      USE flags, ONLY : la2F
      USE filenames, ONLY : filmass
      USE matrices, ONLY : q_out, w_out,wr,rcrys,                   &
     &                     ifc, gamma_r,mass
      implicit none

      integer iq,isig, ir, inu, jnu

      open(unit=10, file="ifc.dat", status="unknown")
      write(10,*)
      write(10,*)nnu
      write(10,*)nr
      write(10,'(15i5)')nint(wr(:))
      DO ir=1,nr
       DO inu=1,nnu
        DO jnu=1,nnu
         write(10,'(5i5,2e20.12)')nint(rcrys(:,ir)),jnu,inu,             &
     &                           ifc(jnu,inu,ir)
        ENDDO
       ENDDO
      ENDDO
      close(10)

      IF(la2F.eq.1)THEN
       open(unit=11, file="gamma_r.dat", status="unknown")
       write(11,*)
       write(11,*)nnu
       write(11,*)nr
       write(11,'(15i5)')nint(wr(:))
       DO isig=1,nsig
        write(11,*)"isig=",isig
        DO ir=1,nr
         DO inu=1,nnu
          DO jnu=1,nnu
           write(11,'(5i5,2e20.12)')nint(rcrys(:,ir)),jnu,inu,           &
     &                             gamma_r(jnu,inu,isig,ir)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       close(11)

      ENDIF
      open(unit=12,file=filmass,status="unknown")
      write(12,'(6f20.10)')mass(:)
      close(12)
      end subroutine
