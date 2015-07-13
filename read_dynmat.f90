      subroutine read_dynmat()
      USE filenames, ONLY : file_hr, fildyn,filelph
      USE constants, ONLY: nq1, nq2, nq3, nqBZ, nqtot, nnu, nstar,pi
      USE flags,     ONLY: la2F
      USE matrices, ONLY: dynmat, qcart, mass
      implicit none

      !! LOCAL
      character(6), external :: int_to_char
      integer iq, jq
      integer nat, ntyp
      integer iqtot, ii, jj, imat, iat, jat, inu, jnu
      character(256) filin
      character(256) dummy1, dummy2,dummy3
      real(8) dyn_re(3,3), dyn_im(3,3)
      real(8),allocatable :: mass_read(:)
      real(8),allocatable :: typ(:)



      !! check
      real(8), allocatable :: w_test(:)
      complex(8), allocatable :: vector(:,:)

      iq=0
      filin= "dyndir/" // TRIM(fildyn) // TRIM(int_to_char(iq))
      open(unit=10,file=filin, status="unknown")
      read(10,*)nq1,nq2,nq3
      read(10,*)nqBZ
      close(10)
      nqtot=nq1*nq2*nq3
      allocate(qcart(3,nqtot))
      allocate(nstar(nqBZ))

      iq=1
      filin= "dyndir/" // TRIM(fildyn) // TRIM(int_to_char(iq))
      open(unit=10,file=filin, status="unknown")
      read(10,*)
      read(10,*)
      read(10,*)ntyp, nat
      nnu=3*nat
      allocate(dynmat(nnu,nnu,nqtot))
      allocate(mass(nnu))
      allocate(mass_read(ntyp))
      allocate(typ(nat))
      DO ii=1,ntyp
       read(10,*)dummy1,dummy2, mass_read(ii)
       write(*,*)"dummy1,dummy2,mass_read",dummy1,dummy2,               &
     &             mass_read(ii)
      ENDDO
      DO ii=1,nat
       read(10,*)dummy1, typ(ii)
      ENDDO
      close(10)
      DO ii=1,nat
       mass(3*(ii-1)+1)=mass_read(typ(ii))
       mass(3*(ii-1)+2)=mass_read(typ(ii))
       mass(3*(ii-1)+3)=mass_read(typ(ii))
      ENDDO
      write(*,*)"mass=",mass(:)
      deallocate(mass_read,typ)

      open(unit=11, file="qstar.dat",status="unknown")
      DO iq=1, nqBZ
       read(11,*)nstar(iq)
      ENDDO
      close(11)

      iqtot=0
      DO iq=1, nqBZ
       filin= "dyndir/" // TRIM(fildyn) // TRIM(int_to_char(iq))
       open(unit=10,file=filin, status="unknown")
       read(10,*)
       read(10,*)
       read(10,*)
       DO ii=1,ntyp
        read(10,*)
       ENDDO
       DO ii=1,nat
        read(10,*)
       ENDDO
       DO jq=1,nstar(iq)
        iqtot=iqtot+1
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)dummy1, dummy2, dummy3, qcart(:,iqtot)
        read(10,*)
        DO imat=1, nat*nat
         read(10,*)iat,jat
         DO ii=1,3
          read(10,*)(dyn_re(ii,jj),dyn_im(ii,jj),jj=1,3)
         ENDDO
         DO ii=1,3
          DO jj=1,3
          dynmat((iat-1)*3+ii,(jat-1)*3+jj,iqtot)                       &
     &                 =dcmplx(dyn_re(ii,jj),dyn_im(ii,jj))
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       close(10)
      ENDDO
      IF(iqtot.ne.nqtot)THEN
       write(*,*)"ERROR: iqtot-nqtot does not match."
       stop
      ENDIF
      qcart(:,:)=2d0*pi*qcart(:,:)
!! check
!      allocate(w_test(nnu))
!      allocate(vector(nnu,nnu))
!      DO inu=1,nnu
!       DO jnu=1,nnu
!        dynmat(inu,jnu,:)=dynmat(inu,jnu,:)/dsqrt(mass(inu)*mass(jnu))
!       ENDDO
!      ENDDO
!      CALL diag_herm_sub(nnu,dynmat(:,:,1),w_test(:),         &
!     &                   vector(:,:) )
!      w_test(:)=dsqrt(w_test(:))
!      write(*,*)"w_test=",w_test(:)
!      deallocate(w_test)
!      deallocate(vector)
!!      stop
!!/check
      IF(la2F.eq.1)THEN
       !! generate mode vector
       CALL gen_modevec_in()
       !!/generate mode vector
       CALL read_lambda()
       !! transform gamma_nu --> gamma_ab
       CALL gamma_nu2Cart()
       !!/transform gamma_nu --> gamma_ab
       
      ENDIF
      end subroutine


      subroutine gamma_nu2Cart()
      USE constants, ONLY: nnu, nqtot,nsig
      USE matrices, ONLY:  modevec_in, mass, gamma_nu, gamma_ab

      complex(8), allocatable :: tmp(:,:)
      complex(8), allocatable :: gamma_nu_l(:,:,:,:)
      real(8), allocatable :: mass_factor(:)
      complex(8), allocatable :: avec(:,:,:)
      complex(8), allocatable :: avec_c(:,:,:)

      integer iq, inu, jnu,isig

      !! local matrices
      allocate(tmp(nnu,nnu))
      allocate(gamma_nu_l(nnu,nnu,nsig,nqtot))
      allocate(mass_factor(nnu))
      allocate(avec(nnu,nnu,nqtot))
      allocate(avec_c(nnu,nnu,nqtot))
      !!/local matrices
      !! global matrices
      allocate(gamma_ab(nnu,nnu,nsig,nqtot))
      !!/global matrices

      DO iq=1, nqtot
       DO isig=1, nsig
        DO jnu=1, nnu
         DO inu=1, nnu
          IF(inu.eq.jnu)THEN 
           gamma_nu_l(inu,jnu,isig,iq)=dcmplx(gamma_nu(inu,isig,iq),0d0)
          ELSE
           gamma_nu_l(inu,jnu,isig,iq)=(0d0,0d0)
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      !! test
      write(*,*)"modevec example"
      write(*,'(6e20.10)')(modevec_in(inu,5,24),inu=1,nnu)
      write(*,*)"gamma_nu"
      write(*,'(3e20.10)')(gamma_nu(inu,1,24),inu=1,nnu)
      !!/test

      !! test
      write(*,*)"gamma_nu_l"
      DO inu=1,3
       write(*,'(6e20.10)')(gamma_nu_l(jnu,inu,1,2),jnu=1,3)
      ENDDO
      !!/test
      

      DO inu =1, nnu
       mass_factor(inu)= dsqrt(mass(inu))
      ENDDO
      DO iq  =1, nqtot
       DO inu =1, nnu
        avec(inu,:,iq)=modevec_in(inu,:,iq)*dcmplx(mass_factor(inu),0d0)
       ENDDO
      ENDDO
      DO iq  =1, nqtot
       DO jnu =1, nnu
        DO inu =1, nnu
         avec_c(inu,jnu,iq)=dconjg(avec(jnu,inu,iq))
        ENDDO
       ENDDO
      ENDDO

      DO iq=1, nqtot
       DO isig=1, nsig
        CALL zgemm('n','n',nnu,nnu,nnu,(1d0,0d0),                       &
     &             gamma_nu_l(:,:,isig,iq),                             &
     &      nnu,avec_c(:,:,iq),nnu,(0d0,0d0),tmp(:,:),nnu)
        CALL zgemm('n','n',nnu,nnu,nnu,(1d0,0d0),avec(:,:,iq),          &
     &      nnu,tmp(:,:),nnu,(0d0,0d0),gamma_ab(:,:,isig,iq),nnu)
       ENDDO
      ENDDO

      !! local matrices
      deallocate(tmp,gamma_nu_l)
      deallocate(mass_factor,avec,avec_c)
      !!/local matrices
      !! test
      write(*,*)"gamma_ab"
      DO inu=1,nnu
       write(*,*)
       write(*,'(6e20.10)')(gamma_ab(jnu,inu,1,24),jnu=1,nnu)
      ENDDO
      !!/test

      end subroutine

      subroutine gen_modevec_in()
      USE constants, ONLY: nnu, nqtot
      USE matrices, ONLY: dynmat, modevec_in, mass
      implicit none

      real(8), allocatable :: w_test(:)
      complex(8), allocatable :: vector(:,:)
      complex(8), allocatable :: dynmat_loc(:,:,:)

      integer inu,jnu, iq

      !! local arrays
      allocate(dynmat_loc(nnu,nnu,nqtot))
      allocate(w_test(nnu))
      allocate(vector(nnu,nnu))
      !!/local arrays
      !! global arrays
      allocate(modevec_in(nnu,nnu,nqtot))
      !!/global arrays
      DO inu=1,nnu
       DO jnu=1,nnu
        dynmat_loc(inu,jnu,:)=dynmat(inu,jnu,:)                         &
     &           /dsqrt(mass(inu)*mass(jnu))
       ENDDO
      ENDDO
      DO iq=1, nqtot
       CALL diag_herm_sub(nnu,dynmat_loc(:,:,iq),w_test(:),         &
     &                   vector(:,:) )
       modevec_in(:,:,iq)=vector(:,:)
      ENDDO
      
      !! local arrays
      deallocate(dynmat_loc)
      !!/local arrays


      end subroutine

      subroutine read_lambda()
      USE filenames, ONLY : file_hr, fildyn, filelph
      USE constants, ONLY: nq1, nq2, nq3, nqBZ, nqtot, nnu, nstar,pi,   &
     &                     nsig
      USE matrices, ONLY: gamma_nu
      implicit none

      !! LOCAL
      character(6), external :: int_to_char
      integer iq, jq, isig
      integer iqtot,  inu
      character(256) filin
      character(256) dummy1, dummy2,dummy3
      real(8), allocatable ::  DOS(:)
      real(8), allocatable ::  omega2(:,:)
      real(8), allocatable ::  lambda(:,:,:)


      !! local arrays
      allocate(DOS(nsig))
      allocate(omega2(nnu,nqtot))
      allocate(lambda(nnu, nsig, nqtot))

      !!/local arrays

      !! global arrays
      allocate(gamma_nu(nnu,nsig,nqtot))
      !!/global arrays


      iqtot=0
      DO iq=1, nqBZ
       iqtot=iqtot+1
       filin= "elph/" // TRIM(filelph) // TRIM(int_to_char(iq))//".dat"
       open(unit=10,file=filin, status="unknown")
       read(10,*)
       read(10,*)
       read(10,*)
       read(10,*)
       read(10,*)
       read(10,*)
       DO isig=1, nsig
       DO inu=1, nnu
        read(10,*)omega2(inu,iqtot),lambda(inu,isig,iqtot)
       ENDDO
       ENDDO
       DO jq=2,nstar(iq)
        lambda(:,:,iqtot+1)=lambda(:,:,iqtot)
        omega2(:,iqtot+1)=omega2(:,iqtot)
        iqtot=iqtot+1
       ENDDO
      ENDDO
      DO iq=1,nqtot
       DO inu=1,nnu
        omega2(inu,iq)=omega2(inu,iq)*omega2(inu,iq)
       ENDDO
      ENDDO
      IF(iqtot.ne.nqtot)THEN
       write(*,*)"ERROR: iqtot-nqtot does not match."
       stop
      ENDIF
      !! read DOS (/Ry/spin/u.c.)
      open(unit=11,file="DOS.dat",status="unknown")
      read(11,*)DOS(:)
      close(11)
      !!/read DOS
      !! change units to . . .
      !!/change units to . . .

      !! calculate gamma_nu from lambda_nu
      DO iqtot=1, nqtot
       DO isig=1, nsig
        DO inu=1, nnu
         gamma_nu(inu,isig,iqtot)=pi*DOS(isig)                          &
     &                            *omega2(inu,iqtot)                    &
     &                            *lambda(inu,isig,iqtot)
        ENDDO
       ENDDO
      ENDDO
      DO isig=1, nsig
       write(*,*)"DOS(",isig,")=",DOS(isig)
      ENDDO
!      DO isig=1, nsig
!       DO inu=1, nnu
!        write(*,'(a10,2i6,a3,f18.8)')"gamma_nu(",inu,isig,")=",         &
!     &                               gamma_nu(inu,isig,80)              &
!     &                               *13.6058d0*241.08d0*1000d0
!       ENDDO
!      ENDDO
      
      !stop
      !! local arrays
      deallocate(DOS,omega2, lambda)
      !!/local arrays
      !!/calculate gamma_nu from lambda_nu
      end subroutine
