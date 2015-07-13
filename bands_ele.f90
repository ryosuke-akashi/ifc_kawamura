      program branch_phonon
      implicit none
      character(5)  gomi1,gomi2,gomi3
      !input
      character(20) :: inputfile="freq.dat"
      !/input
      integer  nbnd, nks   
      integer  nk_dummy
      character(10) char_dummy
      integer  nk,i,ibnd,i_ene
      integer  jbnd
      
      real(8),allocatable :: branch_q(:,:),abs_q(:),branch_w(:,:)
      real(8),allocatable :: gam(:,:,:)
      real(8),allocatable :: broad(:)
      real(8) ene
      real(8)  del_ene 
      real(8) temp
      !
      !read phonon branches
      !
      open(unit=100,file=trim(inputfile),status="unknown")
      read(100,*)gomi1,gomi2,nbnd,gomi3,nks

      write(*,*)"gomi1=",gomi1,"  gomi2=",gomi2,"  gomi3=",gomi3
      write(*,*)"nbnd=",nbnd,"nks=",nks

      allocate(branch_q(nks,3))
      allocate(branch_w(nks,nbnd))


      DO nk = 1,nks 
        read(100,*) (branch_q(nk,i),i=1,3)
        read(100,*) (branch_w(nk,ibnd),ibnd=1,nbnd)
      ENDDO

      !write(*,*)"branch_q(2,*)=",(branch_q(2,i),i=1,3)
      !write(*,*)"branch_w(2,*)=",(branch_w(2,ibnd),ibnd=1,nbnd)

      close(100)

      DO nk=1, nks
       DO ibnd=1, nbnd
        DO jbnd=1, nbnd-ibnd
         IF(branch_w(nk,jbnd).gt.branch_w(nk,jbnd+1))THEN
          temp=branch_w(nk,jbnd+1)
          branch_w(nk,jbnd+1)=branch_w(nk,jbnd)
          branch_w(nk,jbnd)=temp
         ENDIF
        ENDDO
       ENDDO
      ENDDO
      !
      !read gammas
      !
    
      allocate(abs_q(nks))

      abs_q(1) = 0d0
      DO nk=2 , nks
        abs_q(nk)=abs_q(nk-1) &
                   + sqrt( ( branch_q(nk,1)-branch_q(nk-1,1) )**2d0 &
           &              +( branch_q(nk,2)-branch_q(nk-1,2) )**2d0 &
           &              +( branch_q(nk,3)-branch_q(nk-1,3) )**2d0 )   
      ENDDO


      open(unit=300,file="out.dat",status="unknown")
      DO nk=1,nks
        write(300,'(f16.8,100f16.8)')abs_q(nk),                         &
     &                        (branch_w(nk,ibnd),ibnd=1,nbnd)  
      ENDDO

      close(300)

      deallocate(branch_q)
      deallocate(branch_w)
      deallocate(abs_q)

      end program
