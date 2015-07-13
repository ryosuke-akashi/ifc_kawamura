      subroutine read_Rp()
      USE filenames, ONLY : file_hr
      USE constants, ONLY: nr,nnu,weight
      USE matrices,  ONLY: wr, rcrys, rcart,ifc
      implicit none

      integer nbase
      integer ir
      integer ibase,jbase
     ! read input nr, nbase
      open(unit=11, file=file_hr,status="unknown")
      rewind(11)
      read(11,*)
      read(11,*)nbase
      read(11,*)nr
      write(*,*)"nbase, nr=",nbase,nr
      !/read input nr, nbase

      allocate(wr(nr))
      allocate(rcrys(3,nr))
      allocate(rcart(3,nr))
      allocate(ifc(nnu,nnu,nr))
      rewind(11)
      read(11,*)
      read(11,*)
      read(11,*)
      read(11,*)wr(:)
      DO ir = 1, nr
       DO jbase = 1, nbase   
       DO ibase = 1, nbase   
        read(11,*)rcrys(1,ir), rcrys(2,ir), rcrys(3,ir)
       ENDDO
       ENDDO
      ENDDO
      close(11)
      weight=0d0
      DO ir=1, nr
       weight=weight+1d0/dble(wr(ir))
      ENDDO

      end subroutine
