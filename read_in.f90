      subroutine read_in()
      USE filenames, ONLY : file_hr, fildyn,filelph,filmass
      USE flags, ONLY : la2F
      USE constants, ONLY: avec1, avec2, avec3, bvec1, bvec2, bvec3,    &
     &                     nsig
      implicit none

      !! filenames
      read(5,*)file_hr
      read(5,*)fildyn
      read(5,*)filelph
      read(5,*)filmass
      read(5,*)la2F
      read(5,*)nsig
      !read(5,*)file_q_out
      
      !! lattice
      read(5,*)avec1(:)
      read(5,*)avec2(:)
      read(5,*)avec3(:)
      read(5,*)bvec1(:)
      read(5,*)bvec2(:)
      read(5,*)bvec3(:)
      !!/lattice
      end subroutine
