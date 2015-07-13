      program main
      implicit none

      CALL read_in()
      ! read input.txt
      CALL read_dynmat()


      ! initialize constant module
      CALL read_Rp()
      !/initialize constant module

      ! convert input data unit, coordinates, etc.
      CALL convert_in()


      CALL ph_InvFT_D()
      write(*,*)"ph inverse FT end"
      CALL ph_ifc_out()
      end program
