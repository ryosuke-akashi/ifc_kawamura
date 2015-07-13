      subroutine convert_in()
      USE constants,    ONLY : nr
      USE matrices,  ONLY : rcrys, rcart
      implicit none

      ! convert input k, q, r, rp to Cartesian (angstrom) coordinate
      CALL cryst_to_cart(nr,rcrys, rcart,-1)

      end subroutine
