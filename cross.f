! --------------------------------------------------------------------
      subroutine cross (r1, r2, r1xr2, magx)
!
!  This subroutine calculates the cross product of 2 three-dimensional
!  vectors, r1 and r2, and outputs the answer in variable r1xr2.
!
! -----------------------------------------------------------------
      implicit none
!
      real, dimension(3), intent(in)  :: r1, r2
      real, dimension(3), intent(out) :: r1xr2
      real, intent(out)               :: magx
! -----------------------
!  input variables
!        r1,r2: the input vectors r1 and r2
! ----------------------- 
!  output variables
!        r1xr2: r1 x r2
!        magx: magnitude of the cross product
! -----------------------------------
!  start:
      r1xr2(1) =   r1(2)*r2(3) - r1(3)*r2(2)
      r1xr2(2) = -(r1(1)*r2(3) - r1(3)*r2(1))
      r1xr2(3) =   r1(1)*r2(2) - r1(2)*r2(1)
      
      magx = sqrt( r1xr2(1)*r1xr2(1) + r1xr2(2)*r1xr2(2) + r1xr2(3)*r1xr2(3) )
      
      return
      end
! --------------------------------------------------------------------
