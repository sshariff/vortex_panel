      module interfish4 
! ----------------------------------------------------------
      interface
!
      subroutine cross (r1, r2, r1xr2, magx)
!
      real, dimension(3), intent(in)  :: r1, r2
      real, dimension(3), intent(out) :: r1xr2
      real,               intent(out) :: magx
!
      end subroutine cross
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine kine (bpnum,bnodex,oldbndx,dt,vref,qinf,nmax)
!
      use panel_type
!
      integer,                         intent(in)  :: bpnum
      type (nodes), dimension(nmax),   intent(in)  :: bnodex
      type (nodes), dimension(nmax),   intent(in)  :: oldbndx
      real,                            intent(in)  :: dt
      real,         dimension(nmax,3), intent(out) :: vref
      real,                            intent(in)  :: qinf
      integer,                         intent(in)  :: nmax
!
      end subroutine kine
!
      end interface
! ----------------------------------------------------------
      end module interfish4
