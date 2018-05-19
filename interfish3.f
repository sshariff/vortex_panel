      module interfish3
! ----------------------------------------------------------
      interface
!
      subroutine freelab (label,use,stak,nmax)
!
      use panel_type
!
      integer,      intent(out)   :: label
      logical,      intent(inout) :: use(nmax)
      type (stack), intent(inout) :: stak
      integer,      intent(in)    :: nmax
!
      end subroutine freelab
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine lnvortx (xp,x1,x2,u)
!
      use interfish4
!
      real, dimension(3), intent(in)  :: xp,x1,x2
      real, dimension(3), intent(out) :: u
!
      end subroutine lnvortx
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine nulabel (label,max,use,stak,nmax)
!
      use panel_type
!
      integer,      intent(out)   :: label
      integer,      intent(inout) :: max
      logical,      intent(inout) :: use(nmax)
      type (stack), intent(inout) :: stak
      integer,      intent(in)    :: nmax
!
      end subroutine nulabel
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine panvel (bpnum,body,bnodex,oldbndx,right,dt,vref,qinf,nmax)
!
      use panel_type; use interfish4
!
      integer,                         intent(in)  :: bpnum
      type (panel), dimension(nmax),   intent(in)  :: body
      type (nodes), dimension(nmax),   intent(in)  :: bnodex
      type (nodes), dimension(nmax),   intent(in)  :: oldbndx
      real,         dimension(nmax),   intent(out) :: right
      real,                            intent(in)  :: dt
      real,         dimension(nmax,3), intent(out) :: vref
      real,                            intent(in)  :: qinf
      integer,                         intent(in)  :: nmax
!
      end subroutine panvel
!
      end interface
! ----------------------------------------------------------
      end module interfish3
