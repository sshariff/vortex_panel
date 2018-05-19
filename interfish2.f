      module interfish2
! ----------------------------------------------------------
      interface
!
      subroutine blchar (bpnum,body,cp,bldelt,seplist,nseplines,nmax)
!
      use panel_type
!
      integer,                       intent(in)    :: bpnum
      type (panel), dimension(nmax), intent(in)    :: body
      real,         dimension(nmax), intent(in)    :: cp
      real,         dimension(nmax), intent(out)   :: bldelt
      type (splst), dimension(25),   intent(inout) :: seplist
      integer,                       intent(inout) :: nseplines
      integer,                       intent(in)    :: nmax
!
      end subroutine blchar
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine diffuse (wake,wpnum,puse,nmax)
!
      use panel_type
!
      type (panel), dimension(nmax), intent(inout) :: wake
      integer, intent(in)                          :: wpnum
      logical , dimension(nmax), intent(in)        :: puse
      integer, intent(in)                          :: nmax
!
      end subroutine diffuse
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine gaussj (a,n,np,b,m,mp)
!
      real, dimension(np,np), intent(inout) :: a
      integer,                  intent(in)    :: n,np
      real, dimension(np,mp), intent(inout) :: b
      integer,                  intent(in)    :: m,mp
!
      end subroutine gaussj
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine gauss90 (a,b)
!
      double precision, dimension(:,:), intent(inout) :: a,b
!
      end subroutine gauss90
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine lubksb (a,n,np,indx,b)
!
      integer,                intent(in)    :: n,np,indx(n)
      real, dimension(np,np), intent(in)    :: a
      real, dimension(n),     intent(inout) :: b
!
      end subroutine lubksb
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine lubksb90 (a,indx,b)
!
      double precision,    dimension(:,:), intent(in)    :: a
      integer, dimension(:),   intent(in)    :: indx
      double precision,    dimension(:),   intent(inout) :: b
!
      end subroutine lubksb90
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine ludcmp (a,n,np,indx,d)
!
      integer,                intent(in)    :: n,np
      integer, dimension(n),  intent(out)   :: indx
      real,                   intent(out)   :: d      
      real, dimension(np,np), intent(inout) :: a
 !
      end subroutine ludcmp
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine lump (wake,wpnum,wednum,wedges,nseplines,lumplist,  &
                     & lumplen,farlist,farlen,nrlist,euse,nuse,puse, &
                     & estack,nstack,pstack,nmax)
!
      use panel_type; use interfish3
!
      type (panel), dimension(nmax), intent(inout) :: wake
      integer,                       intent(inout) :: wpnum
      integer,                       intent(inout) :: wednum
      type (edges), dimension(nmax), intent(inout) :: wedges
      integer,                       intent(in)    :: nseplines
      integer, dimension(25,100),    intent(inout) :: lumplist
      integer, dimension(25),        intent(in)    :: lumplen
      integer, dimension(25,100),    intent(out)   :: farlist
      integer, dimension(25),        intent(out)   :: farlen
      integer, dimension(25,100),    intent(inout) :: nrlist
      logical, dimension(nmax),      intent(inout) :: euse, nuse, puse
      type (stack),                  intent(inout) :: estack,nstack,pstack
      integer,                       intent(in)    :: nmax
!
      end subroutine lump
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine nodeave (ednum,surface,surfedges,cp,wind,vref,ncp,nq, &
                        & ncirc,bodtest,euse,nmax)
!
      use panel_type
!
      integer,                       intent(in)  :: ednum
      type (panel), dimension(nmax), intent(in)  :: surface
      type (edges), dimension(nmax), intent(in)  :: surfedges
      real,         dimension(nmax), intent(in)  :: cp
      real,       dimension(nmax,3), intent(in)  :: wind
      real,       dimension(nmax,3), intent(in)  :: vref
      real,         dimension(nmax), intent(out) :: ncp
      real,       dimension(nmax,3), intent(out) :: nq
      real,         dimension(nmax), intent(out) :: ncirc
      logical,                       intent(in)  :: bodtest
      logical,      dimension(nmax), intent(in)  :: euse
      integer,                       intent(in)  :: nmax
!
      end subroutine nodeave
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine norms (bnodes,bnodex,body,bodyout,bpnum,nmax)
!
      use panel_type; use interfish4
!
      integer,                       intent(in)    :: bnodes
      type (nodes), dimension(nmax), intent(in)    :: bnodex
      type (panel), dimension(nmax), intent(inout) :: body
      logical,      dimension(nmax), intent(in)    :: bodyout
      integer,                       intent(in)    :: bpnum
      integer,                       intent(in)    :: nmax
!
      end subroutine norms
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine pool (wake,wedges,nseplines,lumplist,lumplen,euse, &
                     & nuse,puse,estack,nstack,pstack,step,nmax)
!
      use panel_type; use interfish3
!
      type (panel), dimension(nmax), intent(inout) :: wake
      type (edges), dimension(nmax), intent(in)    :: wedges
      integer,                       intent(in)    :: nseplines
      integer, dimension(25,100),    intent(inout) :: lumplist
      integer, dimension(25),        intent(in)    :: lumplen
      logical, dimension(nmax),      intent(inout) :: euse, nuse, puse
      type (stack),                  intent(inout) :: estack,nstack,pstack
      integer,                       intent(in)    :: step
      integer,                       intent(in)    :: nmax
!
      end subroutine pool
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine rhs (bpnum,wpnum,body,wake,bnodex,oldbndx,wnodex,wedges, &
                    & right,dt,vref,qinf,puse,nmax)
!
      use panel_type; use interfish3
!
      integer, intent(in)                        :: bpnum, wpnum
      type (panel), dimension(nmax), intent(in)  :: body, wake
      type (nodes), dimension(nmax), intent(in)  :: bnodex,oldbndx,wnodex
      type (edges), dimension(nmax), intent(in)  :: wedges
      real,         dimension(nmax), intent(out) :: right
      real,                          intent(in)  :: dt
      real,       dimension(nmax,3), intent(out) :: vref
      real,                          intent(in)  :: qinf
      logical,      dimension(nmax), intent(in)  :: puse
      integer,                       intent(in)  :: nmax
!
      end subroutine rhs
!
      end interface
! ----------------------------------------------------------
!      interface
!
!      subroutine wnorms (wnodes,wnodex,wake,wpnum,puse,nmax)
!
!      use panel_type; use interfish4
!
!      integer,                       intent(in)    :: wnodes
!      type (nodes), dimension(nmax), intent(in)    :: wnodex
!      type (panel), dimension(nmax), intent(inout) :: wake
!      integer,                       intent(in)    :: wpnum
!      logical,      dimension(nmax), intent(in)    :: puse
!      integer,                       intent(in)    :: nmax
!
!      end subroutine wnorms
!
!      end interface
! ----------------------------------------------------------
      end module interfish2
