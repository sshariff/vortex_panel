      module interfish1 
! ----------------------------------------------------------
      interface
!
      subroutine aij (bpnum,wpnum,body,a,bedges,wedges,bnodex,wnodex,wake, &
                    & oldbndx,right,dt,vref,qinf,puse,nmax)
!
      use panel_type; use interfish2; use interfish3
!
      integer, intent(in)                        :: bpnum, wpnum
      type (panel), dimension(nmax), intent(in)  :: body, wake
      real,         dimension(:,:),  intent(out) :: a
      real,                          intent(in)  :: dt
      type (edges), dimension(nmax), intent(in)  :: bedges, wedges
      type (nodes), dimension(nmax), intent(in)  :: bnodex, wnodex
      type (nodes), dimension(nmax), intent(in)  :: oldbndx
      real,         dimension(nmax), intent(out) :: right
      real,       dimension(nmax,3), intent(out) :: vref
      real,                          intent(in)  :: qinf
      logical,      dimension(nmax), intent(in)  :: puse
      integer,                       intent(in)  :: nmax
!
      end subroutine aij
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine aij1 (bpnum,wpnum,body,a,dt,bedges,wedges,bnodex,wnodex, &
                     & wake,oldbndx,right,vref,qinf,nseplines,       &
		     & seplist,nmax)
!
      use panel_type; use interfish3
!
      integer,                            intent(in)  :: bpnum, wpnum
      type (panel), dimension(nmax),      intent(in)  :: body, wake
      real,         dimension(:,:),       intent(out) :: a
      real,                               intent(in)  :: dt
      type (edges), dimension(nmax),      intent(in)  :: bedges, wedges
      type (nodes), dimension(nmax),      intent(in)  :: bnodex, wnodex
      type (nodes), dimension(nmax),      intent(in)  :: oldbndx
      real,         dimension(nmax),      intent(out) :: right    
      real,         dimension(nmax,3),    intent(out) :: vref
      real,                               intent(in)  :: qinf
      integer,                            intent(in)  :: nseplines
      type (splst), dimension(25),        intent(in)  :: seplist
      integer,                            intent(in)  :: nmax
!
      end subroutine aij1
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine bmotion (dt,qinf,body,bnodes,time,step,bnodex,oldbndx, &
                        & bodyout,blcoll,bldelt,bpnum,tau,tratio,nmax)                   
!
      use panel_type; use interfish2
!
      real,                          intent(out)   :: dt
      real,                          intent(out)   :: qinf
      type (panel), dimension(nmax), intent(inout) :: body
      integer,                       intent(in)    :: bnodes
      real,                          intent(inout) :: time
      integer,                       intent(in)    :: step
      type (nodes), dimension(nmax), intent(inout) :: bnodex
      type (nodes), dimension(nmax), intent(out)   :: oldbndx
      logical,      dimension(nmax), intent(in)    :: bodyout
      type (nodes), dimension(nmax), intent(out)   :: blcoll
      real,         dimension(nmax), intent(inout) :: bldelt
      integer,                       intent(in)    :: bpnum
      real,                          intent(in)    :: tau
      real,                          intent(in)    :: tratio
      integer,                       intent(in)    :: nmax
!
      end subroutine bmotion
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine dump (bednum,wednum,bpnum,wpnum,bedges,wedges,     &
                     & bnodes,wnodes,bnodex,wnodex,maxstep,dumpinc, &
                     & oname,cp,wind,vref,lift,drag,euse,nuse,puse, &
                     & nseplines,seplist,body,wake,time,step,nmax)
!
      use panel_type; use interfish2; use interfish4
!
      integer,                       intent(in) :: bednum, wednum
      integer,                       intent(in) :: bpnum,  wpnum
      type (edges), dimension(nmax), intent(in) :: bedges, wedges
      integer,                       intent(in) :: bnodes, wnodes
      type (nodes), dimension(nmax), intent(in) :: bnodex, wnodex
      integer,                       intent(in) :: maxstep
      integer,                       intent(in) :: dumpinc
      character*30,                  intent(in) :: oname
      real,         dimension(nmax), intent(in) :: cp
      real,       dimension(nmax,3), intent(in) :: wind
      real,       dimension(nmax,3), intent(in) :: vref
      real,                          intent(in) :: lift
      real,                          intent(in) :: drag
      logical ,     dimension(nmax), intent(in) :: euse, nuse, puse
      integer,                       intent(in) :: nseplines
      type (splst), dimension(25),   intent(in) :: seplist
      type (panel), dimension(nmax), intent(in) :: body, wake
      real ,                         intent(in) :: time
      integer,                       intent(in) :: step
      integer,                       intent(in) :: nmax
!
      end subroutine dump
!
      end interface
! ----------------------------------------------------------
      interface
!      
      subroutine kutta (bnodex,wnodes,wnodex,wake,wpnum,wednum, &
                      & wedges,nseplines,oldbndx,lumplist,lumplen,nrlist1,  &
                      & nrlist2,euse,nuse,puse,estack,nstack,pstack,        &
                      & seplist,nmax)
!
      use panel_type; use interfish2; use interfish3
!
      type (nodes), dimension(nmax),   intent(in)    :: bnodex
      integer,                         intent(out)   :: wnodes
      type (nodes), dimension(nmax),   intent(out)   :: wnodex
      type (panel), dimension(nmax),   intent(out)   :: wake
      integer,                         intent(out)   :: wpnum
      integer,                         intent(out)   :: wednum
      type (edges), dimension(nmax),   intent(out)   :: wedges
      integer,                         intent(in)    :: nseplines
      type (nodes), dimension(nmax),   intent(in)    :: oldbndx
      integer,      dimension(25,100), intent(out)   :: lumplist
      integer,      dimension(25),     intent(out)   :: lumplen
      integer,      dimension(25,100), intent(out)   :: nrlist1, nrlist2
      logical,      dimension(nmax),   intent(inout) :: euse,nuse,puse
      type (stack),                    intent(inout) :: estack,nstack,pstack
      type (splst), dimension(25),     intent(inout) :: seplist
      integer,                         intent(in)    :: nmax
!
      end subroutine kutta
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine newpan (bnodes,bnodex,wnodes,wnodex,body,wake,wpnum, &
                       & wednum,wedges,bedges,nseplines,lumplist1,    &
                       & lumplist2,lumplist3,lumplen1,lumplen2,       &
                       & lumplen3,nrlist1,nrlist2,euse,nuse,puse,     &
                       & estack,nstack,pstack,seplist,wlen1,wlen2,    &
                       & wlen3,step,nmax)
!
      use panel_type; use interfish2; use interfish3
!
      type (panel), dimension(nmax), intent(in)    :: body
      type (panel), dimension(nmax), intent(inout) :: wake
      integer, intent(inout)                       :: wednum
      integer, intent(inout)                       :: wpnum
      type (edges), dimension(nmax), intent(in)    :: bedges
      type (edges), dimension(nmax), intent(inout) :: wedges
      integer, intent(in)                          :: bnodes
      integer, intent(inout)                       :: wnodes
      type (nodes), dimension(nmax), intent(in)    :: bnodex
      type (nodes), dimension(nmax), intent(inout) :: wnodex
      integer, intent(in)                          :: nmax
      integer , intent(in)                         :: step
      integer, intent(in)                          :: nseplines
      integer, dimension(25,100), intent(inout)    :: lumplist1
      integer, dimension(25,100), intent(inout)    :: lumplist2
      integer, dimension(25,100), intent(inout)    :: lumplist3
      integer, dimension(25), intent(in)           :: lumplen1
      integer, dimension(25), intent(inout)        :: lumplen2
      integer, dimension(25), intent(inout)        :: lumplen3
      integer, intent(in)                          :: wlen1,wlen2,wlen3
      integer, dimension(25,100), intent(inout)    :: nrlist1, nrlist2
      logical, dimension(nmax), intent(inout)      :: euse, nuse, puse
      type (stack), intent(inout)                  :: estack,nstack
      type (stack), intent(inout)                  :: pstack
      type (splst), dimension(25), intent(inout)   :: seplist
!
      end subroutine newpan
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine presvel (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,  &
                        & puse,cp,wind,cl,cd,bldelt,seplist,nseplines, &
                        & oldcirc,body,wake,vref,step,nmax)
!
      use panel_type; use interfish2; use interfish3; use interfish4
!
      real,                            intent(in)    :: dt
      integer,                         intent(in)    :: bpnum, wpnum
      type (nodes), dimension(nmax),   intent(in)    :: bnodex, wnodex
      type (edges), dimension(nmax),   intent(in)    :: bedges, wedges
      logical ,     dimension(nmax),   intent(in)    :: puse
      real,         dimension(nmax),   intent(out)   :: cp
      real,         dimension(nmax,3), intent(out)   :: wind
      real,                            intent(out)   :: cl, cd
      real,         dimension(nmax),   intent(inout) :: bldelt
      type (splst), dimension(25),     intent(inout) :: seplist
      integer,                         intent(inout) :: nseplines
!      type (panel), dimension(nmax),   intent(in)    :: oldbody
      real,         dimension(nmax),   intent(in)    :: oldcirc
      type (panel), dimension(nmax),   intent(in)    :: body, wake
      real,         dimension(nmax,3), intent(in)    :: vref
      integer,                         intent(in)    :: step
      integer,                         intent(in)    :: nmax
!
      end subroutine presvel
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine rollup (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,wnodes, &
                       & nuse,puse,body,wake,nmax)
!
      use panel_type; use interfish2; use interfish3
!
      real, intent(in)                             :: dt
      integer, intent(in)                          :: bpnum, wpnum
      type (nodes), dimension(nmax), intent(in)    :: bnodex
      type (nodes), dimension(nmax), intent(inout) :: wnodex
      type (edges), dimension(nmax), intent(in)    :: bedges, wedges
      integer, intent(in)                          :: wnodes
      logical , dimension(nmax), intent(in)        :: nuse,puse
      type (panel), dimension(nmax), intent(in)    :: body
      type (panel), dimension(nmax), intent(inout) :: wake
      integer, intent(in)                          :: nmax
!
      end subroutine rollup
!
      end interface
! ----------------------------------------------------------
      interface
!
      subroutine solve (adim,bpnum,kuttas,a,right,body,wake,wedges, &
                      & nseplines,oldcirc,seplist,nmax,bedges,      &
		      & bnodex,wnodex)
! 
      use panel_type; use interfish2; use interfish3
!
! -----------------------------------------------------------------
      integer,                            intent(in)    :: adim
      integer,                            intent(in)    :: bpnum
      integer,                            intent(in)    :: kuttas
      real,         dimension(:,:),       intent(in)    :: a
      real,         dimension(nmax),      intent(in)    :: right
      type (panel), dimension(nmax),      intent(inout) :: body, wake
      type (edges), dimension(nmax),      intent(in)    :: wedges
      integer,                            intent(in)    :: nseplines
      real,         dimension(nmax),      intent(inout) :: oldcirc
      type (splst), dimension(25),        intent(in)    :: seplist
      integer,                            intent(in)    :: nmax
      type (edges), dimension(nmax),      intent(in)    :: bedges
      type (nodes), dimension(nmax),      intent(in)    :: bnodex, wnodex
!
      end subroutine solve
!
      end interface
! ----------------------------------------------------------
      end module interfish1
