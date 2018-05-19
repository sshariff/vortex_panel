! ---------------------------------------------------------------
! ---------------------------------------------------------------
       module triface
! ---------------------------------------------------------------
       interface
       subroutine trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
! ------------------------------
       use panel_type
! ------------------------------
       real,    dimension(3), intent(in) :: x1,x2,x3,x4
       integer,               intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2
       integer,               intent(in) :: nmax
       logical,      dimension(nmax), intent(inout) :: bodyin
       type (nodes), dimension(nmax), intent(inout) :: bnodex
       type (panel), dimension(nmax), intent(inout) :: body
       type (edges), dimension(nmax), intent(inout) :: bedges
! ------------------------------
       end subroutine trir
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine tril (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
! ------------------------------
       use panel_type
! ------------------------------
       real,    dimension(3), intent(in) :: x1,x2,x3,x4
       integer,               intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2
       integer,               intent(in) :: nmax
       logical,      dimension(nmax), intent(inout) :: bodyin
       type (nodes), dimension(nmax), intent(inout) :: bnodex
       type (panel), dimension(nmax), intent(inout) :: body
       type (edges), dimension(nmax), intent(inout) :: bedges
! ------------------------------
       end subroutine tril
       end interface
! ---------------------------------------------------------------
       end module triface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       module wraps
! ---------------------------------------------------------------
       interface
       subroutine rotate (n1,n2,n3,n4,e1,e2,e3,e4,x1,x2,x3,x4)
! ------------------------------
       implicit none
!
       integer,            intent(inout) :: n1,n2,n3,n4,e1,e2,e3,e4
       real, dimension(3), intent(inout) :: x1,x2,x3,x4
!
       end subroutine rotate
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine wrap1 (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	               & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		       & nextnode,nextedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
!
       end subroutine wrap1
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine wrap (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	              & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		      & nextnode,nextedge,newnode,newedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       integer,      dimension(nmax),       intent(inout) :: newnode,newedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
!
       end subroutine wrap
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine wrapmax (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	                 & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		         & nextnode,nextedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
!
       end subroutine wrapmax
       end interface
! ---------------------------------------------------------------
       end module wraps
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       module packblock
! ---------------------------------------------------------------
       interface
       subroutine caudalize (i,ireal,jmaxi,totjmax,x,y,z,xx,yy,zz, &
                           & kkutta,ccaudal,nsmall)
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: i
       integer,                           intent(in)    :: ireal
       integer,                           intent(in)    :: jmaxi
       integer, dimension(nsmall),        intent(in)    :: totjmax
       real,    dimension(nsmall,nsmall), intent(in)    :: x,y,z
       real,    dimension(nsmall,nsmall), intent(inout) :: xx,yy,zz
       logical, dimension(nsmall,nsmall), intent(inout) :: kkutta,ccaudal
       integer,                           intent(in)    :: nsmall  
!
       end subroutine caudalize
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine combine (i,ireal,jmaxi,totjmax,x,y,z,xx,yy,zz,kutta, &
                         & kkutta,nsmall,jjmax,strip)
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: i
       integer,                           intent(in)    :: ireal
       integer,                           intent(in)    :: jmaxi
       integer, dimension(nsmall),        intent(inout) :: totjmax
       real,    dimension(nsmall,nsmall), intent(in)    :: x,y,z
       real,    dimension(nsmall,nsmall), intent(inout) :: xx,yy,zz
       logical, dimension(nsmall,nsmall), intent(in)    :: kutta
       logical, dimension(nsmall,nsmall), intent(inout) :: kkutta
       integer,                           intent(in)    :: nsmall  
       integer, dimension(nsmall,3),      intent(inout) :: jjmax
       integer,                           intent(inout) :: strip
!
       end subroutine combine
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       end module packblock
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       module genblock
! ---------------------------------------------------------------
       interface
       subroutine nextstrip (ireal,imax,jmax,bnodex,body,bedges,bodyin,  &
                           & kutta,caudal,xx,n,e,p,nmax,nsmall,oddy,     &
			   & newnode,newedge,nextnode,nextedge,djl,knode)
! ------------------------------
       use quadface; use panel_type; use wraps
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: ireal,imax
       integer, dimension(nsmall,3),        intent(in)    :: jmax
       logical, dimension(nmax),            intent(inout) :: bodyin
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta,caudal
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nsmall,nsmall),   intent(in)    :: oddy
       integer, dimension(nmax),            intent(inout) :: newnode,newedge
       integer, dimension(nmax),            intent(inout) :: nextnode,nextedge
       integer,                             intent(in)    :: djl
       logical,      dimension(nmax),       intent(inout) :: knode
!
       end subroutine nextstrip
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine notchstuff (i,j,outrn,re,nextnode,nextedge,bnodex, &
                            & body,bedges,bodyin,xx,n,e,p,nmax, &
                            & nsmall,oddy,jstart,jstop)
! ------------------------------
       use panel_type       
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: i,j
       integer,      dimension(nsmall),     intent(inout) :: outrn
       integer,      dimension(nsmall),     intent(inout) :: re
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical,      dimension(nmax),       intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       logical, dimension(nsmall,nsmall),   intent(in)    :: oddy
       integer,      dimension(nmax),       intent(inout) :: jstart,jstop
!
       end subroutine notchstuff
       end interface
! ---------------------------------------------------------------
       end module genblock
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       module sunblock
! ---------------------------------------------------------------
       interface
       subroutine pack (imax,jmax,kutta,caudal,x,y,z,nmax,nsmall,oddy)
! ------------------------------
       use quadface; use panel_type; use packblock
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: imax
       integer, dimension(nsmall,3),      intent(inout) :: jmax
       logical, dimension(nsmall,nsmall), intent(inout) :: kutta,caudal
       real,    dimension(nsmall,nsmall), intent(inout) :: x,y,z
       integer,                           intent(in)    :: nmax,nsmall
       logical, dimension(nsmall,nsmall), intent(out)   :: oddy
!
       end subroutine pack
       end interface
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       interface
       subroutine gen (imax,jmax,bnodex,body,bedges,bodyin,bottom,kutta, &
                     & caudal,xx,n,e,p,nmax,nsmall,le,re,inle,inre,      &
		     & outln,outrn,inln,inrn,oddy,knode)
! ------------------------------
       use quadface; use panel_type; use genblock; use wraps
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: imax
       integer, dimension(nsmall,3),        intent(in)    :: jmax
       logical, dimension(nmax),            intent(inout) :: bodyin
       logical,                             intent(in)    :: bottom
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta,caudal
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(out)   :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       logical, dimension(nsmall,nsmall),   intent(in)    :: oddy
       logical,      dimension(nmax),       intent(inout) :: knode
!
       end subroutine gen
       end interface
! ---------------------------------------------------------------
       end module sunblock
! ---------------------------------------------------------------
! ---------------------------------------------------------------
