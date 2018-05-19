! -------------------------------------------------------------------
      subroutine norms (bnodes,bnodex,body,bodyout,bpnum,nmax)
!
!  This subroutine calculates the unit normal vector and
!  centroid for the body panels
!--------------------------------------------------------------
!
      use panel_type; use interfish4
!
      implicit none
!
      integer,                       intent(in)    :: bnodes
      type (nodes), dimension(nmax), intent(in)    :: bnodex
      type (panel), dimension(nmax), intent(inout) :: body
      logical,      dimension(nmax), intent(in)    :: bodyout
      integer,                       intent(in)    :: bpnum
      integer,                       intent(in)    :: nmax

      real    :: magx
      integer :: n,npanel
      integer :: node1,node2,node3
      real    :: r1(3),r2(3)
      real    :: r1xr2(3)
      real    :: sgn
! ---------------------------------------------------------------      
! input variables
!        bnodes : total number of nodes for body
!        bnodex : node coordinate array (x,y,z for each node)
!        bodyout: logical array used to make body normals point out
!        bpnum: total number of panels for body 
!        nmax: maximum number of edges
! ---------------------------------------------------------------
! output variables
!        body: body panels, contains strength of panel, lists
!              of nodes & edges, "polarity", normal direction
!              and centroid      
! ---------------------------------------------------------------
! local variables
!        magx: magnitude of panel normal vector
!        n,npanel: indices for stepping
!        node1,node2,node3: node numbers for the panel
!        r1,r2: components of two of the panel edges
!        r1xr2: components of panel normal vector
!        sgn: used to switch sign of the normal for panels with in=true
! --------------------------------------------------------------
!  start:
!
      do 150 npanel = 1, bpnum
!
!  locate nodes associated with this panel
        node1 = body(npanel)%node(1)
	node2 = body(npanel)%node(2)
	node3 = body(npanel)%node(3)
!
!  compute vectors for two edges
!
!  vector pointing from the first node of the
!  panel to the second node of the panel:
	r1 = bnodex(node2)%x - bnodex(node1)%x
!
!  vector pointing from the first node of the
!  panel to the third node of the panel:
	r2 = bnodex(node3)%x - bnodex(node1)%x
!
!  compute cross product of two edges
        call cross (r1, r2, r1xr2, magx)
!
!  normal vector -> unit normal vector, pointing out of the body.
!    but we have to switch the sign for panels with in=.true.
        if (bodyout(npanel)) then
	   sgn =  1.0
	else
	   sgn = -1.0
        end if
!
        body(npanel)%normal = sgn*r1xr2/magx
!
!  compute centroid
        body(npanel)%centr = &
	    & ( bnodex(node1)%x + bnodex(node2)%x + bnodex(node3)%x )/3.0
!
 150  continue
!            
      return
      end
! --------------------------------------------------------------------
