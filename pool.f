! -------------------------------------------------------------------
      subroutine pool (wake,wedges,nseplines,lumplist,lumplen,euse, &
                     & nuse,puse,estack,nstack,pstack,step,nmax)
!
!  This subroutine lumps the last two panels into one.
!
!  Here's how it's numbered:
!
!
!   n6----n4-----n2  
!   |    / |    / |   --eh4-----eh2---      -------------
!   | /    | /    |   |     / |    / |      | p4/ | p2/ |
!   n5----n3-----n1  ev3 ed2 ev2 ed1 ev1    |  /  |  /  |
!                     |  /    | /    |      | / p3| / p1|
!                     --eh3-----eh1---      -------------
!
!
!   n6-----n4   ---eh4---       --------
!    |     /|   |       /|      | p4  /|
!    |   /  |  ev3  ed2  |      |    / |
!    |  /   |   |  /    ev2     |  /   |
!    |/     |   |/       |      |/ p3  |
!   n5------n3   ---eh3---      ---------
!  --------------------------------------------------------------
      use panel_type; use interfish3
!
      implicit none
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

      integer :: i,line,edg
      integer :: ev1,ev2,ev3
      integer :: eh1,eh2,eh3,eh4
      integer :: ed1,ed2
      integer :: np1,np2,neh1,neh2,ned
      integer :: n1,n2,n3,n4,n5,n6
      integer :: p1,p2,p3,p4
! --------------------------------------------------------------
! input variables
!        lumplen: number of edges at lumping sites
!        nmax: maximum number of edges 
!        nseplines: number of separation lines
!        step: current step number
!        wedges: edge-node array, lists nodes for each edge                       
! --------------------------------------------------------------
! output variables
!        estack,nstack,pstack: stacks of element numbers 
!                               available for the wake
!        euse,nuse,puse: lookup tables, is element in use?
!        lumplist: list of the edges which will be lumped next.
!                  lumplist (i,j) is the jth lumping edge for 
!                  the ith separation line
!        wake: wake panels, contains strength of 
!              panel, lists of nodes & edges, "polarity",
!              normal direction and centroid
! --------------------------------------------------------------
! local variables
!        i,line,edg: stepping indices
!        element labels: ev1 - vertical edge 1
!                         eh2 - horizontal edge 2
!                         ed1 - diagonal edge 1
!                         n1  - node 1
!                         p1  - panel 1 
! --------------------------------------------------------------
! start:
!            write (12,*)
!	    write (12,*) 'in lump....'
!	    write (12,*)
!
!
      do 9999 line = 1,nseplines
!
         do 8888 edg = 1,lumplen(line)
!
! get info
            ev1 = lumplist(line,edg)
	    n1  = wedges(ev1)%node1
	    n2  = wedges(ev1)%node2
	    p1  = wedges(ev1)%panel2
	    eh1 = wake(p1)%edge(2)
	    ed1 = wake(p1)%edge(1)
	    n3  = wedges(eh1)%node1
	    p2  = wedges(ed1)%panel1
 	    ev2 = wake(p2)%edge(1)
            eh2 = wake(p2)%edge(2)
	    n4  = wedges(ev2)%node2
            p3  = wedges(ev2)%panel2
            ed2 = wake(p3)%edge(1)
	    eh3 = wake(p3)%edge(2)
	    p4  = wedges(ed2)%panel1
	    ev3 = wake(p4)%edge(1)
	    eh4 = wake(p4)%edge(2)
	    n5  = wedges(ev3)%node1
	    n6  = wedges(ev3)%node2
!
! use average of panel strengths
! we must be careful when considering the small panels who have
! an opposite normal to the new big panel
!
! here np1 will consist of p1 and parts of p2 and p3
            wake(p3)%circ = &
	      & ( 2.0*wake(p1)%circ - wake(p2)%circ + wake(p3)%circ )/4.0
! here np2 will consist of p4 and parts of p2 and p3  
            wake(p4)%circ = & 
	      & ( 2.0*wake(p4)%circ - wake(p3)%circ + wake(p2)%circ )/4.0
!
! free up old elements
            call freelab(p1,puse,pstack,nmax)
	    call freelab(p2,puse,pstack,nmax)
	    call freelab(ev1,euse,estack,nmax)
	    call freelab(eh1,euse,estack,nmax)
	    call freelab(ed1,euse,estack,nmax)
	    call freelab(n1,nuse,nstack,nmax)
!
! update lumplist
            lumplist(line,edg) = ev2
!
 8888    continue
!
! for last edge, free up last two nodes and edge
         call freelab(eh2,euse,estack,nmax)
         call freelab(n2,nuse,nstack,nmax) 
!
 9999 continue
      return
      end
! --------------------------------------------------------------------
