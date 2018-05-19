! -------------------------------------------------------------------
      subroutine lump (wake,wpnum,wednum,wedges,nseplines,lumplist,  &
                     & lumplen,farlist,farlen,nrlist,euse,nuse,puse, &
                     & estack,nstack,pstack,nmax)
!
!  This subroutine lumps 8 vortex panels into 2.
!  This is done fairly crudely: the outer corners are used
!  to form two large panels, and the circulation is set
!  to the average of the circulations of the panels which
!  were lumped.
!
!  Here's how it's numbered:
!
!
!   n9----n6-----n3   --eh6-----eh2---      ------------
!   |      |      |   |     / |    / |      | p8/ | p3/ |
!   |      |      |  ev6 ed4 ev3 ed2 ev2    |  /  |  /  |
!   n8----n5-----n2   |  /    | /    |      | / p7| / p2|
!   |      |      |   --eh5-----eh3---      -------------
!   |      |      |   |     / |    / |      | p6/ | p4/ |
!   n7----n4-----n1  ev5 ed3 ev3 ed1 ev1    |  /  |  /  |
!                     |  /    | /    |      | / p5| / p1|
!                     --eh4-----eh1---      -------------
!
!
!   n9-----------n3   -------neh2-----      -------------
!    |           /|   |             /|      |          /|
!    |         /  |   |           /  |      |        /  |
!    |       /    |   |         /    |      |  np2 /    |
!    |     /      |  nev2   ned1   nev1     |     /     |
!    |    /       |   |    /         |      |    / np1  |         
!    |  /         |   |  /           |      |  /        |
!    |/           |   |/             |      |/          |
!   n7-----------n1   -------neh1-----      -------------
!  --------------------------------------------------------------
      use panel_type; use interfish3
!
      implicit none
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
      integer :: i,line,edg
      integer :: ev1,ev2,ev3,ev4,ev5,ev6
      integer :: eh1,eh2,eh3,eh4,eh5,eh6
      integer :: ed1,ed2,ed3,ed4
      integer :: np1,np2,neh1,neh2,nev1,nev2,ned
      integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9
      integer :: p1,p2,p3,p4,p5,p6,p7,p8
      integer :: oldedge
! --------------------------------------------------------------
! input variables
!        nmax: maximum number of edges 
!        nseplines: number of separation lines   
! --------------------------------------------------------------
! output variables
!        estack,nstack,pstack: stacks of element numbers 
!                               available for the wake
!        euse,nuse,puse: lookup tables, is element in use?
!        farlen: number of far edges generated during lumping
!        farlist: list of far edges generated during lumping
!        lumplen: number of edges at lumping sites        
!        lumplist: list of the edges which will be lumped next.
!                  lumplist (i,j) is the jth lumping edge for 
!                  the ith separation line
!        nrlist: list of far edges generated during lumping, used 
!                to reconnect newly lumped panels to old ones
!        wake: wake panels, contains strength of 
!              panel, lists of nodes & edges, "polarity",
!              normal direction and centroid
!        wedges: edge-node array, lists nodes for each edge                       
!        wednum: total number of edges for  wake
!        wpnum: total number of panels for wake
! --------------------------------------------------------------
! local variables
!        i,line,edg: stepping indices
!        element labels: ev1 - vertical edge 1
!                        eh2 - horizontal edge 2
!                        ed1 - diagonal edge 1
!                        n1  - node 1
!                        p1  - panel 1 
!        oldedge: used to connect to previous panels
! --------------------------------------------------------------
! start:
!            write (12,*)
!	    write (12,*) 'in lump....'
!	    write (12,*)
!
!
      do 9999 line = 1,nseplines
!
         farlen(line) = 0
!
         do 8888 edg = 1,lumplen(line),2
!
	    farlen(line) = farlen(line) + 1
! get info
            ev1 = lumplist(line,edg)
	    n1  = wedges(ev1)%node1
	    n2  = wedges(ev1)%node2
	    p1  = wedges(ev1)%panel2
	    ev2 = lumplist(line,edg+1)
	    n3  = wedges(ev2)%node2
	    p2  = wedges(ev2)%panel2
	    ed1 = wake(p1)%edge(1)
	    eh1 = wake(p1)%edge(2)
	    n4  = wedges(ed1)%node2
	    p4  = wedges(ed1)%panel1
	    ev3 = wake(p4)%edge(1)
	    eh3 = wake(p4)%edge(2)
	    n5  = wedges(ev3)%node2
            ed2 = wake(p2)%edge(1)
            p3  = wedges(ed2)%panel1
            ev4 = wake(p3)%edge(1)
            eh2 = wake(p3)%edge(2)
	    n6  = wedges(eh2)%node1
            p5  = wedges(ev3)%panel2
	    ed3 = wake(p5)%edge(1)
	    eh4 = wake(p5)%edge(2)
	    n7  = wedges(eh4)%node1
	    p6  = wedges(ed3)%panel1
	    eh5 = wake(p6)%edge(2)
	    ev5 = wake(p6)%edge(1)
	    n8  = wedges(eh5)%node1
	    p7  = wedges(ev4)%panel2
	    ed4 = wake(p7)%edge(1)
	    p8  = wedges(ed4)%panel1
	    ev6 = wake(p8)%edge(1) 
	    eh6 = wake(p8)%edge(2)
	    n9  = wedges(eh6)%node1
!
!            if ((step.eq.36).and.(nrlist(1,1).eq.0)) then
!	       write (*,*) 'lumplist(1,1) =',lumplist(1,1)
!	       write (*,*) 'lumplist(1,2) =',lumplist(1,2)
!            write (12,*) 'p1,p2,p3,p4 =',p1,p2,p3,p4
!	    write (12,*) 'p5,p6,p7,p8 =',p5,p6,p7,p8
!	    write (12,*)
!	    write (12,*) 'n1,n2,n3 =',n1,n2,n3
!	    write (12,*) 'n4,n5,n6 =',n4,n5,n6
!	    write (12,*) 'n7,n8,n9 =',n7,n8,n9
!	    write (12,*)
!	    write (12,*) 'eh1,eh2,eh3 = ',eh1,eh2,eh3
!	    write (12,*) 'eh4,eh5,eh6 = ',eh4,eh5,eh6
!	    write (12,*)
!	    write (12,*) 'ev1,ev2,ev3 = ',ev1,ev2,ev3
!	    write (12,*) 'ev4,ev5,ev6 = ',ev4,ev5,ev6
!	    write (12,*)
!            write (12,*) 'ed1,ed2,ed3,ed4 =',ed1,ed2,ed3,ed4
!   	    write (12,*)
!	    endif
!
!
! set up new edges
!
! check if we already have the first edge
            if (nrlist(line,farlen(line)).eq.0) then
	       call nulabel(nev1,wednum,euse,estack,nmax)
	    else
	       nev1 = nrlist(line,farlen(line))
	    endif
            wedges(nev1)%node1 = n1
	    wedges(nev1)%node2 = n3
	    call nulabel(nev2,wednum,euse,estack,nmax)
            wedges(nev2)%node1 = n7
	    wedges(nev2)%node2 = n9
	    if (edg.eq.1) then
	       call nulabel(neh1,wednum,euse,estack,nmax)
	       wedges(neh1)%node1 = n7
	       wedges(neh1)%node2 = n1
	    else
	       neh1 = oldedge
	    endif
	    call nulabel(neh2,wednum,euse,estack,nmax)
	    wedges(neh2)%node1 = n9
	    wedges(neh2)%node2 = n3
	    call nulabel(ned,wednum,euse,estack,nmax)
	    wedges(ned)%node1 = n3
	    wedges(ned)%node2 = n7
! set up new panels 
            call nulabel(np1,wpnum,puse,pstack,nmax)
            wake(np1)%node(1) = n3
	    wake(np1)%node(2) = n7
	    wake(np1)%node(3) = n1
	    wake(np1)%edge(1) = ned
	    wake(np1)%edge(2) = neh1
	    wake(np1)%edge(3) = nev1
	    wedges(nev1)%panel2 = np1
	    wedges(ned)%panel2 = np1
	    wedges(neh1)%panel1 = np1
	    if (edg.eq.1) then
	       wedges(neh1)%panel2 = 0
	    endif
	    call nulabel(np2,wpnum,puse,pstack,nmax)
	    wake(np2)%node(1) = n7
	    wake(np2)%node(2) = n9
	    wake(np2)%node(3) = n3
	    wake(np2)%edge(1) = nev2
	    wake(np2)%edge(2) = neh2
	    wake(np2)%edge(3) = ned	    
	    wedges(neh2)%panel2 = np2
	    wedges(ned)%panel1 = np2
	    wedges(nev2)%panel1 = np2
	    wedges(nev2)%panel2 = 0
!
! use average of panel strengths
! we must be careful when considering the small panels who have
! an opposite normal to the new big panel
!            wake(np1)%circ = (  wake(p2)%circ - wake(p4)%circ &
!                            & + wake(p1)%circ + wake(p5)%circ )/4.0 
            wake(np1)%circ = (  wake(p1)%circ + wake(p2)%circ &
                            & + wake(p5)%circ )/3.0 
! 
!            wake(np2)%circ = (  wake(p3)%circ - wake(p7)%circ &
!                            & + wake(p6)%circ + wake(p8)%circ )/4.0
            wake(np2)%circ = (  wake(p3)%circ + wake(p6)%circ &
                            & + wake(p8)%circ)/3.0
!
! 
!
! reconnect small wake panels to large ones
            wedges(ev5)%panel1 = np2
	    wedges(ev6)%panel1 = np2
!
! free up old elements
            if (p1.eq.0.or.p2.eq.0.or.p3.eq.0.or.p4.eq.0.or.p5.eq.0 &
              & .or.p6.eq.0.or.p7.eq.0.or.p8.eq.0) then
               write (*,*) 'panels...'
               write (*,*) p1,p2,p3,p4,p5,p6,p7,p8
	    endif
            call freelab(p1,puse,pstack,nmax)
	    call freelab(p2,puse,pstack,nmax)
            call freelab(p3,puse,pstack,nmax)
 	    call freelab(p4,puse,pstack,nmax)
            call freelab(p5,puse,pstack,nmax)
	    call freelab(p6,puse,pstack,nmax)
            call freelab(p7,puse,pstack,nmax)
	    call freelab(p8,puse,pstack,nmax)
            if (ev1.eq.0.or.ev2.eq.0.or.ev3.eq.0.or.ev4.eq.0)then
               write (*,*) 'evs...'
               write (*,*) ev1,ev2,ev3,ev4
	    endif
	    call freelab(ev1,euse,estack,nmax)
	    call freelab(ev2,euse,estack,nmax)
	    call freelab(ev3,euse,estack,nmax)
	    call freelab(ev4,euse,estack,nmax)
! we can't free ev5 and ev6; they are still used by the small panels
	    call freelab(eh1,euse,estack,nmax)
!
! we can't free eh2, eh6, or n6; they are still in use
!	    call freelab(eh2,euse,estack,nmax)
            if (eh1.eq.0.or.eh3.eq.0.or.eh4.eq.0.or.eh5.eq.0)then
               write (*,*) 'ehs...'
               write (*,*) eh1,eh3,eh4,eh5
	    endif
	    call freelab(eh3,euse,estack,nmax)
	    call freelab(eh4,euse,estack,nmax)
	    call freelab(eh5,euse,estack,nmax)
!	    call freelab(eh6,euse,estack,nmax)
            if (ed1.eq.0.or.ed2.eq.0.or.ed3.eq.0.or.ed4.eq.0)then
               write (*,*) 'eds...'
               write (*,*) ed1,ed2,ed3,ed4
	    endif
	    call freelab(ed1,euse,estack,nmax)
	    call freelab(ed2,euse,estack,nmax)
	    call freelab(ed3,euse,estack,nmax)
	    call freelab(ed4,euse,estack,nmax)
            if (n2.eq.0.or.n4.eq.0.or.n5.eq.0)then
               write (*,*) 'ns...'
               write (*,*) n2,n4,n5
	    endif
	    call freelab(n2,nuse,nstack,nmax)
	    call freelab(n4,nuse,nstack,nmax)
	    call freelab(n5,nuse,nstack,nmax)
!	    call freelab(n6,nuse,nstack,nmax)
!
! keep track of bottom edge and panel for next lumping
            oldedge = neh2
!
! update lumplists
            lumplist(line,edg) = ev5
	    lumplist(line,edg+1) = ev6
	    farlist(line,farlen(line)) = nev1
	    nrlist(line,farlen(line)) = nev2
!	    
!	    write (*,*) 'oldedge =',oldedge
 8888    continue
!
! for last edge, free up last two nodes and edge
         if (eh2.eq.0.or.eh6.eq.0.or.n6.eq.0)then
	    write (*,*) 'last stuff'
	    write (*,*) eh2,eh6,n6
	 endif
         call freelab(eh2,euse,estack,nmax)
         call freelab(eh6,euse,estack,nmax)
         call freelab(n6,nuse,nstack,nmax) 
!
 9999 continue
      return
      end
! --------------------------------------------------------------------

