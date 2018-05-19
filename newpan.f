! -------------------------------------------------------------------
      subroutine newpan (bnodes,bnodex,wnodes,wnodex,body,wake,wpnum, &
                       & wednum,wedges,bedges,nseplines,lumplist1,    &
                       & lumplist2,lumplist3,lumplen1,lumplen2,       &
                       & lumplen3,nrlist1,nrlist2,euse,nuse,puse,     &
                       & estack,nstack,pstack,seplist,wlen1,wlen2,    &
                       & wlen3,step,nmax)
!
!  This subroutine connects the nodes to create the edges and 
!  elements for the new timestep. New panels are created at
!  separation locations, and some vortices are lumped away from 
!  the body. Also calculates the normal direction for the
!  wake panels.
!
!  At each separation edge, there are two coincident edges
!  (as well as corresponding nodes), one associated with the body
!  and one associated with the wake. These are referred to here
!  as the 'body-side edge' and the 'wake-side edge', respectively.
!  At each timestep, the wake-side nodes and edges are convected
!  according to the induced velocity. This is now referred
!  to as the old wake-side edge. A new wake-side edge is then
!  created at the same position as the body-side edge. These
!  edges are then connected with two new panels.
!
!  Here's how it's numbered:
!
!     nrnode1     nrnode2 
!     nwsnd2      wsnode2
!        --------------      -----nredge-----      -------------
!       ||           /|     ||             /|     ||          /|
!  sep  ||         /  |     ||           /  |     || nearpan/  |
!  edge ||       /    |     ||         /    |     ||      /    |
!       ||     /      |    nwside  diedge  wside  ||     /     |
!       ||    /       |     ||    /         |     ||    /      |         
!       ||  /         |     ||  /           |     ||  / farpan |
!       ||/           |     ||/             |     ||/          |
!        --------------      -----fredge-----      -------------
!     nwsnd1      wsnode1
!     frnode1     frnode2
!  --------------------------------------------------------------
      use panel_type; use interfish2; use interfish3
!
      implicit none
!
      type (panel), dimension(nmax), intent(in)    :: body
      type (panel), dimension(nmax), intent(inout) :: wake
      integer,                       intent(inout) :: wednum
      integer,                       intent(inout) :: wpnum
      type (edges), dimension(nmax), intent(in)    :: bedges
      type (edges), dimension(nmax), intent(inout) :: wedges
      integer,                       intent(in)    :: bnodes
      integer,                       intent(inout) :: wnodes
      type (nodes), dimension(nmax), intent(in)    :: bnodex
      type (nodes), dimension(nmax), intent(inout) :: wnodex
      integer,                       intent(in)    :: nmax
      integer,                       intent(in)    :: step
      integer,                       intent(in)    :: nseplines
      integer, dimension(25,100),    intent(inout) :: lumplist1
      integer, dimension(25,100),    intent(inout) :: lumplist2
      integer, dimension(25,100),    intent(inout) :: lumplist3
      integer, dimension(25),        intent(in)    :: lumplen1
      integer, dimension(25),        intent(inout) :: lumplen2
      integer, dimension(25),        intent(inout) :: lumplen3
      integer,                       intent(in)    :: wlen1,wlen2,wlen3
      integer, dimension(25,100),    intent(inout) :: nrlist1, nrlist2
      logical, dimension(nmax),      intent(inout) :: euse, nuse, puse
      type (stack),                  intent(inout) :: estack,nstack,pstack
      type (splst), dimension(25),   intent(inout) :: seplist

      integer :: bnext1, bnext2
      integer :: bnode1,bnode2
      integer :: diedge
      integer :: dinode1,dinode2
      integer :: edg
      integer :: fredge
      integer :: frnode1, frnode2
      integer :: frpan
      real    :: kcirc1,kcirc2
      integer :: kedge
      integer :: kpan1,kpan2
      integer :: lastpan
      integer :: line
      integer :: farlen(25)
      integer :: farlist(25,100)
      integer :: n
      logical :: near
      integer :: next1, next2, nexte
      integer :: nredge
      integer :: nrnode1, nrnode2
      integer :: nrpan
      integer :: nwside
      integer :: nwsnd1, nwsnd2
      integer :: oldgam
      integer :: seplen
      integer :: wside
      integer :: wsnode1, wsnode2
!  ---------------------------------
!  input variables
!        bedges: edge-node array, lists nodes for each edge                       
!        bnodes: total number of nodes for body
!        bnodex: node coordinate array (x,y,z for each node) for body
!        body: body panels, contains strength of 
!              panel, lists of nodes & edges, "polarity",
!              normal direction and centroid
!        nseplines: number of separation lines
!        nmax: maximum number of edges 
!        seplist: list of separation edges 
!        step: current step number
!        wlen1,wlen2,wlen3: lengths of three wake zones, in panels. 
!              panels are lumped at the ends of these zones, where
!              panels double in size. at the end of the last zone, 
!              the vorticity enters the "pool" which is the final       
!              destination of all vorticity
! --------------------------------------------------------------
! output variables
!        estack,nstack,pstack: stacks of element numbers 
!                               available for the wake
!        euse,nuse,puse: lookup tables, is element in use?
!        lumplen1,2,3: number of edges at lumping sites
!        lumplist1,2,3: list of the edges which will be lumped  next
!                       lumplist (i,j) is the jth lumping edge for 
!                       the ith separation line
!        nrlist1,nrlist2: list of near edges generated during 
!                         lumping, used to reconnect newly lumped
!                         panels to old ones
!        wake: wake panels, contains strength of 
!              panel, lists of nodes & edges, "polarity",
!              normal direction and centroid
!        wedges: edge-node array, lists nodes for each edge                       
!        wednum: total number of edges for  wake
!        wnodes: total number of nodes used to describe wake
!        wnodex: node coordinate array (x,y,z for each node) for wake
!        wpnum: total number of panels for wake
! --------------------------------------------------------------
! local variables
!        bnext1,bnext2: next two body nodes of the separation edge
!        bnode1,bnode2: body nodes associated with Kutta edge
!        diedge: diagonal edge
!        dinode1,dinode2: diagonal edge nodes
!        edg: current separation edge
!        fredge: edge of far panel connecting old & new wake-side edges
!        frnode1,frnode2: nodes of fredge (above)
!        frpan: far panel
!        kcirc1,kcirc2: circulation of body nodes at Kutta edges
!        kedge: body edge associated with Kutta edge
!        kpan1,kpan2: body panels corresponding to Kutta edge
!        lastpan: far panel from the last iteration
!        line: looping index for separation lines
!        farlen: number of near edges generated during lumping
!        farlist: list of near edges generated during lumping
!        n: looping index
!        near: do we already have nredge from the last step?
!              (otherwise we already have fredge) 
!        next1,next2,nexte: elements to use for next panel's nodes, 
!                           and panel #
!        nredge: edge of near panel connecting old and new wake-side edges
!        nrnode1,nrnode2: nodes of nredge (above)
!        nrpan: near panel
!        nwside: new wake-side edge
!        nwsnd1,nwsnd2: new wake-side nodes
!        oldgam: circulation of old wake-side panel
!        seplen: current farlen
!        wside: wake-side edge from the shed panel
!        wsnode1, wsnode2: wake-side nodes from the shed panel
! --------------------------------------------------------------
! start:
! connect shed panels to separation line
!
! loop through separation lines
      do 6000 line = 1,nseplines
!  
         near = .false.
!	 if (step.eq.2) then
!	    write (*,*) 'beginning of step 2'
!	    write (*,*) 'wpnum =',wpnum
!	    write (*,*)
!	 endif
! loop through the edge list for each separation line
         do 5900 edg = 1,seplist(line)%nedges
!	    if(step.eq.2)then
!	    write (*,*) 'edge ',edg
!	    write (*,*) 'wpnum =',wpnum
!	    write (*,*)
!            endif
!
! store information from the old wake-side edge
            wside = seplist(line)%sepedge(edg)%wsedge
	    wsnode1 = seplist(line)%sepedge(edg)%wsnode(1)
	    wsnode2 = seplist(line)%sepedge(edg)%wsnode(2)
!
! create new wake-side edge and nodes
            call nulabel(nwside,wednum,euse,estack,nmax)
! for the first panel, we have to make both nodes
	    if (edg.eq.1) then
	       call nulabel(nwsnd1,wnodes,nuse,nstack,nmax)
	       call nulabel(nwsnd2,wnodes,nuse,nstack,nmax)	
! do we already have the nredge info?
            elseif (near) then
	       call nulabel(nwsnd1,wnodes,nuse,nstack,nmax)
	       nwsnd2 = next1
	    else
	       nwsnd1 = next1
	       call nulabel(nwsnd2,wnodes,nuse,nstack,nmax)
	    endif
!	    call nulabel(nwsnd1,wnodes,nuse,nstack,nmax)
!	    call nulabel(nwsnd2,wnodes,nuse,nstack,nmax)
!
! fill in edge connectivity
            wedges(nwside)%node1 = nwsnd1
	    wedges(nwside)%node2 = nwsnd2
	    wedges(nwside)%panel1 = 0
	    wedges(nwside)%panel2 = 0
!
! set new node coordinates
            bnode1 = seplist(line)%sepedge(edg)%bsnode(1)
	    bnode2 = seplist(line)%sepedge(edg)%bsnode(2)
            wnodex(nwsnd1)%x = bnodex(bnode1)%x
	    wnodex(nwsnd2)%x = bnodex(bnode2)%x
!
! replace wake-side edge
            seplist(line)%sepedge(edg)%wsedge = nwside
	    seplist(line)%sepedge(edg)%wsnode(1) = nwsnd1
	    seplist(line)%sepedge(edg)%wsnode(2) = nwsnd2
!
! now we construct two panels between the new wake-side edge and
! the old wake-side edge by connecting the ends of these two edges,
! with another edge making a diagonal. this forms the near panel,
! connected to the new wake-side edge, and the far panel, connected
! to the old wake-side edge
!
! connect ends of wake-side panels with new edges
! near edge: connects node 2 of the new wake-side edge
!            to node 2 of the old wake-side edge.
!            the near panel is formed by the new wake-side
!            edge, the near edge, and the diagonal.
!
            if (near.and.(edg.ne.1)) then
               nredge = nexte
	    else
	       call nulabel(nredge,wednum,euse,estack,nmax)
	    endif
	    nrnode1 = nwsnd2
	    nrnode2 = wsnode2 
!
! order is important!
            wedges(nredge)%node1 = nrnode1
	    wedges(nredge)%node2 = nrnode2
	    wedges(nredge)%panel1 = 0
	    wedges(nredge)%panel2 = 0
!
! far edge: connects node 1 of the new wake-side edge
!           to node 1 of the old wake-side edge.
!            the far panel is formed by the old wake-side
!            edge, the diagonal, and the far edge.
!
            if (near.or.(edg.eq.1)) then
               call nulabel(fredge,wednum,euse,estack,nmax)
	    else
	       fredge = nexte
	    endif
	    frnode1 = nwsnd1
	    frnode2 = wsnode1
!
! order is important!
            wedges(fredge)%node1 = frnode1
	    wedges(fredge)%node2 = frnode2
	    wedges(fredge)%panel1 = 0
	    wedges(fredge)%panel2 = 0   
!
! make diagonal edge
	    call nulabel(diedge,wednum,euse,estack,nmax)
! order is important!
	    dinode1 = wsnode2
	    dinode2 = nwsnd1 
!
            wedges(diedge)%node1 = dinode1
	    wedges(diedge)%node2 = dinode2
	    wedges(diedge)%panel1 = 0
	    wedges(diedge)%panel2 = 0
!
! construct near panel
	    call nulabel(nrpan,wpnum,puse,pstack,nmax)
! order is important!
	    wake(nrpan)%node(1) = nwsnd1
	    wake(nrpan)%node(2) = nwsnd2
	    wake(nrpan)%node(3) = wsnode2
	    wake(nrpan)%edge(1) = nwside
	    wake(nrpan)%edge(2) = nredge
	    wake(nrpan)%edge(3) = diedge
!
! construct far panel
	    call nulabel(frpan,wpnum,puse,pstack,nmax)
! order is important!
	    wake(frpan)%node(1) = wsnode2
	    wake(frpan)%node(2) = nwsnd1
	    wake(frpan)%node(3) = wsnode1
	    wake(frpan)%edge(1) = diedge
	    wake(frpan)%edge(2) = fredge
	    wake(frpan)%edge(3) = wside
!
! put panel data into edge structure
            wedges(nwside)%panel1 = nrpan
	    wedges(diedge)%panel1 = nrpan
	    wedges(diedge)%panel2 = frpan
	    wedges(wside)%panel2  = frpan
!	    if(step.eq.2)then
!	    write (*,*)
!	    write (*,*)
!	    write (*,*) 'connecting ',frpan,' to panel2 of ',wside
!	    write (*,*)
!	    endif
!
            if (edg.eq.1) then
	       wedges(nredge)%panel1 = nrpan
	       wedges(fredge)%panel1 = frpan
	    elseif (near) then
	       wedges(nredge)%panel2 = nrpan
	       wedges(fredge)%panel1 = frpan
	    else
	       wedges(nredge)%panel1 = nrpan
	       wedges(fredge)%panel2 = frpan
	    endif
!
! keep track of which node to use for the next panel:
! look at the body nodes corresponding to the next edge
            if (edg.ne.seplist(line)%nedges) then
               bnext1 = seplist(line)%sepedge(edg+1)%bsnode(1)
               bnext2 = seplist(line)%sepedge(edg+1)%bsnode(2)
	    else
	       bnext1 = 0
	       bnext2 = 0
	    endif
! we want to know whether to keep nredge or fredge for the next step
	    if ((bnext1.eq.bnode2).or.(bnext2.eq.bnode2)) then	       
! we want nredge (associated with new and old wake-side nodes 2)
               next1 = nwsnd2
	       next2 = wsnode2
	       nexte = nredge
	       if (bnext1.eq.bnode2) then
	          near = .false.
	       else
	          near = .true.
	       endif
	    else
! we want fredge (associated with new and old wake-side nodes 1)
               next1 = nwsnd1
               next2 = wsnode1
	       nexte = fredge
	       if (bnext1.eq.bnode1) then
	          near = .false.
	       else
	          near = .true.
	       endif	          
	    endif 
!
! set strengths of new panels:
! enforce Kutta condition by setting sum of vorticities at
! trailing edge to be zero.
! set signs of new panels to be opposite, which will give
! them the same normal direction, and the diagonal vorticities
! will cancel, effectively leaving a vortex ring
!
! get edge number
            kedge = seplist(line)%sepedge(edg)%bsedge
! get panel numbers
            kpan1 = bedges(kedge)%panel1
	    kpan2 = bedges(kedge)%panel2
	    kcirc1 = body(kpan1)%circ
	    if (kpan2.ne.0) then
	       kcirc2 = body(kpan2)%circ
	    else
	       kcirc2 = 0.0
	    endif
            wake(nrpan)%circ =  - (kcirc1 + kcirc2) 
	    wake(frpan)%circ =    (kcirc1 + kcirc2)
!	    
5900    continue
 6000 continue
!
      if (wednum.gt.nmax) then
         write (*,*) 'too many wake edges (newpan.f)'
	 write (*,*) 'wednum =',wednum
	 write (*,*) 'nmax =',nmax
	 write (*,*) 'change nmax in fish.f'
	 stop
      endif
!
! vortex lumping
!
! first level of lumping
      if ((step.ge.wlen1).and.(mod(step,2).eq.0)) then
 	 call lump (wake,wpnum,wednum,wedges,nseplines,lumplist1,  &
                 & lumplen1,farlist,farlen,nrlist1,euse,nuse,puse, &
                 & estack,nstack,pstack,nmax)
!
! when entering a new lumping zone, create new lumplist
         if  (step.eq.wlen1) then
            write (*,*)
            write (*,*) 'starting 1st level of lumping....'
	    write (*,*)
	    do 6600 line = 1,nseplines
	       seplen = farlen(line)
	       lumplist2(line,1:seplen) = farlist(line,1:seplen)                        
               lumplen2(line) = seplen
!	       write (12,*) 'seplen =',seplen
!	       write (12,*) 'lumplen = ',lumplen2(line)
!	       do 27 n=1,seplen
!	         write (12,*) 'lumplist(',n,') =',
!     &                         lumplist2(line,n)
! 27            continue
 6600       continue
	 endif
!
! second level of lumping	 
         if ((step.ge.wlen2).and.(mod(step,4).eq.0)) then
!	 
 	    call lump (wake,wpnum,wednum,wedges,nseplines,lumplist2,   &
                     & lumplen2,farlist,farlen,nrlist2,euse,nuse,puse, &
                     & estack,nstack,pstack,nmax)
!
            if (step.eq.wlen2) then
	       write (*,*) 'starting 2nd level of lumping....'
	       do 6700 line = 1,nseplines
	          seplen = farlen(line)
	          lumplist3(line,1:seplen) = farlist(line,1:seplen)
                  lumplen3(line) = seplen
 6700          continue
	    endif    
!
! the pool
!
            if (step.ge.wlen3) then
	       if (step.eq.wlen3) write (*,*) 'starting pool...'
!
	       call pool (wake,wedges,nseplines,lumplist3,lumplen3,euse, &
                        & nuse,puse,estack,nstack,pstack,step,nmax) 
	    endif
!
	 endif
!
      endif
!
!  we don't need wake normals or centroid
!      call wnorms (wnodes,wnodex,wake,wpnum,puse,nmax)

      write (12,*) 
      write (12,*) 'wake panels, node, edges = ',wpnum,wnodes,wednum
      write (12,*)

      return
      end
! --------------------------------------------------------------------
