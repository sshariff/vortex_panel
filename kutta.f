! -------------------------------------------------------------------
      subroutine kutta (bnodex,wnodes,wnodex,wake,wpnum,wednum, &
                      & wedges,nseplines,oldbndx,lumplist,lumplen,nrlist1,  &
                      & nrlist2,euse,nuse,puse,estack,nstack,pstack,        &
                      & seplist,nmax)                
!
!  This subroutine creates the wake panels at the Kutta edges
!  for the first timestep. Also calculates the normal direction 
!  for the wake panels.
!  
!  This follows the same algorithm as the panel creation in
!  subroutine 'newpan.f', but now we must also create the
!  'old wake-side edge'. This is placed behind the Kutta edge
!  at a distance which corresponds to the distance the body
!  moves in the first time-step, which is calculated using
!  the value qinf, which must be estimated by bmotion.
!  --------------------------------------------------------------
      use panel_type; use interfish2; use interfish3
!
      implicit none
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

      integer :: bnext1, bnext2
      integer :: bnode1,bnode2
      integer :: diedge
      integer :: dinode1,dinode2
      integer :: edg
      integer :: fredge
      integer :: frnode1, frnode2
      integer :: frpan
      integer :: lastpan
      integer :: line
      logical :: near
      integer :: next1, next2, nexte
      integer :: nredge
      integer :: nrnode1, nrnode2
      integer :: nrpan
      integer :: nwside
      integer :: nwsnd1, nwsnd2
      integer :: wside
      integer :: wsnode1, wsnode2
!  ---------------------------------
!  input variables
!        bnodex: body node coordinate array (x,y,z for each node)
!        nmax: maximum number of edges
!        nseplines: number of separation lines
!        oldbndx: node coords from last timestep, used for rhs
!                 to determine the velocity of each panel
!        seplist: list of separation edges 
! --------------------------------------------------------------
! local variables
!        bnext1,bnext2: the next two body nodes of the separation edge
!        bnode1,bnode2: body nodes corresponding to Kutta edge
!        diedge: diagonal edge
!        dinode1,dinode2: diagonal edge nodes
!        edg: current separation edge
!        edge of far panel connecting old and new wake-side edges
!        frnode1, frnode2: nodes of fredge 
!        frpan: far panel
!        lastpan: far panel from the last iteration
!        line: looping index for separation lines
!        near: do we already have nredge from the last step?
!              (otherwise we already have fredge) 
!        next1, next2, nexte: elements to use for next panel's nodes,
!                             and panel #
!        nredge: edge of near panel connecting old and new 
!                wake-side edges
!        nrnode1, nrnode2: nodes of nredge
!        nrpan: near panel
!        nwside: new wake-side edge
!        nwsnd1, nwsnd2: new wake-side nodes
!        wside: wake-side edge from the shed panel
!        wsnode1, wsnode2: wake-side nodes from the shed panel
! --------------------------------------------------------------
! output variables
!        estack,nstack,pstack: stacks of element numbers 
!                               available for the wake
!        euse,nuse,puse: lookup tables, is element in use?
!        nrlist1,nrlist2: list of near edges generated during 
!                         lumping, used to reconnect newly lumped 
!                         panels to old ones
!        lumplen: number of edges at lumping site
!        lumplist: list of the edges which will be lumped next
!                  lumplist (i,j) is the jth lumping edge for 
!                  the ith separation line
!        wedges: wake edge-node array, lists nodes that correspond
!                to each edge
!        wednum: total number of wake edges
!        wnodes: total number of nodes used to describe wake
!        wake: wake panels, contains strength of panel
!              and lists of nodes & edges, normal direction, etc.
!        wnodex: node coordinate array (x,y,z for each node)
!        wpnum: total number of panels used to describe wake
! --------------------------------------------------------------
! start:
!
      wnodes = 0
      wpnum  = 0
      wednum = 0
! connect shed panels to separation line
!
! loop through separation lines
      do 6000 line = 1,nseplines
!
         near = .false.
!
! initialize some data for lumplists
         lumplen(line) = seplist(line)%nedges
         nrlist1(line,1:100) = 0
	 nrlist2(line,1:100) = 0
!
! loop through the edge list for each separation line
         do 5900 edg = 1,seplist(line)%nedges
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
!
! fill in edge connectivity
            wedges(nwside)%node1 = nwsnd1
	    wedges(nwside)%node2 = nwsnd2
	    wedges(nwside)%panel1 = 0
	    wedges(nwside)%panel2 = 0
!
! set new node coordinates
!            write(*,*)
!	    write(*,*) 'separation edge = ',
!     &                  seplist(line)%sepedge(edg)%bsedge
            bnode1 = seplist(line)%sepedge(edg)%bsnode(1)
	    bnode2 = seplist(line)%sepedge(edg)%bsnode(2)
! 	    write(12,*) 'sepnode1 = ',bnode1
!	    write(12,*) 'sepnode2 = ',bnode2
            wnodex(nwsnd1)%x = bnodex(bnode1)%x
	    wnodex(nwsnd2)%x = bnodex(bnode2)%x
!
!
! store this information as the wake-side edge
            seplist(line)%sepedge(edg)%wsedge = nwside
	    seplist(line)%sepedge(edg)%wsnode(1) = nwsnd1
	    seplist(line)%sepedge(edg)%wsnode(2) = nwsnd2
!
! create far wake points
            call nulabel(wside,wednum,euse,estack,nmax)
! for the first panel, we have to make both nodes            
            if (edg.eq.1) then
	       call nulabel(wsnode1,wnodes,nuse,nstack,nmax)
	       call nulabel(wsnode2,wnodes,nuse,nstack,nmax)	    
! do we already have the nredge info?
            elseif (near) then
	       call nulabel(wsnode1,wnodes,nuse,nstack,nmax)
	       wsnode2 = next2
	    else
	       wsnode1 = next2
	       call nulabel(wsnode2,wnodes,nuse,nstack,nmax)
	    endif
!
! store this edge in the lump list
            lumplist(line,edg) = wside
!
! fill in edge connectivity
            wedges(wside)%node1 = wsnode1
	    wedges(wside)%node2 = wsnode2
	    wedges(wside)%panel1 = 0
	    wedges(wside)%panel2 = 0
!
! set new node coordinates for kutta edge from old coordinates 
            wnodex(wsnode1)%x(1) = oldbndx(bnode1)%x(1)           
            wnodex(wsnode1)%x(2) =  bnodex(bnode1)%x(2)
	    wnodex(wsnode1)%x(3) =  bnodex(bnode1)%x(3)
!	    
            wnodex(wsnode2)%x(1) = oldbndx(bnode2)%x(1)           
            wnodex(wsnode2)%x(2) =  bnodex(bnode2)%x(2)
	    wnodex(wsnode2)%x(3) =  bnodex(bnode2)%x(3)
!
! now we construct two panels between the new wake-side edge and
! the old wake-side edge by connecting the ends of these two edges,
! with another edge making a diagonal. this forms the near panel,
! connected to the new wake-side edge, and the far panel, connected
! to the old wake-side edge.
!
! we must pay attention to the order of the nodes in the edges data
! structure, so that the direction of
! the vorticity of each edge is consistent (see comments below).
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
!               write (*,*) 'bnext1,bnext2 =',bnext1,bnext2
!	       write (*,*) 'bnode1,bnode2 =',bnode1,bnode2
	       
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
 5900    continue
 6000 continue
!
! we don't need this wake information
!      call wnorms (wnodes,wnodex,wake,wpnum,puse,nmax)
!
      return
      end
! --------------------------------------------------------------------

