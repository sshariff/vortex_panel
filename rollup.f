! -------------------------------------------------------------------
      subroutine rollup (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges, &
                       & wnodes,nuse,puse,body,wake,nmax)
!
!  This subroutine computes the rollup of the wake vortex panels
!  and vortex points. 
! -----------------------------------------
      use panel_type; use interfish2; use interfish3
!
      implicit none
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

      real, dimension(3)           :: coll,x1,x2
      real, dimension(3)           :: dq
      integer                      :: i,j,n,side
      real                         :: gamma
      real, dimension(nmax,3)      :: q
! ----------------------------------------
!  input variables
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        bedges, wedges: edge-node arrays, lists nodes that 
!                        correspond to each edge
!        bnodex,wnodex: node coordinate arrays (x,y,z for each node)
!                       for body and wake
!        bpnum,  wpnum: total number of panels for body, wake
!        dt: size of current timestep
!        nmax: maximum number of edges 
!        nuse,puse: lookup tables: is element in use?
!        wnodes: number of nodes used to describe wake
! -----------------------------------------
!  output variables
!        wnodex: node coordinate arrays (x,y,z) for wake
! -----------------------------------------
!  local variables
!        coll,x1,x2: coordinates of collocation point and nodes
!        dq: induced velocity at from an edge
!        gamma: strength of current panel
!        i,j,n,side: stepping indices for panels, nodes, & edges
!        q: matrix of total induced velocities at each wake node
! ------------------------------------------
!  start:
!
!  calculate induced velocity ww for each wake node
      do 4900 i = 1,wnodes
!
!        check if this wake node is in current use
!         -----------------------------------------------
!        (start of IF structure 1)
!
         if (nuse(i)) then
!
	    coll = wnodex(i)%x
!  reset induced velocity
            q(i,1:3) = 0.0
! 
!  compute contribution of body panels
            do 2800 j = 1,bpnum
!  panel strength
	       gamma = body(j)%circ
!  loop through panel edges
	       do 2700 side = 1,3
	          x1 = bnodex(bedges(body(j)%edge(side))%node1)%x
	          x2 = bnodex(bedges(body(j)%edge(side))%node2)%x
                  call lnvortx (coll,x1,x2,dq)
	          q(i,1:3) = q(i,1:3) + gamma*dq
 2700          continue
 2800       continue	 
!
!  compute contribution of wake panels 
            do 3800 j = 1,wpnum
!
!              check if this wake panel is in use
!              -------------------------
!              (start of IF structure 2)
!          
               if (puse(j)) then
!  panel strength
	          gamma = wake(j)%circ
!  loop through panel edges
	          do 3700 side = 1,3
	             x1 = wnodex(wedges(wake(j)%edge(side))%node1)%x
	             x2 = wnodex(wedges(wake(j)%edge(side))%node2)%x
                     call lnvortx (coll,x1,x2,dq)
	             q(i,1:3) = q(i,1:3) + gamma*dq
 3700             continue
!
               endif
!             (end of IF structure 2)
!             -------------------------
!
 3800       continue	 
!
         endif
!         -------------------------------------------------
!        (end of IF structure 1)
!
 4900 continue
!
!  convect node according to local velocity
!      write (*,*) 'dt =',dt
      do 7788 i = 1,wnodes
         if (nuse(i)) wnodex(i)%x = wnodex(i)%x + q(i,1:3)*dt
!	 write (*,'(a,i4,3e16.6)') 'i,q(i) = ',i,q(i,1:3)
 7788 continue
!
!  allow wake vorticity to diffuse
      call diffuse(wake,wpnum,puse,nmax)
!
      return
      end
! --------------------------------------------------------------------

