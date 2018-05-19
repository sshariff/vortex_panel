! -------------------------------------------------------------------------
      subroutine rhs (bpnum,wpnum,body,wake,bnodex,oldbndx,wnodex,wedges, &
                    & right,dt,vref,qinf,puse,nmax)
!
!  This subroutine computes the RHS for the matrix inversion. The
!  RHS is the normal velocity of collocation point, which
!  is the superposition of the panel velocity (which includes the
!  freestream) and the influence of the wake panels.
!  
!  RHS_i = - (U + u_w, V + v_w, W + w_w)_i dot n_i
!          (based on Eq. 13.117, p. 472, Katz)
!
! ----------------------------------------------------
      use panel_type; use interfish3
! 
      implicit none 
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
      
      integer            :: bpan,wpan,side
      real, dimension(3) :: coll
      real               :: gamma
      real, dimension(3) :: q
      real, dimension(3) :: x1,x2
! ---------------------------------------
!   input variables
!        bnodex: node coordinate array (x,y,z for each node)
!        bpnum, wpnum: total number of panels for body, wake
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        dt: size of current timestep
!        nmax: maximum number of any element (for memory allocation)
!        oldbndx: node coords from last timestep, used for rhs
!        puse: lookup tables: is element in use?
!        qinf: magnitude of freestream velocity
! -----------------------------------------
!   output variable
!        right: rhs for matrix inversion
!        vref: reference velocity at the collocation point of 
!              each panel. vref is the kinematic velocity of
!              the panel, consisting of the freestream plus the
!              velocity due to translation and rotation      
! -----------------------------------------
!   local variables
!        bpan,wpan,side: stepping indices for elements
!        coll : coordinates of collocation point
!        gamma: strength of current wake panel
!        q: induced velocity
!        x1, x2: coordinates of nodes at each end of the segment
!                whose influence is being computed
! ------------------------------------------
!  start:
!
      right = 0.0
!
!  calculate velocity of panels
!      call panvel (bpnum,body,bnodex,oldbndx,right,dt,vref,qinf,nmax)
      call panvel (bpnum,body,bnodex,oldbndx,right,dt,vref,qinf,nmax)
!
!  calculate effect of wake panels on each body panel 
!  (wake-induced velocity)
      do 5600 bpan = 1,bpnum
	 coll = body(bpan)%centr
!
!  loop through wake panels
         do 5500 wpan = 1,wpnum
!	 
	    if (puse(wpan)) then 
!
!               right(bpan) = 0.
!  determine wake panel strength for velocity calculation
	       gamma = wake(wpan)%circ
!  compute normal velocity induced at body panel by the 
!  edges of wake panel 'wpan'
               do 5405 side = 1,3
        	  x1 = wnodex(wedges(wake(wpan)%edge(side))%node1)%x
        	  x2 = wnodex(wedges(wake(wpan)%edge(side))%node2)%x
        	  call lnvortx (coll,x1,x2,q)
		  right(bpan) = right(bpan) &
		              & - gamma*dot_product(q,body(bpan)%normal)
 5405          continue
            endif
!
 5500    continue
 5600 continue
!
      return
      end
! --------------------------------------------------------------------

