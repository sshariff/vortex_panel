! -------------------------------------------------------------------
      subroutine presvel (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,  &
                        & puse,cp,wind,cl,cd,bldelt,seplist,nseplines, &
                        & oldcirc,body,wake,vref,step,nmax)
!
!  This subroutine computes the pressure and velocity at the
!  collocation points using e, the edge/panel
!  influence matrix, which is calculated in subroutine aij
! -----------------------------------------
      use panel_type; use interfish2; use interfish3; use interfish4
!
      implicit none
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

      real                    :: area
      real, dimension(3)      :: coll
      integer                 :: pan,j,n,side
      real                    :: dgdt
      real, dimension(3)      :: df
      real                    :: drag
      real                    :: lift
      real                    :: gamma
      real                    :: magx
      real, dimension(3)      :: q,qtot,qn,qt
      real                    :: qsq,vsq
      real, dimension(3)      :: r1,r2
      real, dimension(3)      :: r1xr2
      real                    :: totalarea
      real, dimension(nmax,3) :: wp
      real, dimension(3)      :: x1,x2
      
      real :: dumq
!      real :: lastdumq

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
!        oldbody: body panels from last time step. this is needed
!                 for the calculation of some time derivatives
!        puse: lookup tables: is element in use?
!        step: current step number
!        vref: reference velocity at the collocation point of each 
!              panel. vref is the kinematic velocity of each panel,
!              consisting of the freestream plus the velocity due to
!              translation and rotation
! -----------------------------------------
!  output variables
!        bldelt: boundary layer displacement thickness of panels
!        cl, cd: lift and drag coefficients
!        cp: pressure coefficient for each panel
!        nseplines: number of separation lines
!        seplist: list of separation edges
!        wind: induced velocity at the collocation points
! -----------------------------------------
!  local variables
!        area: body panel area
!        coll : coordinates of collocation point
!        df: differential vector force
!        dgdt: time derivative of panel's circulation
!        drag: drag
!        gamma: strength of current panel
!        lift: lift
!        magx: magnitude of r1xr2 (cross product of two panel edges)
!        pan,j,n,side: stepping indices for panels, nodes, and edges
!        q: velocity induced by a single edge
!        qn: normal velocity
!        qt: tangential velocity
!        qtot: total velocity
!        qsq, vsq: square of total velocity and reference velocity
!        r1,r2: components of two of the panel edges
!        r1xr2: cross product of two panel edges
!        totalarea: total body area
!        wp: matrix of total induced velocities at each body panel
!        x1, x2: coordinates of nodes at each end of the segment
!                whose influence is being computed
! -------------------------------------------
! start:
!
      lift = 0.0
      drag = 0.0
      totalarea = 0.0
!  calculate induced velocity ww for each body panel
!      write (12,*)
      write (12,*) 'bpnum =',bpnum
      write (12,*)
      write (12,*) 'vref = ',vref(1,1:3)
      write (12,*)
      
!      write (12,*) 'values in body...'
!      do 1267 side = 1,bpnum
!        write (12,'(i3,e22.12,a,3e13.3)') side,body(side)%circ, &
!	     & '      norm =',body(side)%normal
! 1267 continue
!      write (12,*)
      
      do 4900 pan = 1,bpnum
!
!  reset induced velocity
         wind(pan,:) = 0.0
	 coll = body(pan)%centr
         dumq = 0.0
!
!	 write (12,*)
!         write (12,*)
!	 write (12,*) '---------------------------------------------------'
!         write (12,*) 'looking at effect on body panel ',pan
!	 write (12,*) '       coll = ',coll
!	 write (12,'(a,3g12.4)') '                         n = ',body(pan)%normal
!	 write (12,*)
! 
!  compute contribution of body panels
         do 2800 j = 1,bpnum
!  panel strength
	    gamma = body(j)%circ
!	 write (12,*)
!         write (12,*) '     looking at effect from !BODY! panel ',j
!	 write (12,*) '             which has strength ',gamma
!	    lastdumq = dumq
!  loop through panel edges
	    do 2700 side = 1,3
	       x1 = bnodex(bedges(body(j)%edge(side))%node1)%x
	       x2 = bnodex(bedges(body(j)%edge(side))%node2)%x
!	 write (12,*) '            from side ',side
!         write (12,'(a,3g12.4)') '                   x1 = ',x1
!	 write (12,'(a,3g12.4)') '                   x2 = ',x2
               call lnvortx (coll,x1,x2,q)
!	 write (12,*)
!	 write (12,'(a,3g12.4)') '                         q = ',gamma*q
!	 write (12,*)
!	 write (12,'(a,g12.4)')  '                         q dot n = ', &
!	             & gamma*dot_product(q,body(pan)%normal)
!	 write (12,*)
	       wind(pan,:) = wind(pan,:) + gamma*q
	       dumq = dumq + gamma*dot_product(q,body(pan)%normal)
 2700       continue
!	    write (12,9978) 'At panel ',pan,', effect from body ', &
!	                  & j,'  (',gamma,'), = ',dumq - lastdumq
 2800    continue	 
!
!  compute contribution of wake panels 
         do 3800 j = 1,wpnum
!
	    if (puse(j)) then
!  panel strength
	       gamma = wake(j)%circ
!	 write (12,*)
!         write (12,*) '     looking at effect from !WAKE! panel ',j
!	 write (12,*) '             which has strength ',gamma
!	       lastdumq = dumq
!  loop through panel edges
	       do 3700 side = 1,3
		  x1 = wnodex(wedges(wake(j)%edge(side))%node1)%x
		  x2 = wnodex(wedges(wake(j)%edge(side))%node2)%x
!	 write (12,*) '            from side ',side
!         write (12,'(a,3g12.4)') '                   x1 = ',x1
!	 write (12,'(a,3g12.4)') '                   x2 = ',x2
        	  call lnvortx (coll,x1,x2,q)
!	 write (12,*)
!	 write (12,'(a,3g12.4)') '                         q = ',gamma*q
!	 write (12,*)
!	 write (12,'(a,g12.4)')  '                         q dot n = ', &
!	             & gamma*dot_product(q,body(pan)%normal)
!	 write (12,*)
		  wind(pan,:) = wind(pan,:) + gamma*q
                  dumq = dumq + gamma*dot_product(q,body(pan)%normal)
 3700          continue
!
!	       write (12,9978) 'At panel ',pan,', effect from wake ', &
!	                   & j,'  (',gamma,'), = ',dumq - lastdumq
 9978          format(a,i4,a,i4,a,f9.3,a,f9.3)
            endif
!
 3800    continue	 
!
!  compute total velocity Q
!         qtot = wind(pan,:) + vref(pan,:)
         qtot = 2.0*(wind(pan,:) + vref(pan,:))
	 qsq = dot_product(qtot,qtot)
	 vsq = dot_product(vref(pan,:),vref(pan,:))
!
         qn = body(pan)%normal*(dot_product(qtot,body(pan)%normal))
!	 magx  = sqrt(dot_product(qn,qn))
!	 qn = qn/magx
	 qt = qtot - qn
!        write (12,'(a,i4,a,3f9.2,a,3f9.2,a)') &
!	 & 'panel = ',pan,', qn = (',qn,'),   qt = (',qt,')'
!	 write (12,'(a,3f9.2,a)') 'normalized    qn = (',qn/magx,')'
!	 write (12,*) '     magx = ',magx
!
!
!
!         qsq = 0.0
!	 vsq = 0.0
!         do 3905 n = 1,3
!	    qn = wind(n,:) + vref(pan,n)
!	    qsq = qsq + qn*qn
!	    vsq = vsq + vref(pan,n)*vref(pan,n)
! 3905    continue
!         
!         dgdt = ( body(pan)%circ - oldbody(pan)%circ ) / dt
         dgdt = ( body(pan)%circ - oldcirc(pan) ) / dt
!
! pressure computation:
!   cp = (p - pref)/(1/2 rho vref^2) 
!      = 1 - Q^2/vref^2 - 2/vref^2 (dphi/dt)  (Eq 13.28, p 430, Katz)
!
!   where Q is the total velocity, and vref is the kinematic velocity
!   (the velocity due to translation and rotation of the body panels)
!
!   For the vortex ring model, dphi/dt = 1/2 dgamma/dt
!   (see Eq 13.149, p 489, Katz), so
!   cp = 1 - (Q/vref)^2 - (dgamma/dt)/vref^2
!
!

!
! don't use the dgdt for the first step
!
          if (step.eq.1) then
	     cp(pan) = 1.0 - qsq/vsq
          else
	     cp(pan) = 1.0 - qsq/vsq - dgdt/vsq
          end if
	  
!          write (12,*)
!          write (12,*) 'panel ',pan
!	  write (12,*) 'qsq/vsq = ',qsq/vsq
!	  write (12,*) 'dgdt/vsq = ',dgdt/vsq
!
!
!          cp(pan) = dot_product(wind(pan,:),wind(pan,:))/vsq
!          cp(pan) = qsq/vsq
!
!          write (12,*)
!          write (12,*) 'panel ',pan
!	  write (12,*) '                  normal = ',body(pan)%normal
!	  write (12,*) 'qtot.qtot = ',qsq,qtot/sqrt(qsq)
!	  write (12,*) '|qtot|,qtot = ',sqrt(qsq),'(',qtot,')'
!	  write (12,*) '|vref|,vref = ',sqrt(vsq),vref(pan,:)
!	  write (12,*) '|wind|,wind = ', &
!	              & sqrt(dot_product(wind(pan,:),wind(pan,:))), &
!		      & '(',wind(pan,:),')'
!                      & /sqrt(qsq)
!		      & sqrt(dot_product(wind(pan,:),wind(pan,:)))
!          write (12,*) '   qsq/vsq  = ',cp(pan)
!          write (12,*) 'pan, qsq/vsq = ',pan,cp(pan)
!	  write (12,*)

!
!  compute panel area for force calculation
	 r1 = bnodex(body(pan)%node(2))%x - bnodex(body(pan)%node(1))%x
	 r2 = bnodex(body(pan)%node(3))%x - bnodex(body(pan)%node(1))%x
         call cross (r1, r2, r1xr2, magx)
	 area = magx/2.0
	 totalarea = totalarea + area
	 df = -cp(pan)*area*body(pan)%normal
!	 write (*,'(a,i4,a,3e16.6)') '  df (',pan,') =',df(1),df(2),df(3)
!	 if (body(pan)%normal(2).lt.0) df = -df
	 lift = lift + df(2)
	 drag = drag + df(1)
!	 if (df(2).lt.0) then
!	    write (*,*) 'df =',df
!	    write (*,*) 'pan =',pan
!	    write (*,*) 'cp = ',cp(pan)
!	    write (*,*) 'normal =',body(pan)%normal
!	 endif
!
 4900 continue
         cl = lift/totalarea
	 cd = drag/totalarea
!
! compute new boundary layer characteristics (bl displacement, 
! separation locations, etc)
      call blchar(bpnum,body,cp,bldelt,seplist,nseplines,nmax)
!
      return
      end
! --------------------------------------------------------------------

