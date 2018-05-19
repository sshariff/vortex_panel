! --------------------------------------------------------------------
      subroutine aij1 (bpnum,wpnum,body,a,dt,bedges,wedges,bnodex,wnodex, & 
                     & wake,oldbndx,right,vref,qinf,nseplines,            &
		     & seplist,nmax)
!
!  This subroutine computes the matrix of influence coefficients, a,
!  for the first timestep. This includes the influence of the first
!  set of wake panels. All of the inital wake panels are assumed
!  to have the same strength, and a version of the Kelvin condition 
!  (total circulation = 0) is used to calculate the strength.
!  The body is moving forward in the x-direction initially,
!  the y-direction is toward the top of the body, and the 
!  z-direction toward the side. The Kutta panels are assumed to 
!  lie in the xy-plane, where the Kutta edge was at the beginning 
!  of the timestep.
!
!  The Kutta condition is enforced by requiring that the vorticity 
!  at the trailing edge is zero (see code below).
!
!  The effect of each panels is computed using Biot-Savart,
!  i.e. by summing the effect of each edge. Each edge is associated
!  with 2 panel strengths, but we want to avoid doing the (dl)x(r)
!  calculation twice, so we construct the matrix of edge/panel
!  influence e (which contains (dl)x(r)), and use this to compute
!  a. e is passed back to the main program, for use later for 
!  velocity and pressure computations.
! -----------------------------------------
      use panel_type; use interfish3
!
      implicit none
!--------------------      
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
      integer  ::   bfile, pfile, efile, ofile, keybd, scrn
      common/files/ bfile, pfile, efile, ofile, keybd, scrn
!
      real               :: aij
      real, dimension(3) :: coll
      integer            :: edg,edgnum,bed,wed
      integer            :: i,j,n,side,line,wpan,bpan
      integer            :: newpans
      integer            :: p1,p2,wp
      real, dimension(3) :: q
      integer            :: row
      real, dimension(3) :: x1,x2
      
      character*8           :: dayt
      character*10          :: tim
      character*5           :: zone
      integer, dimension(8) :: values
      real                  :: sec,oldsec
! ------------------------------------------
! common variables
!        bfile: unit number of body motion data file
!        ofile: unit number of output data file
!        keybd: unit number for keyboard input
!        scrn:  unit number for screen output
! ----------------------------------------
! input variables
!        bedges, wedges: edge-node arrays, lists nodes that 
!                        correspond to each edge
!        bnodex, wnodex: node coordinate array (x,y,z for each node)
!        bpnum, wpnum: total number of panels for body, wake
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        dt: size of current timestep
!        nmax: maximum number of any element (for memory allocation)
!        nseplines: number of separation lines
!        oldbndx: node coords from last timestep, used for rhs
!        qinf: magnitude of freestream velocity
!        seplist: list of separation edges
!  -----------------------------------------
!  output variables
!        a: matrix of influence coefficients. a(i,j) is the effect
!           at collocation point of panel i from panel j
!        right: rhs for matrix inversion
!        vref: reference velocity at the collocation point of 
!              each panel. vref is the kinematic velocity of
!              the panel, consisting of the freestream plus the
!              velocity due to translation and rotation      
! -----------------------------------------
!  local variables
!        aij: local variable for calculation of a(i,j)
!        coll : coordinates of collocation point
!        edg,edgnum,bed,wed: indices for edges
!        i,j,n,side,line,wpan,bpan: indices for stepping
!        newpans: number of initial wake panels
!        p1,p2,wp: panel numbers associated with current edge
!        q: induced velocity
!        row: row number, used for wake panels
!        x1, x2: coordinates of nodes at each end of the segment
!                whose influence is being computed
! -------------------------------------------------------------
!  start:
!  
!  determine the number of wake panels (2 panels created per edge)
      newpans = 0
      do 900 line = 1,nseplines
         newpans = newpans + 2*seplist(line)%nedges
 900  continue
      if (newpans.ne.wpnum) then
         write (*,*) 'Discrepancy in number of wake panels'
	 write (*,*) 'newpans =',newpans
	 write (*,*) 'wpnum =',wpnum
!	 stop
      endif
!
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
! 
!  compute a:
!     part of a corresponding to body panels
      do 4000 i = 1,bpnum
!	 write (12,*)
!         write (12,*)
!	 write (12,*) '---------------------------------------------------'
!         write (12,*) 'looking at effect on body panel ',i
	 coll = body(i)%centr
!	 write (12,*) '       coll = ',coll
!	 write (12,'(a,3g12.4)') '                         n = ',body(i)%normal
!	 write (12,*)
!
         do 3000 j = 1,bpnum
!	 write (12,*)
!         write (12,*) '     looking at effect from body panel ',j
            aij = 0.
	    do 2900 side = 1,3
               x1 = bnodex(bedges(body(j)%edge(side))%node1)%x
               x2 = bnodex(bedges(body(j)%edge(side))%node2)%x
!	 write (12,*) '            from side ',side
!         write (12,'(a,3g12.4)') '                   x1 = ',x1
!	 write (12,'(a,3g12.4)') '                   x2 = ',x2
               call lnvortx (coll,x1,x2,q)
	       aij = aij + dot_product(q,body(i)%normal)
!	 write (12,*)
!	 write (12,'(a,3g12.4)') '                         q = ',q
!	 write (12,*)
!	 write (12,'(a,g12.4)')  '                         q dot n = ', &
!	             & dot_product(q,body(i)%normal)
!	 write (12,*)
 2900       continue
 
!	    write (12,*) '      a = ',aij
!	    write (12,*)
            a(i,j) = aij
!
 3000    continue
!	 write (12,*) '---------------------------------------------------'
 4000 continue
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
      write (*,'(a,f10.3,a)') '     computing ''a'' from body panels took ', &
			      & sec-oldsec,' seconds'
!
!  calculate effect of new wake panels on each body panel 
      do 5600 bpan = 1,bpnum
!	 write (12,*)
!         write (12,*)
!         write (12,*) 'looking at effect on body panel ',i
!	 write (12,*)
	 coll = body(bpan)%centr
!
!  loop through wake panels
         do 5500 wpan = 1,wpnum
!	    write (12,*)
!            write (12,*) 'looking at effect from wake panel ',wpan
!  compute normal velocity induced at body panel by the 
!  edges of wake panel 'wpnum'
            aij = 0.
!  loop through edges
            do 5405 side = 1,3
               x1 = wnodex(wedges(wake(wpan)%edge(side))%node1)%x
               x2 = wnodex(wedges(wake(wpan)%edge(side))%node2)%x
               call lnvortx (coll,x1,x2,q)
!	 write (12,*) '      from side ',side
!	 write (12,'(a,3f9.3)') '                       q = ',q
!	 write (12,'(a,f9.3)') '       q dot n = ', &
!	             & dot_product(q,body(i)%normal)
	       aij = aij + dot_product(q,body(bpan)%normal)
!
 5405       continue
            a(bpan,bpnum+wpan) = aij
!
 5500    continue
 5600 continue
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
      write (*,'(a,f10.3,a)') '     computing ''a'' from wake panels took ', &
			      & sec-oldsec,' seconds'
!
!
!  before there are any wake panels, right=vref
      right = 0.0
!
!  calculate velocity of panels 
      call panvel (bpnum,body,bnodex,oldbndx,right,dt,vref,qinf,nmax)       
!
!  add entries for Kutta panels:
!
!  initialize new rows of a(i,j) for wake
      do 888 i = bpnum+1,bpnum+wpnum
	 a(i,1:wpnum) = 0.
 888  continue
!
      row = bpnum + 1
      do 1888 line = 1,nseplines
!  loop through separation edges
         do 1887 edg = 1,seplist(line)%nedges
!	    write(ofile,*)'row =',row
!  find panels associated with this edge
            bed = seplist(line)%sepedge(edg)%bsedge
	    wed = seplist(line)%sepedge(edg)%wsedge
            p1 = bedges(bed)%panel1
	    p2 = bedges(bed)%panel2
	    wp = wedges(wed)%panel1
!  enforce Kutta condition by setting the vorticity at the trailing
!  edge to zero. the vorticity is the sum of the contributions from
!  the two body panels (or one body panel at a single panel Kutta edge)
!  and the new wake panel, and circulation is proportional to 
!  vorticity for these constant-strength vortex lines
	    a(row,p1)  = 1.0
	    if (p2.ne.0) a(row,p2) = 1.0
	    a(row,row) = 1.0
! set signs of new panels to be opposite, which will give
! them the same normal direction, and the diagonal vorticities
! will cancel, leaving horseshoe vortices only
            a(row+1,row)   =  1.0
	    a(row+1,row+1) =  1.0
	    row = row + 2
 1887    continue
 1888 continue
!
!      write (ofile,*) 'a''s...'
!      write (ofile,'(3x,22i11)') ((j),j=1,bpnum+wpnum)
!      do 667 i = 1,bpnum+wpnum
!         write (ofile,'(i3,22f11.3)') i,((a(i,j)), j=1,bpnum+wpnum)
! 667  continue
!      write (ofile,*)
!      write (ofile,*) 'rhs...'
!      do 669 i = 1,bpnum+wpnum
!         write (ofile,'(i3,f16.6)') i,right(i)
! 669  continue
!       
      return
      end
! --------------------------------------------------------------------

