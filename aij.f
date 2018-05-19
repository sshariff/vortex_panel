! --------------------------------------------------------------------------
      subroutine aij (bpnum,wpnum,body,a,bedges,wedges,bnodex,wnodex,wake, &
                    & oldbndx,right,dt,vref,qinf,puse,nmax)
!
!  This subroutine computes the matrix of influence coefficients, a.
!  The effect of each panels is computed using Biot-Savart,
!  i.e. by summing the effect of each edge. Each edge is associated
!  with 2 panel strengths, but we want to avoid doing the (dl)x(r)
!  calculation twice, so we construct the matrix of edge/panel
!  influence e (which contains (dl)x(r)), and use this to compute
!  a. e is passed back to the main program, for use later for 
!  velocity and pressure computations.
! -----------------------------------------
      use panel_type; use interfish2; use interfish3
! -------------------------------------------------------
      implicit none
!
      integer,                       intent(in)  :: bpnum, wpnum
      type (panel), dimension(nmax), intent(in)  :: body, wake
      real,         dimension(:,:),  intent(out) :: a
      real,                          intent(in)  :: dt
      type (edges), dimension(nmax), intent(in)  :: bedges, wedges
      type (nodes), dimension(nmax), intent(in)  :: bnodex, wnodex
      type (nodes), dimension(nmax), intent(in)  :: oldbndx
      real,         dimension(nmax), intent(out) :: right
      real,       dimension(nmax,3), intent(out) :: vref
      real,                          intent(in)  :: qinf
      logical,      dimension(nmax), intent(in)  :: puse
      integer,                       intent(in)  :: nmax


      real, dimension(3) :: coll
      integer            :: edg,edgnum,bed,wed
      real               :: aa
      integer            :: i,j,n,side
      integer            :: newpans
      integer            :: p1,p2
      real, dimension(3) :: q
      real, dimension(3) :: x1,x2
      
      character*8           :: dayt
      character*10          :: tim
      character*5           :: zone
      integer, dimension(8) :: values
      real                  :: sec,oldsec
! ----------------------------------------
!  input variables
!        bedges, wedges: edge-node arrays, lists nodes that 
!                        correspond to each edge
!        bnodex, wnodex: node coordinate array (x,y,z for each node)
!        bpnum, wpnum: total number of panels for body, wake
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        dt: size of current timestep
!        puse: lookup tables: is element in use?
!        nmax: maximum number of any element (for memory allocation)
!        nseplines: number of separation lines
!        oldbndx: node coords from last timestep, used for rhs
!        qinf: magnitude of freestream velocity
! -----------------------------------------
!  output variables
!        a: matrix of influence coefficients. a(i,j) is the effect
!           at collocation point of panel i from panel j
!        right: rhs for matrix inversion
!        vref: reference velocity at the collocation point of 
!              each panel. vref is the kinematic velocity of
!              the panel, consisting of the freestream plus the
!              velocity due to translation and rotation      
! ------------------------------------------
!  local variables
!        aa: local variable for calculation of a(i,j)
!        coll : coordinates of collocation point
!        edg,edgnum,bed,wed: indices for edges
!        i,j,n,side,line,wpan,bpan: indices for stepping
!        newpans: number of initial wake panels
!        p1,p2,wp: panel numbers associated with current edge
!        q: induced velocity
!        row: row number, used for wake panels
!        x1, x2: coordinates of nodes at each end of the segment
!                whose influence is being computed
!  -----------------------------------------
!  start:
!  
      go to 5555

      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
! 
!  compute a
      do 4000 i = 1,bpnum
!        effect on panel i
	 coll = body(i)%centr
!   
         do 3000 j = 1,bpnum
!           from panel j
            aa = 0.
	    do 2900 side = 1,3 
               x1 = bnodex(bedges(body(j)%edge(side))%node1)%x
               x2 = bnodex(bedges(body(j)%edge(side))%node2)%x
               call lnvortx (coll,x1,x2,q)
!  use only normal contribution
	       aa = aa + dot_product(q,body(i)%normal)            
 2900       continue
            a(i,j) = aa
 3000    continue
 4000 continue
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
      write (*,'(a,f10.3,a)') '     computing ''a'' took ', &
			      & sec-oldsec,' seconds'
!
 5555 continue
! 
!  compute right hand side
      call rhs (bpnum,wpnum,body,wake,bnodex,oldbndx,wnodex,wedges, &
              & right,dt,vref,qinf,puse,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
      write (*,'(a,f10.3,a)') '     computing ''rhs'' took ', &
			      & sec-oldsec,' seconds'
! 
      return
      end
! --------------------------------------------------------------------

