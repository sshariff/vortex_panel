       program whale
c
       use panel_type
c
       implicit none
       
       integer, parameter :: nmax=20000, nsmall=1000

       type (panel)   body(nmax)
c         body and wake panels, contains strength of panel
c         and lists of nodes & edges, normal direction
c         and centroid
cc
       type (edges)   bedges(nmax)
c         edge-node arrays, lists nodes that correspond to
c         each edge
c
       type (nodes)   bnodex(nmax)
c         node coordinate array, gives x,y,z for each node
c
       real,    dimension(nsmall,nsmall,3) :: xt
       integer, dimension(nsmall)          :: topbkn,topbke
       integer, dimension(nsmall)          :: topln,tople
       integer, dimension(nsmall)          :: toprn,topre
       real,    dimension(nsmall)          :: x,y,xa,ya
       logical, dimension(nmax)            :: bodyin

       integer :: pan,i,j,node,imax,jmax,nose,diff
       integer :: n,e,p,edge
       real    :: ar,s,b,c,aoa,pi,dy,lc,amp
       
       
       
       integer :: cedge1,cedge2,cedge3,cedge4
       
       real    :: db,z,theta,fac,xi
       
       
c ---------------------------------------------------------------
       interface
          subroutine topgen (imax,jmax,bnodex,body,bedges,bodyin,
     &                       topbkn,topbke,
!     &                       topln,tople,
!     &                       toprn,topre,
     &                       xx,n,e,p,nmax,nsmall)
c
            use quadface; use panel_type
            integer,                         intent(in) :: imax,jmax
            type (nodes), dimension(nmax),   intent(inout) :: bnodex
            type (panel), dimension(nmax),   intent(inout) :: body
            type (edges), dimension(nmax),   intent(inout) :: bedges
            logical,      dimension(nmax),   intent(inout) :: bodyin
            integer,  dimension(nsmall), intent(out)   :: topbkn,topbke
!            integer,  dimension(nsmall),  intent(out)   :: topln,tople
!            integer,  dimension(nsmall),  intent(out)   :: toprn,topre
            real,  dimension(nsmall,nsmall,3), intent(in)  :: xx
            integer,                      intent(inout) :: n,e,p
	    integer,                      intent(in)    :: nmax,nsmall
c
	 end subroutine topgen
       end interface
c ---------------------------------------------------------------
       interface
          subroutine polint (xa,ya,n,x,y,dy)
c
            real, dimension(n), intent(in)  :: xa,ya
	    integer,            intent(in)  :: n
	    real,               intent(in)  :: x
	    real,               intent(out) :: y,dy
c
	 end subroutine polint
       end interface
c ---------------------------------------------------------------
       interface
          subroutine ratint (xa,ya,n,x,y,dy)
c
            real, dimension(n), intent(in)  :: xa,ya
	    integer,            intent(in)  :: n
	    real,               intent(in)  :: x
	    real,               intent(out) :: y,dy
c
	 end subroutine ratint
       end interface
c ------------------------------
!
!   imax = 2*m, where m>=2 is an integer
!   jmax = 3 + 2*n, where n is an integer
!
       jmax = 43


       ar = 3.05
       c = 0.1815
       b = ar*c
       s = b*c
       
       amp = 0.0085
!       amp = 0.0
       
       write (*,*) 'ar =',ar
       write (*,*) 'c =',c
       write (*,*) 'b =',b
       

       pi = acos(-1.0)
       aoa = 10.0*pi/180.0
 
 
!       
       xa(1)  = -1.0 + 1.0
       xa(2)  = -1.0 + 0.95
       xa(3)  = -1.0 + 0.90
       xa(4)  = -1.0 + 0.85
       xa(5)  = -1.0 + 0.80
       xa(6)  = -1.0 + 0.75
       xa(7)  = -1.0 + 0.70
       xa(8)  = -1.0 + 0.65
       xa(9)  = -1.0 + 0.60
       xa(10) = -1.0 + 0.55
       xa(11) = -1.0 + 0.50
       xa(12) = -1.0 + 0.45
       xa(13) = -1.0 + 0.40
       xa(14) = -1.0 + 0.35
       xa(15) = -1.0 + 0.30
       xa(16) = -1.0 + 0.25
       xa(17) = -1.0 + 0.20
       xa(18) = -1.0 + 0.15
       xa(19) = -1.0 + 0.10
       xa(20) = -1.0 + 0.075
       xa(21) = -1.0 + 0.05
       xa(22) = -1.0 + 0.025
       xa(23) = -1.0 + 0.0125
       xa(24) = -1.0 + 7.5e-3
       xa(25) = -1.0 + 5.0e-3
       xa(26) = -1.0 + 0.0
       
!       x(27) = -1.0 + 0.0375
!       x(28) = -1.0 + 2.50e-2
!       x(29) = -1.0 + 1.28e-2
!       x(30) = -1.0 + 6.40e-3
!       x(31) = -1.0 + 3.20e-3
!       x(32) = -1.0 + 1.60e-3
!       x(33) = -1.0 + 8.00e-4
!       x(34) = -1.0 + 4.00e-4
!       x(35) = -1.0 + 2.00e-4
!       x(36) = -1.0 + 1.00e-4
!       x(37) = -1.0 + 0.0
       
       ya(1)  = 0.0
       ya(2)  = 0.00392
       ya(3)  = 0.0113
       ya(4)  = 0.02021
       ya(5)  = 0.03054
       ya(6)  = 0.04160
       ya(7)  = 0.05290
       ya(8)  = 0.06396
       ya(9)  = 0.07441
       ya(10) = 0.08390
       ya(11) = 0.09206
       ya(12) = 0.09854
       ya(13) = 0.10298
       ya(14) = 0.10500
       ya(15) = 0.10412
       ya(16) = 0.10053
       ya(17) = 0.09410
       ya(18) = 0.08441
       ya(19) = 0.07080
       ya(20) = 0.06182
       ya(21) = 0.05065
       ya(22) = 0.03577
       ya(23) = 0.02527
       ya(24) = 0.01937
       ya(25) = 0.01583
       ya(26) = 0.0
 
       go to 60
!       go to 66
       
!       go to 40
 
!       go to 56      
!       go to 52

       x = xa
       y = ya

       imax = 50
       nose = imax/2 + 1

       x(24) = -1.0 + 5.00e-3
       x(25) = -1.0 + 1.25e-3

       go to 666
       

            
   40  continue
 
       imax = 40
       nose = imax/2 + 1

       x(1)  = -1.0 + 1.0
       x(2)  = -1.0 + 0.89
       x(3)  = -1.0 + 0.81
       x(4)  = -1.0 + 0.72
       x(5)  = -1.0 + 0.64
       x(6)  = -1.0 + 0.57
       x(7)  = -1.0 + 0.50
       x(8)  = -1.0 + 0.43
       x(9)  = -1.0 + 0.36
       x(10) = -1.0 + 0.29
       x(11) = -1.0 + 0.23
       x(12) = -1.0 + 0.18
       x(13) = -1.0 + 0.14
       x(14) = -1.0 + 0.10
       x(15) = -1.0 + 0.06
       x(16) = -1.0 + 3.0e-2
       x(17) = -1.0 + 1.0e-3
       x(18) = -1.0 + 5.0e-3
       x(19) = -1.0 + 2.5e-3
       x(20) = -1.0 + 1.0e-3
       x(21) = -1.0 + 0.0
 
 
 
!       go to 54
 
 
 66    continue
 
       imax = 66
       nose = imax/2 + 1

       x(1)  = -1.0 + 1.0
       x(2)  = -1.0 + 0.93
       x(3)  = -1.0 + 0.87
       x(4)  = -1.0 + 0.81
       x(5)  = -1.0 + 0.75
       x(6)  = -1.0 + 0.69
       x(7)  = -1.0 + 0.63
       x(8)  = -1.0 + 0.57
       x(9)  = -1.0 + 0.52
       x(10) = -1.0 + 0.46
       x(11) = -1.0 + 0.40
       x(12) = -1.0 + 0.35
       x(13) = -1.0 + 0.30
       x(14) = -1.0 + 0.25
       x(15) = -1.0 + 0.21
       x(16) = -1.0 + 0.17
       x(17) = -1.0 + 0.13
       x(18) = -1.0 + 0.105
       x(19) = -1.0 + 0.085
       x(20) = -1.0 + 0.075
       x(21) = -1.0 + 6.50e-2
       x(22) = -1.0 + 5.30e-2
       x(23) = -1.0 + 4.30e-2
       x(24) = -1.0 + 3.60e-2
       x(25) = -1.0 + 3.00e-2
       x(26) = -1.0 + 2.50e-2
       x(27) = -1.0 + 2.00e-2
       x(28) = -1.0 + 1.50e-2
       x(29) = -1.0 + 9.55e-3
       x(30) = -1.0 + 5.50e-3
       x(31) = -1.0 + 2.75e-3
       x(32) = -1.0 + 1.50e-3
       x(33) = -1.0 + 5.00e-4
       x(34) = -1.0 + 0.0

       go to 666
 
  60   continue
 
       imax = 60
       nose = imax/2 + 1

       x(1)  = -1.0 + 1.0
       x(2)  = -1.0 + 0.93
       x(3)  = -1.0 + 0.88
       x(4)  = -1.0 + 0.82
       x(5)  = -1.0 + 0.76
       x(6)  = -1.0 + 0.70
       x(7)  = -1.0 + 0.63
       x(8)  = -1.0 + 0.55
       x(9)  = -1.0 + 0.48
       x(10) = -1.0 + 0.41
       x(11) = -1.0 + 0.34
       x(12) = -1.0 + 0.28
       x(13) = -1.0 + 0.23
       x(14) = -1.0 + 0.17
       x(15) = -1.0 + 0.13
       x(16) = -1.0 + 0.102
       x(17) = -1.0 + 0.086
       x(18) = -1.0 + 7.00e-2
       x(19) = -1.0 + 5.15e-2
       x(20) = -1.0 + 4.20e-2
       x(21) = -1.0 + 3.60e-2
       x(22) = -1.0 + 3.00e-2
       x(23) = -1.0 + 2.50e-2
       x(24) = -1.0 + 2.00e-2
       x(25) = -1.0 + 1.50e-2
       x(26) = -1.0 + 9.55e-3
       x(27) = -1.0 + 5.50e-3
       x(28) = -1.0 + 2.75e-3
       x(29) = -1.0 + 1.50e-3
       x(30) = -1.0 + 5.00e-4
       x(31) = -1.0 + 0.0

       go to 666
 

 
 56    continue
 
       imax = 56
       nose = imax/2 + 1

       x(1)  = -1.0 + 1.0
       x(2)  = -1.0 + 0.93
       x(3)  = -1.0 + 0.86
       x(4)  = -1.0 + 0.78
       x(5)  = -1.0 + 0.71
       x(6)  = -1.0 + 0.63
       x(7)  = -1.0 + 0.54
       x(8)  = -1.0 + 0.46
       x(9)  = -1.0 + 0.40
       x(10) = -1.0 + 0.35
       x(11) = -1.0 + 0.31
       x(12) = -1.0 + 0.27
       x(13) = -1.0 + 0.24
       x(14) = -1.0 + 0.21
       x(15) = -1.0 + 0.18
       x(16) = -1.0 + 0.15
       x(17) = -1.0 + 0.13
       x(18) = -1.0 + 0.11
       x(19) = -1.0 + 8.0e-2
       x(20) = -1.0 + 5.5e-2
       x(21) = -1.0 + 4.00e-2
       x(22) = -1.0 + 2.75e-2
       x(23) = -1.0 + 1.50e-2
       x(24) = -1.0 + 8.00e-3
       x(25) = -1.0 + 5.0e-3
       x(26) = -1.0 + 2.75e-3
       x(27) = -1.0 + 1.3e-3
       x(28) = -1.0 + 3.00e-4
       x(29) = -1.0 + 0.0

       go to 666

  52   continue
       
       imax = 52
       nose = imax/2 + 1

       x(1)  = -1.0 + 1.00
       x(2)  = -1.0 + 0.92
       x(3)  = -1.0 + 0.84
       x(4)  = -1.0 + 0.77
       x(5)  = -1.0 + 0.70
       x(6)  = -1.0 + 0.63
       x(7)  = -1.0 + 0.57
       x(8)  = -1.0 + 0.51
       x(9)  = -1.0 + 0.45
       x(10) = -1.0 + 0.41
       x(11) = -1.0 + 0.37
       x(12) = -1.0 + 0.33
       x(13) = -1.0 + 0.30
       x(14) = -1.0 + 0.27
       x(15) = -1.0 + 0.24
       x(16) = -1.0 + 0.21
       x(17) = -1.0 + 0.18
       x(18) = -1.0 + 0.15
       x(19) = -1.0 + 0.125
       x(20) = -1.0 + 7.5e-2
       x(21) = -1.0 + 4.00e-2
       x(22) = -1.0 + 2.00e-2
       x(23) = -1.0 + 1.0e-2
       x(24) = -1.0 + 5.0e-3
       x(25) = -1.0 + 2.5e-3
       x(26) = -1.0 + 1.00e-3
       x(27) = -1.0 + 0.0

!       go to 4056t

  666  continue
    
!       fac = -1.7
       do 223 i = 1, nose
!          xi = float(i-1)/float(nose-1)
!	  if (i.gt.nose/2) then
!!             x(i) = -fac*xi*xi*xi + 1.5*fac*xi*xi - (2.0+fac)*xi/2.0! - 1.0
!             x(i) = -1.0 + (1.0-xi)**2.0
!          else
!	     x(i) = -1.5*xi
!	  endif
	  
          call ratint (xa(1:26),ya(1:26),26,x(i),y(i),dy)
	  
!	  write (*,*) 
!	  write (*,*) 'interpolating at x =',x(i)
!	  write (*,*) '  f(x) = ',y(i)
!	  write (*,*) '  dy =',dy
	  
	  
 223   continue
 
 4056  continue
 
 
 

       do 41 i = nose+1, imax+1
          diff = i - nose
	  x(i) = x(nose-diff)
	  y(i) = y(nose-diff)
!	  write (*,*) 'i =',i,',   x =',x(i)
   41  continue
   
   
!       do 990 j = 1,jmax
          x(nose) = -1.0 + 2.5e-4
!  990  continue
   
!      go to 223
   
   
   
!       do 1 i = 1, imax + 1
!          if (x(i).ne.0.0) then
!	     ya(i) = 0.12/0.20*(0.29690*sqrt(1.0+x(i)) - 
!     &              0.126*(1.0+x(i)) - 0.3516*(1.0+x(i))*(1.0+x(i)) 
!     &            + 0.2843*(1.0+x(i))**3.0 - 0.1015*(1.0+x(i))**4.0)
!          else
!	     ya(i) = 0.0
!	  endif
!	  write (*,*) 'i =',i,',   y =',y(i)
! 1     continue

!       x = xa
!       y = ya
       
 
       
!   
!       y(2)    = y(1) + (y(3)-y(1))*(x(2)-x(1))/(x(3)-x(1))
!       y(imax) = y(imax+1)+(y(imax-1)-y(imax+1))*
!     &          (x(imax)-x(imax+1))/(x(imax-1)-x(imax+1))
!
! since wing is pointing in the -x direction, rotate about -aoa       
       aoa = -aoa
       
       do 78 j = 1,jmax

           z = b/2.0 - b*(real(j-1)/real(jmax-1))
       
           lc = c + amp*sin(2.0*pi*z/0.073)

!           if (j.eq.1) then
!	      z =  b/2.0
!	   else if (j.eq.jmax) then
!	      z = -b/2.0
!	   else 
!              z = 1.1*(b/2.0 - b*(real(j-1)/real(jmax-1)))
!	   endif
!           z = b/2.0*sin(pi/2.0-real(j-1)/real(jmax-1)*pi)
!           theta = pi/2.0-real(j-1)/real(jmax-1)*pi
!           z = b/4*(sin(theta)+sin(theta/2)/sin(pi/4.0))
!
          do 77 i = 1,imax + 1

             if (j.eq.1.or.j.eq.jmax) then
	        xt(i,j,1) = lc*x(i)*cos(aoa) 
                xt(i,j,2) = lc*x(i)*sin(aoa)
	     else if (j.eq.jmax) then
	        xt(i,j,1) = lc*x(i)*cos(aoa) 
                xt(i,j,2) = lc*x(i)*sin(aoa)
	     elseif (i.le.nose) then
!             if (i.le.nose) then
	        xt(i,j,1) = lc*(x(i)*cos(aoa) - y(i)*sin(aoa))
                xt(i,j,2) = lc*(x(i)*sin(aoa) + y(i)*cos(aoa))
             else
	        xt(i,j,1) = lc*(x(i)*cos(aoa) + y(i)*sin(aoa))
                xt(i,j,2) = lc*(x(i)*sin(aoa) - y(i)*cos(aoa))

             end if
	         
             xt(i,j,3) =  z
	     
!             if (j.eq.1) then
!             if (i.eq.1.or.i.eq.2.or.i.eq.imax.or.i.eq.imax+1) then
!	        write (*,*) 'i =',i,',   xt =',xt(i,j,:)	        
!	     end if
!	     end if

	     
!	     
 77       continue
 78    continue



       go to 7364

       
       xt(1,1,1) = (xt(1,1,1)+xt(2,2,1))/2.0
       xt(1,1,3) = (xt(1,1,3)+xt(2,2,3))/2.0
       xt(1,jmax,1) = (xt(2,jmax-1,1)+xt(1,jmax,1))/2.0
       xt(1,jmax,3) = (xt(2,jmax-1,3)+xt(1,jmax,3))/2.0
       xt(nose,1,1) = (xt(nose,1,1)+xt(nose-1,2,1))/2.0
       xt(nose,1,3) = (xt(nose,1,3)+xt(nose-1,2,3))/2.0
       xt(nose,jmax,1) = (xt(nose-1,jmax-1,1)+xt(nose,jmax,1))/2.0
       xt(nose,jmax,3) = (xt(nose-1,jmax-1,3)+xt(nose,jmax,3))/2.0
       xt(imax+1,1,1) = (xt(imax+1,1,1)+xt(imax,2,1))/2.0
       xt(imax+1,1,3) = (xt(imax+1,1,3)+xt(imax,2,3))/2.0
       xt(imax+1,jmax,1) = (xt(imax+1,jmax,1)+xt(imax,jmax-1,1))/2.0
       xt(imax+1,jmax,3) = (xt(imax+1,jmax,3)+xt(imax,jmax-1,3))/2.0

 7364 continue

c       
       p = 0
       e = 0
       n = 0
c
       call topgen (imax,jmax,bnodex,body,bedges,bodyin,
     &              topbkn,topbke,
!     &              topln,tople,
!     &              toprn,topre,
     &              xt,n,e,p,nmax,nsmall)
c

       open (3,file='bfish')
       open (4,file='efish')
       open (5,file='pfish')
       open (8,file='fish.tec')

       write (3,'(f12.6)') 1.0/float(imax-1)
       write (3,*) n
       write (3,*) 1
       do 333 node = 1,n
          write (3,'(i4,3e16.6)') node,bnodex(node)%x
 333   continue
       close (3)
       
       write (4,*) e
       do 334 edge = 1,e
          write (4,'(i8)') edge
	  write (4,'(2i8)') bedges(edge)%node1,bedges(edge)%node2
	  write (4,'(2i8)') bedges(edge)%panel1,bedges(edge)%panel2 
 334   continue
c
       write (4,*) 1
       write (4,*) jmax-1
       do 335 j = 1,jmax-1
          write (4,*) topbke(j)
 335   continue
       
       write (5,*) p
       do 336 pan = 1,p
          write (5,'(i8)') pan
	  write (5,'(3i8)') body(pan)%node
	  write (5,'(3i8)') body(pan)%edge
	  if (bodyin(pan)) then
	     write (5,'(i6)') 1
	  else
	     write (5,'(i6)') 0
	  endif
 336   continue

       write(8,*) 'title="whale"'
       write(8,*) 'variables = x,y,z'
       write(8,909)  'zone t=body, n=',n,', e=',p,
     &                   ', f=fepoint, et=triangle'
 909   format(a,i6,a,i6,a)
       do 5 i=1,n
          write(8,6) bnodex(i)%x
  6       format(3e16.6)
  5    continue
       do 8 i=1,p
          write(8,'(3i6)') body(i)%node
  8    continue

       close(3)
       close(4)
       close(5)
       close(8)
       
       
       
       end
c -------------------------------------------------------------------
c -------------------------------------------------------------------
       subroutine topgen (imax,jmax,bnodex,body,bedges,bodyin,
     &                    topbkn,topbke,
!     &                    topln,tople,
!     &                    toprn,topre,
     &                    xx,n,e,p,nmax,nsmall)
c ------------------------------
       use quadface; use panel_type
c ------------------------------
       implicit none
            integer,                       intent(in)    :: imax,jmax
            type (nodes), dimension(nmax), intent(inout) :: bnodex
            type (panel), dimension(nmax), intent(inout) :: body
            type (edges), dimension(nmax), intent(inout) :: bedges
 	    integer,    dimension(nsmall), intent(out)   :: topbkn,topbke
            logical,      dimension(nmax), intent(inout) :: bodyin
            real, dimension(nsmall,nsmall,3), intent(in) :: xx
            integer,                       intent(inout) :: n,e,p
	    integer,                       intent(in)    :: nmax,nsmall
c ------------------------------
!            integer,    dimension(nsmall), intent(out) :: topln,tople
!            integer,    dimension(nsmall), intent(out) :: toprn,topre

!       integer, dimension(nsmall) :: topbkn,topbke
       integer, dimension(nsmall) :: topln,tople
       integer, dimension(nsmall) :: toprn,topre

       integer, dimension(nmax) :: newnode,newedge
       real,    dimension(3)    :: x1,x2,x3,x4
       integer :: i,j
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2
       logical :: top
c ------------------------------
c
!       write (*,*) 'in topgen...'
c
       top = .true.
c      
       newnode  = 0
       newedge  = 0
c
       i = 1
       j = 1
c
       n  = n + 4
       n1 = n - 3
       n2 = n - 2
       n3 = n - 1
       n4 = n
       e  = e + 5
       e1 = e - 4
       e2 = e - 3
       e3 = e - 2
       e4 = e - 1
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p
!
       x1 = xx(i  ,j  ,1:3)
       x2 = xx(i  ,j+1,1:3)
       x3 = xx(i+1,j  ,1:3)
       x4 = xx(i+1,j+1,1:3)
!
       call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
       newnode(j)   = n3
       newnode(j+1) = n4
       newedge(j)   = e2
c
       topln(i)   = n1
       topln(i+1) = n3
       tople(i)   = e1
c
       topbkn(j)   = n1
       topbkn(j+1) = n2
       topbke(j)   = e4
c	  
       n  = n + 2
       n3 = n2
       n1 = n - 1
       n2 = n 
!      n4 = n4
       e  = e + 4
       e1 = e - 3
       e2 = e3
       e3 = e - 2
       e4 = e - 1
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p
c       
       x1 = xx(i  ,j+2,1:3)
       x2 = xx(i+1,j+2,1:3)
       x3 = xx(i  ,j+1,1:3)
       x4 = xx(i+1,j+1,1:3)
c
       call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
       newnode(j+1) = n4
       newnode(j+2) = n2
       newedge(j+1) = e3
c
       topbkn(j+1) = n3
       topbkn(j+2) = n1
       topbke(j+1) = e1
c
       toprn(i)   = n1
       toprn(i+1) = n2
       topre(i)   = e4
c
       n  = n + 2
       n1 = n - 1
       n2 = newnode(j)
       n3 = n
       n4 = newnode(j+1)
       e  = e + 4
       e1 = e - 3
       e2 = e - 2
       e3 = newedge(j)
       e4 = e - 1
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p
c       
       x1 = xx(i+2,j  ,1:3)
       x2 = xx(i+1,j  ,1:3)
       x3 = xx(i+2,j+1,1:3)
       x4 = xx(i+1,j+1,1:3)
c
       call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
       newnode(j)   = n1
       newnode(j+1) = n3
       newedge(j)   = e1
c
       topln(i+1) = n2
       topln(i+2) = n1
       tople(i+1) = e4
c
       n  = n + 1
       n1 = n 
       n2 = n3
       n3 = newnode(j+2)
!       n4 = n4
       e  = e + 3
       e1 = e - 2
       e3 = e2
       e2 = newedge(j+1)
       e4 = e - 1
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p
c       
       x1 = xx(i+2,j+2,1:3)
       x2 = xx(i+2,j+1,1:3)
       x3 = xx(i+1,j+2,1:3)
       x4 = xx(i+1,j+1,1:3)
c
       call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
       newnode(j+1) = n2
       newnode(j+2) = n1
       newedge(j+1) = e4
c
       toprn(i+1) = n3
       toprn(i+2) = n1
       topre(i+1) = e1
c
       do 6798 j = 3, jmax-2, 2 
       
c
	  n  = n + 2
	  n1 = toprn(i)
	  n2 = n - 1
	  n3 = toprn(i+1)
	  n4 = n
	  e  = e + 4
	  e1 = topre(i)
	  e2 = e - 3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j  ,1:3)
	  x2 = xx(i  ,j+1,1:3)
	  x3 = xx(i+1,j  ,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j)   = n3
	  newnode(j+1) = n4
	  newedge(j)   = e2
c
          topbkn(j)   = n1
	  topbkn(j+1) = n2
	  topbke(j)   = e4
c	  
	  n  = n + 2
	  n1 = n - 1
	  n3 = n2
	  n2 = n 
!         n4 = n4
	  e  = e + 4
	  e1 = e - 3
	  e2 = e3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j+2,1:3)
	  x2 = xx(i+1,j+2,1:3)
	  x3 = xx(i  ,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
	  call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
	  newnode(j+1) = n4
	  newnode(j+2) = n2
	  newedge(j+1) = e3
c
	  topbkn(j+1) = n3
	  topbkn(j+2) = n1
	  topbke(j+1) = e1
c
	  toprn(i)   = n1
	  toprn(i+1) = n2
	  topre(i)   = e4
c
	  n  = n + 1
	  n1 = toprn(i+2)
	  n2 = newnode(j) 
	  n3 = n
	  n4 = newnode(j+1)
	  e  = e + 3
	  e1 = e - 2 
	  e2 = e - 1 
	  e3 = newedge(j)
	  e4 = topre(i+1)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i+2,j  ,1:3)
	  x2 = xx(i+1,j  ,1:3)
	  x3 = xx(i+2,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
	  call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
	  newnode(j)   = n1
	  newnode(j+1) = n3
	  newedge(j)   = e1
c
	  n  = n + 1
	  n1 = n 
	  n2 = n3
	  n3 = newnode(j+2)
!         n4 = n4
	  e  = e + 3
	  e1 = e - 2
	  e3 = e2
	  e2 = newedge(j+1)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i+2,j+2,1:3)
	  x2 = xx(i+2,j+1,1:3)
	  x3 = xx(i+1,j+2,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
	  call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &            e1,e2,e3,e4,e5,p1,p2,
     &            bodyin,bnodex,body,bedges,top)
c
	  newnode(j+1) = n2
	  newnode(j+2) = n1
	  newedge(j+1) = e4
c
	  toprn(i+1) = n3
	  toprn(i+2) = n1
	  topre(i+1) = e1
c
 6798  continue
c

       do 6998 i = 3, imax-3, 2
c
	  j  = 1
	  
	  n  = n + 2
	  n1 = newnode(j)
	  n2 = newnode(j+1)
	  n3 = n - 1
	  n4 = n
	  e  = e + 4
	  e1 = e - 3
	  e2 = e - 2 
	  e3 = e - 1
	  e4 = newedge(j)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j  ,1:3)
	  x2 = xx(i  ,j+1,1:3)
	  x3 = xx(i+1,j  ,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j)   = n3
	  newnode(j+1) = n4
	  newedge(j)   = e2
c
          topln(i)   = n1
          topln(i+1) = n3
          tople(i)   = e1
c
	  n  = n + 1
	  n1 = newnode(j+2)
	  n3 = n2
	  n2 = n 
!	  n4 = n4
	  e  = e + 3
	  e1 = newedge(j+1)
	  e2 = e3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j+2,1:3)
	  x2 = xx(i+1,j+2,1:3)
	  x3 = xx(i  ,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j+1) = n4
	  newnode(j+2) = n2
	  newedge(j+1) = e3
c
          toprn(i)   = n1
	  toprn(i+1) = n2
	  topre(i)   = e4
c
	  n  = n + 2
	  n1 = n - 1
	  n2 = newnode(j)
	  n3 = n
	  n4 = newnode(j+1)
	  e  = e + 4
	  e1 = e - 3
	  e2 = e - 2
	  e3 = newedge(j)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c       
	  x1 = xx(i+2,j  ,1:3)
	  x2 = xx(i+1,j  ,1:3)
	  x3 = xx(i+2,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j)   = n1
	  newnode(j+1) = n3
	  newedge(j)   = e1
c
          topln(i+1) = n2
	  topln(i+2) = n1
	  tople(i+1) = e4
c
	  n  = n + 1
	  n1 = n 
	  n2 = n3
	  n3 = newnode(j+2)
!         n4 = n4
	  e  = e + 3
	  e1 = e - 2
	  e3 = e2
	  e2 = newedge(j+1)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c       
	  x1 = xx(i+2,j+2,1:3)
	  x2 = xx(i+2,j+1,1:3)
	  x3 = xx(i+1,j+2,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j+1) = n2
	  newnode(j+2) = n1
	  newedge(j+1) = e4
c
          toprn(i+1) = n3
          toprn(i+2) = n1
          topre(i+1) = e1
c
          do 6898 j = 3, jmax-2, 2 
c
	     n  = n + 1
	     n1 = toprn(i)
	     n2 = newnode(j+1)
	     n3 = toprn(i+1)
	     n4 = n
	     e  = e + 3
	     e1 = topre(i)
	     e2 = e - 2
	     e3 = e - 1
	     e4 = newedge(j)
	     e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c       
	     x1 = xx(i  ,j  ,1:3)
	     x2 = xx(i  ,j+1,1:3)
	     x3 = xx(i+1,j  ,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j)   = n3
	     newnode(j+1) = n4
	     newedge(j)   = e2
c
	     n  = n + 1
	     n1 = newnode(j+2)
	     n3 = n2
	     n2 = n 
!	     n4 = n4
	     e  = e + 3
	     e1 = newedge(j+1)
	     e2 = e3
	     e3 = e - 2
	     e4 = e - 1
	     e5 = e
	     p  = p + 2 
	     p1 = p - 1
	     p2 = p
c
	     x1 = xx(i  ,j+2,1:3)
	     x2 = xx(i+1,j+2,1:3)
	     x3 = xx(i  ,j+1,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j+1) = n4
	     newnode(j+2) = n2
	     newedge(j+1) = e3
c
             toprn(i)   = n1
	     toprn(i+1) = n2
	     topre(i)   = e4
c
	     n  = n + 1
	     n1 = toprn(i+2)
	     n2 = newnode(j) 
	     n3 = n
	     n4 = newnode(j+1)
	     e  = e + 3
	     e1 = e - 2 
	     e2 = e - 1 
	     e3 = newedge(j)
	     e4 = topre(i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
c       
	     x1 = xx(i+2,j  ,1:3)
	     x2 = xx(i+1,j  ,1:3)
	     x3 = xx(i+2,j+1,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j)   = n1
	     newnode(j+1) = n3
	     newedge(j)   = e1
c
	     n  = n + 1
	     n1 = n 
	     n2 = n3
	     n3 = newnode(j+2)
!            n4 = n4
	     e  = e + 3
	     e1 = e - 2
	     e3 = e2
	     e2 = newedge(j+1)
	     e4 = e - 1
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
c       
	     x1 = xx(i+2,j+2,1:3)
	     x2 = xx(i+2,j+1,1:3)
	     x3 = xx(i+1,j+2,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j+1) = n2
	     newnode(j+2) = n1
	     newedge(j+1) = e4
c
             toprn(i+1) = n3
	     toprn(i+2) = n1
	     topre(i+1) = e1

 6898     continue
c
 6998  continue
c
       i = imax - 1
c
	  j  = 1
	  
	  n  = n + 2
	  n1 = newnode(j)
	  n2 = newnode(j+1)
	  n3 = n - 1
	  n4 = n
	  e  = e + 4
	  e1 = e - 3
	  e2 = e - 2 
	  e3 = e - 1
	  e4 = newedge(j)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j  ,1:3)
	  x2 = xx(i  ,j+1,1:3)
	  x3 = xx(i+1,j  ,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j)   = n3
	  newnode(j+1) = n4
	  newedge(j)   = e2
c
          topln(i)   = n1
          topln(i+1) = n3
          tople(i)   = e1
c
	  n  = n + 1
	  n1 = newnode(j+2)
	  n3 = n2
	  n2 = n 
!	  n4 = n4
	  e  = e + 3
	  e1 = newedge(j+1)
	  e2 = e3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i  ,j+2,1:3)
	  x2 = xx(i+1,j+2,1:3)
	  x3 = xx(i  ,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j+1) = n4
	  newnode(j+2) = n2
	  newedge(j+1) = e3
c
          toprn(i)   = n1
	  toprn(i+1) = n2
	  topre(i)   = e4
c
!	  n  = n 
	  n1 = topbkn(j)
	  n2 = newnode(j)
	  n3 = topbkn(j+1)
	  n4 = newnode(j+1)
	  e  = e + 3
	  e1 = topbke(j)
	  e2 = e - 2
	  e3 = newedge(j)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i+2,j  ,1:3)
	  x2 = xx(i+1,j  ,1:3)
	  x3 = xx(i+2,j+1,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j)   = n1
	  newnode(j+1) = n3
	  newedge(j)   = e1
c
          topln(i+1) = n2
	  topln(i+2) = n1
	  tople(i+1) = e4
c
!	  n  = n 
	  n1 = topbkn(j+2)
	  n2 = n3
	  n3 = newnode(j+2)
!         n4 = n4
	  e  = e + 2
	  e1 = e - 1
	  e3 = e2
	  e2 = newedge(j+1)
	  e4 = topbke(j+1)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c
	  x1 = xx(i+2,j+2,1:3)
	  x2 = xx(i+2,j+1,1:3)
	  x3 = xx(i+1,j+2,1:3)
	  x4 = xx(i+1,j+1,1:3)
c
          call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &               e1,e2,e3,e4,e5,p1,p2,
     &               bodyin,bnodex,body,bedges,top)
c
          newnode(j+1) = n2
	  newnode(j+2) = n1
	  newedge(j+1) = e4
c
          toprn(i+1) = n3
          toprn(i+2) = n1
          topre(i+1) = e1
c
          do 7898 j = 3, jmax-2, 2 
c
	     n  = n + 1
	     n1 = toprn(i)
	     n2 = newnode(j+1)
	     n3 = toprn(i+1)
	     n4 = n
	     e  = e + 3
	     e1 = topre(i)
	     e2 = e - 2
	     e3 = e - 1
	     e4 = newedge(j)
	     e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
c       
	     x1 = xx(i  ,j  ,1:3)
	     x2 = xx(i  ,j+1,1:3)
	     x3 = xx(i+1,j  ,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j)   = n3
	     newnode(j+1) = n4
	     newedge(j)   = e2
c
	     n  = n + 1
	     n1 = newnode(j+2)
	     n3 = n2
	     n2 = n 
!	     n4 = n4
	     e  = e + 3
	     e1 = newedge(j+1)
	     e2 = e3
	     e3 = e - 2
	     e4 = e - 1
	     e5 = e
	     p  = p + 2 
	     p1 = p - 1
	     p2 = p
c
	     x1 = xx(i  ,j+2,1:3)
	     x2 = xx(i+1,j+2,1:3)
	     x3 = xx(i  ,j+1,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j+1) = n4
	     newnode(j+2) = n2
	     newedge(j+1) = e3
c
             toprn(i)   = n1
	     toprn(i+1) = n2
	     topre(i)   = e4
c
!            n  = n
	     n1 = topbkn(j)
	     n2 = newnode(j)
	     n3 = topbkn(j+1)
	     n4 = newnode(j+1)
	     e  = e + 2
	     e1 = topbke(j)
	     e2 = e - 1
	     e3 = newedge(j)
	     e4 = topre(i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
c       
	     x1 = xx(i+2,j  ,1:3)
	     x2 = xx(i+1,j  ,1:3)
	     x3 = xx(i+2,j+1,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j)   = n1
	     newnode(j+1) = n3
	     newedge(j)   = e1
c
!	     n  = n
	     n1 = topbkn(j+2) 
	     n2 = n3
	     n3 = newnode(j+2)
!            n4 = n4
	     e  = e + 2
	     e1 = e - 1
	     e3 = e2
	     e2 = newedge(j+1)
	     e4 = topbke(j+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
c       
	     x1 = xx(i+2,j+2,1:3)
	     x2 = xx(i+2,j+1,1:3)
	     x3 = xx(i+1,j+2,1:3)
	     x4 = xx(i+1,j+1,1:3)
c
             call quad (x1,x2,x3,x4,n1,n2,n3,n4,
     &                  e1,e2,e3,e4,e5,p1,p2,
     &                  bodyin,bnodex,body,bedges,top)
c
             newnode(j)   = n2
	     newnode(j+1) = n1
	     newedge(j)   = e4
c
             toprn(i+1) = n3
	     toprn(i+2) = n1
	     topre(i+1) = e1

 7898     continue
c
 7998  continue
c
c
       end subroutine topgen
c -------------------------------------------------------------------
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=40)
      INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
c -------------------------------------------------------------------
      SUBROUTINE ratint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=40,TINY=1.e-25)
      INTEGER i,m,ns
      REAL dd,h,hh,t,w,c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.)pause 'failure in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
c -------------------------------------------------------------------
