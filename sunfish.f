! ---------------------------------------------------------------
       program sunfish
!
       use panel_type; use sunblock
!
       implicit none
       
       integer, parameter :: nmax=20000, nsmall=1000

       type (panel)   body(nmax)
!         body and wake panels, contains strength of panel
!         and lists of nodes & edges, normal direction
!         and centroid
!!
       type (edges)   bedges(nmax)
!         edge-node arrays, lists nodes that correspond to
!         each edge
!
       type (nodes)   bnodex(nmax)
!         node coordinate array, gives x,y,z for each node
!
       real,    dimension(nsmall,nsmall,3) :: xt,xb
       integer, dimension(nsmall)          :: topln,toprn
       integer, dimension(nsmall)          :: tople,topre
       integer, dimension(nsmall)          :: botle,botre
       integer, dimension(nsmall)          :: botln,botrn
       real,    dimension(nsmall,nsmall)   :: x,y,z
       logical, dimension(nmax)            :: bodyin
       logical, dimension(nsmall,nsmall)   :: kuttat,kuttab
       logical, dimension(nmax)            :: kedge,knode
       logical, dimension(nsmall,nsmall)   :: caudalt,caudalb
       logical, dimension(nsmall,nsmall)   :: oddyt,oddyb
       
       character*40 :: inname,outname,bname,pname,ename,tname

       integer                      :: pan,i,j,jj,node,imax,nose,diff
       integer, dimension(nsmall,3) :: jmaxt,jmaxb
       integer                      :: jubt,jubb,last,mid,jdum
       integer                      :: n,e,p,edge,kutt,caud
       real                         :: ar,s,b,c,aoa,pi,lc,amp
       logical, dimension(nsmall)   :: kdum
       real,    dimension(nsmall)   :: xdum
       integer                      :: dj,dj2,di
       integer :: ptl,ptr,pbl,pbr
       logical :: bottom
       real    :: db,theta,fac,xi
       
       integer :: infile = 9
       integer :: bfile  = 12
       integer :: pfile  = 13
       integer :: efile  = 14
       integer :: tfile  = 3
       integer :: keybd  = 5
       integer :: scrn   = 6

! ---------------------------------------------------------------
!
! start:
!
!     
      p = 0
      e = 0
      n = 0

      kuttat  = .false.
      caudalt = .false.
      kuttab  = .false.
      caudalb = .false.

      kedge = .false.
      knode = .false.
      
      pi = acos(-1.0)
      aoa = 0.0*pi/180.0
!
! since wing is pointing in the -x direction, rotate about -aoa       
      aoa = -aoa
 
      write (scrn,*) 'Name of input file?'
      read  (keybd,'(a)') inname
      write (scrn,*)
      write (scrn,*) 'Name of new node file?'
      read  (keybd,'(a)') bname
      write (scrn,*)
      write (scrn,*) 'Name of new panel file?'
      read  (keybd,'(a)') pname
      write (scrn,*)
      write (scrn,*) 'Name of new edge file?'
      read  (keybd,'(a)') ename
      write (scrn,*)
      write (scrn,*) 'Name of Tecplot file?'
      read  (keybd,'(a)') tname
      
      open (infile,file=inname, status='old')
      open (bfile,file=bname)
      open (pfile,file=pname)
      open (efile,file=ename)
      open (tfile,file=tname)
      
      
!      go to 7877
      
      read (infile,*) imax,jubt
      
      write (*,*) imax,jubt

      write (*,*)
      write (*,*) 'Generating top...'
      
      do i = 1,imax 

         j = 1
         read (infile,*) jdum,jmaxt(i,1),x(i,j),y(i,j),z(i,j),kutt,caud
	 if (kutt/=0) kuttat(i,j)  = .true.
	 if (caud/=0) caudalt(i,j) = .true.
!         write (*,*) j,jmax(i,1),x(i,j),y(i,j),z(i,j),kutta(i,j),caudal(i,j)
      
         do jj = 2,jmaxt(i,1)

            read (infile,*) jdum,jdum,x(i,jj),y(i,jj),z(i,jj),kutt,caud
            if (kutt/=0) kuttat(i,jj)  = .true.
	    if (caud/=0) caudalt(i,jj) = .true.
!         write (*,*) j,jmax(i,1),x(i,j),y(i,j),z(i,j),kutta(i,j),caudal(i,j)
         end do
! 
      end do
      
!      write (*,*)
!      write (*,*) 'imax =',imax
!      do i = 1,imax
!         write (*,*) 'jmaxt(',i,',1) =',jmaxt(i,1)
!      end do
!      write (*,*)
!      write (*,*) 'calling pack....' 
      
! pack strips together in proper order & deal with caudal points
      call pack (imax,jmaxt,kuttat,caudalt,x,y,z,nmax,nsmall,oddyt)	     

      do i = 1,imax
         mid = jmaxt(i,1)
	 last = 2*mid - 1
!
         xt (i, 1:mid, 1) = x (i, 1:mid)
         xt (i, 1:mid, 2) = y (i, 1:mid)
         xt (i, 1:mid, 3) = z (i, 1:mid)

! the third number in the subscript triplets is the stride. this is needed
! for the subarray with decreasing indices	 
         xt (i, mid+1:last:1, 1) =  x (i, mid-1:1:-1)
         xt (i, mid+1:last:1, 2) =  y (i, mid-1:1:-1)
         xt (i, mid+1:last:1, 3) = -z (i, mid-1:1:-1)

	 strip2: if (jmaxt(i,2) /= 0) then
	 
	    dj = last
	    mid = jmaxt(i,2)
	    dj2 = jmaxt(i,1)
	    last = 2*mid - 1
	   
	    xt (i, 1+dj:mid+dj, 1) = x (i, 1+dj2:mid+dj2)
	    xt (i, 1+dj:mid+dj, 2) = y (i, 1+dj2:mid+dj2)
	    xt (i, 1+dj:mid+dj, 3) = z (i, 1+dj2:mid+dj2)
	   
            xt (i, mid+1+dj:last+dj:1, 1) =  x (i, mid-1+dj2:1+dj2:-1)
            xt (i, mid+1+dj:last+dj:1, 2) =  y (i, mid-1+dj2:1+dj2:-1)
            xt (i, mid+1+dj:last+dj:1, 3) = -z (i, mid-1+dj2:1+dj2:-1)
	    
	    strip3: if (jmaxt(i,3) /= 0) then
	    
	       dj = last + dj
	       mid = jmaxt(i,3)
	       dj2 = jmaxt(i,1)+jmaxt(i,2)
	       last = 2*mid - 1

	       xt (i, 1+dj:mid+dj, 1) = x (i, 1+dj2:mid+dj2)
	       xt (i, 1+dj:mid+dj, 2) = y (i, 1+dj2:mid+dj2)
	       xt (i, 1+dj:mid+dj, 3) = z (i, 1+dj2:mid+dj2)

               xt (i, mid+1+dj:last+dj:1, 1) =  x (i, mid-1+dj2:1+dj2:-1)
               xt (i, mid+1+dj:last+dj:1, 2) =  y (i, mid-1+dj2:1+dj2:-1)
               xt (i, mid+1+dj:last+dj:1, 3) = -z (i, mid-1+dj2:1+dj2:-1)

	    end if strip3
	 
	 end if strip2
	 	 
      end do
!
      bottom = .false.
!      
      botle = 0
      botre = 0
      
!      
! note: botle & botre are dummies here. Those arguments are used in the
!       next call to gen
!

      call gen (imax,jmaxt,bnodex,body,bedges,bodyin,bottom,kuttat,   &
              & caudalt,xt,n,e,p,nmax,nsmall,tople,topre,botle,botre, &
	      & topln,toprn,botln,botrn,oddyt,knode)
!
      write (*,*)
      write (*,*) '...top finished.'
      write (*,*)

!      go to 7878

 7877 continue
      
      
      write (*,*) 'Generating bottom...'
!  
      read (infile,*) imax,jubb
      
      write (*,*) imax,jubb

      di = 0

      do i = 1,imax 

         j = 1
         read (infile,*) jdum,jmaxb(i,1),x(i,j),y(i,j),z(i,j),kutt,caud
	 if (kutt/=0) kuttab(i,j)  = .true.
	 if (caud/=0) caudalb(i,j) = .true.
!         write (*,*) j,jmax(i),x(i,j),y(i,j),z(i,j),kutta(i,j),caudal(i,j)
      
         do jj = 2,jmaxb(i,1)

           read (infile,*) jdum,jmaxb(i,1),x(i,jj),y(i,jj),z(i,jj),kutt,caud
	   if (kutt/=0) kuttab(i,jj)  = .true.
	   if (caud/=0) caudalb(i,jj) = .true.
!         write (*,*) j,jmax(i),x(i,j),y(i,j),z(i,j),kutta(i,j),caudal(i,j)
         end do
! 
! apparently, we have to reverse the order of the j's on the bottom 
         xdum (1:jmaxb(i,1):1) = x (i,1:jmaxb(i,1):1)
	 x  (i,1:jmaxb(i,1):1) = xdum (jmaxb(i,1):1:-1)
	 
         xdum (1:jmaxb(i,1):1) = y (i,1:jmaxb(i,1):1)
	 y  (i,1:jmaxb(i,1):1) = xdum  (jmaxb(i,1):1:-1)
	 
         xdum (1:jmaxb(i,1):1) = z (i,1:jmaxb(i,1):1)
	 z  (i,1:jmaxb(i,1):1) = xdum  (jmaxb(i,1):1:-1)
	 
         kdum     (1:jmaxb(i,1):1) = kuttab (i,1:jmaxb(i,1):1)
	 kuttab (i,1:jmaxb(i,1):1) = kdum       (jmaxb(i,1):1:-1)
	 
         kdum      (1:jmaxb(i,1):1) = caudalb (i,1:jmaxb(i,1):1)
	 caudalb (i,1:jmaxb(i,1):1) = kdum        (jmaxb(i,1):1:-1)
	 
! apparently, we also have to shift the data and copy the center-line data 
! from the top half (for strips bordering the top of the fish)

         if (y(i,1)==0.0) then

	    xdum (1:jubt) = x (i,1:jubt)
	    x (i,2:jubt+1) = xdum (1:jubt)
	    x (i,1) = xt (i-di,1,1)

	    xdum (1:jubt) = y (i,1:jubt)
	    y (i,2:jubt+1) = xdum (1:jubt)
	    y (i,1) = xt (i-di,1,2)

	    xdum (1:jubt) = z (i,1:jubt)
	    z (i,2:jubt+1) = xdum (1:jubt)
	    z (i,1) = xt (i-di,1,3)

	    kdum (1:jubt) = kuttab (i,1:jubt)
	    kuttab (i,2:jubt+1) = kdum (1:jubt)
	    kuttab (1,1) = kuttat (i-di,1)

	    kdum (1:jubt) = caudalb (i,1:jubt)
	    caudalb (i,2:jubt+1) = kdum (1:jubt)
	    caudalb (1,1) = caudalt (i-di,1)

	    jmaxb (i,1) = jmaxb (i,1) + 1
	 
	 else
	    
	    di = di + 1
	    
	 end if

      end do

!
      call pack (imax,jmaxb,kuttab,caudalb,x,y,z,nmax,nsmall,oddyb)	     
!
!	
      do i = 1,imax
         mid = jmaxb(i,1)
	 last = 2*mid - 1
!
         xb (i, 1:mid, 1) = x (i, 1:mid)
         xb (i, 1:mid, 2) = y (i, 1:mid)
         xb (i, 1:mid, 3) = z (i, 1:mid)

! the third number in the subscript triplets is the stride. this is needed
! for the subarray with decreasing indices	 
         xb (i, mid+1:last:1, 1) =  x (i, mid-1:1:-1)
         xb (i, mid+1:last:1, 2) =  y (i, mid-1:1:-1)
         xb (i, mid+1:last:1, 3) = -z (i, mid-1:1:-1)

	 strip2: if (jmaxb(i,2) /= 0) then
	 
	    dj = last
	    mid = jmaxb(i,2)
	    dj2 = jmaxb(i,1)
	    last = 2*mid - 1
	   
	    xb (i, 1+dj:mid+dj, 1) = x (i, 1+dj2:mid+dj2)
	    xb (i, 1+dj:mid+dj, 2) = y (i, 1+dj2:mid+dj2)
	    xb (i, 1+dj:mid+dj, 3) = z (i, 1+dj2:mid+dj2)
	   
            xb (i, mid+1+dj:last+dj:1, 1) =  x (i, mid-1+dj2:1+dj2:-1)
            xb (i, mid+1+dj:last+dj:1, 2) =  y (i, mid-1+dj2:1+dj2:-1)
            xb (i, mid+1+dj:last+dj:1, 3) = -z (i, mid-1+dj2:1+dj2:-1)

	    strip3: if (jmaxb(i,3) /= 0) then
	    
	       dj = last + dj
	       mid = jmaxb(i,3)
	       dj2 = jmaxb(i,1)+jmaxb(i,2)
	       last = 2*mid - 1

	       xb (i, 1+dj:mid+dj, 1) = x (i, 1+dj2:mid+dj2)
	       xb (i, 1+dj:mid+dj, 2) = y (i, 1+dj2:mid+dj2)
	       xb (i, 1+dj:mid+dj, 3) = z (i, 1+dj2:mid+dj2)

               xb (i, mid+1+dj:last+dj:1, 1) =  x (i, mid-1+dj2:1+dj2:-1)
               xb (i, mid+1+dj:last+dj:1, 2) =  y (i, mid-1+dj2:1+dj2:-1)
               xb (i, mid+1+dj:last+dj:1, 3) = -z (i, mid-1+dj2:1+dj2:-1)

	    end if strip3
	 
	 end if strip2
	 	 
      end do
!
      bottom = .true.
      
      call gen (imax,jmaxb,bnodex,body,bedges,bodyin,bottom,kuttab, &
              & caudalb,xb,n,e,p,nmax,nsmall,botle,botre,tople,topre, &
	      & botln,botrn,topln,toprn,oddyb,knode)
!
      write (*,*)
      write (*,*) '...bottom finished.'

!      do i = 1, imax
!         ptl = bedges(tople (i))%panel1
!         ptr = bedges(topre (i))%panel2
!         pbl = bedges(botle (i))%panel1
!         pbr = bedges(botre (i))%panel2
!         bedges(tople (i))%panel2 = ptl
!         bedges(topre (i))%panel1 = ptr
!         bedges(botle (i))%panel2 = pbl
!         bedges(botre (i))%panel1 = pbr
!      end do

 7878 continue      
       
      kutt = 0
! sort out kutta edges
      write (*,*) 'kutta edges.......'
      do edge = 1,e
         if (knode(bedges(edge)%node1).and.knode(bedges(edge)%node2)) then
	    kedge(edge) = .true.
	    kutt = kutt + 1
	    write (*,*) kutt,edge
	 end if
      end do

       
      write (bfile,'(f12.6)') 1.0/float(imax-1)
      write (bfile,*) n
      write (bfile,*) 1
      do 333 node = 1,n
         write (bfile,'(i4,3e16.6)') node,bnodex(node)%x
 333  continue
       
      write (efile,*) e
      do 334 edge = 1,e
         write (efile,'(i8)') edge
	 write (efile,'(2i8)') bedges(edge)%node1,bedges(edge)%node2
	 write (efile,'(2i8)') bedges(edge)%panel1,bedges(edge)%panel2 
 334  continue

      write (pfile,*) p
      do 336 pan = 1,p
         write (pfile,'(i8)') pan
	 write (pfile,'(3i8)') body(pan)%node
	 write (pfile,'(3i8)') body(pan)%edge
	 if (bodyin(pan)) then
	    write (pfile,'(i6)') 1
	 else
	    write (pfile,'(i6)') 0
	 endif
 336  continue

      write(tfile,*) 'title="whale"'
      write(tfile,*) 'variables = x,y,z'
      write(tfile,909)  'zone t=body, n=',n,', e=',p, &
                       & ', f=fepoint, et=triangle'
 909  format(a,i6,a,i6,a)
      do 5 i=1,n
         write(tfile,6) bnodex(i)%x
  6      format(3e16.6)
  5   continue
      do 8 i=1,p
         write(tfile,'(3i6)') body(i)%node
  8   continue

      close(bfile)
      close(efile)
      close(pfile)
      close(tfile)
       
      end
! ---------------------------------------------------------------
! ---------------------------------------------------------------
       subroutine pack (imax,jmax,kutta,caudal,x,y,z,nmax,nsmall,oddy)	     
! ------------------------------
       use quadface; use panel_type; use packblock
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: imax
       integer, dimension(nsmall,3),      intent(inout) :: jmax
       logical, dimension(nsmall,nsmall), intent(inout) :: kutta,caudal
       real,    dimension(nsmall,nsmall), intent(inout) :: x,y,z
       integer,                           intent(in)    :: nmax,nsmall
       logical, dimension(nsmall,nsmall), intent(out)   :: oddy
! ------------------------------
       integer, dimension(nsmall,3)      :: jjmax
       integer, dimension(nsmall)        :: totjmax
       logical, dimension(nsmall,nsmall) :: kkutta,ccaudal
       real,    dimension(nsmall,nsmall) :: xx,yy,zz
       integer                           :: i,ireal,j,jj,jjj,strip
       integer                           :: jlb,jub,ii
       integer, dimension(6)             :: jm
! ---------------------------------------------------------------
!  This redistributes the node orders. It is hardwired for a maximum of
!  three different strips at a given x-location.
!
       kkutta  = .false.
       ccaudal = .false.
       i = 0
       ireal = 0
       totjmax = 0
       jjmax = 0
       xx = 0.
       yy = 0.
       zz = 0.
       
 999   continue
 
       strip = 1
       i = i + 1
       ireal = ireal + 1

!       write (*,*) 'i =',i
!       write (*,*) 'ireal =',ireal
!       write (*,*) 'strip =',strip

       if ((i+5)<=imax) then       
          jm(1:6) = jmax(i:i+5,1)
       else
          jm(1:imax-i+1) = jmax(i:imax,1)
       end if
       
!       write (*,*)
!       write (*,'(a,6i3)') 'jm(1:6) =',jm(1:6)
       
       totjmax (ireal) = totjmax (ireal) + jm(1)
       jjmax (ireal,1) = jm(1)

       xx (ireal,1:jm(1)) = x (i,1:jm(1))
       yy (ireal,1:jm(1)) = y (i,1:jm(1))
       zz (ireal,1:jm(1)) = z (i,1:jm(1))
       kkutta (ireal,1:jm(1)) = kutta (i,1:jm(1))

!       do ii = 1,jm(1)
!         write (*,*) 'xx (ireal,',ii,') = ',xx (ireal,ii)
!	  write (*,*) 'yy (ireal,',ii,') = ',yy (ireal,ii)
!	  write (*,*) 'zz (ireal,',ii,') = ',zz (ireal,ii)
!	  write (*,*)
!       end do 


       if (i==imax) go to 1000

! this next bit of nastiness covers all combinations of strips which are
! possible with up to 3 strips at each x-location

       c: if (caudal(i+1,1)) then
! this strip is caudal; replace nodes to the left of them w/ these values
          call caudalize (i,ireal,jm(2),totjmax,x,y,z,xx,yy,zz,kkutta, &
	                & ccaudal,nsmall)  
          c_s: if (x(i+1,1)==xx(ireal,1)) then
! the next strip is the same x-value (i.e., multiple strips for this ireal)
             call combine (i,ireal,jm(3),totjmax,x,y,z,xx,yy,zz,kutta, &
                         & kkutta,nsmall,jjmax,strip)
             cs_c: if (caudal(i+1,1)) then
	        call caudalize (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz, &
		              & kkutta,ccaudal,nsmall)  
	        csc_s: if (x(i+1,1)==xx(ireal,1)) then
		   call combine (i,ireal,jm(5),totjmax,x,y,z,xx,yy,zz,kutta, &
                               & kkutta,nsmall,jjmax,strip)
                   cscs_c: if (caudal(i+1,1)) then
	              call caudalize (i,ireal,jm(6),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall) 
                   end if cscs_c
		else csc_s
		   csc_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(5),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall) 
                   end if csc_c
                end if csc_s
             else cs_c
	        cs_s: if (x(i+1,1)==xx(ireal,1)) then
		   call combine (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz,kutta, &
                               & kkutta,nsmall,jjmax,strip)
	           css_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(5),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall) 
                   end if css_c
	        end if cs_s
	     end if cs_c
	  else c_s
             c_c: if (caudal(i+1,1)) then
                call caudalize (i,ireal,jm(3),totjmax,x,y,z,xx,yy,zz, &
		              & kkutta,ccaudal,nsmall)
	        cc_s: if (x(i+1,1)==xx(ireal,1)) then
		   call combine (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz,kutta, &
                               & kkutta,nsmall,jjmax,strip)
                   ccs_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(5),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall)
                   end if ccs_c
                else cc_s
		   cc_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall)
                   end if cc_c
                end if cc_s
             end if c_c
          end if c_s
       else c
          s: if (x(i+1,1)==xx(ireal,1).or.x(i+1,2)==xx(ireal,1)) then
             call combine (i,ireal,jm(2),totjmax,x,y,z,xx,yy,zz,kutta, &
                         & kkutta,nsmall,jjmax,strip)
             s_c: if (caudal(i+1,1)) then
                call caudalize (i,ireal,jm(3),totjmax,x,y,z,xx,yy,zz, &
		              & kkutta,ccaudal,nsmall)
                sc_s: if (x(i+1,1)==xx(ireal,1)) then
		   call combine (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz,kutta, &
                               & kkutta,nsmall,jjmax,strip)
		   scs_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(5),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall)
                   end if scs_c
                else sc_s
		   sc_c: if (caudal(i+1,1)) then
		      call caudalize (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall)
		   end if sc_c
                end if sc_s
             else s_c
	        s_s: if (x(i+1,1)==xx(ireal,1)) then
		   call combine (i,ireal,jm(3),totjmax,x,y,z,xx,yy,zz,kutta, &
                               & kkutta,nsmall,jjmax,strip)
	           ss_c: if (caudal(i+1,1)) then
                      call caudalize (i,ireal,jm(4),totjmax,x,y,z,xx,yy,zz, &
		                    & kkutta,ccaudal,nsmall)
	           end if ss_c
                end if s_s
             end if s_c
          end if s
       end if c
      
       if (i==imax) then
          go to 1000
       else	
          go to 999
       end if
       
 1000  continue  
 
 
       imax = ireal
       jmax = jjmax 
       kutta  = kkutta
       caudal = ccaudal
       x = xx
       y = yy
       z = zz
        
       oddy = .false.
       do i = 1,imax
          do j = 1,jmax (i)
	     if (not(caudal(i,j)).and.(z(i,j)==0)) oddy (i,j) = .true.	
          end do
       end do
	
!       write (*,*)
!       write (*,*) 'at end of pack'
!       write (*,*) 'imax =',imax
!       write (*,*) 
!       do i = 1,imax
!         write (*,*) 'jmax(',i,',1) =',jmax(i,1)
!	 if (jmax(i,2)>0) write (*,*) 'jmax(',i,',2) =',jmax(i,2)
!      end do

       
       return
!
       end subroutine pack
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine caudalize (i,ireal,jmaxi,totjmax,x,y,z,xx,yy,zz, &
                           & kkutta,ccaudal,nsmall)
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: i
       integer,                           intent(in)    :: ireal
       integer,                           intent(in)    :: jmaxi
       integer, dimension(nsmall),        intent(in)    :: totjmax
       real,    dimension(nsmall,nsmall), intent(in)    :: x,y,z
       real,    dimension(nsmall,nsmall), intent(inout) :: xx,yy,zz
       logical, dimension(nsmall,nsmall), intent(inout) :: kkutta,ccaudal
       integer,                           intent(in)    :: nsmall  
! ------------------------------
       integer :: j,jj
! ------------------------------
!
       i = i+1
       
       do j = 1,jmaxi
         do jj = 1,totjmax(ireal)

	    if (y(i,j)==yy(ireal,jj)) then
	       xx (ireal,jj) = x (i,j)
	       yy (ireal,jj) = y (i,j)
	       zz (ireal,jj) = z (i,j)
	       ccaudal (ireal,jj) = .true.
	       kkutta (ireal,jj) = kkutta (i,j)
	       
!	       write (*,*)
!	       write (*,*) 'caudalize...'
!	       write (*,*) 'replacing (',ireal,',',jj,') with(',i,',',j 
!	       write (*,*)
            end if
	    
          end do
       end do
       
       return
!
       end subroutine caudalize
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine combine (i,ireal,jmaxi,totjmax,x,y,z,xx,yy,zz,kutta, &
                         & kkutta,nsmall,jjmax,strip)
! ------------------------------
       implicit none
!
       integer,                           intent(inout) :: i
       integer,                           intent(in)    :: ireal
       integer,                           intent(in)    :: jmaxi
       integer, dimension(nsmall),        intent(inout) :: totjmax
       real,    dimension(nsmall,nsmall), intent(in)    :: x,y,z
       real,    dimension(nsmall,nsmall), intent(inout) :: xx,yy,zz
       logical, dimension(nsmall,nsmall), intent(in)    :: kutta
       logical, dimension(nsmall,nsmall), intent(inout) :: kkutta
       integer,                           intent(in)    :: nsmall  
       integer, dimension(nsmall,3),      intent(inout) :: jjmax
       integer,                           intent(inout) :: strip
! ------------------------------
       integer :: jlb,jub,ii
! ------------------------------
!
       i = i + 1
       strip = strip + 1
       
       jlb = totjmax(ireal) + 1
       jub = totjmax(ireal) + jmaxi

       xx (ireal,jlb:jub) = x (i,1:jmaxi)
       yy (ireal,jlb:jub) = y (i,1:jmaxi)
       zz (ireal,jlb:jub) = z (i,1:jmaxi)
       kkutta (ireal,jlb:jub) = kutta (i,1:jmaxi)
       
       totjmax (ireal) = totjmax (ireal) + jmaxi
       jjmax (ireal,strip) = jmaxi

!       write (*,*)
!       write (*,*) 'combine...'
!       write (*,*)
!       write (*,*) 'adding strip ',i,' as strip 2 of ',ireal
!       write (*,*)
       
!       do ii = jlb,jub
!          write (*,*) 'xx (ireal,',ii,') = ',xx (ireal,ii)
!	  write (*,*) 'yy (ireal,',ii,') = ',yy (ireal,ii)
!	  write (*,*) 'zz (ireal,',ii,') = ',zz (ireal,ii)
!	  write (*,*)
!       end do 

       return
!
       end subroutine combine
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine gen (imax,jmax,bnodex,body,bedges,bodyin,bottom,kutta, &
                     & caudal,xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
		     & outln,outrn,inln,inrn,oddy,knode)
! ------------------------------
       use quadface; use panel_type; use genblock; use wraps
! ------------------------------
       implicit none
!       
       integer,                             intent(in)    :: imax
       integer, dimension(nsmall,3),        intent(in)    :: jmax
       logical, dimension(nmax),            intent(inout) :: bodyin
       logical,                             intent(in)    :: bottom
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta,caudal
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(out)   :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       logical, dimension(nsmall,nsmall),   intent(in)    :: oddy
       logical,      dimension(nmax),       intent(inout) :: knode
! ------------------------------
       integer, dimension(nmax) :: newnode,newedge,nextnode,nextedge
       integer, dimension(nmax) :: jstart,jstop
       real,    dimension(3)    :: x1,x2,x3,x4
       integer :: i,j,jend,offset,jj,jmaxl,jmaxr,dum,L,R
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2,strip,diff
       integer :: djl,djr,djl2
       logical :: top,nextnotch,nownotch,notched
! ------------------------------
!
!      note!!!
!
!      certain trivial cases are NOT considered here.
!      there are separate loops for i=1,2,  i=3,4...n-3,n-2, and i=n-1,n
!      it is assumed that the fish is more than 2 or 3 cells long,
!      (i.e. there is no provision for caudal points in the loop for
!      the first cells). 

       write (*,*)
       write (*,*) 'in gen....'
       write (*,*)
       
       top = .true.
       nextnotch = .false.
       nownotch = .false.
       notched = .false.
!      
       jstart = 0
       jstop = 0
       
       newnode  = 0
       newedge  = 0
       nextnode = 0
       nextedge = 0
!
       i = 1
       j = 1
       offset = 0
       djl = 0
       djr = 0

       jmaxl = jmax(i,1)
       jmaxr = jmax(i+1,1) 
       jend  = 2*min(jmaxl,jmaxr) - 1
!
       if (j==jmaxl) offset = jmaxr - jmaxl
       
       wrapone: if (offset>=1) then
!          write (*,*)
!	  write (*,*) 'entering wrapone...'

          call wrap1 (i,j,offset,bnodex,body,bedges,bodyin,xx,n,e,p,     &
	            & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		    & nextnode,nextedge,kutta,knode)
       else wrapone

       L = j

       if (bottom) then
          n  = n + 2
	  n1 = inln (i)
	  n2 = n - 1
	  n3 = inln (i+1)
          e  = e + 4
	  e1 = inle (i)
       else
          n  = n + 4
          n1 = n - 3
	  n2 = n - 2
          n3 = n - 1
	  e  = e + 5
	  e1 = e - 4
       end if
       n4 = n
       e2 = e - 3
       e3 = e - 2
       e4 = e - 1
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
       
 324   format (2(a,i3))
 325   format (3(a,i4))
 326   format (4(a,i4))
 327   format (2(a,i4))
 328   format (5(a,i4))

       x1 = xx (i  , L  , 1:3)
       x2 = xx (i  , L+1, 1:3)
       x3 = xx (i+1, L  , 1:3)
       x4 = xx (i+1, L+1, 1:3)

       if (kutta (i  , L   )) knode (n1) = .true.
       if (kutta (i  , L+1 )) knode (n2) = .true.
       if (kutta (i+1, L   )) knode (n3) = .true.
       if (kutta (i+1, L+1 )) knode (n4) = .true.
       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad....'

       call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                & e1,e2,e3,e4,e5,p1,p2,    &
                & bodyin,bnodex,body,bedges,top)
!
! newnode is used here to store left side nodes for use after wrapping
! we only need to do this for i=1
       newnode (L)   = n1
       newnode (L+1) = n2
       newedge (L)   = e4   
 
       nextnode (L)   = n3
       nextnode (L+1) = n4
       nextedge (L)   = e2
!
       outln (i) = n1
       outln (i+1) = n3
       le (i) = e1       
!       
       outrn (i) = n2
       outrn (i+1) = n4
       re (i) = e3      
!	  
       if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr)) then

          offset = jmaxr - jmaxl
          if (offset>0) then
	     djr =  2*offset
	  else
	     djl = -2*offset
	  end if
	  
!       write (*,*)
!       write (*,*) 'entering wrap....'
          
          call wrap (i,j+1,offset,bnodex,body,bedges,bodyin, &
                   & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
		   & outln,outrn,inln,inrn,nextnode,nextedge,&
		   & newnode,newedge,kutta,knode)
       end if
  
       L = j + djl + 1
       R = j + djr + 1

       if (bottom.and.(j==jend-2)) then
          n1 = inrn (i)
	  n2 = inrn (i+1)
	  e  = e + 3
	  e1 = e - 2
	  e3 = e - 1
	  e4 = inre (i) 
       else if (offset==0) then
          n  = n + 2
          n1 = n - 1
          e  = e + 4
          e1 = e - 3
	  e3 = e - 2
	  e4 = e - 1
       else
          n  = n + 1
	  n1 = newnode (L+1) 
	  e  = e + 3
	  e1 = newedge (L)
	  e3 = e - 2
	  e4 = e - 1
       end if  	  
       n2 = n 
       n3 = outrn (i)
       n4 = outrn (i+1)
       e2 = re (i)
       e5 = e
       p  = p + 2
       p1 = p - 1
       p2 = p
       
!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
       
       x1 = xx (i  , L+1, 1:3)
       x2 = xx (i+1, R+1, 1:3)
       x3 = xx (i  , L  , 1:3)
       x4 = xx (i+1, R  , 1:3)

       if (kutta (i  , L+1 )) knode (n1) = .true.
       if (kutta (i+1, R+1 )) knode (n2) = .true.
       if (kutta (i  , L   )) knode (n3) = .true.
       if (kutta (i+1, R   )) knode (n4) = .true.
       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #2 ....'

       call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                & e1,e2,e3,e4,e5,p1,p2,    &
                & bodyin,bnodex,body,bedges,top)
!
       newnode (L)   = n3
       newnode (L+1) = n1
       newedge (L)   = e1
!
       nextnode (R)   = n4
       nextnode (R+1) = n2
       nextedge (R)   = e3
!
       outrn (i) = n1
       outrn (i+1) = n2
       re (i) = e4
 
       
       i1j3: do j = 3, jend-2, 2
       
          if (((j)==jmaxl.or.(j)==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
!	  
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if

!       write (*,*)
!       write (*,*) 'entering wrap #2 ....'

             call wrap (i,j,offset,bnodex,body,bedges,bodyin,    &
                      & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
		      & outln,outrn,inln,inrn,nextnode,nextedge, &
		      & newnode,newedge,kutta,knode)
		      
	  end if

          L = j + djl
	  R = j + djr

          if (offset==0.and.j<jmaxl) then
	     n  = n + 2
	     n2 = n - 1
             e  = e + 4
	     e2 = e - 3
	     e3 = e - 2
	     e4 = e - 1
          else
	     diff = j+1-jmaxl 
	     n  = n + 1
	     n2 = newnode (jmaxl-diff)
	     e  = e + 3
	     e2 = e - 2
	     e3 = e - 1
	     e4 = newedge (jmaxl-diff) 
	  end if	     
          n1 = outrn (i)
          n3 = outrn (i+1)
	  n4 = n
          e1 = re (i)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!
!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i  , L  , 1:3)
	  x2 = xx (i  , L+1, 1:3)
	  x3 = xx (i+1, R  , 1:3)
	  x4 = xx (i+1, R+1, 1:3)

          if (kutta (i  , L   )) knode (n1) = .true.
	  if (kutta (i  , L+1 )) knode (n2) = .true.
	  if (kutta (i+1, R   )) knode (n3) = .true.
	  if (kutta (i+1, R+1 )) knode (n4) = .true.
	   
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #3....'

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          newnode (L)   = n1
          newnode (L+1) = n2
          newedge (L)   = e4   

          nextnode (R)   = n3
	  nextnode (R+1) = n4
	  nextedge (R)   = e2
!       
          outrn (i) = n2
          outrn (i+1) = n4
          re (i) = e3      
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                        & offset==0) then
!	  
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if
	     
!       write (*,*)
!       write (*,*) 'entering wrap #3....'

             call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                      & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
		      & outln,outrn,inln,inrn,nextnode,nextedge, &
		      & newnode,newedge,kutta,knode)
          end if
 
          L = j + djl + 1
	  R = j + djr + 1

          if (bottom.and.(j==jend-2)) then
             n1 = inrn (i)
	     n2 = inrn (i+1)
	     e  = e + 3
	     e1 = e - 2
	     e3 = e - 1
	     e4 = inre (i)
	  else if (offset==0.and.j+1<jmaxl) then	  
             n  = n + 2
	     n1 = n - 1
	     n2 = n
	     e  = e + 4
	     e1 = e - 3
	     e3 = e - 2
	     e4 = e - 1
	  else
	     diff = j+2-jmaxl
	     n  = n + 1
	     n1 = newnode (jmaxl-diff)
	     n2 = n
	     e  = e + 3
	     e1 = e - 2
	     e3 = e - 1
	     e4 = newedge (jmaxl-diff)  
	  end if
	  n3 = outrn (i)
          n4 = outrn (i+1)
	  e2 = re (i)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i  , L+1, 1:3)
	  x2 = xx (i+1, R+1, 1:3)
	  x3 = xx (i  , L, 1:3)
	  x4 = xx (i+1, R, 1:3)

	  if (kutta (i  , L+1 )) knode (n1) = .true.
	  if (kutta (i+1, R+1 )) knode (n2) = .true.
	  if (kutta (i  , L   )) knode (n3) = .true.
	  if (kutta (i+1, R   )) knode (n4) = .true.
       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #4....'

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          newnode (L)   = n3
          newnode (L+1) = n1
          newedge (L)   = e1

	  nextnode (R)   = n4
	  nextnode (R+1) = n2
	  nextedge (R)   = e3
!
	  outrn (i)   = n1
	  outrn (i+1) = n2
	  re (i) = e4
!       
       end do i1j3
!
       end if wrapone

!       write (*,*)
!       write (*,'(a,i3,a,3i5)') 'after slice ',i,',  n,p,e =',n,p,e
!       write (*,*) 


       j = 1
       newnode = nextnode
       newedge = nextedge
       offset = 0
       djl = 0
       djr = 0
       
       jmaxl = jmax(i+1,1)
       jmaxr = jmax(i+2,1)    
       jend  = 2*min(jmaxl,jmaxr) - 1
!
       wraptwo: if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and. &
                                              & offset==0) then
!
          offset = jmaxr - jmaxl
          if (offset>0) then
             djr =  2*offset
	  else
	     djl = -2*offset
	  end if

!       write (*,*)
!       write (*,*) 'entering wrap #4....'

          call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                   & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                   & outln,outrn,inln,inrn,nextnode,nextedge, &
		   & newnode,newedge,kutta,knode)
			   
       else wraptwo
 
          L = j
	        		       
	  if (bottom) then
             n  = n + 1
	     n1 = inln (i+2)
	     e  = e + 3
	     e4 = inle (i+1) 
	  else
             n  = n + 2
             n1 = n - 1
             e  = e + 4
             e4 = e - 1
	  end  if
	  n2 = newnode (L)
	  n3 = n
	  n4 = newnode (L+1)
	  e1 = e - 3
	  e2 = e - 2
	  e3 = newedge (L)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i+2, L  , 1:3)
	  x2 = xx (i+1, L  , 1:3)
	  x3 = xx (i+2, L+1, 1:3)
	  x4 = xx (i+1, L+1, 1:3)

	  if (kutta (i+2, L   )) knode (n1) = .true.
	  if (kutta (i+1, L   )) knode (n2) = .true.
	  if (kutta (i+2, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L+1 )) knode (n4) = .true.

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #5....'

	  call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
	  nextnode (L)   = n1
	  nextnode (L+1) = n3
	  nextedge (L)   = e1
!
	  outln (i+1) = n2
	  outln (i+2) = n1
	  le (i+1) = e4
!       
	  outrn (i+1) = n4
	  outrn (i+2) = n3
	  re (i+1) = e2      
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and. &
	    & (jmaxl/=jmaxr).and.(offset==0)) then
!
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if

!       write (*,*)
!       write (*,*) 'entering wrap #5....'

	     call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                      & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	              & outln,outrn,inln,inrn,nextnode,nextedge,  &
		      & newnode,newedge,kutta,knode)
          end if

          L = j + djl + 1
	  R = j + djr + 1

          if (bottom.and.(j==jend-2)) then
	     n1 = inrn (i+2)
	     e  = e + 2
	     e1 = inre (i+1)
	  else
	     n  = n + 1
	     n1 = n 
	     e  = e + 3
	     e1 = e - 2
	  end if
	  n2 = outrn (i+2)
	  n3 = newnode (L+1)
	  n4 = outrn (i+1)
	  e2 = newedge (L)
	  e3 = re (i+1)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i+2, R+1, 1:3)
	  x2 = xx (i+2, R  , 1:3)
	  x3 = xx (i+1, L+1, 1:3)
	  x4 = xx (i+1, L  , 1:3)

          if (kutta (i+2, R+1 )) knode (n1) = .true.
	  if (kutta (i+2, R   )) knode (n2) = .true.
	  if (kutta (i+1, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L   )) knode (n4) = .true.

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #6....'

	  call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
	  nextnode (R)   = n2
	  nextnode (R+1) = n1
	  nextedge (R)   = e4
!
	  outrn (i+1) = n3
	  outrn (i+2) = n1
	  re (i+1) = e1
!
       i2j3: do j = 3, jend-2, 2 

          if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
!	  
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if

!       write (*,*)
!       write (*,*) 'entering wrap #6....'

             call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                      & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                      & outln,outrn,inln,inrn,nextnode,nextedge, &
		      & newnode,newedge,kutta,knode)
		      
	  end if
          	  
	  L = j + djl
	  R = j + djr  

	  n  = n + 1
          n1 = outrn (i+2)
          n2 = outrn (i+1)
	  n3 = n
	  n4 = newnode (L+1)
	  e  = e + 3
	  e4 = re (i+1)
	  e1 = e - 2 
	  e2 = e - 1 
	  e3 = newedge (L)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i+2, R  , 1:3)
	  x2 = xx (i+1, L  , 1:3)
	  x3 = xx (i+2, R+1, 1:3)
	  x4 = xx (i+1, L+1, 1:3)

          if (kutta (i+2, R   )) knode (n1) = .true.
	  if (kutta (i+1, L   )) knode (n2) = .true.
	  if (kutta (i+2, R+1 )) knode (n3) = .true.
	  if (kutta (i+1, L+1 )) knode (n4) = .true.
       
!          write (*,*) 'x1 =',x1
!          write (*,*) 'x2 =',x2
!          write (*,*) 'x3 =',x3
!          write (*,*) 'x4 =',x4


!       write (*,*)
!       write (*,*) 'entering quad #7....'


          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
	  nextnode (R)   = n1
	  nextnode (R+1) = n3
	  nextedge (R)   = e1
!       
          outrn (i+1) = n4
          outrn (i+2) = n3
          re (i+1) = e2      
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and.&
	                    & offset==0) then
!
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if

!       write (*,*)
!       write (*,*) 'entering wrap #7....'

             call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                      & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	              & outln,outrn,inln,inrn,nextnode,nextedge,  &
		      & newnode,newedge,kutta,knode)
	  end if

          L = j + djl + 1
	  R = j + djr + 1

          if (bottom.and.(j==jend-2)) then
             n1 = inrn (i+2)
	     e  = e + 2
	     e1 = inre (i+1)
	  else 	  
             n  = n + 1
	     n1 = n 
	     e  = e + 3
	     e1 = e - 2
	  end if
	  n2 = outrn (i+2)
	  n3 = newnode (L+1)
          n4 = outrn (i+1)
	  e2 = newedge (L)
	  e3 = re (i+1)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

	  x1 = xx (i+2, R+1, 1:3)
	  x2 = xx (i+2, R  , 1:3)
	  x3 = xx (i+1, L+1, 1:3)
	  x4 = xx (i+1, L  , 1:3)

          if (kutta (i+2, R+1 )) knode (n1) = .true.
	  if (kutta (i+2, R   )) knode (n2) = .true.
	  if (kutta (i+1, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L   )) knode (n4) = .true.
       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4
!

!       write (*,*)
!       write (*,*) 'entering quad #8....'


          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
	  nextnode (R)   = n2
	  nextnode (R+1) = n1
	  nextedge (R)   = e4
!
	  outrn (i+1) = n3
	  outrn (i+2) = n1
	  re (i+1) = e1
!
       end do i2j3       
!
       end if wraptwo

!       write (*,*)
!       write (*,'(a,i3,a,3i5)') 'after slice ',i+1,',  n,p,e =',n,p,e
!       write (*,*) 

       i3: do i = 3, imax-1, 2

	  j = 1
	  newnode = nextnode
	  newedge = nextedge
	  offset = 0
	  djl = 0
	  djr = 0
!
	  jmaxl = jmax(i,1)
	  jmaxr = jmax(i+1,1)
	  if (nownotch) then
!	     write (*,*) 'busting a notch at i=',i
	     jmaxl = jstart(i) 
	  end if
          jend  = 2*min(jmaxl,jmaxr) - 1 
	  	  
          if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr)) then
	  
             offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if
	     if (i==imax-1.and.jmaxr==1) then
	     
!       write (*,*)
!       write (*,*) 'entering wrapmax....'
	     
	        call wrapmax (i,j,offset,bnodex,body,bedges,bodyin,   &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                            & outln,outrn,inln,inrn,newnode,newedge,  &
			    & kutta,knode)
	     else
	     
!       write (*,*)
!       write (*,*) 'entering wrap #8....'     
	     
	        call wrap (i,j,offset,bnodex,body,bedges,bodyin,    &
                         & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                         & outln,outrn,inln,inrn,nextnode,nextedge, &
			 & newnode,newedge,kutta,knode)
	     end if
	  end if

          L = j + djl
	  R = j + djr           

          if (bottom) then
             n  = n + 1
	     n3 = inln (i+1)
	     e  = e + 3
	     e1 = inle (i)
	  else
             n  = n + 2
             n3 = n - 1
	     e  = e + 4
	     e1 = e - 3
	  end if
	  n1 = newnode (L)
	  n2 = newnode (L+1)
	  n4 = n
	  e2 = e - 2 
	  e3 = e - 1
	  e4 = newedge (L)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!
	  x1 = xx (i  , L  , 1:3)
	  x2 = xx (i  , L+1, 1:3)
	  x3 = xx (i+1, R  , 1:3)
	  x4 = xx (i+1, R+1, 1:3)

          if (kutta (i  , L   )) knode (n1) = .true.
	  if (kutta (i  , L+1 )) knode (n2) = .true.
	  if (kutta (i+1, R   )) knode (n3) = .true.
	  if (kutta (i+1, R+1 )) knode (n4) = .true.
	   
          if ((x4(3)==0.0).and.(j>=jmaxr)) then
	     diff = j+1-jmaxr
	     if(not(bottom)) then
                n  = n - 1
	        n4 = nextnode (jmaxr-diff)
	        if (x3(3)==0.0) then
	           e  = e - 1
		   e2 = nextedge (jmaxr-diff)
		   e3 = e - 1
		   e5 = e
	        end if 
	     else
	        e  = e - 1
		e2 = nextedge (jmaxr-diff)
		e3 = e - 1
		e5 = e
	     end if
	  end if

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*) 
!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #9....'


          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n3
	  nextnode (R+1) = n4
	  nextedge (R)   = e2
!
          outln (i) = n1
          outln (i+1) = n3
          le (i) = e1
!       
          outrn (i) = n2
          outrn (i+1) = n4
          re (i) = e3      
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr)) then
	     
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if
	     if (i==imax-1.and.jmaxr==1) then
	     
!       write (*,*)
!       write (*,*) 'entering wrapmax #2....'
	          
	        call wrapmax (i,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                            & outln,outrn,inln,inrn,newnode,newedge,  &
			    & kutta,knode)
              else

!       write (*,*)
!       write (*,*) 'entering wrap #9....'

	         call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                          & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
	                  & outln,outrn,inln,inrn,nextnode,nextedge, &
			  & newnode,newedge,kutta,knode)
	      end if
	  end if

          L = j + djl + 1
	  R = j + djr + 1

          if (bottom.and.(j==jend-2)) then
	     n2 = inrn (i+1)
	     e  = e + 2
	     e3 = e - 1
	     e4 = inre (i)
	  else
	     n  = n + 1
	     n2 = n 
	     e  = e + 3
	     e3 = e - 2
	     e4 = e - 1
	  end if
	  n1 = newnode (L+1)
	  n3 = outrn (i)
          n4 = outrn (i+1)
	  e1 = newedge (L)
	  e2 = re (i)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!
	  x1 = xx (i  , L+1, 1:3)
	  x2 = xx (i+1, R+1, 1:3)
	  x3 = xx (i  , L  , 1:3)
	  x4 = xx (i+1, R  , 1:3)

          if (kutta (i  , L+1 )) knode (n1) = .true.
	  if (kutta (i+1, R+1 )) knode (n2) = .true.
	  if (kutta (i  , L   )) knode (n3) = .true.
	  if (kutta (i+1, R   )) knode (n4) = .true.

          if ((x2(3)==0.0).and.(j+1>=jmaxr).and. &
	     & not(bottom.and.(j==jend-2))) then
! if we're connecting to the top, then the numbers of nodes and edges
! is already correct	     
             n = n - 1
	     diff = j+2-jmaxr
	     n2 = nextnode (jmaxr-diff)
	     if (x4(3)==0.0) then
		e  = e - 1
		e3 = nextedge (jmaxr-diff)
		e4 = e - 1
		e5 = e
	     end if 
	  end if

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #10....'


          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n4
	  nextnode (R+1) = n2
	  nextedge (R)   = e3
!
          outrn (i)   = n1
	  outrn (i+1) = n2
	  re (i) = e4
!
          i3j3: do j = 3, jend-2, 2

             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	        
		offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr =  2*offset
		else
		   djl = -2*offset
		end if
	        if (i==imax-1.and.jmaxr==1) then
		
!       write (*,*)
!       write (*,*) 'entering wrapmax #3....'
		
	           call wrapmax (i,j,offset,bnodex,body,bedges,bodyin,   &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else
		
!       write (*,*)
!       write (*,*) 'entering wrap #10....'
		
	           call wrap (i,j,offset,bnodex,body,bedges,bodyin,    &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if
!
             if (nownotch.and.(j==jmaxl)) then
	        djl = djl + 2*(jmax(i,1)-jmaxl)
!                write (*,*) 'j =',j
!		write (*,*) 'adding ',djl,' to djl'
	     end if

             L = j + djl
	     R = j + djr

	     n  = n + 1
	     e  = e + 3
             n1 = outrn (i)
	     n2 = newnode (L+1)
             n3 = outrn (i+1)
	     n4 = n
             e1 = re (i)
	     e2 = e - 2
	     e3 = e - 1
	     e4 = newedge (L)
	     e5 = e
             p  = p + 2
	     p1 = p - 1
	     p2 = p

	     x1 = xx (i  , L  , 1:3)
	     x2 = xx (i  , L+1, 1:3)
	     x3 = xx (i+1, R  , 1:3)
	     x4 = xx (i+1, R+1, 1:3)

             if (kutta (i  , L   )) knode (n1) = .true.
	     if (kutta (i  , L+1 )) knode (n2) = .true.
	     if (kutta (i+1, R   )) knode (n3) = .true.
	     if (kutta (i+1, R+1 )) knode (n4) = .true.

             if ((x4(3)==0.0).and.(j>=jmaxr)) then
! we are at the back side of the last row of cells, and want to reuse
! the Kutta nodes and edges     
                n = n - 1
		diff = j+1-jmaxr
		n4 = nextnode (jmaxr-diff)
		if (x3(3)==0.0) then
		   e  = e - 1
		   e2 = nextedge (jmaxr-diff)
		   e3 = e - 1
		   e5 = e
	        end if 
	     end if

	     if (oddy(i+1,j).and.(j/=jmaxr)) nextnotch = .true.
! this is the "notch" between the fin and the body, i.e. the local minimum in
! fish length. we must reach this before any wrapping has occurred, so djr=0    

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #11....'

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)

             nextnode (R)   = n3
	     nextnode (R+1) = n4
	     nextedge (R)   = e2
!       
             outrn (i)   = n2
             outrn (i+1) = n4
             re (i) = e3      
	     
	     if (nextnotch.and.not(notched)) then

!		write (*,*) 'entering notchstuff....'

! we need to remember where the notch starts	      
	        jstart(i+1) = j
! we also need to know where the fin starts (i.e., the notch can consist
! of more than one oddy point
                jstop (i+1) = j
   672          continue
                if (oddy(i+1,jstop(i+1)+1)) then
		   jstop(i+1) = jstop(i+1) + 1
		   go to 672
		end if 
		
	        notched = .true.

	     end if
!
             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then
	     
	        offset = jmaxr - jmaxl
                if (offset>0) then
	           djr =  2*offset
	        else
		   djl = -2*offset
	        end if
	        if (i==imax-1.and.jmaxr==1) then

!       write (*,*)
!       write (*,*) 'entering wrapmax #4....'

	           call wrapmax (i,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
                 else

!       write (*,*)
!       write (*,*) 'entering wrap #11....'

	            call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                             & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
	                     & outln,outrn,inln,inrn,nextnode,nextedge, &
			     & newnode,newedge,kutta,knode)
	         end if
	     end if
!
             if (nownotch.and.(j+1==jmaxl)) then
	        djl = djl + 2*(jmax(i,1)-jmaxl)
!                write (*,*) 'j =',j+1
!		write (*,*) 'adding ',djl,' to djl'
	     end if

             L = j + djl + 1
	     R = j + djr + 1

             if (bottom.and.(j==jend-2)) then
		n2 = inrn (i+1)
	        e  = e + 2
		e3 = e - 1
	        e4 = inre (i)
             else
                n  = n + 1
	        n2 = n 
                e  = e + 3
		e3 = e - 2
		e4 = e - 1
	     end if
             n1 = newnode (L+1) 
	     n3 = outrn (i)
             n4 = outrn (i+1)
             e1 = newedge (L)
	     e2 = re (i)
	     e5 = e
	     p  = p + 2 
	     p1 = p - 1
	     p2 = p
!
	     x1 = xx (i  , L+1, 1:3)
	     x2 = xx (i+1, R+1, 1:3)
	     x3 = xx (i  , L  , 1:3)
	     x4 = xx (i+1, R  , 1:3)

             if (kutta (i  , L+1 )) knode (n1) = .true.
	     if (kutta (i+1, R+1 )) knode (n2) = .true.
	     if (kutta (i  , L   )) knode (n3) = .true.
	     if (kutta (i+1, R   )) knode (n4) = .true.

             if ((x2(3)==0.0).and.(j+1>=jmaxr).and. &
	       & not(bottom.and.(j==jend-2))) then
		diff = j+2-jmaxr
                n  = n - 1
		n2 = nextnode (jmaxr-diff)
		if (x4(3)==0.0) then
		   e  = e - 1
		   e3 = nextedge (jmaxr-diff)
		   e4 = e - 1
		   e5 = e
	        end if 
	     end if

	     if (oddy(i+1,j+djr+1).and.(j+djr+1/=jmaxr)) then
                nextnotch = .true.
! this is a special case of the notch which requires a different quadrilateral	     
!	        call rotate (n1,n2,n3,n4,e1,e2,e3,e4,x1,x2,x3,x4)
             end if
!
!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #12....'

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n4
	     nextnode (R+1) = n2
	     nextedge (R)   = e3
!
             outrn (i) = n1
	     outrn (i+1) = n2
	     re (i) = e4
	     
!             if (bottom.and.i==33.and.j==9) return

!
	     if (nextnotch.and.not(notched)) then

!		write (*,*) 'entering notchstuff #2....'

	        jstart(i+1) = j+1
                jstop (i+1) = j+1
   673          continue
                if (oddy(i+1,jstop(i+1)+1)) then
		   jstop(i+1) = jstop(i+1) + 1
		   go to 673
		end if 

	        notched = .true.

             end if


          end do i3j3


          if (jmax(i+1,2)>0) then
	     if (nownotch) then
		djl2 = jstop(i)-1
		nownotch = .false.
	     else
		djl2 = 2*jmax(i,1)-1
	     end if
	     
!	     write (*,*)
!	     write (*,*) 'entering nextstrip'
!	     write (*,*) 'djl2 =',djl2
!	     write (*,*)
	     
	     call nextstrip (i,imax,jmax,bnodex,body,bedges,bodyin,kutta, &
                           & caudal,xx,n,e,p,nmax,nsmall,oddy,newnode,    &
			   & newedge,nextnode,nextedge,djl2,knode)
	           
	  end if			
	  
	  if (notched) nownotch = .true.
	  nextnotch = .false.
	  notched = .false.	       

!       write (*,*)
!       write (*,'(a,i3,a,3i5)') 'after slice ',i,',  n,p,e =',n,p,e
!       write (*,*) 

          if (i==(imax-1)) exit


          j = 1
	  
	  newnode = nextnode
	  newedge = nextedge
	  offset = 0
	  djl = 0
	  djr = 0
!	  
	  jmaxl = jmax(i+1,1)
	  jmaxr = jmax(i+2,1)
	  if (nownotch) jmaxl = jstart(i+1)
          jend  = 2*min(jmaxl,jmaxr) - 1
	  	  
!	  write (*,*) 'jend =',jend
!	  write (*,*) 'jmaxl =',jmaxl
!	  write (*,*) 'jmaxr =',jmaxr
!	  write (*,*)

          if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	     
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if
	     if (i+1==imax-1.and.jmaxr==1) then
	     
!       write (*,*)
!       write (*,*) 'entering wrapmax #5....'
     
	        call wrapmax (i+1,j,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                            & outln,outrn,inln,inrn,newnode,newedge,  &
			    & kutta,knode)
	     else

!       write (*,*)
!       write (*,*) 'entering wrap #12....'

	        call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                         & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                         & outln,outrn,inln,inrn,nextnode,nextedge, &
			 & newnode,newedge,kutta,knode)
	     end if
	  end if

          L = j + djl
	  R = j + djr

          if (bottom) then
             n  = n + 1
             n1 = inln (i+2)
	     e  = e + 3
	     e4 = inle (i+1) 
          else
             n  = n + 2
             n1 = n - 1
             e  = e + 4
             e4 = e - 1
          end  if
          n2 = newnode (L)
          n3 = n
          n4 = newnode (L+1)
          e1 = e - 3
          e2 = e - 2
          e3 = newedge (L)
          e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!       
	  x1 = xx (i+2, R  , 1:3)
	  x2 = xx (i+1, L  , 1:3)
	  x3 = xx (i+2, R+1, 1:3)
	  x4 = xx (i+1, L+1, 1:3)

          if (kutta (i+2, R   )) knode (n1) = .true.
	  if (kutta (i+1, L   )) knode (n2) = .true.
	  if (kutta (i+2, R+1 )) knode (n3) = .true.
	  if (kutta (i+1, L+1 )) knode (n4) = .true.

          if ((x3(3)==0.0).and.(j>=jmaxr)) then
	     diff = j+1-jmaxr
	     if (not(bottom)) then
        	n  = n - 1
		n3 = nextnode (jmaxr-diff)
		if (x1(3)==0.0) then
		   e  = e - 1
		   e1 = nextedge (jmaxr-diff)
		   e2 = e - 2
		   e4 = e - 1
		   e5 = e 
		end if 
	     else
		e  = e - 1
		e1 = nextedge (jmaxr-diff)
		e2 = e - 1
		e5 = e 
	     end if
	  end if

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #13....'

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
		   
!          write (*,*)
!	  write (*,'(a,3i5)') 'i4j1 .........  n,p,e =',n,p,e
		   
          nextnode (R)   = n1
	  nextnode (R+1) = n3
	  nextedge (R)   = e1
!
          outln (i+1) = n2
	  outln (i+2) = n1
	  le (i+1) = e4
!       
          outrn (i+1) = n4
          outrn (i+2) = n3
          re (i+1) = e2  
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                        & offset==0) then
	  
	     offset = jmaxr - jmaxl
	     if (offset>0) then
		djr =  2*offset
	     else
		djl = -2*offset
	     end if
	     if (i+1==imax-1.and.jmaxr==1) then

!       write (*,*)
!       write (*,*) 'entering wrapmax #6....'

	        call wrapmax (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
                            & outln,outrn,inln,inrn,newnode,newedge,    &
			    & kutta,knode)
             else

!       write (*,*)
!       write (*,*) 'entering wrap #13....'

	        call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                         & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	                 & outln,outrn,inln,inrn,nextnode,nextedge,  &
			 & newnode,newedge,kutta,knode)
	     end if
	  end if

          L = j + djl + 1
	  R = j + djr + 1

          if (bottom.and.(j==jend-2)) then
	     n1 = inrn (i+2)
	     e  = e + 2
	     e1 = inre (i+1)
	  else
	     n  = n + 1
	     n1 = n 
	     e  = e + 3
	     e1 = e - 2
	  end if
	  n2 = outrn (i+2)
	  n3 = newnode (L+1)
          n4 = outrn (i+1)
	  e2 = newedge (L)
	  e3 = re (i+1)
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!       
	  x1 = xx (i+2, R+1, 1:3)
	  x2 = xx (i+2, R  , 1:3)
	  x3 = xx (i+1, L+1, 1:3)
	  x4 = xx (i+1, L  , 1:3)

	  if (kutta (i+2, R+1 )) knode (n1) = .true.
	  if (kutta (i+2, R   )) knode (n2) = .true.
	  if (kutta (i+1, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L   )) knode (n4) = .true.

          if ((x1(3)==0.0).and.(j+1>=jmaxr).and. &
	    & not(bottom.and.(j==jend-2))) then
	     diff = j+2-jmaxr
             n  = n - 1
	     n1 = nextnode (jmaxr-diff)
	     if (x2(3)==0) then
		e  = e - 1
		e4 = nextedge (jmaxr-diff)
		e5 = e
	     end if 
	  end if

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #14....'

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
!          write (*,*)
!	  write (*,'(a,3i5)') 'i4j2 .........  n,p,e =',n,p,e
		   
          nextnode (R)   = n2
	  nextnode (R+1) = n1
	  nextedge (R)   = e4
!
          outrn (i+1) = n3
          outrn (i+2) = n1
          re (i+1) = e1
!
          i4j3: do j = 3, jend-2, 2 
!
             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	     
!	        write (*,*)
!		write (*,*) 'jmaxl =',jmaxl
!		write (*,*) 'jmaxr =',jmaxr
!		write (*,*) 'offset =',offset
     
	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr =  2*offset
		else
		   djl = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then

!       write (*,*)
!       write (*,*) 'entering wrapmax #7....'

	           call wrapmax (i+1,j,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else

!       write (*,*)
!       write (*,*) 'entering wrap #14....'

	           call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             if (nownotch.and.(j==jmaxl)) then
	        djl = djl + 2*(jmax(i+1,1)-jmaxl)
!                write (*,*) 'j =',j
!		write (*,*) 'adding ',djl,' to djl'
	     end if

             L = j + djl
	     R = j + djr

	     n  = n + 1
             n1 = outrn (i+2)
             n2 = outrn (i+1)
	     n3 = n
	     n4 = newnode (L+1)
	     e  = e + 3
	     e1 = e - 2 
	     e2 = e - 1 
	     e3 = newedge (L)
             e4 = re (i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!       
	     x1 = xx (i+2, R  , 1:3)
	     x2 = xx (i+1, L  , 1:3)
	     x3 = xx (i+2, R+1, 1:3)
	     x4 = xx (i+1, L+1, 1:3)

             if (kutta (i+2, R   )) knode (n1) = .true.
	     if (kutta (i+1, L   )) knode (n2) = .true.
	     if (kutta (i+2, R+1 )) knode (n3) = .true.
	     if (kutta (i+1, L+1 )) knode (n4) = .true.

             if ((x3(3)==0.0).and.(j>=jmaxr)) then
		diff = j+1-jmaxr
                n  = n - 1
		n3 = nextnode (jmaxr-diff)
		if (x1(3)==0.0) then
		   e  = e - 1
		   e1 = nextedge (jmaxr-diff)
		   e2 = e - 1
		   e5 = e
	        end if 
	     end if

	     if (oddy(i+2,j+djr).and.(j+djr/=jmaxr)) then
	        nextnotch = .true.
!	        call rotate (n1,n2,n3,n4,e1,e2,e3,e4,x1,x2,x3,x4)
             end if


!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #15....'

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!          write (*,*)
!	  write (*,'(a,3i5)') 'i4j3+ .........  n,p,e =',n,p,e	   
!
             nextnode (R)   = n1
	     nextnode (R+1) = n3
	     nextedge (R)   = e1
!       
             outrn (i+1) = n4
             outrn (i+2) = n3
             re (i+1) = e2      
!
             if (nextnotch.and.not(notched)) then
!		write (*,*)
!		write (*,*) 'entering notchstuff #3....'
!		write (*,*)
		
	        jstart(i+2) = j
                jstop (i+2) = j
   674          continue
                if (oddy(i+2,jstop(i+2)+1)) then
		   jstop(i+2) = jstop(i+2) + 1
		   go to 674
		end if 

	        notched = .true.

	     end if
	     
             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then

	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr =  2*offset
		else
		   djl = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then

!       write (*,*)
!       write (*,*) 'entering wrapmax #8....'

	           call wrapmax (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
                               & outln,outrn,inln,inrn,newnode,newedge,   &
			       & kutta,knode)
                else
		
!       write (*,*)
!       write (*,*) 'entering wrap #15....'
		
	           call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	                    & outln,outrn,inln,inrn,nextnode,nextedge,  & 
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             if (nownotch.and.(j+1==jmaxl)) then
	        djl = djl + 2*(jmax(i+1,1)-jmaxl)
!                write (*,*) 'j =',j+1
!		write (*,*) 'adding ',djl,' to djl'
	     end if

             L = j + djl + 1
	     R = j + djr + 1

             if (bottom.and.(j==jend-2)) then
                n1 = inrn (i+2)
	        e  = e + 2
	        e1 = inre (i+1)
             else
                n  = n + 1
	        n1 = n 
	        e  = e + 3
	        e1 = e - 2
	     end if
	     n2 = outrn (i+2)
             n3 = newnode (L+1)
             n4 = outrn (i+1)
	     e2 = newedge (L)
	     e3 = re (i+1)
	     e4 = e - 1
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!       
	     x1 = xx (i+2, R+1, 1:3)
	     x2 = xx (i+2, R  , 1:3)
	     x3 = xx (i+1, L+1, 1:3)
	     x4 = xx (i+1, L  , 1:3)

             if (kutta (i+2, R+1 )) knode (n1) = .true.
	     if (kutta (i+2, R   )) knode (n2) = .true.
	     if (kutta (i+1, L+1 )) knode (n3) = .true.
	     if (kutta (i+1, L   )) knode (n4) = .true.

             if ((x1(3)==0.0).and.(j+1>=jmaxr).and. &
	       & not(bottom.and.(j==jend-2))) then
		diff = j+2-jmaxr
                n  = n - 1
		n1 = nextnode (jmaxr-diff)
		if (x2(3)==0) then
		   e  = e - 1
		   e4 = nextedge (jmaxr-diff)
		   e5 = e
	        end if 
	     end if
	     
	     if (oddy(i+2,j+djr+1).and.(j+djr+1/=jmaxr)) nextnotch = .true.

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)

!	  write (*,*) 'x1 =',x1
!	  write (*,*) 'x2 =',x2
!	  write (*,*) 'x3 =',x3
!	  write (*,*) 'x4 =',x4

!       write (*,*)
!       write (*,*) 'entering quad #16....'

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!          write (*,*)
!	  write (*,'(a,3i5)') 'i4j4+ .........  n,p,e =',n,p,e
		   
             nextnode (R)   = n2
	     nextnode (R+1) = n1
	     nextedge (R)   = e4
!
             outrn (i+1) = n3
	     outrn (i+2) = n1
	     re (i+1) = e1
!
	     if (nextnotch.and.not(notched)) then
!		write (*,*)
!		write (*,*) 'entering notchstuff #4.....'
!		write (*,*)

	        jstart(i+2) = j+1
                jstop (i+2) = j+1
   675          continue
                if (oddy(i+2,jstop(i+2)+1)) then
		   jstop(i+2) = jstop(i+2) + 1
		   go to 675
		end if 

	        notched = .true.

	     end if
	     
          end do i4j3


          if (jmax(i+2,2)>0) then
	     if (nownotch) then
		djl2 = jstop(i+1)-1
		nownotch = .false.
	     else
		djl2 = 2*jmax(i+1,1)-1
	     end if

	     call nextstrip (i+1,imax,jmax,bnodex,body,bedges,bodyin,kutta, &
	                   & caudal,xx,n,e,p,nmax,nsmall,oddy,newnode,      &
			   & newedge,nextnode,nextedge,djl2,knode) 

	  end if	
	  
	  if (notched) nownotch = .true.
	  nextnotch = .false.
	  notched = .false.
	  			       
!       write (*,*)
!       write (*,'(a,i3,a,3i5)') 'after slice ',i+1,',  n,p,e =',n,p,e
!       write (*,*) 

       end do i3

!       write (*,*)
!       write (*,*) 'n =',n
!       write (*,*) 'p =',p
!       write (*,*) 'e =',e

       return
!
       end subroutine gen
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine wrap1 (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	               & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		       & nextnode,nextedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
! ------------------------------
       real,    dimension(3)    :: x1,x2,x3,x4
       integer :: j,jfar
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2
! ------------------------------
       jfar = jcorn
!
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
       x1 = xx (i  , jcorn  , 1:3)
       x2 = xx (i+1, jcorn  , 1:3)
       x3 = xx (i+1, jcorn+1, 1:3)
       x4 = xx (i+1, jcorn+2, 1:3)

       if (kutta (i  , jcorn   )) knode (n1) = .true.
       if (kutta (i+1, jcorn   )) knode (n2) = .true.
       if (kutta (i+1, jcorn+1 )) knode (n3) = .true.
       if (kutta (i+1, jcorn+2 )) knode (n4) = .true.
       
       call trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                & e1,e2,e3,e4,e5,p1,p2,     &
                & bodyin,bnodex,body,bedges,nmax)
!       
       nextnode (jfar)   = n2
       nextnode (jfar+1) = n3
       nextnode (jfar+2) = n4
       nextedge (jfar)   = e2
       nextedge (jfar+1) = e4
!       
       outln (i)   = n1
       outln (i+1) = n2
       le (i) = e1
!
       outrn (i)   = n1
       outrn (i+1) = n4
       re (i) = e5
!       
       do j = 2,offset
!
          jfar = jcorn + 2*(j-1)
!	  
	  n  = n + 2
	  n1 = outrn (i)
	  n2 = outrn (i+1)
	  n3 = n - 1
	  n4 = n
	  e  = e + 4
	  e1 = re (i)
	  e2 = e - 3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!      
	  x1 = xx (i,   jcorn , 1:3)
	  x2 = xx (i+1, jfar  , 1:3)
	  x3 = xx (i+1, jfar+1, 1:3)
	  x4 = xx (i+1, jfar+2, 1:3)

	  if (kutta (i  , jcorn  )) knode (n1) = .true.
	  if (kutta (i+1, jfar   )) knode (n2) = .true.
	  if (kutta (i+1, jfar+1 )) knode (n3) = .true.
	  if (kutta (i+1, jfar+2 )) knode (n4) = .true.
       
	  call trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                   & e1,e2,e3,e4,e5,p1,p2,     &
                   & bodyin,bnodex,body,bedges,nmax)
!       
	  nextnode (jfar)   = n2
	  nextnode (jfar+1) = n3
	  nextnode (jfar+2) = n4
	  nextedge (jfar)   = e2
	  nextedge (jfar+1) = e4
!       
	  outrn (i)   = n1
	  outrn (i+1) = n4
	  re (i) = e5
!
       end do
!       
       return
!       
       end subroutine wrap1
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine wrap (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	              & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		      & nextnode,nextedge,newnode,newedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       integer,      dimension(nmax),       intent(inout) :: newnode,newedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
! ------------------------------
       real, dimension(3) :: x1,x2,x3,x4
       integer :: j,jfar
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2
! ------------------------------
!
       if (offset>0) then
	  
	  do j = 1,offset
!
	     jfar = jcorn + 2*(j-1)
!
	     n  = n + 2
	     n1 = outrn (i)
	     n2 = outrn (i+1)
	     n3 = n - 1
	     n4 = n
	     e  = e + 4
	     e1 = re (i)
	     e2 = e - 3
	     e3 = e - 2
	     e4 = e - 1
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!      
	     x1 = xx (i,   jcorn , 1:3)
	     x2 = xx (i+1, jfar  , 1:3)
	     x3 = xx (i+1, jfar+1, 1:3)
	     x4 = xx (i+1, jfar+2, 1:3)
       
	     if (kutta (i  , jcorn  )) knode (n1) = .true.
	     if (kutta (i+1, jfar   )) knode (n2) = .true.
	     if (kutta (i+1, jfar+1 )) knode (n3) = .true.
	     if (kutta (i+1, jfar+2 )) knode (n4) = .true.

	     call trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
!       
	     nextnode (jfar)   = n2
	     nextnode (jfar+1) = n3
	     nextnode (jfar+2) = n4
	     nextedge (jfar)   = e2
	     nextedge (jfar+1) = e4
!       
	     outrn (i)   = n1
	     outrn (i+1) = n4
	     re (i) = e5
!
	  end do
!
       else    

	  do j = 1,-offset
!
	     jfar = jcorn + 2*(j-1)
!
	     n1 = outrn (i+1)
	     n2 = outrn (i)
	     n3 = newnode (jfar+1)
	     n4 = newnode (jfar+2)
	     e  = e + 2
	     e1 = re (i)
	     e2 = newedge (jfar)
	     e3 = e - 1
	     e4 = newedge (jfar+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!      
	     x1 = xx (i+1, jcorn , 1:3)
	     x2 = xx (i,   jfar  , 1:3)
	     x3 = xx (i,   jfar+1, 1:3)
	     x4 = xx (i,   jfar+2, 1:3)

       	     if (kutta (i+1, jcorn  )) knode (n1) = .true.
	     if (kutta (i  , jfar   )) knode (n2) = .true.
	     if (kutta (i  , jfar+1 )) knode (n3) = .true.
	     if (kutta (i  , jfar+2 )) knode (n4) = .true.
       
	     call trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
!       
	     outrn (i)   = n4
	     outrn (i+1) = n1
	     re (i) = e5
!
	  end do
!	  
       end if
       
       return
!       
       end subroutine wrap
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine wrapmax (i,jcorn,offset,bnodex,body,bedges,bodyin,xx,n,e,p, &
	                 & nmax,nsmall,le,re,inle,inre,outln,outrn,inln,inrn, &
		         & nextnode,nextedge,kutta,knode)
! ------------------------------
       use triface; use panel_type
! ------------------------------
       implicit none
!
       integer,                             intent(in)    :: i
       integer,                             intent(in)    :: jcorn
       integer,                             intent(in)    :: offset
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nmax),            intent(inout) :: bodyin
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       integer,      dimension(nsmall),     intent(inout) :: le,re
       integer,      dimension(nsmall),     intent(in)    :: inle,inre
       integer,      dimension(nsmall),     intent(inout) :: outln,outrn
       integer,      dimension(nsmall),     intent(in)    :: inln,inrn
       integer,      dimension(nmax),       intent(inout) :: nextnode,nextedge
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta
       logical,      dimension(nmax),       intent(inout) :: knode
! ------------------------------
       real,    dimension(3)    :: x1,x2,x3,x4
       integer :: j,jfar
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2
! ------------------------------
!
       do j = 1,-offset
!
	  jfar = jcorn + 2*(j-1)
!
	  n  = n + 1
	  n1 = n
	  n2 = outrn (i+1)
	  n3 = nextnode (jfar+1)
	  n4 = nextnode (jfar+2)
	  e  = e + 3
	  e1 = e - 2
	  e2 = nextedge (jfar)
	  e3 = e - 1
	  e4 = nextedge (jfar+1)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!      
	  x1 = xx (i+1, jcorn , 1:3)
	  x2 = xx (i  , jfar  , 1:3)
	  x3 = xx (i  , jfar+1, 1:3)
	  x4 = xx (i  , jfar+2, 1:3)

       	  if (kutta (i+1, jcorn  )) knode (n1) = .true.
	  if (kutta (i  , jfar   )) knode (n2) = .true.
	  if (kutta (i  , jfar+1 )) knode (n3) = .true.
	  if (kutta (i  , jfar+2 )) knode (n4) = .true.
       
          call trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                   & e1,e2,e3,e4,e5,p1,p2,     &
                   & bodyin,bnodex,body,bedges,nmax)
!       
          outln (i)   = n2
	  outln (i+1) = n1
	  le (i) = e1
!	     
	  outrn (i)   = n4
	  outrn (i+1) = n1
	  re (i) = e5
!
       end do
!	  
       return
!       
       end subroutine wrapmax
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine trir (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
! ------------------------------
       use panel_type
! ------------------------------
       implicit none
       
       real,    dimension(3), intent(in) :: x1,x2,x3,x4
       integer,               intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2
       integer,               intent(in) :: nmax
       logical,      dimension(nmax), intent(inout) :: bodyin
       type (nodes), dimension(nmax), intent(inout) :: bnodex
       type (panel), dimension(nmax), intent(inout) :: body
       type (edges), dimension(nmax), intent(inout) :: bedges
! ------------------------------
       bnodex(n1)%x = x1
       bnodex(n2)%x = x2
       bnodex(n3)%x = x3
       bnodex(n4)%x = x4
       
       bedges(e1)%node1 = n1
       bedges(e1)%node2 = n2
       bedges(e2)%node1 = n2
       bedges(e2)%node2 = n3
       bedges(e3)%node1 = n3
       bedges(e3)%node2 = n1
       bedges(e4)%node1 = n4
       bedges(e4)%node2 = n3
       bedges(e5)%node1 = n1
       bedges(e5)%node2 = n4
       
       bedges(e1)%panel1 = p1
       bedges(e2)%panel1 = p1
       bedges(e3)%panel2 = p1
       bedges(e3)%panel1 = p2
       bedges(e4)%panel2 = p2
       bedges(e5)%panel2 = p2
       
       body(p1)%node(1) = n1
       body(p1)%node(2) = n2
       body(p1)%node(3) = n3
       body(p1)%edge(1) = e1
       body(p1)%edge(2) = e2
       body(p1)%edge(3) = e3
       
       body(p2)%node(1) = n1
       body(p2)%node(2) = n4
       body(p2)%node(3) = n3
       body(p2)%edge(1) = e5
       body(p2)%edge(2) = e4
       body(p2)%edge(3) = e3
       
       bodyin(p1) = .true.
       bodyin(p2) = .false.
!
       return
!       
       end subroutine trir
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine tril (x1,x2,x3,x4,n1,n2,n3,n4,  &
                      & e1,e2,e3,e4,e5,p1,p2,     &
                      & bodyin,bnodex,body,bedges,nmax)
! ------------------------------
       use panel_type
! ------------------------------
       implicit none
       
       real,    dimension(3), intent(in) :: x1,x2,x3,x4
       integer,               intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2
       integer,               intent(in) :: nmax
       logical,      dimension(nmax), intent(inout) :: bodyin
       type (nodes), dimension(nmax), intent(inout) :: bnodex
       type (panel), dimension(nmax), intent(inout) :: body
       type (edges), dimension(nmax), intent(inout) :: bedges
! ------------------------------
       bnodex(n1)%x = x1
       bnodex(n2)%x = x2
       bnodex(n3)%x = x3
       bnodex(n4)%x = x4
       
       bedges(e1)%node1 = n2
       bedges(e1)%node2 = n1
       bedges(e2)%node1 = n3
       bedges(e2)%node2 = n2
       bedges(e3)%node1 = n1
       bedges(e3)%node2 = n3
       bedges(e4)%node1 = n3
       bedges(e4)%node2 = n4
       bedges(e5)%node1 = n4
       bedges(e5)%node2 = n1
       
       bedges(e1)%panel1 = p1
       bedges(e2)%panel1 = p1
       bedges(e3)%panel2 = p1
       bedges(e3)%panel1 = p2
       bedges(e4)%panel2 = p2
       bedges(e5)%panel2 = p2
       
       body(p1)%node(1) = n1
       body(p1)%node(2) = n3
       body(p1)%node(3) = n2
       body(p1)%edge(1) = e1
       body(p1)%edge(2) = e2
       body(p1)%edge(3) = e3
       
       body(p2)%node(1) = n1
       body(p2)%node(2) = n3
       body(p2)%node(3) = n4
       body(p2)%edge(1) = e5
       body(p2)%edge(2) = e4
       body(p2)%edge(3) = e3
       
       bodyin(p1) = .false.
       bodyin(p2) = .true.
!
       return
!       
       end subroutine tril
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine rotate (n1,n2,n3,n4,e1,e2,e3,e4,x1,x2,x3,x4)
! ------------------------------
       implicit none
       
       integer,            intent(inout) :: n1,n2,n3,n4,e1,e2,e3,e4
       real, dimension(3), intent(inout) :: x1,x2,x3,x4
       
       integer            :: dumn1,dumn2,dumn3,dumn4,dume1,dume2,dume3,dume4
       real, dimension(3) :: dumx1,dumx2,dumx3,dumx4
! ------------------------------
! this subroutine rotates the orientation of a quadrilateral   
   
       dumn1 = n1
       dumn2 = n2
       dumn3 = n3
       dumn4 = n4
       dume1 = e1
       dume2 = e2
       dume3 = e3
       dume4 = e4
       dumx1 = x1
       dumx2 = x2
       dumx3 = x3
       dumx4 = x4
       
       n1 = dumn2
       n2 = dumn4
       n3 = dumn1 
       n4 = dumn3
       e1 = dume4
       e2 = dume1
       e3 = dume2
       e4 = dume3
       x1 = dumx2
       x2 = dumx4
       x3 = dumx1
       x4 = dumx3
 
       return
       
       end subroutine rotate
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine nextstrip (ireal,imax,jmax,bnodex,body,bedges,bodyin,  &
                           & kutta,caudal,xx,n,e,p,nmax,nsmall,oddy,     &
			   & newnode,newedge,nextnode,nextedge,djl,knode)
! ------------------------------
       use quadface; use panel_type; use wraps
! ------------------------------
       implicit none
!       
       integer,                             intent(in)    :: ireal,imax
       integer, dimension(nsmall,3),        intent(in)    :: jmax
       logical, dimension(nmax),            intent(inout) :: bodyin
       logical, dimension(nsmall,nsmall),   intent(in)    :: kutta,caudal
       real,    dimension(nsmall,nsmall,3), intent(in)    :: xx
       integer,                             intent(inout) :: n,e,p
       integer,                             intent(in)    :: nmax,nsmall
       type (nodes), dimension(nmax),       intent(inout) :: bnodex
       type (panel), dimension(nmax),       intent(inout) :: body
       type (edges), dimension(nmax),       intent(inout) :: bedges
       logical, dimension(nsmall,nsmall),   intent(in)    :: oddy
       integer, dimension(nmax),            intent(inout) :: newnode,newedge
       integer, dimension(nmax),            intent(inout) :: nextnode,nextedge
       integer,                             intent(in)    :: djl
       logical,      dimension(nmax),       intent(inout) :: knode
! ------------------------------
       integer, dimension(nsmall) :: le,re
       integer, dimension(nsmall) :: inle,inre
       integer, dimension(nsmall) :: outln,outrn
       integer, dimension(nsmall) :: inln,inrn
       real,    dimension(3)    :: x1,x2,x3,x4,xdum
       integer :: i,j,jend,offset,jj,jmaxl,jmaxr,L,R
       integer :: e1,e2,e3,e4,e5,n1,n2,n3,n4,p1,p2,diff,djr
       integer :: djl2,djr2,drl,jend1,jend2,dum,jmaxlo,jhi
       logical :: top
! ------------------------------

!       write (*,*)
!       write (*,*) 'in nextstrip!!!'
!       write (*,*)
!       write (*,*) 'ireal =',ireal
!       write (*,*) 'djl   =',djl

!	   write (*,*)
!           write (*,*) 'newnodes...'
!           do dum = 1,24
!	   write (*,*) dum,newnode (dum)
!	   end do
!           write (*,*) 'nextnodes...'
!           do dum = 1,24
!	   write (*,*) dum,nextnode (dum)
!	   end do


       top = .true.

       offset = 0
       
       djl2 = 0
       djr2 = 0
!
       oddi: if (mod(ireal,2)==1) then
! this strip is located at an odd i-value        
       
          i = ireal
	  
          oddioddj: if (mod(djl,2)==0) then
! this strip starts at an odd j-value	  

          j = djl + 1

! compare sizes of strips (using same offset for left and right)
          if (djl>=jmax(i,1)) then
	     jmaxl = jmax(i,2) + djl
	  else
	     jmaxl = jmax(i,1)
	  end if
	  jmaxr = jmax(i+1,2) + djl
          jend1 = 2*(jmaxl-djl) - 1 + djl
	  jend2 = 2*(jmaxr-djl) - 1 + djl
	  jend  = min(jend1,jend2)
	  jmaxlo = min(jmaxl,jmaxr)
! find the offset for the right side 
          djr = 2*jmax(i+1,1) - (djl+1)

          L = j
	  R = j + djr

          n  = n + 2
	  n1 = newnode (L)
	  n2 = newnode (L+1)
!	  write (*,*)
!	  write (*,*) 'pulling ',n3,' from newnode ...',j
!	  write (*,*) 'pulling ',n2,' from newnode ...',j
          n3 = n - 1
	  n4 = n
	  e  = e + 4
	  e1 = e - 3
	  e2 = e - 2 
	  e3 = e - 1
	  e4 = newedge (L)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!
	  x1 = xx (i  , L  , 1:3)
	  x2 = xx (i  , L+1, 1:3)
	  x3 = xx (i+1, R  , 1:3)
	  x4 = xx (i+1, R+1, 1:3)

          if (kutta (i  , L   )) knode (n1) = .true.
	  if (kutta (i  , L+1 )) knode (n2) = .true.
	  if (kutta (i+1, R   )) knode (n3) = .true.
	  if (kutta (i+1, R+1 )) knode (n4) = .true.
	   
!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
       
 324   format (2(a,i3))
 325   format (3(a,i4))
 326   format (4(a,i4))
 327   format (2(a,i4))
 328   format (5(a,i4))

!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4
	  
          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
		   
          nextnode (R)   = n3
	  nextnode (R+1) = n4
	  nextedge (R)   = e2

!	  write (*,*)
!	  write (*,*) 'putting ',n3,' into nextnode ...',j+djr
!	  write (*,*) 'putting ',n2,' into nextnode ...',j+djr+1
!
          outln (i) = n1
          outln (i+1) = n3
          le (i) = e1
!       
          outrn (i) = n2
          outrn (i+1) = n4
          re (i) = e3      

          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr)) then
	     
	     offset = jmaxr - jmaxl
             if (offset>0) then
		djr2 =  2*offset
	     else
		djl2 = -2*offset
	     end if
	     if (i==imax-1.and.jmaxr==1) then
	        call wrapmax (i,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                            & outln,outrn,inln,inrn,newnode,newedge,  &
			    & kutta,knode)
              else
	         call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                          & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
	                  & outln,outrn,inln,inrn,nextnode,nextedge, &
			  & newnode,newedge,kutta,knode)
	      end if
	  end if

          L = j + djl2 + 1
	  R = j + djr2 + djr + 1

	  if (j+2==jend) then
	     n2 = outln (i+1)
	     e  = e + 2
	     e3 = e - 1
	     e4 = le (i+1)
	  else	     
             n  = n + 1
	     n2 = n 
             e  = e + 3
	     e3 = e - 2
	     e4 = e - 1
	  end if
	  n1 = newnode (L+1) 
!	  write (*,*)
!	  write (*,*) 'pulling ',n1,' from newnode ...',l+1
          n3 = outrn (i)
          n4 = outrn (i+1)
          e1 = newedge (L)
	  e2 = re (i)
	  e5 = e
	  p  = p + 2 
	  p1 = p - 1
	  p2 = p

	  x1 = xx (i  , L+1, 1:3)
	  x2 = xx (i+1, R+1, 1:3)
	  x3 = xx (i  , L  , 1:3)
	  x4 = xx (i+1, R  , 1:3)

          if (kutta (i  , L+1 )) knode (n1) = .true.
	  if (kutta (i+1, R+1 )) knode (n2) = .true.
	  if (kutta (i  , L   )) knode (n3) = .true.
	  if (kutta (i+1, R   )) knode (n4) = .true.

          if ((x2(3)==0.0).and.(j+1>=jmaxlo)) then
	     jhi = 2*jmax(i+1,1) - 1 + jmax(i+1,2)
	     diff = R - jhi + 1
	     if (j+2/=jend) then
        	n = n - 1
		n2 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n2,' from nextnode ...',jhi-diff
!	  write (*,*)
!	   write (*,*)
!           write (*,*) 'newnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,newnode (dum)
!	   end do
!           write (*,*) 'nextnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,nextnode (dum)
!	   end do
	  
		if (x4(3)==0.0) then
		   e  = e - 1
		   e3 = nextedge (jhi-diff)
		   e4 = e - 1
		   e5 = e
		end if 
	     else
	        e = e - 1
		e3 = nextedge (jhi-diff)
	        e5 = e
	     end if
	  end if

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4


          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n4
	  nextnode (R+1) = n2
	  nextedge (R)   = e3
!	  write (*,*)
!	  write (*,*) 'putting ',n4,' into nextnode ...',j+djr2+djr+1
!	  write (*,*) 'putting ',n2,' into nextnode ...',j+djr2+djr+2
!
          outrn (i)   = n1
	  outrn (i+1) = n2
	  re (i) = e4
	  
          do j = 3+djl, jend-2, 2
!	  
             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	        
		offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i==imax-1.and.jmaxr==1) then
	           call wrapmax (i,j,offset,bnodex,body,bedges,bodyin,   &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else
!		   write (*,*) 'entering wrap (#1 of next)'
	           call wrap (i,j,offset,bnodex,body,bedges,bodyin,    &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2
	     R = j + djr2 + djr

	     n  = n + 1
             n1 = outrn (i)
	     n2 = newnode (L+1)
!	  write (*,*)
!	  write (*,*) 'pulling ',n2,' from newnode ...',L+1
             n3 = outrn (i+1)
	     n4 = n
	     
	     e  = e + 3
             e1 = re (i)
	     e2 = e - 2
	     e3 = e - 1
	     e4 = newedge (L)
	     e5 = e
             p  = p + 2
	     p1 = p - 1
	     p2 = p

	     x1 = xx (i  , L  , 1:3)
	     x2 = xx (i  , L+1, 1:3)
	     x3 = xx (i+1, R  , 1:3)
	     x4 = xx (i+1, R+1, 1:3)

             if (kutta (i  , L   )) knode (n1) = .true.
	     if (kutta (i  , L+1 )) knode (n2) = .true.
	     if (kutta (i+1, R   )) knode (n3) = .true.
	     if (kutta (i+1, R+1 )) knode (n4) = .true.

             if ((x4(3)==0.0).and.(j>=jmaxlo)) then
                jhi = 2*jmax(i+1,1) - 1 + jmax(i+1,2)
		diff = R - jhi + 1
                n  = n - 1
		n4 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n4,' from nextnode ...',jhi-diff
!	  write (*,*)
!	   write (*,*)
!           write (*,*) 'newnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,newnode (dum)
!	   end do
!           write (*,*) 'nextnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,nextnode (dum)
!	   end do
		if (x3(3)==0.0) then
		   e = e - 1
		   e2 = nextedge (jhi-diff)
		   e3 = e - 1
		   e5 = e
	        end if 
	     end if
	     
!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4


             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n3
	     nextnode (R+1) = n4
	     nextedge (R)   = e2
!	  write (*,*)
!	  write (*,*) 'putting ',n3,' into nextnode ...',j+djr2+djr
!	  write (*,*) 'putting ',n2,' into nextnode ...',j+djr2+djr+1
!       
             outrn (i)   = n2
             outrn (i+1) = n4
             re (i) = e3      
	     
             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then
	     
	        offset = jmaxr - jmaxl
                if (offset>0) then
	           djr2 =  2*offset
	        else
		   djl2 = -2*offset
	        end if
	        if (i==imax-1.and.jmaxr==1) then
	           call wrapmax (i,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
                else
	            call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                             & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
	                     & outln,outrn,inln,inrn,nextnode,nextedge, &
			     & newnode,newedge,kutta,knode)
!       write (*,*)
!       write (*,324) ' i =',i
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
			     
			     
	         end if
	     end if

             L = j + djl2 + 1
	     R = j + djr2 + djr + 1
	     
	     if (j+2==jend) then
	        n2 = outln (i+1)
		e  = e + 2
		e3 = e - 1
		e4 = le (i)
	     else	     
        	n  = n + 1
		n2 = n 
        	e  = e + 3
		e3 = e - 2
		e4 = e - 1
	     end if
	     n1 = newnode (L+1) 
!	  write (*,*)
!	  write (*,*) 'pulling ',n1,' from newnode ...',L+1
             n3 = outrn (i)
             n4 = outrn (i+1)
             e1 = newedge (L)
	     e2 = re (i)
	     e5 = e
	     p  = p + 2 
	     p1 = p - 1
	     p2 = p
	     
	     x1 = xx (i  , L+1, 1:3)
	     x2 = xx (i+1, R+1, 1:3)
	     x3 = xx (i  , L  , 1:3)
	     x4 = xx (i+1, R  , 1:3)

             if (kutta (i  , L+1 )) knode (n1) = .true.
	     if (kutta (i+1, R+1 )) knode (n2) = .true.
	     if (kutta (i  , L   )) knode (n3) = .true.
	     if (kutta (i+1, R   )) knode (n4) = .true.

             if ((x2(3)==0.0).and.(j+1>=jmaxr)) then
                jhi = 2*jmax(i+1,1) - 1 + jmax(i+1,2)
		diff = R - jhi + 1
	        if (j+2/=jend) then
        	   n = n - 1
		   n2 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n2,' from nextnode ...',jhi-diff
!	  write (*,*)
!	   write (*,*)
!           write (*,*) 'newnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,newnode (dum)
!	   end do
!           write (*,*) 'nextnodes...'
!           do dum = 1,jhi+1
!	   write (*,*) dum,nextnode (dum)
!	   end do
		   if (x4(3)==0.0) then
		      e  = e - 1
		      e3 = nextedge (jhi-diff)
		      e4 = e - 1
		      e5 = e
		   end if 
	        else
		   e  = e - 1
		   e3 = nextedge (jhi-diff)
		   e5 = e
		end if
	     end if

!       write (*,*)
!       write (*,324) ' i =',i,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4


             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n4
	     nextnode (R+1) = n2
	     nextedge (R)   = e3
!	  write (*,*)
!	  write (*,*) 'putting ',n4,' into nextnode ...',R
!	  write (*,*) 'putting ',n2,' into nextnode ...',R+1
!
             outrn (i)   = n1
	     outrn (i+1) = n2
	     re (i) = e4

          end do

       else oddioddj
! we are at an even j-value. in order to use the same loops as we did in 
! gen (which are based on odd values of j incremented by 2), we set j to
! one less than its true value (cf above)

          j = djl       

! compare sizes of strips (using same offset for left and right)
          if (djl>=jmax(i,1)) then
	     jmaxl = jmax(i,2) + djl
	  else
	     jmaxl = jmax(i,1)
	  end if
	  jmaxr = jmax(i+1,2) + djl
          jend1 = 2*(jmaxl-djl) - 1 + djl
	  jend2 = 2*(jmaxr-djl) - 1 + djl
	  jend  = min(jend1,jend2)
	  jmaxlo = min(jmaxl,jmaxr)
! now, shift right offset appropriately
          djr = 2*jmax(i+1,1) - (djl+1)

          L = j + 1
	  R = j + djr + 1

	  n  = n + 2
	  n1 = newnode (L+1)
	  n2 = n - 1
	  n3 = outrn (i)
          n4 = n
	  e  = e + 4
	  e1 = newedge (L)
	  e2 = e - 3
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

	  x1 = xx (i  , L+1, 1:3)
	  x2 = xx (i+1, R+1, 1:3)
	  x3 = xx (i  , L  , 1:3)
	  x4 = xx (i+1, R  , 1:3)

          if (kutta (i  , L+1 )) knode (n1) = .true.
	  if (kutta (i+1, R+1 )) knode (n2) = .true.
	  if (kutta (i  , L   )) knode (n3) = .true.
	  if (kutta (i+1, R   )) knode (n4) = .true.

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n4
	  nextnode (R+1) = n2
	  nextedge (R)   = e3
!
          outrn (i)   = n1
	  outrn (i+1) = n2
	  re (i) = e4
!
          do j = 2+djl, jend, 2
!	  
             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	        
		offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i==imax-1.and.jmaxr==1) then
	           call wrapmax (i,j,offset,bnodex,body,bedges,bodyin,   &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else
	           call wrap (i,j,offset,bnodex,body,bedges,bodyin,    &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2
	     R = j + djr2 + djr

             if (j+1==jend) then
	        n4 = outln (i+1)
		e  = e + 2
		e2 = e - 1
		e3 = le (i)
	     else
        	n  = n + 1 
		n4 = n
		e  = e + 3
		e2 = e - 2
		e3 = e - 1
	     end if 
             n1 = outrn (i)
	     n2 = newnode (L+1)
             n3 = outrn (i+1)
             e1 = re (i)
	     e4 = newedge (L)
	     e5 = e
             p  = p + 2
	     p1 = p - 1
	     p2 = p
	     
	     x1 = xx (i  , L  , 1:3)
	     x2 = xx (i  , L+1, 1:3)
	     x3 = xx (i+1, R  , 1:3)
	     x4 = xx (i+1, R+1, 1:3)

             if (kutta (i  , L   )) knode (n1) = .true.
	     if (kutta (i  , L+1 )) knode (n2) = .true.
	     if (kutta (i+1, R   )) knode (n3) = .true.
	     if (kutta (i+1, R+1 )) knode (n4) = .true.

             if ((x4(3)==0.0).and.(j>=jmaxr)) then
                jhi = 2*jmax(i+1,1) - 1 + jmax(i+1,2)
		diff = R - jhi + 1
	        if (j+1/=jend) then
                   n = n - 1
		   n4 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n4,' from nextnode ...',jhi-diff
!	  write (*,*)
		   if (x3(3)==0.0) then
		      e  = e - 1
		      e2 = nextedge (jhi-diff)
		      e3 = e - 1
		      e5 = e
	           end if 
		else
		   e  = e - 1
                   e2 = nextedge (jhi-diff)
		   e5 = e
	        end if
	     end if


             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n3
	     nextnode (R+1) = n4
	     nextedge (R)   = e2
!       
             outrn (i)   = n2
             outrn (i+1) = n4
             re (i) = e3      
	     
	     if (j==jend-1) exit
	     
             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then
	     
	        offset = jmaxr - jmaxl
                if (offset>0) then
	           djr2 =  2*offset
	        else
		   djl2 = -2*offset
	        end if
	        if (i==imax-1.and.jmaxr==1) then
	           call wrapmax (i,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
                 else
	            call wrap (i,j+1,offset,bnodex,body,bedges,bodyin,  &
                             & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
	                     & outln,outrn,inln,inrn,nextnode,nextedge, &
			     & newnode,newedge,kutta,knode)
	         end if
	     end if
	     
	     L = j + djl2 + 1
	     R = j + djr2 + djr + 1
	  
             n  = n + 1
	     n2 = n 
	     n1 = newnode (L+1) 
             n3 = outrn (i)
             n4 = outrn (i+1)
             e  = e + 3
             e1 = newedge (L)
	     e2 = re (i)
	     e3 = e - 2
	     e4 = e - 1
	     e5 = e
	     p  = p + 2 
	     p1 = p - 1
	     p2 = p
!
	     x1 = xx (i  , L+1, 1:3)
	     x2 = xx (i+1, R+1, 1:3)
	     x3 = xx (i  , L  , 1:3)
	     x4 = xx (i+1, R  , 1:3)

             if (kutta (i  , L+1 )) knode (n1) = .true.
	     if (kutta (i+1, R+1 )) knode (n2) = .true.
	     if (kutta (i  , L   )) knode (n3) = .true.
	     if (kutta (i+1, R   )) knode (n4) = .true.

             if ((x2(3)==0.0).and.(j>=jmaxlo)) then
                jhi = 2*jmax(i+1,1) - 1 + jmax(i+1,2)
		diff = R - jhi + 1 
                n = n - 1
		n2 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n2,' from nextnode ...',jhi-diff
!	  write (*,*)
		if (x4(3)==0.0) then
		   e  = e - 1
		   e3 = nextedge (jhi-diff)
		   e4 = e - 1
		   e5 = e
	        end if 
	     end if

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n4
	     nextnode (R+1) = n2
	     nextedge (R)   = e3
!
             outrn (i)   = n1
	     outrn (i+1) = n2
	     re (i) = e4

          end do
          
	  end if oddioddj
	  
       else oddi
! this strip is located at an even i-value. i is shifted by one (see comment
! at oddioddj loop above   
	  
          i = ireal - 1
	  
	  evenioddj: if (mod(djl,2)==0) then
! this strip starts at an odd j-value	  

          j = 1 + djl 
       
! compare sizes of strips (using same offset for left and right)
          if (djl>=jmax(i+1,1)) then
	     jmaxl = jmax(i+1,2) + djl
	  else
	     jmaxl = jmax(i+1,1)
	  end if
	  jmaxr = jmax(i+2,2) + djl
          jend1 = 2*(jmaxl-djl) - 1 + djl
	  jend2 = 2*(jmaxr-djl) - 1 + djl
	  jend  = min(jend1,jend2)
	  jmaxlo = min(jmaxl,jmaxr)
! now, shift right offset appropriately
          djr = 2*jmax(i+2,1) - (djl+1)

          L = j
	  R = j + djr

          n  = n + 2
          n1 = n - 1
          n2 = newnode (L)
          n3 = n
          n4 = newnode (L+1)
          e  = e + 4
          e1 = e - 3
          e2 = e - 2
          e3 = newedge (L)
          e4 = e - 1
          e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!       
	  x1 = xx (i+2, R  , 1:3)
	  x2 = xx (i+1, L  , 1:3)
	  x3 = xx (i+2, R+1, 1:3)
	  x4 = xx (i+1, L+1, 1:3)

          if (kutta (i+2, R   )) knode (n1) = .true.
	  if (kutta (i+1, L   )) knode (n2) = .true.
	  if (kutta (i+2, R+1 )) knode (n3) = .true.
	  if (kutta (i+1, L+1 )) knode (n4) = .true.

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n1
	  nextnode (R+1) = n3
	  nextedge (R)   = e1
!
          outln (i+1) = n2
	  outln (i+2) = n1
	  le (i+1) = e4
!       
          outrn (i+1) = n4
          outrn (i+2) = n3
          re (i+1) = e2  
!
          if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                        & offset==0) then
	  
	     offset = jmaxr - jmaxl
	     if (offset>0) then
		djr2 =  2*offset
	     else
		djl2 = -2*offset
	     end if
	     if (i+1==imax-1.and.jmaxr==1) then
	        call wrapmax (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
                            & outln,outrn,inln,inrn,newnode,newedge,    &
			    & kutta,knode)
             else
	        call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                         & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	                 & outln,outrn,inln,inrn,nextnode,nextedge,  &
			 & newnode,newedge,kutta,knode)
	     end if
	  end if

          L = j + djl2 + 1
	  R = j + djr2 + djr + 1

          if (j+2==jend) then
	     n1 = outln (i+2)
	     e  = e + 2
	     e1 = le (i+2)
             e4 = e - 1
	  else
             n  = n + 1
	     n1 = n 
	     e  = e + 3
	     e1 = e - 2
	     e4 = e - 1
	  end if
	  n2 = outrn (i+2)
          n3 = newnode (L+1)
          n4 = outrn (i+1)
	  e2 = newedge (L)
	  e3 = re (i+1)
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p

	  x1 = xx (i+2, R+1, 1:3)
	  x2 = xx (i+2, R  , 1:3)
	  x3 = xx (i+1, L+1, 1:3)
	  x4 = xx (i+1, L  , 1:3)

          if (kutta (i+2, R+1 )) knode (n1) = .true.
	  if (kutta (i+2, R   )) knode (n2) = .true.
	  if (kutta (i+1, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L   )) knode (n4) = .true.

          if ((x1(3)==0.0).and.(j+1>=jmaxlo)) then
             jhi = 2*jmax(i+2,1) - 1 + jmax(i+2,2)
	     diff = R - jhi + 1
	     if (j+2/=jend) then
        	n = n - 1
		n1 = nextnode (jhi-diff)
		if (x2(3)==0.0) then
		   e  = e - 1
		   e4 = nextedge (jhi-diff)
		   e5 = e
		end if 
	     else
	        e  = e - 1
		e4 = nextedge (jhi-diff)
		e5 = e
	     end if
	  end if

          call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                   & e1,e2,e3,e4,e5,p1,p2,    &
                   & bodyin,bnodex,body,bedges,top)
!

          nextnode (R)   = n2
	  nextnode (R+1) = n1
	  nextedge (R)   = e4
!
          outrn (i+1) = n3
          outrn (i+2) = n1
          re (i+1) = e1
!
          do j = 3+djl, jend-2, 2 
!
             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	     
	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then
	           call wrapmax (i+1,j,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else
	           call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2
	     R = j + djr2 + djr

	     n  = n + 1
             n1 = outrn (i+2)
             n2 = outrn (i+1)
	     n3 = n
	     n4 = newnode (L+1)
	     e  = e + 3
	     e1 = e - 2 
	     e2 = e - 1 
	     e3 = newedge (L)
             e4 = re (i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!       
	     x1 = xx (i+2, R  , 1:3)
	     x2 = xx (i+1, L  , 1:3)
	     x3 = xx (i+2, R+1, 1:3)
	     x4 = xx (i+1, L+1, 1:3)

             if (kutta (i+2, R   )) knode (n1) = .true.
	     if (kutta (i+1, L   )) knode (n2) = .true.
	     if (kutta (i+2, R+1 )) knode (n3) = .true.
	     if (kutta (i+1, L+1 )) knode (n4) = .true.

             if ((x3(3)==0.0).and.(j>=jmaxlo)) then
                jhi = 2*jmax(i+2,1) - 1 + jmax(i+2,2)
                diff = R - jhi + 1
                n = n - 1
		n3 = nextnode (jhi-diff)
		if (x1(3)==0.0) then
		   e  = e - 1
		   e1 = nextedge (jhi-diff)
		   e2 = e - 1
		   e5 = e
	        end if 
	     end if

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n1
	     nextnode (R+1) = n3
	     nextedge (R)   = e1
!       
             outrn (i+1) = n4
             outrn (i+2) = n3
             re (i+1) = e2      
	     
             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then

	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then
	           call wrapmax (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
                               & outln,outrn,inln,inrn,newnode,newedge,    &
			       & kutta,knode)
                else
	           call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	                    & outln,outrn,inln,inrn,nextnode,nextedge,  &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2 + 1
	     R = j + djr2 + djr + 1

             if (j+2==jend) then
	        n1 = outln (i+2)
		e  = e + 2
		e1 = le (i+1)
                e4 = e - 1
	     else
                n  = n + 1
	        n1 = n 
	        e  = e + 3
	        e1 = e - 2
	        e4 = e - 1
	     end if
	     n2 = outrn (i+2)
             n3 = newnode (L+1)
             n4 = outrn (i+1)
	     e2 = newedge (L)
	     e3 = re (i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
       
	     x1 = xx (i+2, R+1, 1:3)
	     x2 = xx (i+2, R  , 1:3)
	     x3 = xx (i+1, L+1, 1:3)
	     x4 = xx (i+1, L  , 1:3)

             if (kutta (i+2, R+1 )) knode (n1) = .true.
	     if (kutta (i+2, R   )) knode (n2) = .true.
	     if (kutta (i+1, L+1 )) knode (n3) = .true.
	     if (kutta (i+1, L   )) knode (n4) = .true.

             if ((x1(3)==0.0).and.(j+1>=jmaxlo)) then
                jhi = 2*jmax(i+2,1) - 1 + jmax(i+2,2)
		diff = R - jhi + 1
	        if (j+2/=jend) then
                   n  = n - 1
		   n1 = nextnode (jhi-diff)
		   if (x2(3)==0.0) then
		      e  = e - 1
		      e4 = nextedge (jhi-diff)
		      e5 = e
		   end if 
	        else
		   e = e - 1
		   e4 = nextedge (jhi-diff)
		   e5 = e
	        end if 
	     end if

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n2
	     nextnode (R+1) = n1
	     nextedge (R)   = e4
!
             outrn (i+1) = n3
	     outrn (i+2) = n1
	     re (i+1) = e1

          end do
!
          else evenioddj
! we are at an even j-value

          j = djl
       
! compare sizes of strips (using same offset for left and right)
          if (djl>=jmax(i+1,1)) then
	     jmaxl = jmax(i+1,2) + djl
	  else
	     jmaxl = jmax(i+1,1)
	  end if
	  jmaxr = jmax(i+2,2) + djl
          jend1 = 2*(jmaxl-djl) - 1 + djl
	  jend2 = 2*(jmaxr-djl) - 1 + djl
	  jend  = min(jend1,jend2)
	  jmaxlo = min(jmaxl,jmaxr)
! now, shift right offset appropriately
          djr = 2*jmax(i+2,1) - (djl+1)

!	  write (*,*)
!	  write (*,*) 'djl =',djl
!	  write (*,*) 'djr =',djr
!	  write (*,*) 'drl =',drl
!	  write (*,*) 'jmaxl =',jmaxl
!	  write (*,*) 'jmaxr =',jmaxr
!	  write (*,*) 'jend =',jend

          L = j + 1
	  R = j + djr + 1
	  
	  n  = n + 2
	  n1 = n - 1
	  n2 = n
	  n3 = newnode (L+1)
!	  write (*,*)
!	  write (*,*) 'pulling ',n3,' from newnode ...',j+2
          n4 = outrn (i+1)
	  e  = e + 4
	  e1 = e - 3
	  e2 = newedge (L)
	  e3 = e - 2
	  e4 = e - 1
	  e5 = e
	  p  = p + 2
	  p1 = p - 1
	  p2 = p
!       
	  x1 = xx (i+2, R+1, 1:3)
	  x2 = xx (i+2, R  , 1:3)
	  x3 = xx (i+1, L+1, 1:3)
	  x4 = xx (i+1, L  , 1:3)

          if (kutta (i+2, R+1 )) knode (n1) = .true.
	  if (kutta (i+2, R   )) knode (n2) = .true.
	  if (kutta (i+1, L+1 )) knode (n3) = .true.
	  if (kutta (i+1, L   )) knode (n4) = .true.

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

          call quad (x1,x2,x3,x4,n1,n2,n3,n4,       &
                   & e1,e2,e3,e4,e5,p1,p2,          &
                   & bodyin,bnodex,body,bedges,top)
!
          nextnode (R)   = n2
	  nextnode (R+1) = n1
	  nextedge (R)   = e4
!	  write (*,*)
!	  write (*,*) 'putting ',n2,' into nextnode ...',j+djr+1
!	  write (*,*) 'putting ',n1,' into nextnode ...',j+djr+2
!
          outrn (i+1) = n3
          outrn (i+2) = n1
          re (i+1) = e1
!
          do j = 2+djl, jend, 2 
!
             if ((j==jmaxl.or.j==jmaxr).and.(jmaxl/=jmaxr).and.offset==0) then
	     
	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then
	           call wrapmax (i+1,j,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,   &
                               & outln,outrn,inln,inrn,newnode,newedge,  &
			       & kutta,knode)
	        else
!		   write (*,*) 'wrapping in nextstrip!'
		
	           call wrap (i+1,j,offset,bnodex,body,bedges,bodyin,  &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,    &
                            & outln,outrn,inln,inrn,nextnode,nextedge, &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2
	     R = j + djr2 + djr

	     n  = n + 1
             n1 = outrn (i+2)
             n2 = outrn (i+1)
	     n3 = n
	     n4 = newnode (L+1)
!	  write (*,*)
!	  write (*,*) 'pulling ',n4,' from newnode ...',j+1
	     e  = e + 3
	     e1 = e - 2 
	     e2 = e - 1 
	     e3 = newedge (L)
             e4 = re (i+1)
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!       
	     x1 = xx (i+2, R  , 1:3)
	     x2 = xx (i+1, L  , 1:3)
	     x3 = xx (i+2, R+1, 1:3)
	     x4 = xx (i+1, L+1, 1:3)

             if (kutta (i+2, R   )) knode (n1) = .true.
	     if (kutta (i+1, L   )) knode (n2) = .true.
	     if (kutta (i+2, R+1 )) knode (n3) = .true.
	     if (kutta (i+1, L+1 )) knode (n4) = .true.

             if ((x3(3)==0.0).and.(j>=jmaxr)) then
                jhi = 2*jmax(i+2,1) - 1 + jmax(i+2,2)
		diff = R - jhi + 1
                n  = n - 1
		n3 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n3,' from nextnode ...',jhi-diff
!	  write (*,*)
		if (x1(3)==0.0) then
		   e  = e - 1
		   e1 = nextedge (jhi-diff)
		   e2 = e - 1
		   e5 = e
		end if 
	     end if

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4
       
             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n1
	     nextnode (R+1) = n3
	     nextedge (R)   = e1
!	  write (*,*)
!	  write (*,*) 'putting ',n1,' into nextnode ...',j+djr2+djr
!	  write (*,*) 'putting ',n3,' into nextnode ...',j+djr2+djr+1
!       
             outrn (i+1) = n4
             outrn (i+2) = n3
             re (i+1) = e2      
!	     
             if (j==jend-1) exit

             if (((j+1)==jmaxl.or.(j+1)==jmaxr).and.(jmaxl/=jmaxr).and. &
	                                           & offset==0) then

	        offset = jmaxr - jmaxl
        	if (offset>0) then
		   djr2 =  2*offset
		else
		   djl2 = -2*offset
		end if
	        if (i+1==imax-1.and.jmaxr==1) then
	           call wrapmax (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                               & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
                               & outln,outrn,inln,inrn,newnode,newedge,    &
			       & kutta,knode)
                else
!		   write (*,*) 'wrapping 2 in nextstrip!'
	           call wrap (i+1,j+1,offset,bnodex,body,bedges,bodyin, &
                            & xx,n,e,p,nmax,nsmall,le,re,inle,inre,     &
	                    & outln,outrn,inln,inrn,nextnode,nextedge,  &
			    & newnode,newedge,kutta,knode)
	        end if
	     end if

             L = j + djl2 + 1
	     R = j + djr2 + djr + 1

             n  = n + 1
             n1 = n 
	     n2 = outrn (i+2)
             n3 = newnode (L+1)
!	  write (*,*)
!	  write (*,*) 'pulling ',n3,' from newnode ...',j+2
             n4 = outrn (i+1)
	     e  = e + 3
	     e1 = e - 2
	     e2 = newedge (L)
	     e3 = re (i+1)
	     e4 = e - 1
	     e5 = e
	     p  = p + 2
	     p1 = p - 1
	     p2 = p
!       
	     x1 = xx (i+2, R+1, 1:3)
	     x2 = xx (i+2, R  , 1:3)
	     x3 = xx (i+1, L+1, 1:3)
	     x4 = xx (i+1, L  , 1:3)

             if (kutta (i+2, R+1 )) knode (n1) = .true.
	     if (kutta (i+2, R   )) knode (n2) = .true.
	     if (kutta (i+1, L+1 )) knode (n3) = .true.
	     if (kutta (i+1, L   )) knode (n4) = .true.

             if ((x1(3)==0.0).and.(j+1>=jmaxr)) then
                jhi = 2*jmax(i+2,1) - 1 + jmax(i+2,2)
		diff = R - jhi + 1
!		write (*,*) 'jhi=',jhi
!		write (*,*) 'diff =',diff
                n  = n - 1
		n1 = nextnode (jhi-diff)
!	  write (*,*) 'resetting ',n1,' from nextnode ...',jhi-diff
!	  write (*,*)
		if (x2(3)==0) then
		   e  = e - 1
		   e4 = nextedge (jhi-diff)
		   e5 = e
	        end if 
	     end if

!       write (*,*)
!       write (*,324) ' i =',i+1,'  ...   j =',j+1
!       write (*,*) '__________________________________________________________'
!       write (*,*) 
!       write (*,325) ' n =',n ,' ...   p =',p ,' ...   e =',e
!       write (*,*) 
!       write (*,326) 'n1 =',n1,' ...  n2 =',n2,' ...  n3 =',n3,' ...  n4 =',n4
!       write (*,*) 
!       write (*,327) 'p1 =',p1,' ...  p2 =',p2
!       write (*,*) 
!       write (*,328) 'e1 =',e1,' ...  e2 =',e2,' ...  e3 =',e3,' ...  e4 =',e4, &
!                             & ' ...  e5 =',e5
!       write (*,*)
!       
!       write (*,*) 'x1 =',x1
!       write (*,*) 'x2 =',x2
!       write (*,*) 'x3 =',x3
!       write (*,*) 'x4 =',x4

             call quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
             nextnode (R)   = n2
	     nextnode (R+1) = n1
	     nextedge (R)   = e4
!	  write (*,*)
!	  write (*,*) 'putting ',n3,' into nextnode ...',j+djr2+djr+1
!	  write (*,*) 'putting ',n2,' into nextnode ...',j+djr2+djr+2
             outrn (i+1) = n3
	     outrn (i+2) = n1
	     re (i+1) = e1
!
          end do
	         
          end if evenioddj
       
       end if oddi

      open (8,file='bnode')
      write(8,*) 'title="whale"'
      write(8,*) 'variables = x,y,z'
      write(8,909)  'zone t=body, n=',n,', e=',p, &
                       & ', f=fepoint, et=triangle'
 909  format(a,i6,a,i6,a)
      do 5 dum=1,n
         write(8,6) bnodex(dum)%x
  6      format(3e16.6)
  5   continue
      do 81 dum=1,p
         write(8,'(3i6)') body(dum)%node
  81  continue

      return
       
      end subroutine nextstrip
! ---------------------------------------------------------------
