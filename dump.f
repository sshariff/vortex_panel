! -------------------------------------------------------------------
      subroutine dump (bednum,wednum,bpnum,wpnum,bedges,wedges,     &
                     & bnodes,wnodes,bnodex,wnodex,maxstep,dumpinc, &
                     & oname,cp,wind,vref,lift,drag,euse,nuse,puse, &
		     & nseplines,seplist,body,wake,time,step,nmax) 
!
!  This subroutine dumps data.
!  --------------------------------------------------------------
      use panel_type; use interfish2; use interfish4
!
      implicit none
!
      integer,                       intent(in) :: bednum, wednum
      integer,                       intent(in) :: bpnum,  wpnum
      type (edges), dimension(nmax), intent(in) :: bedges, wedges
      integer,                       intent(in) :: bnodes, wnodes
      type (nodes), dimension(nmax), intent(in) :: bnodex, wnodex
      integer,                       intent(in) :: maxstep
      integer,                       intent(in) :: dumpinc
      character*30,                  intent(in) :: oname
      real,         dimension(nmax), intent(in) :: cp
      real,       dimension(nmax,3), intent(in) :: wind
      real,       dimension(nmax,3), intent(in) :: vref
      real,                          intent(in) :: lift
      real,                          intent(in) :: drag
      logical ,     dimension(nmax), intent(in) :: euse, nuse, puse
      integer,                       intent(in) :: nseplines
      type (splst), dimension(25),   intent(in) :: seplist
      type (panel), dimension(nmax), intent(in) :: body, wake
      real ,                         intent(in) :: time
      integer,                       intent(in) :: step
      integer,                       intent(in) :: nmax

      character*40             :: name
      integer                  :: dfile
      character*3              :: dumpnum
      integer                  :: length
      integer                  :: nused, pused
      integer, dimension(nmax) :: wakemap
      
      real, dimension(bpnum,bpnum,3) :: lx,lu,lx2
      real, dimension(bpnum,bpnum)   :: lcp,lgamma 
      
      integer :: i,j,n1,n2,e1,e2,p1,p2,p3,e3
      integer :: firste,imax,jmax,lpan,jmid
      
      real :: normsign,cl,cd,area,totalarea,magx
      real :: clmid,cdmid,totalamid
      real, dimension(3) :: qtot,qn,df,r1,r2,r1xr2
      
      integer, dimension(bpnum,bpnum) :: pannum
      
!      real :: cl2,cd2,df2,totalarea2,dcl,dcd

      integer   ::  bfile, pfile, efile, ofile, keybd, scrn
      common/files/ bfile, pfile, efile, ofile, keybd, scrn
! --------------------------------------------------------------
! input variables
!        bednum, wednum: total number of edges for  body, wake
!        bedges, wedges: edge-node arrays, lists nodes that 
!                        correspond to each edge
!        bnodes, wnodes: total number of nodes for body, wake
!        bnodex,wnodex: node coordinate arrays (x,y,z for each node)
!                       for body and wake
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        bpnum,  wpnum: total number of panels for body, wake
!        cp: pressure coefficient for each panel
!        drag: drag
!        dumpinc: number of steps between data dump
!        euse,nuse,puse: lookup table: is wake element being used?
!        lift: lift
!        maxstep: final step number
!        nmax: maximum number of edges 
!        oname: name of output data file
!        step: current step number
!        time: current value of time
!        vref: reference velocity at the collocation point of each 
!              panel. vref is the kinematic velocity of each panel,
!              consisting of the freestream plus the velocity due to
!              translation and rotation
!        wind: induced velocity at the collocation points
! --------------------------------------------------------------
! local variables
!        barea: body panel areas
!        bncirc, wncirc: circulation averaged at a node
!        bncp,wncp: pressure coefficient averaged at a node
!        bnwind, wnwind: induced velocity at the nodes
!        bodtest: are we sending body to nodeave?
!        dfile: file number for the dump file
!        dumpnum: dump number
!        i: stepping index
!        len: length of output filename
!        magx: magnitude of r1xr2 (cross product of two panel edges)
!        name: name of output data file
!        nused, pused: number of wake nodes actually used
!        r1, r2: components of two of the panel edges
!        r1xr2: cross product of two panel edges
!        wakemap: mapping of used wake node numbers
!        warea: wake panel areas
! --------------------------------------------------------------
! common variables
!        bfile: unit number of body motion data file
!        pfile: unit number of body panel data file
!        efile: unit number of body edge-node data
!        ofile: unit number of output data file
!        keybd: unit number for keyboard input
!        scrn:  unit number for screen output
! --------------------------------------------------------------
! start:
!
      dfile = 266

      if (mod(step,dumpinc).eq.0) then
!       
         j = 1
         e1 = 4
         firste = e1

 708     continue
! 
         i = 1
!
         if (mod(j,4).eq.3) then
	    p1 = bedges(e1)%panel1
	 else
            p1 = bedges(e1)%panel2
	 end if
	 
!	 write (*,*)
!	 write (*,'(a,2i3)') 'i,j = ',i,j
!	 write (*,*)
	 
!	 write (*,'(a,i3)') 'p1 = ',p1
!	 write (*,'(a,i3)') 'e1 = ',e1
         n1 = bedges(e1)%node1 
	 n2 = bedges(e1)%node2
!
         if (mod(j,4).eq.3) then
	    pannum(i,j)=p1
	    lx2(i,j  ,:) = bnodex(n2)%x
	    lx2(i,j+1,:) = bnodex(n1)%x
         else
	    pannum(i,j+1)=p1
	    lx2(i,j  ,:) = bnodex(n1)%x
	    lx2(i,j+1,:) = bnodex(n2)%x
	 end if

!	 pannum(i,j) = p1
!	 
!	 write (*,*)
!	 write (*,'(a,2i3)') 'n1,n2        = ',n1,n2
!         lx2(i,j,:) = bnodex(n1)%x
!	 write (*,'(a,2i3,a,3f7.2)') 'lx2 (',i,j,') = ',lx2(i,j,:)
!
         e2 = body(p1)%edge(3)
	 if (mod(j,4).eq.3) then
	    p2 = bedges(e2)%panel1
	 else
	    p2 = bedges(e2)%panel2
	 end if 
	 
!	 write (*,*)
!	 write (*,'(a,2i3)') 'i,j = ',i+1,j
!	 write (*,*)
!	 write (*,'(a,i3)') 'p2 = ',p2
!	 write (*,'(a,i3)') 'e2 = ',e2
!	    
!	 pannum(i+1,j) = p2
!
         if (mod(j,4).eq.3) then
	    pannum(i,j+1) = p2
	 else
	    pannum(i,j) = p2
	 end if
	 
	 e1 = body(p2)%edge(2)
	 if (mod(j,4).eq.3) then
	    p1 = bedges(e1)%panel1
	 else
	    p1 = bedges(e1)%panel2
	 end if
	 
!	 write (*,*)
!	 write (*,'(a,2i3)') 'i,j = ',i+2,j
!	 write (*,*)
!	 write (*,'(a,i3)') 'p1 = ',p1
!	 write (*,'(a,i3)') 'e1 = ',e1
	 
!	 
!         n1 = bedges(e1)%node1 
!	 n2 = bedges(e1)%node2
!	 write (*,*)
!	 write (*,'(a,2i3)') 'n1,n2        = ',n1,n2
!         lx2(i+1,j,:) = bnodex(n1)%x
!	 write (*,'(a,2i3,a,3f7.2)') 'lx2 (',i+1,j,') = ',lx2(i+1,j,:)
!
	 if (p1.eq.0.or.e1.eq.firste) go to 1400

 788	 continue   
 
         i = i + 1
!
         n1 = bedges(e1)%node1 
	 n2 = bedges(e1)%node2

         if ((mod(i,2)+mod(j,4).eq.1).or.&
	   & (mod(i,2)+mod(j,4).eq.4)) then
	    pannum(i,j) = p1
	 else
	    pannum(i,j+1) = p1
	 end if

         if (mod(i,2).eq.1) then
	    if (mod(j,4).eq.3) then
	       lx2(i,j  ,:) = bnodex(n2)%x
	       lx2(i,j+1,:) = bnodex(n1)%x
	    else
	       lx2(i,j  ,:) = bnodex(n1)%x
	       lx2(i,j+1,:) = bnodex(n2)%x
	    end if
	 end if

!	 pannum(i,j) = p1
!	 
!         n1 = bedges(e1)%node1 
!	 n2 = bedges(e1)%node2
!	 write (*,*)
!	 write (*,'(a,2i3)') 'n1,n2        = ',n1,n2
!         lx2(i,j,:) = bnodex(n1)%x
!	 write (*,'(a,2i3,a,3f7.2)') 'lx2 (',i,j,') = ',lx2(i,j,:)
!
	 e2 = body(p1)%edge(3)
	 if (mod(j,4).eq.3) then
	    p2 = bedges(e2)%panel1
	 else
	    p2 = bedges(e2)%panel2
	 end if
!	 
!	 write (*,*)
!	 write (*,'(a,2i3)') 'i,j = ',i+1,j
!	 write (*,*)
!	 write (*,'(a,i3)') 'p2 = ',p2
!	 write (*,'(a,i3)') 'e2 = ',e2

         if ((mod(i,2)+mod(j,4).eq.1).or.&
	   & (mod(i,2)+mod(j,4).eq.4)) then
	    pannum(i,j+1) = p2
	 else
	    pannum(i,j) = p2
	 end if


!	 pannum(i+1,j) = p2
!
         if (mod(i,2).eq.0) then
            e1 = body(p2)%edge(1)
         else
	    e1 = body(p2)%edge(2)
         end if
!
         if (mod(j,4).eq.3) then
	    p1 = bedges(e1)%panel1
         else
	    p1 = bedges(e1)%panel2
	 end if
!
!	 write (*,*)
!	 write (*,'(a,2i3)') 'i,j = ',i+2,j
!	 write (*,*)
!	 write (*,'(a,i3)') 'p1 = ',p1
!	 write (*,'(a,i3)') 'e1 = ',e1
!	 
         n1 = bedges(e1)%node1 
	 n2 = bedges(e1)%node2
!	 write (*,*)
!	 write (*,'(a,2i3)') 'n1,n2        = ',n1,n2
!         lx2(i+1,j,:) = bnodex(n1)%x
!	 write (*,'(a,2i3,a,3f7.2)') 'lx2 (',i+1,j,') = ',lx2(i+1,j,:)
!
         if (mod(i,2).eq.0) then
	    if (mod(j,4).eq.3) then
	       lx2(i,j  ,:) = bnodex(n2)%x
	       lx2(i,j+1,:) = bnodex(n1)%x
	    else
	       lx2(i,j  ,:) = bnodex(n1)%x
	       lx2(i,j+1,:) = bnodex(n2)%x
	    end if
	 end if



	 if (p1.eq.0.or.e1.eq.firste) go to 1400

	 go to 788

 1400    continue
 
!         write (*,*) 'at 1400,     i,j,p1,e1,firste =',i,j,p1,e1,firste
!         if (j.eq.1.and.p1.eq.0)  write (*,*) 'Surface is not closed...'
!	 if (j.eq.1.and.e1.eq.firste)  write (*,*) 'Surface is closed...'
!
         i = i + 1
	 
         n1 = bedges(e1)%node1 
	 n2 = bedges(e1)%node2
!
         if (mod(j,4).eq.3) then
	    pannum(i,j)=p1
	    lx2(i,j  ,:) = bnodex(n2)%x
	    lx2(i,j+1,:) = bnodex(n1)%x
         else
	    pannum(i,j+1)=p1
	    lx2(i,j  ,:) = bnodex(n1)%x
	    lx2(i,j+1,:) = bnodex(n2)%x
	 end if
	 
         e2 = body(p1)%edge(3)
	 if (mod(j,4).eq.3) then
	    p2 = bedges(e2)%panel1
	 else
	    p2 = bedges(e2)%panel2
	 end if 
	 
         if (mod(j,4).eq.3) then
	    pannum(i,j+1) = p2
	 else
	    pannum(i,j) = p2
	 end if
	 
!	 pannum(i,j) = p1
!	 
!         n1 = bedges(e1)%node1 
!	 n2 = bedges(e1)%node2
!	 write (*,*)
!	 write (*,'(a,2i3)') 'n1,n2        = ',n1,n2
	 
!         lx2(i,j,:) = bnodex(n1)%x
!	 write (*,'(a,2i3,a,3f7.2)') 'lx2 (',i,j,') = ',lx2(i,j,:)
!
         imax = i
!	 write (*,*) 'imax =',imax
!	  
         if (mod(j,4).eq.3) then
	    p1 = bedges(e1)%panel1
!	    write (*,*) 'p1 = ',p1
	    e3 = body(p1)%edge(3)
!	    write (*,*) 'e3 = ',e3
	    p3 = bedges(e3)%panel1
!	    write (*,*) 'p3 = ',p3
	    e3 = body(p3)%edge(1)
!	    write (*,*) 'e3 = ',e3
	    if (bedges(e3)%panel1.eq.0) go to 1468
	    if (j.eq.3) then
	       e1 = firste + 13
	    else
	       e1 = firste + 12
	    end if
         else 
	    e1 = firste + 2
	 end if
!
	 write (*,*) 'resetting e1 to ',e1,' for next strip'
	 
	 firste = e1
	 
	 if (mod(j,4).eq.3) then
	    p1 = bedges(e1)%panel2
	 else
	    p1 = bedges(e1)%panel1
	 end if
!
!	 write (*,'(a,i3)') 'p1 = ',p1
!	 write (*,'(a,i3)') 'e1 = ',e1
!        e1 = firste + e2 - 1
!	 write (*,*) 'firste =',firste
!	 write (*,*) 'e2 =',e2
         if (e1.gt.bednum) go to 1468
!
! 1466    continue
!
  1467   continue

         j = j + 2
         go to 708
!	 
 1468    continue
         jmax = j+1
	 write (*,*) 'jmax = ',jmax
! 
         write(dumpnum,'(i3)') int(step/dumpinc)
         do 17 i=1,3
	    if(dumpnum(i:i).eq.' ') dumpnum(i:i)='0'
17	 continue
	 i = 1
9	 if(oname(i:i).ne.' ')then
	    i = i + 1
	    length = i
	    go to 9
	 end if
         length = length - 1
	 name = oname(1:length)//'.d'//dumpnum//'.tec'
	 open(dfile,file=name)

!         do 766 i=1,bpnum
!            write (ofile,765) 'body (',i,') = ',body(i)%circ
! 765        format(a,i4,a,f16.6)
! 766     continue
!
!         write (ofile,*)
!         do 866 i=1,wpnum
!            write (ofile,765) 'wake (',i,') = ',wake(i)%circ
!! 865        format(a,i4,a,f16.6)
! 866     continue
!
         write (ofile,*)
 
         write (ofile,*)
	 write (ofile,'(a,f16.6)') 'lift =',lift
	 write (ofile,'(a,f16.6)') 'drag =',drag
!
         write(dfile,*) 'title="',oname(1:length)//'.d'//dumpnum,'"'
         write(dfile,*) 'variables = x,y,z,u,v,w,cp,gamma'
!         write(dfile,909) 'zone t=body, n=',bnodes,', e=',bpnum, &
!                        & ', f=fepoint, et=triangle'
         write(dfile,909) 'zone t=body, i=',imax,', j=',jmax, &
                        & ', f=point'
 909     format(a,i4,a,i4,a)
!
         do 5 j = 1,jmax
	    do 4 i = 1,imax 
	       lpan = pannum(i,j)
	       if (mod(i,2).eq.0) then
	          normsign = -1.0
	       else
	          normsign =  1.0
	       end if
	       qtot = 2.0*(wind(lpan,:)+vref(lpan,:))
               qn = body(lpan)%normal*(dot_product(qtot,body(lpan)%normal))
!               write (12,'(a,i4,a,3f9.2,a)') 'panel = ',lpan,', qn = (',qn,')'
               write (dfile,6) body(lpan)%centr, &
	     & 2.0*(wind(lpan,:)+vref(lpan,:)),  &
	     & cp(lpan),normsign*body(lpan)%circ
 6             format(8e16.6)
 4          continue
 5       continue
!
         write (dfile,909) 'zone t=mesh, i=',imax,', j=',jmax, &
                        & ', f=point'
         do 15 j = 1,jmax
	    do 14 i = 1,imax
               write (dfile,6) lx2(i,j,:),0.0,0.0,0.0,0.0,0.0
 14          continue
 15       continue
!
!  determine number of wnodes in use, set up wakemap.
!  we have to collapse the unused nodes for tecplot,
!  so wakemap contains the node number counting only
!  used nodes.
!
         nused = 0
         do 9854 i = 1,wnodes
            if(nuse(i)) then
               nused = nused + 1
	       wakemap(i) = nused
            endif
 9854    continue
         pused = 0
         do 9855 i = 1,wpnum
            if(puse(i)) pused = pused + 1
 9855    continue
!
         write(dfile,909) 'zone t=wake, n=',nused,', e=',pused, &
                        & ', f=fepoint, et=triangle'
         do 15 i=1,wnodes
            if(nuse(i)) then
	       write(dfile,6) wnodex(i)%x,0.,0.,0.,0.,0.
            endif
 15      continue
         do 18 i=1,wpnum
            if(puse(i)) then
	       write(dfile,'(3i6)') wakemap(wake(i)%node(1)), &
                                  & wakemap(wake(i)%node(2)), &
                                  & wakemap(wake(i)%node(3))
            endif
 18      continue     
!
1001     close(dfile)
      
! calculate center-line lift      
      totalarea = 0.0
      totalamid = 0.0
!      totalarea2 = 0.0
      cl = 0.0
      cd = 0.0
      clmid = 0.0
      cdmid = 0.0
!      cl2 = 0.0
!      cd2 = 0.0
      jmid = jmax/2
      
      do 7213 j = 1,jmax
!      do 7213 j = jmid-1,jmid+2
         do 7212 i = 1, imax - 1
	 
!	 write (*,*)
!	 write (*,*)
!	 write (*,*) 'i, j =',i,j
!	 write (*,*) 'pannum =',pannum(i,j)
	 
	 r1 = bnodex(body(pannum(i,j))%node(2))%x &
	  & - bnodex(body(pannum(i,j))%node(1))%x
!	  
	 r2 = bnodex(body(pannum(i,j))%node(3))%x &
	  & - bnodex(body(pannum(i,j))%node(1))%x
!
         call cross (r1, r2, r1xr2, magx)
	 area = magx/2.0
	 
!	 write (*,*)
!	 write (*,*) 'node 1  =',body(pannum(i,j))%node(1)
!	 write (*,*) 'node 2  =',body(pannum(i,j))%node(2)
!	 write (*,*) 'node 3  =',body(pannum(i,j))%node(3)
!	 write (*,*) 'x_node1 =',bnodex(body(pannum(i,j))%node(1))%x
!	 write (*,*) 'x_node2 =',bnodex(body(pannum(i,j))%node(2))%x
!	 write (*,*) 'x_node3 =',bnodex(body(pannum(i,j))%node(3))%x
!	 write (*,*)
!	 write (*,*) ' r1   =',r1
!	 write (*,*) '|r1|  =',sqrt(dot_product(r1,r1))
!	 write (*,*) ' r2   =',r2
!	 write (*,*) '|r2|  =',sqrt(dot_product(r2,r2))
!	 write (*,*) 'r1xr2 =',r1xr2
!	 write (*,*) 'magx  =',magx
!	 write (*,*) 'area  =',area
	 
	 totalarea = totalarea + area*abs(body(pannum(i,j))%normal(2))
	 df = -cp(pannum(i,j))*area*body(pannum(i,j))%normal
!	 write (*,'(a,i4,a,3e16.6)') '  df (',pan,') =',df(1),df(2),df(3)
!	 if (body(pannum(i,j))%normal(2).lt.0) df = -df
	 cl = cl + df(2)
	 cd = cd + df(1)
!	 write (*,*)
!	 write (*,*) 'cp =',cp(pannum(i,j))
!	 write (*,*) 'df =',df
!	 write (*,*) 'norm  =',body(pannum(i,j))%normal
!	 write (*,*) 'df(2) =',df(2)
!	 write (*,*) 'df(1) =',df(1)
	 
	 if (j.ge.jmid-1.and.j.le.jmid+2) then
	    totalamid = totalamid + area*abs(body(pannum(i,j))%normal(2))
	    clmid = clmid + df(2)
	    cdmid = cdmid + df(1)
	  end if
	 
 7212    continue
 7213 continue
 
!      write (*,*) '2*totalarea =',totalarea
 
      totalarea = totalarea/2.0
      cl = cl/totalarea
      cd = cd/totalarea
      totalamid = totalamid/2.0
      clmid = clmid/totalamid
      cdmid = cdmid/totalamid
 
!      write (*,*)
!      write (*,*) 'totalarea =',totalarea
!      write (*,*) 'totalarea2 =',totalarea2
!      write (*,*) 'cl =',cl
!      write (*,*) 'cd =',cd
!      write (*,*) 'cl2 =',cl2
!      write (*,*) 'cd2 =',cd2
 
      write (ofile,*)
      write (ofile,'(a,g16.7)') 'cl =',cl
      write (ofile,'(a,g16.7)') 'cd =',cd
      write (ofile,*)
      write (ofile,'(a,g16.7)') 'clmid =',clmid
      write (ofile,'(a,g16.7)') 'cdmid =',cdmid

      end if
      

      return
      end
! --------------------------------------------------------------------

