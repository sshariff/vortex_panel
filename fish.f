! -------------------------------------------------------------------
      program fish
!
! This code uses vortex-panel methods to solve for the flow around
! a fish with specified body motions.
!
      use panel_type; use interfish1
! -----------------------------------------------------------------
      implicit none
! ----------------------------------------------------------
! common declarations
      integer :: bfile= 9
      integer :: pfile=10
      integer :: efile=11
      integer :: ofile=12
      integer :: keybd= 5
      integer :: scrn = 6
      common/files/ bfile, pfile, efile, ofile, keybd, scrn
!
! parameter declarations
      integer, parameter :: nmax=16000       
      real               :: tratio=10.0
      integer            :: wlen1=200, wlen2=200, wlen3=200
!
! allocatable arrays
      real, allocatable, dimension(:,:) :: a
!
! global variables
      integer                           :: adim
      integer                           :: bednum, wednum
      real,        dimension(nmax)      :: bldelt
      integer                           :: bnodes, wnodes
      integer                           :: bpnum,  wpnum
      logical,     dimension(nmax)      :: bodyout
      real,        dimension(nmax)      :: cp
      real                              :: drag
      real                              :: dt
      integer                           :: dumpinc
      logical,     dimension(nmax)      :: euse
      real                              :: lift
      integer,     dimension(25)        :: lumplen1, lumplen2, lumplen3
      integer,     dimension(25,100)    :: lumplist1, lumplist2, lumplist3
      integer                           :: maxstep
      integer,     dimension(25,100)    :: nrlist1, nrlist2
      integer                           :: nseplines
      logical,     dimension(nmax)      :: nuse
      real,        dimension(nmax)      :: oldcirc
      character*30                      :: oname
      logical,     dimension(nmax)      :: puse
      real                              :: qinf
      real,        dimension(nmax)      :: right
      integer                           :: step
      real                              :: tau
      real                              :: thyme
      real,        dimension(nmax,3)    :: vref
      real,        dimension(nmax,3)    :: wind
       
      character*8                       :: dayt
      character*10                      :: tim
      character*5                       :: zone
      integer,     dimension(8)         :: values
      real                              :: sec,oldsec


!
! global panel variables using defined structures
! (structures are defined in panel_type.f)
      type (panel), dimension(nmax) ::  body, wake
!      type (panel), dimension(nmax) ::  oldbody
!  !!!!!!! for oldbody, do we need more than just the last gamma???!!!!!
      type (edges), dimension(nmax) ::  bedges, wedges
      type (nodes), dimension(nmax) ::  bnodex, wnodex
      type (nodes), dimension(nmax) ::  blcoll
      type (nodes), dimension(nmax) ::  oldbndx
      type (splst), dimension(25)   ::  seplist
      type (stack)                  ::  estack,nstack,pstack
!
! local variables (which were in init.f)
      character*40 :: bname, pname, ename
      integer      :: in
      integer      :: kedge
      integer      :: n
      integer      :: n1, n2, n3
      integer      :: e1, e2, e3
      integer      :: p, edg, line
! -----------------------------------------------------------
! common declarations
!        bfile: unit number of body motion data file
!        pfile: unit number of body panel data file
!        efile: unit number of body edge-node data
!        ofile: unit number of output data file
!        keybd: unit number for keyboard input
!        scrn:  unit number for screen output
! -----------------------------------------------------------
! parameter declarations
!        nmax: maximum number of body nodes. For N nodes, 
!              where N is large, there will be about
!              2N panels, and 3N edges
!        tratio: approximate number of time steps per period
!        wlen1, wlen2, wlen3: lengths of three wake zones, in panels.
!               panels are lumped at the ends of these zones, where
!               panels double in size. at the end of the last zone,
!               the vorticity enters the "pool" which is the final                            
!               destination of all vorticity
! ----------------------------------------------------------      
! global variables
!        a: matrix of influence coefficients. a(i,j) is the effect
!           at collocation point of panel i from panel j
!        bednum, wednum: total number of edges for body, wake
!        bldelt: boundary layer disp. thickness for each panel
!        bnodes, wnodes: total number of nodes for body, wake
!        bodyout: logical array used to make body normals point out
!        bpnum,  wpnu: total number of panels for body, wake
!        cp: pressure coefficient for each panel
!        drag: drag
!        dt: size of current timestep
!        dumpinc: number of steps between data dump
!        euse: lookup table: is wake edge being used?
!        lift: lift
!        lumplen1,2,3: number of edges at lumping sites
!        lumplist1,2,3: list of the edges which will be lumped  next
!                       lumplist (i,j) is the jth lumping edge for 
!                       the ith separation line
!        maxstep: maximum timesteps of prescribed body data
!        nrlist1,nrlist2: list of near edges generated during lumping,
!                   used to reconnect newly lumped panels to old ones
!        nseplines: number of separation lines
!        nuse: lookup table: is wake node being used?
!        oname: name of output data file
!        puse: lookup table: is wake panel being used?
!        qinf: magnitude of freestream velocity
!        right: rhs for matrix inversion
!        step: current step number
!        tau: fundamental period of the fish motion
!        time: current value of time
!        vref: reference velocity at the collocation point of each 
!              panel. vref is the kinematic velocity of each panel,
!              consisting of the freestream plus the velocity due to
!              translation and rotation
!        wind: induced velocity at the collocation points
! ----------------------------------------------------------
! global panel variables using defined structures
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
!        oldbody: body panels from last time step. this is needed
!                 for the calculation of some time derivatives
!        bedges, wedges: edge-node arrays, lists nodes that 
!                        correspond to each edge
!        bnodex,wnodex: node coordinate arrays (x,y,z for each node)
!                       for body and wake
!        blcoll: collocation points displaced by boundary layer
!        oldbndx: node coords from last timestep, used by rhs
!                 to determine the velocity of each panel
!        seplist: list of separation edges
!        estack, nstack, pstack: stacks of element numbers 
!                                 available for the wake
! -----------------------------------------------------------
! start:
!
! initializations
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
!
! now we have to perform the steps of init.f in the main program, so that
! we can properly take advantage of allocation
!
!!      call init (maxstep,thyme,step,dumpinc,oname,body,bpnum,bodyout, &
!!               & bnodes,bednum,bedges,wpnum,wnodes,wednum,nseplines,  &
!!	       & oldcirc, &
!!               & seplist,euse,nuse,puse,estack,nstack,pstack,tau,nmax)
      bfile = 9
      pfile = 10
      efile = 11
      ofile = 12
      keybd = 5
      scrn = 6

! initialize some variables
      thyme = 0.
      step = 1
      wpnum = 0
      wnodes = 0
      wednum = 0
!
! open files
      write (scrn,*) 'Body node file?'
      read  (keybd,'(a)') bname
      write (scrn,*)
      write (scrn,*) 'Panel connectivity file?'
      read  (keybd,'(a)') pname
      write (scrn,*)
      write (scrn,*) 'Edge-node file?'
      read  (keybd,'(a)') ename
      write (scrn,*)
      write (scrn,*) 'Output file?'
      read  (keybd,'(a)') oname
      write (scrn,*) 'Step increment between data dumps?'
      read  (keybd,*) dumpinc 
      
      open (bfile,file=bname,status='old')
      open (pfile,file=pname,status='old')
      open (efile,file=ename,status='old')
      open (ofile,file=oname,status='new')
      
! read headers
      read (bfile,*) tau
      read (bfile,*) bnodes
      read (bfile,*) maxstep
      if (bnodes.gt.nmax) then
        write (scrn,*) 'Number of body nodes is too large!'
	write (scrn,*) 'bnodes,nmax=',bnodes,nmax
	write (scrn,*) 'Change nmax parameter...'
	stop
      endif
      
      read (pfile,*) bpnum
      if (bpnum.gt.nmax) then
        write (scrn,*) 'Number of body panels is too large!'
	write (scrn,*) 'bpnum,nmax=',bpnum,nmax
	write (scrn,*) 'Change nmax parameter...'
	stop
      endif
     
      read (efile,*) bednum
      if (bednum.gt.nmax) then
         write (scrn,*) 'Number of edges is too large!'
   	 write (scrn,*) 'bednum,nmax=',bednum,nmax
	 write (scrn,*) 'Change nmax parameter...'
	 stop
      endif
!      
! set all wake elements to non-use
      euse(:) = .false.
      nuse(:) = .false.
      puse(:) = .false.       
!
      estack%bottom = 0
      nstack%bottom = 0
      pstack%bottom = 0
!
! read in panel connectivity information
      do 50 n=1,bpnum
         read (pfile,*) p,n1,n2,n3,e1,e2,e3,in
 	 body(p)%node(1) = n1
 	 body(p)%node(2) = n2
 	 body(p)%node(3) = n3
 	 body(p)%edge(1) = e1
 	 body(p)%edge(2) = e2
 	 body(p)%edge(3) = e3
	 if (in.eq.0) then
	    bodyout(p) = .true.
	 else
	    bodyout(p) = .false.
	 endif
 50   continue
!
! read in node list for edges      
      do 100 n=1,bednum
         read (efile,*) edg, n1, n2, e1, e2
         bedges(edg)%node1 = n1
         bedges(edg)%node2 = n2
         bedges(edg)%panel1 = e1
	 bedges(edg)%panel2 = e2
 100  continue
!
      oldcirc = 0.0    
! 
! get Kutta edges
      read (efile,*) nseplines
!      write (ofile,*) 'nseplines = ',nseplines
      do 99 line = 1,nseplines
         read (efile,*) seplist(line)%nedges
! 	 write (ofile,*) 'nedges = ',seplist(line)%nedges
	 do 88 edg = 1,seplist(line)%nedges
	    read (efile,*) kedge
	    seplist(line)%sepedge(edg)%bsedge = kedge
! 	    write (ofile,*) 'edge = ',seplist(line)%sepedge(edg)%bsedge
	    seplist(line)%sepedge(edg)%bsnode(1) = bedges(kedge)%node1                     
            seplist(line)%sepedge(edg)%bsnode(2) = bedges(kedge)%node2
  88     continue
  99  continue

      close (pfile)
      close (efile)

!      call init (maxstep,thyme,step,dumpinc,oname,body,bpnum,bodyout, &
!               & bnodes,bednum,bedges,wpnum,wnodes,wednum,nseplines,  &
!               & oldcirc, &
!               & seplist,euse,nuse,puse,estack,nstack,pstack,tau,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'init    done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
! set up body nodes
      call bmotion (dt,qinf,body,bnodes,thyme,step,bnodex,oldbndx, &
                  & bodyout,blcoll,bldelt,bpnum,tau,tratio,nmax)             
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'bmotion done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! set up initial wake at Kutta surfaces
      call kutta (bnodex,wnodes,wnodex,wake,wpnum,wednum,  &
                & wedges,nseplines,oldbndx,lumplist1,lumplen1,         &
                & nrlist1,nrlist2,euse,nuse,puse,estack,nstack,pstack, &
                & seplist,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'kutta   done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'

      adim = bpnum + wpnum
      allocate(a(adim,adim))
!
! compute influence coefficients (Kutta edges are included)
!      call aij1 (bpnum,wpnum,body,a,dt,bedges,wedges,bnodex,wnodex,wake, &
!               & oldbndx,right,vref,qinf,puse,nseplines,seplist,nmax)
      call aij1 (bpnum,wpnum,body,a,dt,bedges,wedges,bnodex,wnodex,wake, &
               & oldbndx,right,vref,qinf,nseplines,seplist,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'aij1    done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! solve for circulations
!      call solve (bpnum+wpnum,wpnum,a,right,body,wake,wedges,nseplines, &
!                & seplist,nmax)        
!      call solve (bpnum+wpnum,wpnum,a,right,body,wake,wedges,nseplines, &
!                & oldcirc,seplist,nmax)        
!
!
!
      call solve (adim,bpnum,wpnum,a,right,body,wake,wedges,nseplines, &
                & oldcirc,seplist,nmax,bedges,bnodex,wnodex)        
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'solve   done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! calculate pressures and velocities
      call presvel (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,      &
                  & puse,cp,wind,lift,drag,bldelt,seplist,nseplines, &       
                  & oldcirc,body,wake,vref,step,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'presvel done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! output data
      call dump (bednum,wednum,bpnum,wpnum,bedges,wedges,bnodes,wnodes, &
               & bnodex,wnodex,maxstep,dumpinc,oname,cp,wind,vref,      &
	       & lift,drag,euse,nuse,puse,nseplines,seplist,body,wake,  &
	       & thyme,step,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'dump    done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! compute wake rollup     
      call rollup (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,wnodes, &
                 & nuse,puse,body,wake,nmax)
!
      oldsec = sec
      call date_and_time(dayt,tim,zone,values)
      sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
          & real(values(8))/1000
!
      write (*,*)
      write (*,'(2a,f12.3,a)') 'rollup  done in ....................', &
                             & '.......... ',sec-oldsec,' seconds'
!
! main loop
      do 9999 step = 2,maxstep
!
! move body nodes
        call bmotion (dt,qinf,body,bnodes,thyme,step,bnodex,oldbndx, &
	            & bodyout,blcoll,bldelt,bpnum,tau,tratio,nmax)
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000
!
        write (*,*)
        write (*,'(2a,f12.3,a)') 'bmotion done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
! create new panels, edges, vortices
	call newpan (bnodes,bnodex,wnodes,wnodex,body,wake,wpnum,         &   
                   & wednum,wedges,bedges,nseplines,lumplist1,lumplist2,  &
		   & lumplist3,lumplen1,lumplen2,lumplen3,nrlist1,        &
		   & nrlist2,euse,nuse,puse,estack,nstack,pstack,seplist, &
		   & wlen1,wlen2,wlen3,step,nmax)
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'newpan  done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
! compute influence coefficients
	call aij (bpnum,wpnum,body,a,bedges,wedges,bnodex,wnodex,wake, &
                & oldbndx,right,dt,vref,qinf,puse,nmax)  
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'aij     done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
		
!
! solve for circulations
! the 0 is a dummy argument because we don't have to solve for the
! wake panels
!       call solve (bpnum+wpnum,wpnum,a,right,body,wake,wedges,nseplines, &
!                 & oldcirc,seplist,nmax)        
!
!
!
	call solve (adim,bpnum,0,a,right,body,wake,wedges,nseplines, &
                  & oldcirc,seplist,nmax,bedges,bnodex,wnodex)        
!
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'solve   done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
! calculate pressures and velocities
	call presvel (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,     &
                    & puse,cp,wind,lift,drag,bldelt,seplist,nseplines, &
                    & oldcirc,body,wake,vref,step,nmax)
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'presvel done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
! output data
        call dump (bednum,wednum,bpnum,wpnum,bedges,wedges,bnodes,wnodes, &
	         & bnodex,wnodex,maxstep,dumpinc,oname,cp,wind,vref,lift, &
		 & drag,euse,nuse,puse,nseplines,seplist,body,wake,       &
		 & thyme,step,nmax)
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'dump    done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
! compute wake rollup          
	call rollup (dt,bpnum,wpnum,bnodex,wnodex,bedges,wedges,wnodes, &
	           & nuse,puse,body,wake,nmax)
!
        oldsec = sec
        call date_and_time(dayt,tim,zone,values)
        sec = real(values(5)*3600 + values(6)*60 + values(7)) + &
            & real(values(8))/1000

        write (*,*)
        write (*,'(2a,f12.3,a)') 'rollup  done in ....................', &
                               & '.......... ',sec-oldsec,' seconds'
!
 9999 continue
! end main loop

! close files
      close(bfile)
      close(ofile)
!
      write(scrn,*)
      write(scrn,*)'Simulation complete...'
      stop
      end
! --------------------------------------------------------------------
