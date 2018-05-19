! -------------------------------------------------------------------
      subroutine bmotion (dt,qinf,body,bnodes,time,step,bnodex,oldbndx, &
                        & bodyout,blcoll,bldelt,bpnum,tau,tratio,nmax)                   
!
!  This subroutine reads in the prescribed body-node positions
!  for the current timestep, and the size of the timestep.
!
! --------------------------------------------------------------
      use panel_type; use interfish2
!
      implicit none
!
      real,                          intent(out)   :: dt
      real,                          intent(out)   :: qinf
      type (panel), dimension(nmax), intent(inout) :: body
      integer,                       intent(in)    :: bnodes
      real,                          intent(inout) :: time
      integer,                       intent(in)    :: step
      type (nodes), dimension(nmax), intent(inout) :: bnodex
      type (nodes), dimension(nmax), intent(out)   :: oldbndx
      logical,      dimension(nmax), intent(in)    :: bodyout
      type (nodes), dimension(nmax), intent(out)   :: blcoll
      real,         dimension(nmax), intent(inout) :: bldelt
      integer,                       intent(in)    :: bpnum
      real,                          intent(in)    :: tau
      real,                          intent(in)    :: tratio
      integer,                       intent(in)    :: nmax

      integer :: n,inode,panl
      real    :: amp
      real    :: py
      real    :: xn,yn,zn

      integer ::    bfile, pfile, efile, ofile, keybd, scrn
      common/files/ bfile, pfile, efile, ofile, keybd, scrn
! ---------------------------------------------------------------      
! input variables
!        bldelt: boundary layer displacement thickness for each 
!                panel, based on calculation done during the 
!                previous step
!        bnodes: total number of nodes used to describe body
!        body: panel info structure, contains normals for boundary
!              layer displacements
!        bodyout: logical array used to make body normals point out
!        bpnum: total number of body panels
!        nmax: maximum number of edges 
!        oldbndx: node coordinates from last timestep, used to 
!                 compute panel velocity for rhs
!        step: current step number
!        tau: fundamental period of the fish motion
!        tratio: approximate number of time steps per period
! ---------------------------------------------------------------
! output variables
!        blcoll: collocation points, displaced by boundary layer
!        bnodex: node coordinates
!        dt: size of current timestep
!        qinf: magnitude of freestream velocity, used in kutta.f
!        time: current value of time
! --------------------------------------------------------------- 
! common variables
!        bfile: unit number of body motion data file
!        pfile: unit number of body panel data file
!        efile: unit number of body edge-node data
!        ofile: unit number of output data file
!        keybd: unit number for keyboard input
!        scrn:  unit number for screen output
! ---------------------------------------------------------------
! local variables
!        n, inode, panl: indices for stepping
!        amp: amplitude of vertical oscillation
!        py: pi
!        xn, yn , zn: coordinates read from body node datafile
! --------------------------------------------------------------
! start:
      amp = 0.5    
      qinf = 5.0
!      dt = tau/tratio
      dt = tau*1.0/qinf
!        (assuming c = 1.0)

      time = time + dt
      
      py = acos(-1.0)
!
!      
      if(step.eq.1)then 
!  read initial position. this is put into oldbndx, and it is
!  used to set up the initial kutta panels.
         do 931 n = 1,bnodes
!
            read (bfile,*) inode,xn,yn,zn 
	    bnodex(inode)%x(1) = xn
	    bnodex(inode)%x(2) = yn
	    bnodex(inode)%x(3) = zn
!	 write (*,*) 'inode,xn,yn,zn=',inode,xn,yn,zn
!	 write (*,*) 'bnodex(inode)=',bnodex(inode)      
            if(inode.ne.n) then
	       write(*,*) ' Nodes are out of order...'
            endif
 931     continue
         bldelt = 0.0
      endif
      oldbndx = bnodex
      do 932 n = 1,bnodes
         bnodex(n)%x(1) = bnodex(n)%x(1) - qinf*dt
!	 bnodex(n)%x(2) = bnodex(n)%x(2) +
!     &     amp*(sin(2.0*py*time/tau)-sin(2.0*py*(time-dt)/tau))
 932  continue
! 
!  compute new body normal vectors and centroids
      call norms (bnodes,bnodex,body,bodyout,bpnum,nmax)
!
!  add boundary layer displacement to collocation points
      do 999 panl = 1, bpnum
!
         blcoll(panl)%x = body(panl)%centr + &
	                & bldelt(panl)*body(panl)%normal
 999  continue
   
      write(scrn,*)
      write(scrn,*)
      write(scrn,'(a,i4)') ' Timestep ',step
      write(scrn,*)
!            
      return
      end
! --------------------------------------------------------------------

