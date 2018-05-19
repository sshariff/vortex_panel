! --------------------------------------------------------------------------
      subroutine panvel (bpnum,body,bnodex,oldbndx,right,dt,vref,qinf,nmax)
!
!  This subroutine clears the RHS vector, and sets it to the 
!  local normal velocity of the panel
!
! -----------------------------------------------------------------
      use panel_type; use interfish4
!
      implicit none 
!
      integer,                         intent(in)  :: bpnum
      type (panel), dimension(nmax),   intent(in)  :: body
      type (nodes), dimension(nmax),   intent(in)  :: bnodex
      type (nodes), dimension(nmax),   intent(in)  :: oldbndx
      real,         dimension(nmax),   intent(out) :: right
      real,                            intent(in)  :: dt
      real,         dimension(nmax,3), intent(out) :: vref
      real,                            intent(in)  :: qinf
      integer,                         intent(in)  :: nmax

      integer :: pan, n
! ---------------------------------------
!   input variables
!        body: body panels, contains strength of panel, lists
!              of nodes & edges, "polarity", normal direction
!              and centroid      
!        bpnum: total number of panels for body
!        bnodex: body node coordinate array (x,y,z for each node)
!        dt: size of current timestep
!        nmax: maximum number of edges
!        oldbndx: node coords from last timestep, used for rhs
!                 to determine the velocity of each panel
!        qinf: magnitude of freestream velocity
! -----------------------------------------
!   output variables
!        right: rhs for matrix inversion
!        vref: reference velocity at the collocation point of each 
!              panel. vref is the kinematic velocity of each panel,
!              consisting of the freestream plus the velocity due to
!              translation and rotation
! -----------------------------------------
!   local variables
!        pan,n: indices for stepping
! ------------------------------------------
!  start:
! 
!  get the kinematic velocities (sum of freestream velocity
!  and the velocity of the panel. note here freestream = 0,
!  because we have a ground-fixed reference frame)
! 
      call kine (bpnum,bnodex,oldbndx,dt,vref,qinf,nmax)
!
      do 1600 pan = 1,bpnum
!        take normal component of kinematic velocity 
         right(pan) = dot_product((-vref(pan,1:3)),body(pan)%normal)
!	 write(*,*) 'normal =',body(pan)%normal
!	 write(*,*) 'vref = ',vref(pan,1:3)
 1600 continue
!
      return
      end
! --------------------------------------------------------------------
