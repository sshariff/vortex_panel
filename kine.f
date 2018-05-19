! --------------------------------------------------------------------
      subroutine kine (bpnum,bnodex,oldbndx,dt,vref,qinf,nmax)
!
!  This subroutine computes the kinematic velocity (the velocity
!  of the collocation point due to translation and rotation of the
!  body panel)
!
! -----------------------------------------------------------------
      use panel_type
!
      implicit none 
!
      integer,                         intent(in)  :: bpnum
      type (nodes), dimension(nmax),   intent(in)  :: bnodex
      type (nodes), dimension(nmax),   intent(in)  :: oldbndx
      real,                            intent(in)  :: dt
      real,         dimension(nmax,3), intent(out) :: vref
      real,                            intent(in)  :: qinf
      integer,                         intent(in)  :: nmax

      integer :: pan, n
      real    :: py
      real    :: tau
! ---------------------------------------
!   input variables
!        bpnum: total number of panels for body
!        bnodex: body node coordinate array (x,y,z for each node)
!        dt: size of current timestep
!        nmax: maximum number of edges
!        oldbndx: node coords from last timestep, used for rhs
!                 to determine the velocity of each panel
!        qinf: magnitude of freestream velocity
! -----------------------------------------
!   output variable
!        vref: reference velocity at the collocation point of each 
!              panel. vref is the kinematic velocity of each panel,
!              consisting of the freestream plus the velocity due to
!              translation and rotation
! -----------------------------------------
!   local variables
!        pan,n: indices for stepping
!        py: pi
!        tau: period of oscillation
! -----------------------------------------
!  start:
      py = acos(-1.0)
      tau = 1.0
!      
      do 9999 pan = 1,bpnum
         vref(pan,1) = qinf
!         vref(pan,2) = vref(pan,1)/10.0*2.0*py/tau*cos(2.0*py/tau*time)
         vref(pan,2) = 0.0
         vref(pan,3) = 0.0 
 9999 continue
      return
      end
! --------------------------------------------------------------------
