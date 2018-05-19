! -------------------------------------------------------------------
      subroutine diffuse (wake,wpnum,puse,nmax)
!
!  This subroutine allows the wake panels vorticity diffusion.
!  For the inviscid version, a non-physical dummy diffusion
!  calculation poorly approximates the diffusion by scaling
!  the circulation of the wake panels by a fixed factor f.
!  --------------------------------------------------------------
      use panel_type
!
      implicit none
!
      type (panel), dimension(nmax), intent(inout) :: wake
      integer,                       intent(in)    :: wpnum
      logical ,     dimension(nmax), intent(in)    :: puse
      integer,                       intent(in)    :: nmax
!
      real    :: f=0.9998
      integer :: pan
! --------------------------------------------------------------
! input variables
!        nmax: maximum number of edges 
!        puse: lookup tables: is element in use?
!        wpnum: total number of panels for wake
! --------------------------------------------------------------
! local variables
!        f: diffusion factor, here set to 0.02%
!        pan: panel number
! ---------------------------------------------------------------
! output variables
!         wake: wake panels, contains strength of panel
!               and lists of nodes & edges,  "polarity",
!               normal direction and centroid
! --------------------------------------------------------------
! start:
!
      do 6600 pan = 1,wpnum
         if (puse(pan)) then
	    wake(pan)%circ = f*wake(pan)%circ   
	 endif
 6600 continue
!
      return
      end

