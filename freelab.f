! -------------------------------------------------------------------
      subroutine freelab (label,use,stak,nmax)
!
!  This subroutine frees up an element number that is no longer
!  being used
!  --------------------------------------------------------------
      use panel_type
!
      implicit none
!
      integer,      intent(out)   :: label
      logical,      intent(inout) :: use(nmax)
      type (stack), intent(inout) :: stak
      integer,      intent(in)    :: nmax
! --------------------------------------------------------------
! input variables
!        nmax: maximum number of edges 
! --------------------------------------------------------------
! output variables
!        label: label number for new element
!        stak: stack of free labels
!        use: lookup table: is element in use?
! --------------------------------------------------------------
! start:
!
!      write (12,*)
!      write (12,*) 'freeing element ',label
!      write (12,*) 'stak%bottom =',stak%bottom
!      write (12,*) 'stak%element(1) =',stak%element(1)
!      write (12,*) 'use (label)? ',use(label)
      if (label.eq.0) then
         write (*,*) 'error in freelab.f'
         write (*,*) 'trying to free label 0'
         stop
      endif
      stak%bottom = stak%bottom + 1
      stak%element(stak%bottom) = label
      use(label) = .false.
!
      return
      end
! --------------------------------------------------------------------

