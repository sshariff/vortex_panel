! -------------------------------------------------------------------
      subroutine nulabel (label,max,use,stak,nmax)
!
!  This subroutine determines which label to use for a new element,
!  which can be a panel, node, or edge.  When far panels are lumped
!  together, element numbers are freed and placed in a stack.
!  When new wake panels are created at the separation lines,
!  the stack is checked to see if there are any current element 
!  numbers not in use. If necessary, a new element is created, and 
!  the upper bound of element numbers is updated.
!  --------------------------------------------------------------
      use panel_type
!
      implicit none
!
      integer,      intent(out)   :: label
      integer,      intent(inout) :: max
      logical,      intent(inout) :: use(nmax)
      type (stack), intent(inout) :: stak
      integer,      intent(in)    :: nmax
! --------------------------------------------------------------
! input variables
!        nmax: maximum number of edges 
! --------------------------------------------------------------
! output variables
!        label: label number for new element
!        max: largest element number in use
!        stak: stack of free labels
!        use: lookup table: is element in use?
! --------------------------------------------------------------
! start:
!
! are there available numbers?
      if ((stak%bottom).ne.0) then
! pop last entry from stack
!         write (*,*) 'stak%bottom =',stak%bottom
!	 write (*,*) 'element(bottom) =',stak%element(stak%bottom)
         label = stak%element(stak%bottom)
	 if (label.eq.0) then
	    write (*,*) 'error in nulabel.f'
	    write (*,*) 'pulled label 0'
	    stop
	 endif
	 stak%bottom = stak%bottom - 1
	 use(label) = .true.
!         write (12,*) 'pulling freed label',label
      else 
! make new element
         max = max + 1
	 label = max
	 use(label) = .true.
      endif
      return
      end
! --------------------------------------------------------------------

