! -------------------------------------------------------------------
      subroutine blchar (bpnum,body,cp,bldelt,seplist,nseplines,nmax)
!
!  This subroutine determines the boundary layer characteristics
!  for the current timestep (bl displacement, separation locations,
!  etc)
!  --------------------------------------------------------------
      use panel_type
!
      implicit none
!
      integer,                       intent(in)    :: bpnum
      type (panel), dimension(nmax), intent(in)    :: body
      real,         dimension(nmax), intent(in)    :: cp
      real,         dimension(nmax), intent(out)   :: bldelt
      type (splst), dimension(25),   intent(inout) :: seplist
      integer,                       intent(inout) :: nseplines
      integer,                       intent(in)    :: nmax
!
      integer :: pan
! --------------------------------------------------------------
! input variables
!        body: body panels, contains strength of 
!              panel, lists of nodes & edges, "polarity",
!              normal direction and centroid
!        bpnum: total number of panels used to describe body
!        cp: pressure coefficient for each panel
!        nmax: maximum number of edges 
! --------------------------------------------------------------
! local variables
!        pan: index for panel number
! --------------------------------------------------------------
! output variables
!        bldelt: boundary layer displacement for each panel
!        nseplines: number of separation lines
!        seplist: list of separation edges
! --------------------------------------------------------------
! start:
!
      do 7000 pan = 1,bpnum
         bldelt(pan) = 0.0
 7000 continue
      return
      end
! --------------------------------------------------------------------

