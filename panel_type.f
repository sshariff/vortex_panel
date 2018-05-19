! -------------------------------------------------------------------
      module panel_type
! defined data structures
! --------------------------------------------------------------
! note: for the wake, the edge, panel, and node numbers are 
!       reused. available numbers which have been freed due
!       to vortex lumping are listed in the stacks. the stack
!       data structure has an index, which is the label of the
!       last stack entry used.  
! -------------------------------------------------------------
! stack: label numbers freed by lumping.
!        this should be about 10 times the number of
!        separation lines times the number of edges per
!        separation line (because 9 edges are freed with
!        each lumping)
!    ----------------------
      type stack
        integer      bottom
	integer      element(2500)
      end type stack
! -------------------------------------------------------------
! nodes: this data structure contains the coordinates of the
!        body nodes
!   -----------------------
      type nodes
        real         x(3)
      end type nodes
! -------------------------------------------------------------
! edges: this data structure contains the nodes associated
!        with each edge
!   -----------------------
      type edges
        integer      node1
	integer      node2
	integer      panel1
	integer      panel2
      end type edges
! --------------------------------------------------------------
! panel: this data structure contains all the information for
!        a given panel
!  -------------------------
      type panel
! strength of circulation
	real         circ
! nodes 
	integer      node(3)
! edges
	integer      edge(3)
! normal direction
	real         normal(3)
! centroid coordinates
	real         centr(3)
      end type panel
! -------------------------------------------------------------
! sepedges: this data structure contains the nodes and edges
!          corresponding to a separation edge. these are 
!          described by two coincident edges, the body-side edge
!          and the wake-side edge. the body-side edge is a body
!          panel edge, is bound to the body and does not change.
!          the wake-side edge is a wake panel edge, and is shed
!          from the separation line at each timestep.
      type sepedges
! body-side edge
        integer       bsedge
! body-side nodes
        integer       bsnode(2)
! wake-side edge
        integer       wsedge
! wake-side nodes
        integer       wsnode(2)     
      end type sepedges
!  -------------------------------------------------------------
! splsts: this data structure contains (for each separation line
!         list) the number of edges which make up the separation
!         line, and an ordered list of these separation edges. 
!         seplist is hardwired to handle only 50 separation edges 
!         per line. 
      type splst
! number of separation edges
        integer nedges
! separation edges
        type (sepedges) sepedge(100)	
      end type splst
!  -------------------------------------------------------------
      end module panel_type
! --------------------------------------------------------------------
