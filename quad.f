      module quadface
        interface
          subroutine quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                         & e1,e2,e3,e4,e5,p1,p2,    &
                         & bodyin,bnodex,body,bedges,top)
!
             use panel_type
!
             real,         dimension(3), intent(in)    :: x1,x2,x3,x4
             logical,      dimension(:), intent(inout) :: bodyin
             type (nodes), dimension(:), intent(inout) :: bnodex
             type (panel), dimension(:), intent(inout) :: body
             type (edges), dimension(:), intent(inout) :: bedges
	     logical,                    intent(in)    :: top
             integer, intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2
	     
          end subroutine quad
        end interface
      end module quadface
! -------------------------------------------------------------------
! -------------------------------------------------------------------
       subroutine quad (x1,x2,x3,x4,n1,n2,n3,n4, &
                      & e1,e2,e3,e4,e5,p1,p2,    &
                      & bodyin,bnodex,body,bedges,top)
!
       use panel_type
!
       implicit none

       real,         dimension(3), intent(in)    :: x1,x2,x3,x4
       logical,      dimension(:), intent(inout) :: bodyin
       type (nodes), dimension(:), intent(inout) :: bnodex
       type (panel), dimension(:), intent(inout) :: body
       type (edges), dimension(:), intent(inout) :: bedges
       logical,                    intent(in)    :: top
       integer, intent(in) :: n1,n2,n3,n4,e1,e2,e3,e4,e5,p1,p2

       bnodex(n1)%x = x1
       bnodex(n2)%x = x2
       bnodex(n3)%x = x3
       bnodex(n4)%x = x4
       
       bedges(e1)%node1 = n1
       bedges(e1)%node2 = n3
       bedges(e2)%node1 = n3
       bedges(e2)%node2 = n4
       bedges(e3)%node1 = n2
       bedges(e3)%node2 = n4
       bedges(e4)%node1 = n1
       bedges(e4)%node2 = n2
       bedges(e5)%node1 = n4
       bedges(e5)%node2 = n1
       
       if (bedges(e1)%panel1==0) then
          bedges(e1)%panel1 = p1
       else
          bedges(e1)%panel2 = p1
       end if
       
       if (bedges(e2)%panel1==0) then 
          bedges(e2)%panel1 = p1
       else
          bedges(e2)%panel2 = p1
       end if

       if (bedges(e3)%panel2==0) then
          bedges(e3)%panel2 = p2
       else
          bedges(e3)%panel1 = p2
       end if
       
       if (bedges(e4)%panel2==0) then
          bedges(e4)%panel2 = p2
       else 
          bedges(e4)%panel1 = p1
       endif

       bedges(e5)%panel1 = p2
       bedges(e5)%panel2 = p1
       
       body(p1)%node(1) = n1
       body(p1)%node(2) = n3
       body(p1)%node(3) = n4
       body(p1)%edge(1) = e1
       body(p1)%edge(2) = e2
       body(p1)%edge(3) = e5
       
       body(p2)%node(1) = n1
       body(p2)%node(2) = n2
       body(p2)%node(3) = n4
       body(p2)%edge(1) = e4
       body(p2)%edge(2) = e3
       body(p2)%edge(3) = e5
       
       if (top) then
          bodyin(p1) = .true.
	  bodyin(p2) = .false.
       else
          bodyin(p1) = .false.
	  bodyin(p2) = .true.
       endif
!
!       write (*,*) 'quad...'
       
       end subroutine quad
