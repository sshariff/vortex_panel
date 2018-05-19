! -------------------------------------------------------------------
      subroutine lnvortx (xp,x1,x2,u)
!
!  This subroutine calculates the velocity induced at a point
!  P  by a unit strength vortex line segment from point 1
!  to point 2.  The coordinates of the nodes are in xp, x1, 
!  and x2, respectively. Gamma is circulation, and the induced
!  velocity is output in U. After Katz, 1991.  
!
!--------------------------------------------------------------
!
      use interfish4
!
      implicit none
!
      real, dimension(3), intent(in)  :: xp,x1,x2
      real, dimension(3), intent(out) :: u

      real    :: eps
      real    :: k
      real    :: magr1,magr2,magx
      integer :: n
      real    :: pi
      real    :: r0r1,r0r2
      real    :: r1(3),r2(3)
      real    :: r1xr2(3)
!--------------------------------------------------------------
! input variables
!        x1, x2: coordinates of nodes at each end of the segment
!                whose influence is being computed
!        xp: coordinate of collocation point
! ---------------------------------------------------------------
!  output variables
!        u: vector velocity (u,v,w) induced at point P
! ---------------------------------------------------------------
! local variables
!        eps: assumed radius of the vortex; a tolerance for treating 
!             point P as if it were on the segment 
!        k: coefficient for velocity calculation
!        magr1,magr2,magx: the magnitudes of the vectors r1,r2, & r1xr2
!        n: index for stepping
!        pi: pi
!        r0r1,r0r2: the dot products (r0.r1), and (r0.r2), 
!                   r0 is the length of the element (r0=r2-r1)
!        r1,r2: vectors pointing from the collocation point P to the
!               ends of the vortex line segment (points 1 and 2)
!        r1xr2: the cross product of r1 and r2
! --------------------------------------------------------------
! start:
!
! set tolerance
      eps = 1e-10
      pi = acos(-1.0)
!      
! calculate r1 and r2, and magnitudes
      r1 = x1 - xp
      r2 = x2 - xp
!
      magr1 = sqrt( r1(1)*r1(1) + r1(2)*r1(2) + r1(3)*r1(3) )
      magr2 = sqrt( r2(1)*r2(1) + r2(2)*r2(2) + r2(3)*r2(3) )
!      
! calculate r1 cross r2. magnitude
      call cross (r1, r2, r1xr2, magx)
!
! check if point P lies on or very close to the vortex
      if ((magr1.lt.eps).or.(magr2.lt.eps).or.(magx.lt.eps)) then
! if so, assume no self-induced velocity
!         write (*,*) 'no velocity, man.'
!	 write (*,*) 'magr1,magr2,magx = ',magr1,magr2,magx
!	 write (*,*) 'x1 = ',x1
!	 write (*,*) 'x2 = ',x2
!         write (*,*) 'xp = ',xp
         u = 0.0
      else
! otherwise, calculate dot-products, where r0=r2-r1=x2-x1
         r0r1 = dot_product((x2-x1),r1)
	 r0r2 = dot_product((x2-x1),r2)
!
! calculate the coefficient   
!
         k = -1.0/(4.0*pi*magx*magx)*(r0r1/magr1-r0r2/magr2)
         u = k*r1xr2

!         write (12,*) '                   r1 = ',r1
!	 write (12,*) '                   r2 = ',r2
!         write (12,*) '                   magr1 = ',magr1
!	 write (12,*) '                   magr2 = ',magr2
!	 write (12,*) '                   r1xr2 = ',r1xr2
!         write (12,*) '                   magx  = ',magx
!         write (12,*) '                   r0r1  = ',r0r1
!	 write (12,*) '                   r0r2  = ',r0r2
!	 write (12,*) '                   k     = ',k
!	 write (12,*) '                   u     = ',u
!
      endif
!
      return
      end
! --------------------------------------------------------------------

