! --------------------------------------------------------------------
      subroutine solve (adim,bpnum,kuttas,a,right,body,wake,wedges, &
                      & nseplines,oldcirc,seplist,nmax,bedges,      &
		      & bnodex,wnodex)
!
!  This subroutine solves for the unknown circulations. a(i,j) is
!  inverted, then multiplied by the right hand side.
!
! -----------------------------------------------------------------
      use panel_type; use interfish2; use interfish3
!
      implicit none 
!
!--------------------      
!      interface
!      subroutine svdcmp (a,m,n,mp,np,w,v)
!      INTEGER, intent(in)    :: m,mp,n,np  
!      REAL,    intent(inout) :: a(mp,np)
!      REAL,    intent(out)   :: v(n,n),w(n)
!      end subroutine svdcmp
!      end interface
! -----------------------------------------------------------------
!      interface
!      function pythag (a,b)
!      REAL, intent(in)  :: a,b
!      REAL, intent(out) :: pythag
!      end function pythag
!      end interface
! -----------------------------------------------------------------
!      interface
!      subroutine svbksb (u,w,v,m,n,mp,np,b,x)
!      INTEGER, intent(in)  :: m,mp,n,np
!      REAL,    intent(in)  :: b(m),u(mp,np),v(n,n),w(n)
!      REAL,    intent(out) :: x(n)
!      end subroutine svbksb
!      end interface
! -----------------------------------------------------------------
      integer,                            intent(in)    :: adim
      integer,                            intent(in)    :: bpnum
      integer,                            intent(in)    :: kuttas
      real,         dimension(:,:),       intent(in)    :: a
      real,         dimension(nmax),      intent(in)    :: right
      type (panel), dimension(nmax),      intent(inout) :: body, wake
      type (edges), dimension(nmax),      intent(in)    :: wedges
      integer,                            intent(in)    :: nseplines
      real,         dimension(nmax),      intent(inout) :: oldcirc
      type (splst), dimension(25),        intent(in)    :: seplist
      integer,                            intent(in)    :: nmax
      type (edges), dimension(nmax),      intent(in)    :: bedges
      type (nodes), dimension(nmax),      intent(in)    :: bnodex, wnodex

      real                     :: d
      integer, dimension(adim) :: indx
      integer                  :: panl,row,line,edg,wp,newsize
        
      real, dimension(adim,adim) :: a1     
      real, dimension(adim)      :: localr  
!      real, dimension(bpnum+kuttas)       :: w,svdout
!      real, dimension(bpnum+kuttas,bpnum+kuttas) :: v
!      real :: wmax,wmin

       real :: qdum 
       real,    dimension(3)      :: wind,q,coll,x1,x2,qn
       integer, dimension(kuttas) :: wakemap
       integer                    :: k,l,wakepan
! ---------------------------------------
!   input variables
!        a: matrix of influence coefficients 
!        kuttas: the number of kutta panels being solved for
!        nmax: maximum number of edges (for memory allocation)
!        nseplines: number of separation lines
!        right: rhs for matrix inversion, i.e. Ax=B. right = B,
!               but the solution subroutines replace right with x,
!               which is the matrix of circulations
!        seplist: list of separation edges
!        bpnum+kuttas: size of the matrix being inverted, includes body panels
!              plus the equations for the initial set of wake panels.
!              there are two wake panels for each extra equation
!        wedges: edge-node arrays, lists nodes that correspond to
!                each edge
! ---------------------------------------
!   local variables        
!        d: pivoting statistic from ludcmp
!        gam: panel strength
!        indx: permutation vector returned by ludcmp, used in lubksb
!        panl,row,line,edg,wp: stepping indices for various elements
! ---------------------------------------
!   output variables
!        body, wake: body and wake panels, contains strength of 
!                    panel, lists of nodes & edges, "polarity",
!                    normal direction and centroid
! ---------------------------------------
!   start:
!
        localr(1:adim) = right(1:adim)
        a1 = a
!
!	write (12,*) 'a''s...'

!!        write (12,'(3x,26i14)') ((l),l=1,bpnum+kuttas)
!        write (12,'(i12,25i14)') ((l),l=1,bpnum+kuttas)
!        do 1167 k = 1,bpnum+kuttas
!          write (12,'(i3,26g14.4)') k,((a(k,l)), l=1,bpnum+kuttas)
! 1167   continue
!        write (12,*)

!        go to 1172

!        write (12,'(i12,15i13)') ((l),l=1,16)
!        do 1167 k = 1,bpnum+kuttas
!          write (12,'(i3,16g13.3)') k,((a(k,l)), l=1,16)
! 1167   continue
!        write (12,*)
!        write (12,'(i12,16i13)') ((l),l=17,32)
!        do 1168 k = 1,bpnum+kuttas
!          write (12,'(i3,16g13.3)') k,((a(k,l)), l=17,32)
! 1168   continue
 
! 1172   continue
	
!	write (12,*)
!        write (12,*) 'rhs...'
!        do 6169 k = 1,bpnum+kuttas
!           write (12,'(i3,e22.12)') k,right(k)
! 6169   continue
!	write (12,*)
	
	
!	go to 4782
	
	
!        call svdcmp (a1,bpnum+kuttas,bpnum+kuttas,nmax,nmax,w,v) 
	
	
!	wmax = 0.
!	do 9998 k=1,bpnum+kuttas
!	   if(w(k).gt.wmax) wmax = w(k)
!  9998  continue
!        wmin = wmax*1e-6
!	do 9994 k=1,bpnum+kuttas
!	   if(w(k).lt.wmin) w(k) = 0.0
!  9994  continue
!	call svbksb(a1,w,v,bpnum+kuttas,bpnum+kuttas,nmax,nmax,localr,svdout)
	
! 4782   continue
 
 
        call ludcmp (a1,bpnum+kuttas,adim,indx,d)
	call lubksb (a1,bpnum+kuttas,adim,indx,localr)
	
!	write (12,*) 'soln ...'
!        do 1267 k = 1,bpnum+kuttas
!!          write (12,'(i3,e22.12)') k,svdout(k)
!          write (12,'(i3,e22.12)') k,localr(k)
! 1267   continue
!        write (12,*)
	
!	sing = matmul(a(1:bpnum+kuttas,1:bpnum+kuttas),svdout(1:bpnum+kuttas))
!	sing = matmul(a(1:bpnum+kuttas,1:bpnum+kuttas),localr(1:bpnum+kuttas))
	
!        write (12,*) 'checking ax=b...'
!        do 6269 k = 1,bpnum+kuttas
!           write (12,'(i3,e22.12)') k,sing(k)
! 6269   continue
!	write (12,*)
	
!	write (12,*)
!	
!	goto 2456

!        write(12,*)
!	write(12,*) 'double precision.........'
 
 
! 1113   continue
 
 
!        go to 11113
		
!        write (12,*)
!	write (12,*) 'in solve....'
!	write (12,*)
!	
!	write (12,*) 'a...'
!        write (12,'(3x,22i11)') ((l),l=1,bpnum+kuttas)
!        do 1667 k = 1,bpnum+kuttas
!          write (12,'(i3,22e22.12)') k,((a(k,l)), l=1,bpnum+kuttas)
! 1667   continue
!        write (12,*)
!        write (12,*) 'rhs...'
!        do 6169 k = 1,bpnum+kuttas
!           write (12,'(i3,22e22.12)') k,localr(k)
! 6169   continue
	
	
!
!   copy old values of circulation
       oldcirc(1:bpnum+kuttas) = body(1:bpnum+kuttas)%circ

       do 2600 panl = 1,bpnum
!          body(panl)%circ = right(panl)
!          body(panl)%circ = svdout(panl)       
          body(panl)%circ = localr(panl)
 2600  continue

!         do 766 edg=1,bpnum
!            write (12,765) 'body (',edg,') = ',body(edg)%circ
! 765        format(a,i4,a,f16.6)
! 766     continue


       if (kuttas.ne.0) then 
!
!   we have to "unroll" the wake panels from the separation edge list,
!   as in aij1.f  
          row = bpnum
!	  write(12,*) 'body pans =',row
          do 4000 line = 1,nseplines
	     do 3999 edg = 1,seplist(line)%nedges
	        row = row + 1
!		write(*,*) 'on row ',row
!   get corresponding wake number
		wp = wedges(seplist(line)%sepedge(edg)%wsedge)%panel1
!!!!!

                wakemap(row)   = wp
                wakemap(row+1) = wp+1

!!!!!   set near panel
!		wake(wp)%circ = right(row)
!		wake(wp)%circ = svdout(row)
!
		wake(wp)%circ = localr(row)
!   set far panel
!		wake(wp+1)%circ = right(row+1)
!		wake(wp+1)%circ = svdout(row+1)
!
		wake(wp+1)%circ = localr(row+1)
	        row = row + 1
!
!		write(12,*)'row =',row
!		write(12,*)'near panel =',wp
!		write(12,*)'far panel =',wp+1
!		write(12,*)
 3999        continue
 4000     continue
!
       endif


!!  checking induced velocities
!       do 91 line = 1,bpnum

!          wind = 0.0
!	  coll = body(line)%centr
!!  loop through panels 
!          do 90 edg = 1,bpnum+kuttas
!	     do 89 row = 1,3
!	     
!	        if (edg.le.bpnum) then
!		   x1 = bnodex(bedges(body(edg)%edge(row))%node1)%x
!		   x2 = bnodex(bedges(body(edg)%edge(row))%node2)%x
!	        else
!		   wakepan = wakemap(edg)
!		   x1 = wnodex(wedges(wake(wakepan)%edge(row))%node1)%x
!		   x2 = wnodex(wedges(wake(wakepan)%edge(row))%node2)%x
!		endif
!        	call lnvortx (coll,x1,x2,q)
!	        wind = wind + q*localr(edg)
!                
! 89	     continue
!	  
! 90       continue
!          write (12,905) 'wind at panel ',line,' = (',wind,')'
!          wind(1) = wind(1) + 5.0
!	  write (12,905) 'qtot at panel ',line,' = (',wind,')'
!	  qn = body(line)%normal*(dot_product(wind,body(line)%normal))
!	  write (12,905) 'qn   at panel ',line,' = (',qn,')'
!	  write (12,905) 'qt   at panel ',line,' = (',wind-qn,')'
!          write (12,*) '   qsq/vsq  = ', dot_product(wind,wind)/25.0
!	  write (12,*)

 905      format (a,i3,a,3g16.6,a)

 91    continue


       return
       end
! --------------------------------------------------------------------
      SUBROUTINE svdcmp (a,m,n,mp,np,w,v)
      INTEGER, intent(in)    :: m,mp,n,np  
      REAL,    intent(inout) :: a(mp,np)
      REAL,    intent(out)   :: v(n,n),w(n)
!,NMAX
!      PARAMETER (NMAX=3000)
!     USES pythag
      INTEGER i,its,j,jj,k,l,nm
      REAL anorm,c,f,g,h,s,scayle,x,y,z,rv1(n),pythag
      g=0.0
      scayle=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scayle*g
        g=0.0
        s=0.0
        scayle=0.0
        if(i.le.m)then
          do 11 k=i,m
            scayle=scayle+abs(a(k,i))
11        continue
          if(scayle.ne.0.0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scayle
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scayle*a(k,i)
16          continue
          endif
        endif
        w(i)=scayle *g
        g=0.0
        s=0.0
        scayle=0.0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scayle=scayle+abs(a(i,k))
17        continue
          if(scayle.ne.0.0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scayle
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scayle*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0
33      continue
        if(g.ne.0.0)then
          g=1.0/g
          do 36 j=l,n
            s=0.0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((abs(rv1(l))+anorm).eq.anorm)  goto 2
            if((abs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((abs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=pythag(f,1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0)then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
! --------------------------------------------------------------------
      FUNCTION pythag (a,b)
      REAL, intent(in)  :: a,b
      REAL, intent(out) :: pythag
      REAL absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.
        else
          pythag=absb*sqrt(1.+(absa/absb)**2)
        endif
      endif
      return
      END
! --------------------------------------------------------------------
      SUBROUTINE svbksb (u,w,v,m,n,mp,np,b,x)
      INTEGER, intent(in)  :: m,mp,n,np
      REAL,    intent(in)  :: b(m),u(mp,np),v(n,n),w(n)
      REAL,    intent(out) :: x(n)
!,NMAX
!      PARAMETER (NMAX=3000)
      INTEGER i,j,jj
      REAL s,tmp(n)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
! --------------------------------------------------------------------
