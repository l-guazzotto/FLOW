module interpolating_functions

use constant

contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine splinderiv(xv,yv,yv2,x,yp,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		integer n
		real (kind=dkind) :: x,yp	! point and first derivative
		real (kind=dkind), dimension(1:n) :: xv,yv,yv2	!table of x,y,d2y/dx2
		integer j
		real (kind=dkind) :: a,b

		call locate(xv,n,x,j)

		if ( (j>=n).OR.(j<=1) ) then

			yp = 0.d0

		else

			a = (xv(j+1)-x) / (xv(j+1)-xv(j))
			b = 1.d0-a

			yp = (yv(j+1)-yv(j)) / (xv(j+1)-xv(j)) -  &
					(3.d0*a**2-1.d0)/6.d0 * (xv(j+1)-xv(j)) * yv2(j) +  &
					(3.d0*b**2-1.d0)/6.d0 * (xv(j+1)-xv(j)) * yv2(j+1) 
	
		endif

		return

		end subroutine splinderiv


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine spline(x,y,n,yp1,ypn,y2)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
	  integer, parameter :: nmax=2000
	  integer n
	  real (kind=dkind) x,y,y2,u,yp1,ypn
      dimension x(n),y(n),y2(n),u(nmax)
	  real (kind=dkind) sig,p,qn,un
	  integer i,k

      if (yp1.gt..99e30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue

      return

      end subroutine spline

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	  subroutine splint(xa,ya,y2a,n,x,y)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     
	  integer n
	  real (kind=dkind) xa,ya,y2a,x,y
	  dimension xa(n),ya(n),y2a(n)
	  integer klo,khi,k
	  real (kind=dkind) h,a,b

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+  &
           ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      
	  return
      
	  end subroutine splint

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      subroutine locate(xx,n,x,j)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	  integer n,j
	  real(kind=dkind) x
	  real (kind=dkind),dimension(1:n) :: xx
	  integer jl,ju,jm

      jl=0
      ju=n+1

10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif

      j=jl

      return

      end subroutine locate

	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine integrate(nw,fun,w_1D,answer)
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	! executes 1D numerical integration of the function fun
	! NOTE Jocabian is included in fun!

		use constant, only : dkind

		implicit none

		integer :: nw
		real (kind=dkind), dimension (1:nw) :: fun,w_1D
		real (kind=dkind) :: answer
		integer :: i

		answer = 0.d0

		do i=1,nw

			answer = answer+w_1D(i)*fun(i)

		enddo

		continue

		return

	end subroutine integrate

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine set_weights(nw,xl,xr,w_1D,csi_1D) 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! sets the weights for the numerical quadrature formula

	implicit none

	integer :: nw
	real(kind=dkind), dimension(1:nw) :: w_1D,csi_1D
	integer :: i,j
	real (kind=dkind) :: xl, xr
	real (kind=dkind) :: a,b,AA,BB

	call gauleg(xl,xr,csi_1D,w_1D,nw)

	continue

	return

end subroutine set_weights

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine gauleg(x1,x2,x_1D,w_1D,nw)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
	  implicit none
      
      real(kind=dkind) :: eps=3.d-15
	  integer :: nw
	  real(kind=dkind) :: x1,x2	! limits of integration
	  real(kind=dkind), dimension(1:nw) :: x_1D,w_1D	! location of 
					!integration points and corresponding weights
	  integer :: m,i,j
	  real(kind=dkind) :: xl,xm,z,p1,p2,p3,z1,pp

      m=(nw+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(pi*(i-.25d0)/(nw+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,nw
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=nw*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x_1D(i)=xm-xl*z
        x_1D(nw+1-i)=xm+xl*z
        w_1D(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w_1D(nw+1-i)=w_1D(i)
12    continue
      return
 
end subroutine gauleg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function brent(ax,bx,cx,f,tol,xmin)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	! calculates the minimum of the input function f
	! the minimum point is returned as xmin, the minimum value as output of the function

	  implicit none

      integer, parameter :: itmax=100
	  real(dkind), parameter :: cgold=0.3819660112501051d0
	  real(dkind), parameter :: zeps=1.0d-10
	  real(dkind) :: brent

	  real(dkind) :: a, b, ax, bx, cx, tol, xmin
	  real(dkind) :: v, w, x, e, fx, fv, fw, xm, tol1, tol2, r, q, p, etemp, d, u, fu

	  real(dkind) :: f

	  integer :: iter

	  external f

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,itmax
        xm=0.5*(a+b)
        tol1=tol*abs(x)+zeps
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=cgold*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations.'
3     xmin=x
      brent=fx

      return

end function brent

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine indexx(n,arrin,indx)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	integer :: n
	integer, dimension(1:n) :: indx
	real(kind=dkind), dimension(1:n) :: arrin

	integer :: i, j, l, ir, indxt
	real(kind=dkind) :: q

      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
 
end subroutine indexx

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine interp_setup(npoints,int_ord,xdata,ydata,data3,cscoef)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBSNAK, DBSINT

	implicit none

	integer :: int_ord	! order of the interpolation
	integer :: npoints	! dimension of the input data

	real(kind=dkind) :: xdata(npoints)	! the "x" data
	real(kind=dkind) :: ydata(npoints)	! the "y" data
	real(kind=dkind) :: data3(npoints+int_ord)	! the "x" points in the interpolation
	real(kind=dkind) :: cscoef(npoints)	! the interpolation coefficients

	call DBSNAK(npoints, xdata, int_ord,data3)

	call DBSINT(npoints, xdata, ydata, int_ord, data3, cscoef)

	continue

	return

end subroutine interp_setup



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine lin_interp_2D(thing,xQ,zQ,answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! important: thing grid order is (i,j(, (i+1,j), (i+1,j+1), (i,j+1)

    real(kind=dkind), dimension(1:4) :: thing
    real(kind=dkind), dimension(1:4) :: fi
	real(kind=dkind) :: xQ, zQ, answer

	fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
	fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
	fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
	fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

	answer = thing(1)*fi(1) + thing(2)*fi(2) +  &
					thing(3)*fi(3) + thing(4)*fi(4) 

	continue

	return

end subroutine lin_interp_2D


end module interpolating_functions