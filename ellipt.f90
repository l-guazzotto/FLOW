!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION rf(x,y,z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	REAL(kind=dkind) ::  rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
	PARAMETER (ERRTOL=.0025d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,  &
				C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
	REAL (kind=dkind) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

	if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,  &
				z).gt.BIG)pause 'invalid arguments in rf'

	xt=x
	yt=y
	zt=z

	1     continue

		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		xt=.25d0*(xt+alamb)
		yt=.25d0*(yt+alamb)
		zt=.25d0*(zt+alamb)
		ave=THIRD*(xt+yt+zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave

	if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
	e2=delx*dely-delz**2
	e3=delx*dely*delz
	rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
	return

END FUNCTION rf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION rd(x,y,z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


	use constant, only : dkind

	implicit none

	REAL(kind=dkind) :: rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
	PARAMETER (ERRTOL=.0015d0,TINY=1.d-25,BIG=4.5d21,C1=3.d0/14.d0,  &
			C2=1.d0/6.d0,C3=9.d0/22.d0,C4=3.d0/26.d0,C5=.25d0*C3,C6=1.5d0*C4)
	REAL(kind=dkind) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,  &
					sqrtx,sqrty,sqrtz,sum,xt,yt,zt
	if(min(x,y).lt.0.d0.or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG)  &
				pause 'invalid arguments in rd'

	xt=x
	yt=y
	zt=z
	sum=0.d0
	fac=1.d0

	1     continue

		sqrtx=sqrt(xt)
		sqrty=sqrt(yt)
		sqrtz=sqrt(zt)
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
		sum=sum+fac/(sqrtz*(zt+alamb))
		fac=.25d0*fac
		xt=.25d0*(xt+alamb)
		yt=.25d0*(yt+alamb)
		zt=.25d0*(zt+alamb)
		ave=.2d0*(xt+yt+3.d0*zt)
		delx=(ave-xt)/ave
		dely=(ave-yt)/ave
		delz=(ave-zt)/ave

	if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1

	ea=delx*dely
	eb=delz*delz
	ec=ea-eb
	ed=ea-6.d0*eb
	ee=ed+ec+ec
	rd=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*  &
			ec+delz*C4*ea)))/(ave*sqrt(ave))

	return

END FUNCTION rd


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION elle(phi,ak)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!    USES rd,rf

	use constant, only : dkind

	implicit none

	REAL(kind=dkind) :: elle,ak,phi
	REAL(kind=dkind) :: cc,q,s,rd,rf

	s=sin(phi)
	cc=cos(phi)**2
	q=(1.d0-s*ak)*(1.+s*ak)
	elle=s*(rf(cc,q,1.d0)-((s*ak)**2)*rd(cc,q,1.d0)/3.d0)

	return

END FUNCTION elle

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FUNCTION ellf(phi,ak)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!    USES rf

	use constant, only : dkind

	implicit none

	REAL(kind=dkind) :: ellf,ak,phi
	REAL(kind=dkind) :: s,rf

	s=sin(phi)
	ellf=s*rf(cos(phi)**2,(1.d0-s*ak)*(1.d0+s*ak),1.d0)

	return

END FUNCTION ellf

