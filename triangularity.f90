module triangularity

	use array_dimensions
	use constant
	use interpolating_functions

	implicit none

	integer, parameter :: n_tri_lim = 20
	integer theta_points,theta_points1,theta_points2
	integer, parameter :: theta_points_max_lim = 1250
	! tri_type = 0 -> no triangularity: circular or elliptic boundary
	! tri_type = 1 -> r=Sum[r(k)*cos(k*theta)]
	! tri_type = 2 -> r**2=r**2*{k**2*sin**2(theta)+cos**2[theta+asin(delta)*sin(theta)]}
	! tri_type = 3 -> r**2=a**2/Sqrt{sin**2(theta)/k**2+cos**2(theta)-
	!					(delta/4)*[cos(3*theta)-cos(theta)]}
	! tri_type = 4 -> numerical input
    ! tri_type = 5 -> r=Sum[rc(k)*cos(k*theta)+rs(k)*sin(k*theta)]
	! tri_type = 6 -> R=Sum[Rc(k)*cos(k*theta)+Rs(k)*sin(k*theta)]
	!				  Z=Sum[Zc(k)*cos(k*theta)+Zs(k)*sin(k*theta)]
	!				  r=[(R-R0)**2 + Z**2 ]**0.5
	! tri_type = 8 -> numerical input, IMSL periodic spline
	! tri_type = 18 -> "standard" delta, k input
	! tri_type = 13 -> two surfaces numerical input, for vacuum

	! theta_points1 and theta_points2 are used when two definitions of the boundary are needed
	! e.g. for a dipole equilibrium or for a case with X-points





	! for triangularity = 0:
	! if(x^2/a^2+y^2/b^2) >=1 then psi=0.
	! Where x=R-R_0 and y=Z-Z_0 with R_0 and Z_0 being 
	! the center of your square box.
	real (kind=dkind) :: a_elps	
	real (kind=dkind) :: b_elps	
	real(kind = dkind),dimension(0:n_tri_lim) :: rcoeff_u
	real(kind = dkind),dimension(0:n_tri_lim) :: rcoeff_d
	real(kind = dkind),dimension(0:n_tri_lim) :: rcos_u
	real(kind = dkind),dimension(0:n_tri_lim) :: rcos_d
	real(kind = dkind),dimension(0:n_tri_lim) :: rsin_u
	real(kind = dkind),dimension(0:n_tri_lim) :: rsin_d
	real(kind = dkind),dimension(1:theta_points_max_lim) :: theta_dat
	real(kind = dkind),dimension(1:theta_points_max_lim) :: rminor_dat
	real(kind = dkind),dimension(1:theta_points_max_lim) :: d2rminor_dat
	real (kind=dkind), parameter :: ellipt=2.1d0
	real (kind=dkind) :: asin_d_up 
	real (kind=dkind) :: asin_d_down 
	real (kind=dkind), parameter :: smallr0 = 0.59d0		!0.62d0		!
	real (kind=dkind), parameter :: delta_u_o4 = 0.35d0/4.d0		!0.075d0		!
	real (kind=dkind), parameter :: delta_d_o4 = 0.54d0/4.d0		!0.075d0		!
	real(kind = dkind),dimension(0:n_tri_lim) :: bigR_cos
	real(kind = dkind),dimension(0:n_tri_lim) :: Z_cos
	real(kind = dkind),dimension(0:n_tri_lim) :: bigR_sin
	real(kind = dkind),dimension(0:n_tri_lim) :: Z_sin
	real(kind=dkind), dimension(:,:), allocatable :: r_data
	real(kind=dkind), dimension(:,:), allocatable :: r_cscoef
	integer :: r_ord = 3 ! interpolation order

	real(kind=dkind) :: delta_up, delta_down, k_ellipt

!	data rcoeff / .59d0 , -0.0d0 , 0.d0 , 0.0d0 , 0.0d0 , -0.0d0/
!	data rcoeff_u / .9145d0 , -0.07965d0 , -0.3245d0 , 0.07965d0 , 0.0d0 , -0.0d0/
!	data rcoeff_d / .9145d0 , -0.051625d0 , -0.3245d0 , 0.051625d0 , 0.0d0 , -0.0d0/
!	data rcoeff_u / .59d0 , -0.0d0 , 0.15d0 /
!	data rcoeff_d / .59d0 , -0.0d0 , 0.15d0 /

	real (kind = dkind) :: theta

	integer, dimension(:,:), allocatable :: ind_bound 
	! indeces of the external point and internal point surroundings (point "1" only)
	real(kind=dkind), dimension(:,:), allocatable :: coord_bound
	! internal point coordinates
	real(kind=dkind) :: R_P, z_P
	real(kind=dkind) :: radius_ext ! for tri_type=10
	real(kind=dkind) :: theta_temp
	real(kind=dkind), dimension(:,:), allocatable :: bound_13	! for extrapolation near the boundary (use only one array for simplicity)

	real(kind=dkind), dimension(:,:,:), allocatable :: bound_tri
	real(kind=dkind), dimension(:,:), allocatable :: coord_tri

	real(kind=dkind) :: volume
	real(kind=dkind) :: area
	real(kind=dkind) :: R_edge_min, R_edge_max
	real(kind=dkind) :: z_edge_min, z_edge_max

	real(kind=dkind) :: shape_ellipticity


	contains

	!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine set_triangularity(asin_d_up, asin_d_down,   &
						theta_dat,rminor_dat,bound1,bound2,d2rminor_dat,  &
						rcos_u,rcos_d,rsin_u,rsin_d,  &
						bigR_cos,bigR_sin,Z_cos,Z_sin,fname  )
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) asin_d_up
	real(kind=dkind) asin_d_down
	real (kind=dkind) :: theta_dat(1:theta_points_max),rminor_dat(1:theta_points_max),  &
						d2rminor_dat(1:theta_points_max)
	real(kind=dkind) :: bound1,bound2
	real (kind=dkind), dimension(0:n_tri) :: rcos_u,rcos_d,rsin_u,rsin_d
	real (kind=dkind), dimension(0:n_tri) :: bigR_cos,bigR_sin,Z_cos,Z_sin
	character*7 :: fname

	if(input_EQDSK) return

	theta_points = theta_points_max

	if (tri_type==2) call init_triangularity(asin_d_up, asin_d_down)

	if (tri_type==4) then

		call init_r0(theta_dat,rminor_dat)
		call spline(theta_dat,rminor_dat,theta_points,bound1,bound2,d2rminor_dat)

	elseif (tri_type==5) then

        call init_r0_fourier(rcos_u,rcos_d,rsin_u,rsin_d)

	elseif (tri_type==6) then

        call init_R_Z(bigR_cos,bigR_sin,Z_cos,Z_sin,fname)

!	elseif (tri_type==7) then

	elseif (tri_type==8) then

		call init_r0_IMSL

	elseif (tri_type==11) then

		! LDX (hole) boundary
		call init_2rs

	elseif (tri_type==13) then

		! two boundaries (separatrix solution)
		call init_2rs !use same routine to generate the boundary, 
		! the rest is of course different

	elseif (tri_type==18) then

		call init_r0_standard

	endif

	continue

	end subroutine set_triangularity


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine init_triangularity(up,down)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) up
	real(kind=dkind) down

	up = dasin(.54d0)		!0.d0	!
	down = dasin(.35d0)		!0.d0	!
!	up = .54d0
!	down = .35d0

	end subroutine init_triangularity

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine init_r0(theta,r)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind) :: theta(1:theta_points_max),r(1:theta_points_max)
	integer i
!	real (kind=dkind) rloc

	open(15,file='r0_spline2.dat', & 
		status='old',action='read')

!	rloc=0.

	do i=1,theta_points_max
		read(15,*) theta(i),r(i)
	enddo

	theta_points = i-1

	continue

	end subroutine init_r0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine init_r0_fourier(rc_u,rc_d,rs_u,rs_d)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), dimension(0:n_tri) :: rc_u,rc_d,rs_u,rs_d
    integer i

    open(101,file='r_fourier.dat',status='old',action='read')

    do i=0,n_tri
           read(101,*) rc_u(i)
    enddo

    rs_u(0) = 0.d0

    do i=1,n_tri
          read(101,*) rs_u(i)
    enddo

    do i=0,n_tri
         read(101,*) rc_d(i)
    enddo

    rs_d(0) = 0.d0

    do i=1,n_tri
          read(101,*) rs_d(i)
    enddo

    close(101)

    end subroutine init_r0_fourier

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine init_R_Z(rc,rs,zc,zs,fname)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	 real (kind=dkind), dimension(0:n_tri) :: rc,rs,zc,zs
	 character*7 :: fname
	 integer :: fid = 33
	 integer i,j,k
	 real (kind=dkind), dimension(1:theta_points_max-2) :: th_temp,r_temp
	 real (kind=dkind) :: bigR_loc,Z_loc,smallr_loc,theta_loc,thetaRZ
	 real(kind=dkind) :: minim

! first read the Fourier components of the contour

	 open(fid,file=fname,status='old',action='read')

	 read(fid,*) rc(0)
	 read(fid,*) zc(0)

	 rs(0) = 0.d0
	 zs(0) = 0.d0

	 do i=1,n_tri

		read(fid,*) rc(i)
		read(fid,*) rs(i)
		read(fid,*) zc(i)
		read(fid,*) zs(i)

	 enddo

	 close(fid)

! then build a set of contour points

	 do i = 1,theta_points_max-2

		bigR_loc = 0.d0
		Z_loc = 0.d0
		thetaRZ = (i-1.d0)/(theta_points_max-1.d0)*2.d0*pi

		do k=0,n_tri
			bigR_loc = bigR_loc + rc(k)*dcos(k*thetaRZ)+rs(k)*dsin(k*thetaRZ)
			Z_loc = Z_loc + zc(k)*dcos(k*thetaRZ)+zs(k)*dsin(k*thetaRZ)
		enddo

		smallr_loc = dsqrt( (bigR_loc-rmajor)**2 + Z_loc**2 ) 

		if ((bigR_loc-rmajor)==0.d0) then
			theta_loc = pi/2.d0 * dsign(1.d0,Z_loc)
		else
			theta_loc = datan2(Z_loc,(bigR_loc-rmajor))
		endif

		! keep angles between -3pi/2 and pi/2
		if (theta_loc>pi/2.d0) theta_loc=theta_loc-2.d0*pi

		th_temp(i) = theta_loc
		r_temp(i) = smallr_loc

	 enddo

	 ! now ordinate the data

  	 do i=2,theta_points_max-1

		k=0
		minim=1.d3	 

		do j=1,theta_points_max-2

			if(th_temp(j)<minim) then

				minim = th_temp(j)
				k=j

			endif

		enddo

		theta_dat(i) = th_temp(k)
		rminor_dat(i) = r_temp(k)
		th_temp(k) = 1.d3

	 continue

	 enddo

	 ! finally, set the first (i.e. smallest angle) point to -3pi/2
	 theta_dat(1) = -1.5d0*pi
	 rminor_dat(1) = (rminor_dat(2)+rminor_dat(theta_points_max-1))/2.d0
	 ! and the last (i.e. largest angle) to pi/2
	 theta_dat(theta_points_max) = pi/2.d0
	 rminor_dat(theta_points_max) = rminor_dat(1)

	 !now we can build the spline representation of the contour
	 call spline(theta_dat,rminor_dat,theta_points_max,0.d0,0.d0,d2rminor_dat)


	 return

	end subroutine init_R_Z

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine init_r0_IMSL
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!	use IMSL, only : DCSPER, DCSINT
	use pseudo_IMSL, only : DBSNAK, DBSINT

	integer :: i
	real(kind=dkind), dimension(theta_points_max_lim,2) :: dummy

	i = 0

	open(33,file='roftheta.dat', status='old', action='read')

	do 

		i = i+1
		read(33,*,end=99) dummy(i,1),dummy(i,2)

	enddo

99	continue

	close(33)

	theta_points = i-1

	allocate(r_data(theta_points+r_ord,3))
	allocate(r_cscoef(1,theta_points))

	do i=1,theta_points

		r_data(i,1) = dummy(i,1)
		r_data(i,2) = dummy(i,2)

	enddo


	call DBSNAK(theta_points, r_data(1:theta_points,1),  &
					r_ord,r_data(:,3))

	call DBSINT(theta_points, r_data(1:theta_points,1),  &
		 r_data(1:theta_points,2), r_ord,  &
		 r_data(:,3),  &
		 r_cscoef(1,1:theta_points))


!!	call dcsper(theta_points,r_data(:,1),r_data(:,2),r_data(:,3),r_cscoef)
!	call dcsint(theta_points,r_data(:,1),r_data(:,2),r_data(:,3),r_cscoef)

	continue

	return

  end subroutine init_r0_IMSL

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine plasma_shape_conversion(R_b,Z_b)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this subroutine converts an input of the kind (R,Z) to a (theta,r) form
! the output is the same as from init_r0_IMSL

	use constant
	use pseudo_IMSL, only : DBSNAK, DBSINT, interp_check

	implicit none

	real(kind=dkind), dimension(theta_points) :: R_b, Z_b ! the boundary
	real(kind=dkind) :: rloc, thetaloc
	integer, allocatable, dimension(:) :: sorter, i_repeated
	real(kind=dkind), allocatable, dimension(:,:) :: data_temp
	real(kind=dkind) :: theta_min, theta_max
	real(kind=dkind) :: theta1, theta2, tol
	real(kind=dkind), dimension(3) :: th_search

	integer :: i, j
	integer :: ind1, ind2, i0, iskip

	tri_type = 8

	allocate(data_temp(theta_points+1,2))

	do i=1,theta_points

		thetaloc = atan2(Z_b(i),R_b(i)-rmajor)

		if (thetaloc<0.d0) thetaloc = thetaloc + 2.d0*pi

		rloc = sqrt( (R_b(i)-rmajor)**2 + Z_b(i)**2 )

		data_temp(i,1) = thetaloc
		data_temp(i,2) = rloc

	enddo

	! first check if the data range covers a circle

	theta_min = minval(data_temp(1:theta_points,1))
	theta_max = maxval(data_temp(1:theta_points,1))

	if(theta_max-theta_min<2.d0*pi) then

!!$		i = minloc(data_temp(1:theta_points,1))
!!$
!!$		data_temp(theta_points+1,2) = data_temp(i,2)

		i = 0

		do

			i = i+1

			if(theta_min == data_temp(i,1)) then

				data_temp(theta_points+1,2) = data_temp(i,2)
				exit

			endif

		enddo

		theta_max = theta_min + 2.d0*pi
		data_temp(theta_points+1,1) = theta_max

		theta_points = theta_points + 1

	endif

	! reorder the data

	allocate(sorter(theta_points))

	call indexx(theta_points, data_temp(1:theta_points,1), sorter)

	! now check for repeated entries
	! this can happen if the input has a starting point different from theta=0

	iskip = 0
	allocate(i_repeated(theta_points))

	do i = 2,theta_points

		if(data_temp(sorter(i),1)==data_temp(sorter(i-1),1)) then

			iskip = iskip + 1
			i_repeated(iskip) = i

		endif

	enddo

	theta_points = theta_points - iskip

	allocate(r_data(theta_points+r_ord,3))
	allocate(r_cscoef(1,theta_points))

	j = 1

	do i = 1,theta_points+iskip

		if(i==i_repeated(j)) then
		! skip this point

			j = j+1

		else

			r_data(i-j+1,1) = data_temp(sorter(i),1)
			r_data(i-j+1,2) = data_temp(sorter(i),2)

		endif

	enddo

	deallocate(i_repeated)
	deallocate(sorter)
	deallocate(data_temp)

	call interp_setup(theta_points,r_ord, &
		r_data(1:theta_points,1),r_data(1:theta_points,2), &
		r_data(1:theta_points+r_ord,3),r_cscoef(1,1:theta_points))

	continue

end subroutine plasma_shape_conversion

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine init_2rs
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this subroutine initializes the two boundary curves
  ! for non-simply-connected domains

	use pseudo_IMSL, only : DBSNAK, DBSINT

	integer :: i
	real(kind=dkind) :: dummy

! outer boundary

	i = 0

	open(33,file='r1.dat', status='old', action='read')

	do

		read(33,*,end=11) dummy,dummy
		i = i+1

	enddo

11	theta_points1 = i

	rewind(33)

! inner boundary

	i = 0

	open(34,file='r2.dat', status='old', action='read')

	do

		read(34,*,end=12) dummy,dummy
		i = i+1

	enddo

12	theta_points2 = i

	rewind(34)

	theta_points=max(theta_points1,theta_points2)

	allocate(r_data(theta_points+r_ord,6))
	!1-3: outer; 4-6: inner

	allocate(r_cscoef(2,theta_points))
	!1: outer; 2: inner

	! outer boundary interpolation

	do i=1,theta_points1

		read(33,*) r_data(i,1), r_data(i,2)

	enddo

	close(33)

	call DBSNAK(theta_points1, r_data(1:theta_points1,1),  &
					r_ord,r_data(1:theta_points1+r_ord,3))

	call DBSINT(theta_points1, r_data(1:theta_points1,1),  &
		 r_data(1:theta_points1,2), r_ord,  &
		 r_data(1:theta_points1+r_ord,3),  &
		 r_cscoef(1,1:theta_points1))

	! inner boundary interpolation

	do i=1,theta_points2

		read(34,*) r_data(i,4), r_data(i,5)

	enddo

	close(34)

	call DBSNAK(theta_points2, r_data(1:theta_points2,4),  &
					r_ord,r_data(1:theta_points2+r_ord,6))

	call DBSINT(theta_points2, r_data(1:theta_points2,4),  &
		 r_data(1:theta_points2,5), r_ord,  &
		 r_data(1:theta_points2+r_ord,6),  &
		 r_cscoef(2,1:theta_points2))

	continue

	return

  end subroutine init_2rs

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine init_r0_standard
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this initializes the radius using the standard k and delta stuff

	use pseudo_IMSL, only : DBSNAK, DBSINT

	real(kind=dkind) :: xloc, zloc, rloc, theta_fake, theta_loc, delta
	integer :: i

	theta_points = theta_points_max

	allocate(r_data(theta_points+r_ord,3))
	allocate(r_cscoef(1,theta_points))

	do i = 0,theta_points-1

		theta_fake = 2.d0*i*pi/(theta_points-1.d0)

		if(theta_fake<pi) then
			delta = delta_up
		else
			delta = delta_down
		endif

		xloc = a_elps * cos( theta_fake + delta * sin(theta_fake) )
		! + rmajor, not needed
		zloc = k_ellipt * a_elps * sin(theta_fake)
		rloc = ( xloc**2 + zloc**2 )**0.5d0

		theta_loc = atan2(zloc,xloc)

		if (theta_loc<0.d0) theta_loc = theta_loc + 2.d0*pi

		r_data(i+1,1) = theta_loc
		r_data(i+1,2) = rloc

	enddo

	call DBSNAK(theta_points, r_data(1:theta_points,1),  &
					r_ord,r_data(:,3))

	call DBSINT(theta_points, r_data(1:theta_points,1),  &
		 r_data(1:theta_points,2), r_ord,  &
		 r_data(:,3),  &
		 r_cscoef(1,1:theta_points))

	tri_type = 8

	continue

  end subroutine init_r0_standard


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function psifun(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind
	use pseudo_IMSL, only : DBS2VL
	use exp_data, only : s_ord, xknot_psi, zknot_psi, psi_bscoef,  &
						 nx_FLOW, nz_FLOW

	implicit none

	real(kind=dkind) :: x,z
	real(kind=dkind) :: answer

	answer = DBS2VL(x,z,s_ord,s_ord,xknot_psi,zknot_psi, &
						nx_FLOW,nz_FLOW,psi_bscoef)

	continue

	return

end function psifun





end module triangularity

