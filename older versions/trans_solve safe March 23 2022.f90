module solver

  use constant
  use magnetic
  use p_d_profile
  use flow
  use triangularity
  use exp_data
  use interpolating_functions
!  use IMSL, only : dcsval, dcsder
  use pseudo_IMSL, only : dbsval, dbsder
!  real(kind=dkind) :: dcsval, dcsder
  implicit none
!  private :: dofpsi, dddpsi
!  private :: pofpsi, dpdpsi
!  private :: bzero, dbzerodpsi
  private :: dddpsi
!  private :: dpdpsi
!  private :: dbzerodpsi
  private :: mach_theta, dmach_thetadpsi
  private :: mach_phi, dmach_phidpsi
  private :: sofpsi, dsdpsi, hofpsi, dhdpsi, iofpsi, didpsi
  private :: omegaofpsi, domegadpsi, phiofpsi, dphidpsi
  private :: wofpsi, pparofpsi, pperpofpsi
  private :: sparofpsi,sperpofpsi	
  private :: bpol

  real (kind=dkind), private :: bhut_max,bhut_min
  real (kind=dkind) :: psi_pres

  ! ------------------------------------------------------------------
  ! stuff to speed up functions calculation

  real(kind=dkind), private :: psic_flag = 1.d6
  real(kind=dkind), private :: psi_flag = 1.d9
  real(kind=dkind), private :: psi_flag_dep = 1.d9
  real(kind=dkind), private :: psi_flag_dep2 = 1.d9
  real(kind=dkind), private :: psi_flag_ham = 1.d9

  real(kind=dkind), private :: d_loc, dp_loc, p_loc, pp_loc, b0_loc,  &
							   b0p_loc, mth_loc, mthp_loc, mph_loc,  &
							   mphp_loc
  real(kind=dkind), private :: s_loc, sp_loc, phi_loc, phip_loc, i_loc,  &
							   ip_loc, omega_loc, omegap_loc, h_loc, hp_loc
  real(kind=dkind), private :: ppar_loc, pparp_loc, pperp_loc,pperpp_loc,  &
							   tpar_loc, tparp_loc, tperp_loc,tperpp_loc,  &
							   theta_loc, thetap_loc
  real(kind=dkind), private :: w_loc, wp_loc, spar_loc, sparp_loc, &
							   sperp_loc, sperpp_loc
  real(kind=dkind), private :: dc_loc, dcp_loc, pc_loc, pcp_loc,  &
							   b0c_loc, b0cp_loc
real(kind=dkind) :: Bernmax, fBernmax, delta_Bern, psi_Bern

  ! ------------------------------------------------------------------
  ! The following variables are defined to allow the bernoulli 
  ! function and such to be used and its parameters set by another
  ! function in this module.

  real (kind=dkind), private :: g_Phi, g_r, g_dpsidx, g_dpsidz
  real (kind=dkind), private :: g_dpsidx_L, g_dpsidz_L, g_dpsidx_R, g_dpsidz_R
  real (kind=dkind), private :: g_I, g_Omega, g_S, g_H, g_wofpsi
  real (kind=dkind), private :: g_wperp, g_wpar
  integer, private :: g_indi, g_indj, g_nx, g_nz
  real (kind=dkind), private :: g_dx, g_dz, g_spar, g_sperp, g_mtheta
  real (kind=dkind) :: g_b_polc, g_b_torc, g_bfield, g_delta
  real (kind=dkind), private :: g_tpar, g_theta, g_D, g_Lambda
!!$ ---------------------------------------------------------------------
!!$ NOTES:
!!$ 3/14/2002: Change made to dofpsi and dddpsi, used to be 1 and 0.
!!$ 3/14/2002: Change alpha from 3 to 1
!!$ 3/14/2002: Change beta_center = 5.0d0
!!$ 3/14/2002: Change mach_theta_edge = 0.0d0
!!$ ---------------------------------------------------------------------

  integer, private :: m_bound	!

  logical, private :: inside,outside
  real(kind=dkind), private :: bpol_av = 0.d0
  real(kind=dkind) :: bpol_max, bpol_min
  logical, private :: reduce=.false.
  integer :: i_bmin, j_bmin

  ! stuff for the grid
  real(kind=dkind) :: g_ratio
  integer :: n_temp
  real(kind=dkind), private :: der1, der2 ! for derivative interpolation
  real(kind=dkind), dimension(:,:), allocatable, public :: fmax_2D

	real(kind=dkind), dimension(:,:), allocatable :: rho1, rho2
	logical :: clean_edge_rho, clean_edge_rho_psi

	logical :: initialize_r_for_tri_type_m2 = .false.
	logical :: tri_type_m2_ready = .false.

	real(kind=dkind) :: xT1, xT2, xT3, zT1, zT2, zT3, two_A_tri


contains

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : bc_setup_option

	implicit none

	integer :: n

	if(bc_setup_option>0) then
		call initialize_bc_equations(n)
	elseif(bc_setup_option<0) then
		call initialize_bc_minimum(n)
	endif

end subroutine initialize_bc


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc_equations(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! March 26 2013:
! added geometry stuff for extrapolation near the edge
! to avoid messy nested ifs, everything is saved anyway (in bound_temp),
! just not transfered to global arrays if not needed
! once and for all, "P" is the grid point where Bc's need to be assigned, "Q" is the interpolated point
! "Q2" is the inner interpolated point (used for gradients with triangular functions)

! the array variables:

! coord_bound:
! 1 -> xQ (-1<xQ<1)
! 2 -> zQ (-1<zQ<1)
! 3 -> +-1 (inner/outer) or 2 (point on boundary, it never happens)
! 4 -> distance P - boundary
! 5 -> distance Q - boundary

! ind_bound:
! 1 -> iP
! 2 -> jP
! 3 -> i_point1
! 4 -> j_point1

! bound_13
! 0 -> index of P in cell (1-4 if internal, -1 otherwise)
! 1 -> xP (-1<xQ<1)
! 2 -> zP (-1<zQ<1)
! 3 -> cx (cos of unit vector direction)
! 4 -> cz (sin of unit vector direction)
! 5 -> RQ (grid coordinate)
! 6 -> ZQ (grid coordinate)

! bound_tri(m_bound,1:3,1:2): i,j coordinates of points 1-3 in the triangle enclosing Q2
! triangle does not use external points (for stability)

! coord_tri:
! 1 -> RQ2
! 2 -> ZQ2
! 3 -> distance P-Q2 (for convenience)
! 4 -> triangle area (for convenience)

	use constant, only : dx, dz, x_coord, z_coord, dx_a, dz_a, bc_switch

	integer :: n

	integer, dimension(n*n,4) :: ind
	real(kind=dkind), dimension(n*n,6) :: coord
	real(kind=dkind), dimension(n*n,4) :: coord_tri_temp
	integer, dimension(n*n,1:3,1:2) :: bound_tri_temp
	real(kind=dkind), dimension(n*n,0:6) :: bound_temp !bound_13
	! third number is distance ratio

	real(kind=dkind) :: ex,ez, dist0,dist, dist_unit
	real(kind=dkind) :: xb, zb, th
	real(kind=dkind), dimension(1:2) :: xvec, fvec ! xvec is the position of the point on the boundary
	! solution points and function values from IMSL
	real(kind=dkind), dimension(1:2) :: xguess
	real(kind=dkind), dimension(1:2) :: xscale, fscale
	real(kind=dkind) :: fnorm
	integer :: iparam(6)
	real(kind=dkind) ::  rparam(5)
	external  eq_bound, eq_bound2
	integer :: ntrial = 5000
	real(kind=dkind) :: tolx=1.d-12
	real(kind=dkind) :: tolf = 1.d-9
	integer :: zone
	real(kind=dkind) :: dummy(1:7)
	real(kind=dkind) :: r
	real(kind=dkind) :: RQ2, ZQ2, cx, cz

	integer :: i_vec(1:4), j_vec(1:4)
	integer :: i, j, k

	logical :: bound, inner, truebound, zero_dist
	logical :: newton_check
	integer :: i_diff, j_diff

!	if(n>=65) 	fix_orp = 1.1435d-1	!	fix_orp = fix_orp*7.d0	!	accelerate = .true.	!	
!	if(n>=65) max_it = max_it*2 !

	if( (bc_type<3).or.(bc_type==13).or.(bc_type==23).or.  &
		( ((bc_type==4).or.(bc_type==5).or.(bc_type==7).or.(bc_type==8).or.(bc_type==14).or.(bc_type==24)).and.(n<bc_switch)) ) return

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)
	if(allocated(bound_13)) deallocate(bound_13)
	if(allocated(coord_tri)) deallocate(coord_tri)
	if(allocated(bound_tri)) deallocate(bound_tri)

	clean_edge_rho = .true.	! cleans up leftovers from previous grid
	clean_edge_rho_psi = .true.	! cleans up leftovers from previous grid

	m_bound = 0

	xscale = 1.d0
	fscale = 1.d0
	iparam(1) = 0
	coord = 0.d0

	zero_dist = .false.

	i_vec(1) = 0
	i_vec(2) = 1
	i_vec(3) = 1
	i_vec(4) = 0

	j_vec(1) = 0
	j_vec(2) = 0
	j_vec(3) = 1
	j_vec(4) = 1

	do j=1, n
       do i=1, n

			call check_position(i,j,bound,truebound,n)

			if(truebound) then
			! real boundary point

				m_bound = m_bound+1

				call radius(i,j,n,n,ex,ez,r,dx,dz)

				if(ex**2+ez**2<r**2) then
				! the point is very close to the boundary

					R_P = x_coord(i)
					z_P = z_coord(j)

					xguess(1) = R_P
					xguess(2) = z_P

					xvec(1) = xguess(1)
					xvec(2) = xguess(2)

					call newt(ntrial,xvec,2,newton_check,tolx,tolf)

					if((bc_type==7).or.(bc_type==8)) then

						xb = xvec(1)
						zb = xvec(2)
						th = atan2(zb,xb-rmajor)
						if(th<0.d0) th = th + 2.d0*pi

						coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
										ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )

					endif

					dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

					if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

					else

						! proceed as usual

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

						ind(m_bound,1) = i
						ind(m_bound,2) = j

						inner = .false.

						coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
						coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

						bound_temp(m_bound,3) = 2.d0*(R_P - xvec(1))/(dist+dist0)
						bound_temp(m_bound,4) = 2.d0*(z_P - xvec(2))/(dist+dist0)
						! direction of normal unit vector (note that the norm could be <1:
						! this matters for internal points, and thus needs to be rescaled)

						dist_unit = sqrt(bound_temp(m_bound,3)**2+bound_temp(m_bound,4)**2)

						if(dist_unit<1.d0) then
						! need to rescale unit vector

							bound_temp(m_bound,3) = bound_temp(m_bound,3)/dist_unit
							bound_temp(m_bound,4) = bound_temp(m_bound,4)/dist_unit

						endif

						bound_temp(m_bound,5) = R_P + bound_temp(m_bound,3)*(dist+dist0)
						bound_temp(m_bound,6) = z_P + bound_temp(m_bound,4)*(dist+dist0)

						call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
									ind(m_bound,3:4) )

						call check_points(ind(m_bound,3),ind(m_bound,4),n,dx,dz,inner)

						if(inner) then

							coord(m_bound,3) = 1.d0

						else

							coord(m_bound,3) = -1.d0

						endif

						coord(m_bound,4) = 0.d0 !distance from point to surface
						coord(m_bound,5) = dist !distance from surface to interpolated point

						if(bc_type==5) then

							call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
											dummy(3),zone,dummy(4),dummy(5),dummy(6))
							! NOTE: this works only because the inner boundary is a circle,
							! in the general case one should use the coordinates of the 
							! point on the boundary

						endif

						! difference between the "P" point coordinates and corner point coordinates
						i_diff = i - ind(m_bound,3)
						j_diff = j - ind(m_bound,4)

						! this determines which point of the cell corresponds to P
						! if the point is internal obviously it's not a cell point
						if((i_diff==0).and.(j_diff==0)) then
							bound_temp(m_bound,0) = 1
							bound_temp(m_bound,1) = -1.d0
							bound_temp(m_bound,2) = -1.d0
						elseif((i_diff==1).and.(j_diff==0)) then
							bound_temp(m_bound,0) = 2
							bound_temp(m_bound,1) = 1.d0
							bound_temp(m_bound,2) = -1.d0
						elseif((i_diff==1).and.(j_diff==1)) then
							bound_temp(m_bound,0) = 3
							bound_temp(m_bound,1) = 1.d0
							bound_temp(m_bound,2) = 1.d0
						elseif((i_diff==0).and.(j_diff==1)) then
							bound_temp(m_bound,0) = 4
							bound_temp(m_bound,1) = -1.d0
							bound_temp(m_bound,2) = 1.d0
						else
							bound_temp(m_bound,0) = -1
						endif

!						dist_pippa = bound_temp(m_bound,3)**2+bound_temp(m_bound,4)**2
						continue

					endif

				else
				! the point REALLY is on the boundary (numerically),
				! this is extremely unlikely to occur and dealt with easily

					zero_dist = .true.

					continue

				endif

			elseif(bound) then
				! (external) boundary point

				m_bound = m_bound+1

				R_P = x_coord(i)
				z_P = z_coord(j)

				xguess(1) = R_P
				xguess(2) = z_P

				xvec(1) = xguess(1)
				xvec(2) = xguess(2)

				call newt(ntrial,xvec,2,newton_check,tolx,tolf)

				if((bc_type==7).or.(bc_type==8)) then

					xb = xvec(1)
					zb = xvec(2)
					th = atan2(zb,xb-rmajor)
					if(th<0.d0) th = th + 2.d0*pi

					coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
									ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )

				endif

				dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

				if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

				else

					if(dist0<0.2d0*sqrt(dx_a(i)**2+dz_a(j)**2)) then

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

					else

						dist = dist0

					endif

					ind(m_bound,1) = i
					ind(m_bound,2) = j

					inner = .false.

					coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
					coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

					bound_temp(m_bound,3) = 2.d0*(R_P - xvec(1))/(dist+dist0)
					bound_temp(m_bound,4) = 2.d0*(z_P - xvec(2))/(dist+dist0)
					! direction of normal unit vector (note that the norm could be <1:
					! this matters for internal points, and thus needs to be rescaled)

					dist_unit = sqrt(bound_temp(m_bound,3)**2+bound_temp(m_bound,4)**2)

					if(dist_unit<1.d0) then
					! need to rescale unit vector

						bound_temp(m_bound,3) = bound_temp(m_bound,3)/dist_unit
						bound_temp(m_bound,4) = bound_temp(m_bound,4)/dist_unit

					endif

					bound_temp(m_bound,5) = R_P + bound_temp(m_bound,3)*(dist+dist0)
					bound_temp(m_bound,6) = z_P + bound_temp(m_bound,4)*(dist+dist0)

					call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
								ind(m_bound,3:4) )

					call check_points(ind(m_bound,3),ind(m_bound,4),n,dx,dz,inner)

					if(inner) then

						coord(m_bound,3) = 1.d0

					else

						coord(m_bound,3) = -1.d0

					endif

					coord(m_bound,4) = dist0 !distance from point to surface
					coord(m_bound,5) = dist !distance from surface to interpolated point

					! difference between the "P" point coordinates and corner point coordinates
					i_diff = i - ind(m_bound,3)
					j_diff = j - ind(m_bound,4)

					! this determines which point of the cell corresponds to P
					! if the point is internal obviously it's not a cell point
					if((i_diff==0).and.(j_diff==0)) then
						bound_temp(m_bound,0) = 1
						bound_temp(m_bound,1) = -1.d0
						bound_temp(m_bound,2) = -1.d0
					elseif((i_diff==1).and.(j_diff==0)) then
						bound_temp(m_bound,0) = 2
						bound_temp(m_bound,1) = 1.d0
						bound_temp(m_bound,2) = -1.d0
					elseif((i_diff==1).and.(j_diff==1)) then
						bound_temp(m_bound,0) = 3
						bound_temp(m_bound,1) = 1.d0
						bound_temp(m_bound,2) = 1.d0
					elseif((i_diff==0).and.(j_diff==1)) then
						bound_temp(m_bound,0) = 4
						bound_temp(m_bound,1) = -1.d0
						bound_temp(m_bound,2) = 1.d0
					else
						bound_temp(m_bound,0) = -1
					endif

					if(bc_type==5) then

						call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
										dummy(3),zone,dummy(4),dummy(5),dummy(6))
						! NOTE: this works only because the inner boundary is a circle,
						! in the general case one should use the coordinates of the 
						! point on the boundary

					endif

					if((bc_type==23).or.(bc_type==24))  then
						call get_Q2
					endif

!					dist_pippa = bound_temp(m_bound,3)**2+bound_temp(m_bound,4)**2
					continue

				endif

				continue

			endif

			if(zero_dist) then
			! take care of the boundary point
			! (what follows in coord(m,1:2) does not matter, it's just to avoid "funny" numbers)

				ind(m_bound,1) = i
				ind(m_bound,2) = j

				coord(m_bound,1) = rmajor
				coord(m_bound,2) = 0.d0
				coord(m_bound,3) = 2.d0

			endif

       end do
    end do

	if(tri_type==11) then

		allocate(coord_bound(m_bound,6))
		allocate(ind_bound(m_bound,5))

	else

		allocate(ind_bound(m_bound,4))
		if((bc_type==7).or.(bc_type==8)) then
			allocate(coord_bound(m_bound,6))
		else
			allocate(coord_bound(m_bound,5))
		endif

		if((bc_type==13).or.(bc_type==14).or.(bc_type==23).or.(bc_type==24)) then
			allocate(bound_13(m_bound,0:6))
			bound_13 = 0.d0
		endif

		if((bc_type==23).or.(bc_type==24)) then
			allocate(coord_tri(m_bound,1:4))
			coord_tri = 0.d0
			allocate(bound_tri(m_bound,1:3,1:2))
			bound_tri = 0
		endif

	endif

	do i=1,m_bound

		do j=1,4
			ind_bound(i,j) = ind(i,j)
		enddo

		coord_bound(i,1) = 2.d0/dx_a(ind(i,3))*(coord(i,1) - x_coord(ind(i,3))) - 1.d0

		coord_bound(i,2) = 2.d0/dz_a(ind(i,4))*(coord(i,2) - z_coord(ind(i,4))) - 1.d0

		coord_bound(i,3) = coord(i,3) ! "in" or "out"

		coord_bound(i,4) = coord(i,4) ! distance from point to surface

		coord_bound(i,5) = coord(i,5) ! distance from surface to interpolated point

		if((bc_type==7).or.(bc_type==8)) 	coord_bound(i,6) = coord(i,6) ! psi on the boundary

		if(tri_type==11) then

			! what matters is the points outside!

			ind_bound(i,5) = sort_grid(ind(i,1),ind(i,2),1)

			if(bc_type==5) coord_bound(i,6) =  gtheta_bound(coord(i,6))

		endif

		if((bc_type==13).or.(bc_type==14).or.(bc_type==23).or.(bc_type==24)) then
			do j = 0, 6
				bound_13(i,j) = bound_temp(i,j)
			enddo
		endif

		if((bc_type==23).or.(bc_type==24)) then

			do j = 1, 4
				coord_tri(i,j) = coord_tri_temp(i,j)
			enddo

			do j = 1, 3
			do k = 1, 2
				bound_tri(i,j,k) = bound_tri_temp(i,j,k)
			enddo
			enddo

		endif

		continue

	enddo

	continue

	contains

	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_Q2
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: d2max, dist_temp, new_dist
	integer :: ii, jj, iii, jjj, i_far, j_far
	logical :: inner2
	real(kind=dkind), dimension(1:2) :: p1, p2, p3

	d2max = -10.d0

	cx = bound_temp(m_bound,3)
	cz = bound_temp(m_bound,4)

	RQ2 = bound_temp(m_bound,5)
	ZQ2 = bound_temp(m_bound,6)
	! if the point is external, this will be modified
	! Q2 will be moved also if internal, to get it as close to the triangle center as possible

	! define inner2: 
	! inner=.true. and inner2=.true -> point P is really inside;
	! inner=.true. and inner2=.false -> point P is outside, but other points in the cell are inside;
	! inner=.false. and inner2=.false -> point P is outside and at least one more point in the cell is outside.

	! first determine if point is inner or outer

	call check_point1(ind(m_bound,3),ind(m_bound,4),inner2)

	! CAREFUL, USE "INNER2", HERE POINT IS INNER ONLY IF THERE ARE NO OUTER POINTS
!	if(not(inner).and.not(inner2)) then
	if(inner.and.inner2) then
	! no outer points in the cell, point is really inner, keep three closest points,
	! but still move Q2 to "center" of triangle

		! find farthest point

		iii = ind(m_bound,3)
		jjj = ind(m_bound,4)

		do jj = 1, 4

			dist_temp = dist2(x_coord(iii+i_vec(jj)),z_coord(jjj+j_vec(jj)),RQ2,ZQ2)
			if(dist_temp>d2max) then
				d2max = dist_temp
				j_far = jj
			endif

		enddo

		! fill in triangle points

		ii = 1

		do jj = 1, 4

			if(jj/=j_far) then

				bound_tri_temp(m_bound,ii,1) = iii + i_vec(jj)
				bound_tri_temp(m_bound,ii,2) = jjj + j_vec(jj)
				ii = ii + 1

			endif

		enddo

		call optimum_d(new_dist)

		RQ2 = R_P - new_dist*cx
		ZQ2 = Z_P - new_dist*cz

		coord_tri_temp(m_bound,1) = RQ2
		coord_tri_temp(m_bound,2) = ZQ2
		coord_tri_temp(m_bound,3) = new_dist

		p1(1) = x_coord(bound_tri_temp(m_bound,1,1))
		p2(1) = x_coord(bound_tri_temp(m_bound,2,1))
		p3(1) = x_coord(bound_tri_temp(m_bound,3,1))

		p1(2) = z_coord(bound_tri_temp(m_bound,1,2))
		p2(2) = z_coord(bound_tri_temp(m_bound,2,2))
		p3(2) = z_coord(bound_tri_temp(m_bound,3,2))

		coord_tri_temp(m_bound,4) = tri_area(p1,p2,p3)

		return

	endif

	!------------------ external point --------------------

	! we need to distinguish two cases:
	! 1) there is only one outer point (that's P)
	! 2) there is more than one outer point
	! actually, no, just move Q2 in the normal direction until the triangle only has internal points

	call find_triangle

	continue

	end subroutine get_Q2

	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine find_triangle
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	! moves Q2 until there is a surrounding triangle that only has internal points

	real(kind=dkind) :: new_dist, dd
	logical :: found
	integer :: iP, jP
	integer :: p1_tri(1:2)

	found = .false.

	dd = 0.1d0*sqrt(dx**2+dz**2)
	iP = ind(m_bound,1)
	jP = ind(m_bound,2)

	new_dist = dist2(R_P,Z_P,RQ2,ZQ2)
	new_dist = sqrt(new_dist)

	do

		call point1(iP,jP,RQ2,ZQ2,p1_tri)
		! get point1 of element around Q2

		call check_4_tri(p1_tri,found)
		! check four triangles in element, and find the optimal one if there is an internal one
		! (bound_tri_temp is filled in check_4_tri)

		if(found) exit

		new_dist = new_dist + dd

		RQ2 = R_P - new_dist*cx
		ZQ2 = Z_P - new_dist*cz

	enddo

	call optimum_d(new_dist)

	RQ2 = R_P - new_dist*cx
	ZQ2 = Z_P - new_dist*cz

	coord_tri_temp(m_bound,1) = RQ2
	coord_tri_temp(m_bound,2) = ZQ2
	coord_tri_temp(m_bound,3) = new_dist

	continue

	end subroutine find_triangle

	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_4_tri(p1_tri,found)
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: p1_tri(1:2)
	logical :: found
	integer, dimension(1:4,1:3) :: i_tri, j_tri
	! these two contain the triangle indexes, which allows to write a cycle instead of repeating intructions
	integer, dimension(1:4) :: both_inside
	integer :: good_triangles, best_triangle
	real(kind=dkind) :: A1, A2, A3, Atot, dmin, dtemp
	real(kind=dkind), dimension(1:2) :: p1, p2, p3, QQ
	integer ::  iP, jP
	integer :: iii, jjj

	good_triangles = 0
	both_inside= 0
	! how many and which triangles satisfy the test

	iP = ind(m_bound,1)
	jP = ind(m_bound,2)

	QQ(1) = RQ2
	QQ(2) = ZQ2

	! --------------- first set triangles ---------------
	! triangle 1 |-
	i_tri(1,1) = 0
	j_tri(1,1) = 0

	i_tri(1,2) = 1
	j_tri(1,2) = 1

	i_tri(1,3) = 0
	j_tri(1,3) = 1

	! triangle 2 _|
	i_tri(2,1) = 0
	j_tri(2,1) = 0

	i_tri(2,2) = 1
	j_tri(2,2) = 1

	i_tri(2,3) = 1
	j_tri(2,3) = 0

	! triangle 3 |_
	i_tri(3,1) = 0
	j_tri(3,1) = 0

	i_tri(3,2) = 1
	j_tri(3,2) = 0

	i_tri(3,3) = 0
	j_tri(3,3) = 1

	! triangle 4 -|
	i_tri(4,1) = 1
	j_tri(4,1) = 0

	i_tri(4,2) = 1
	j_tri(4,2) = 1

	i_tri(4,3) = 0
	j_tri(4,3) = 1

	! --------------- done with triangles ---------------

	! for each triangle, check if point Q2 is inside; if so, check if triangle is internal

	tri_loop: do jjj = 1, 4

		p1(1) = x_coord(iP+i_tri(jjj,1))
		p2(1) = x_coord(iP+i_tri(jjj,2))
		p3(1) = x_coord(iP+i_tri(jjj,3))

		p1(2) = z_coord(jP+j_tri(jjj,1))
		p2(2) = z_coord(jP+j_tri(jjj,2))
		p3(2) = z_coord(jP+j_tri(jjj,3))

		Atot = tri_area(p1,p2,p3)
		A1 = tri_area(QQ,p2,p3)
		A2 = tri_area(p1,QQ,p3)
		A3 = tri_area(p1,p2,QQ)

		if((A1+A2+A3)>Atot) cycle
		! point Q2 is not in triagle

		! if point Q2 is in triangle, check the triangle

		do iii = 1, 3

			if(sort_grid(iP+i_tri(jjj,iii),jP+j_tri(jjj,iii),0) < 0) cycle tri_loop

		enddo

		! if we made it this far, the triangle is OK

		good_triangles = good_triangles + 1
		both_inside(good_triangles) = jjj

	enddo tri_loop

	if(good_triangles>=1) then

		found = .true.

		if(good_triangles==1) then

			best_triangle = both_inside(good_triangles)

		else
		! use the triangle closest to P

			dmin = 1.d9*(dx**2+dz**2)

			do jjj = 1, good_triangles

				dtemp = 0.d0

				do iii = 1,3

					dtemp = dtemp + dist2(R_P,Z_P,x_coord(iP+i_tri(jjj,iii)),z_coord(jP+j_tri(jjj,iii)))

				enddo

				if(dtemp<dmin) then

					dmin = dtemp
					best_triangle = both_inside(jjj)

				endif

			enddo

		endif

		do iii = 1, 3

			bound_tri_temp(m_bound,iii,1) = iP + i_tri(best_triangle,iii)
			bound_tri_temp(m_bound,iii,2) = jP + j_tri(best_triangle,iii)

		enddo

		jjj = best_triangle

		p1(1) = x_coord(iP+i_tri(jjj,1))
		p2(1) = x_coord(iP+i_tri(jjj,2))
		p3(1) = x_coord(iP+i_tri(jjj,3))

		p1(2) = z_coord(jP+j_tri(jjj,1))
		p2(2) = z_coord(jP+j_tri(jjj,2))
		p3(2) = z_coord(jP+j_tri(jjj,3))

		Atot = tri_area(p1,p2,p3)
		coord_tri_temp(m_bound,4) = Atot

	endif

	continue

	end subroutine check_4_tri

	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine optimum_d(d)
	!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: d
	real(kind=dkind) :: R1, R2, R3, Z1, Z2, Z3

	R1 = x_coord(bound_tri_temp(m_bound,1,1))
	Z1 = z_coord(bound_tri_temp(m_bound,1,2))
	R2 = x_coord(bound_tri_temp(m_bound,2,1))
	Z2 = z_coord(bound_tri_temp(m_bound,2,2))
	R3 = x_coord(bound_tri_temp(m_bound,3,1))
	Z3 = z_coord(bound_tri_temp(m_bound,3,2))

	d = abs((cx*(R1 + R2 + R3 - 3.d0*R_P) + cz*(Z1 + Z2 + Z3 - 3.d0*Z_P))/3.d0)
	! due to the direction ambiguity in cx, cz, this could be negative
	! (actually I think it WILL be negative with the definitions in use CHECK)

	continue

	return

	end subroutine optimum_d

end subroutine initialize_bc_equations

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_bc_minimum(n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!	use IMSL, only : DNEQBF, DNEQNF

	use constant, only : dx, dz, x_coord, z_coord, dx_a, dz_a, bc_switch

	integer :: n

	integer, dimension(n*n,4) :: ind
	real(kind=dkind), dimension(n*n,6) :: coord
	! third number is distance ratio

	real(kind=dkind) :: ex,ez, dist0,dist
	real(kind=dkind) :: x,z, RP, zP, R0, xb, zb, th
	real(kind=dkind) ::  xs, zs ! points on the surface
	real(kind=dkind), dimension(1:2) :: xvec, fvec 
	! solution points and function values from IMSL
	real(kind=dkind), dimension(1:2) :: xguess
	real(kind=dkind), dimension(1:2) :: xscale, fscale
	real(kind=dkind) :: fnorm
	integer :: iparam(6)
	real(kind=dkind) ::  rparam(5)
	external  eq_bound, eq_bound2
	integer :: ntrial = 5000
	real(kind=dkind) :: tolx=1.d-12
	real(kind=dkind) :: tolf = 1.d-9
	integer :: zone
	real(kind=dkind) :: dummy(1:7)
	real(kind=dkind) :: r
	real(kind=dkind), dimension(1:3) :: theta_start
	real(kind=dkind) :: thetamin

	integer :: i,j

	logical :: bound, inner, truebound, zero_dist
	logical :: newton_check

!	if(n>=65) 	fix_orp = 1.1435d-1	!	fix_orp = fix_orp*7.d0	!	accelerate = .true.	!	
!	if(n>=65) max_it = max_it*2 !

	if( (bc_type<3).or.(bc_type==13).or.  &
		( ((bc_type==4).or.(bc_type==5).or.(bc_type==7).or.(bc_type==8).or.(bc_type==14)).and.(n<bc_switch)) ) return

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)

	m_bound = 0

	xscale = 1.d0
	fscale = 1.d0
	iparam(1) = 0
	coord = 0.d0

	zero_dist = .false.

	do j=1, n
       do i=1, n

			call check_position(i,j,bound,truebound,n)

			if(truebound) then
			! real boundary point

				m_bound = m_bound+1

				call radius(i,j,n,n,ex,ez,r,dx,dz)

				if(ex**2+ez**2<r**2) then
				! the point is very close to the boundary

					R_P = x_coord(i)
					z_P = z_coord(j)

					call radius_1_3(R_P,z_P,ex,ez,theta_start(2),r,zone,dummy(1),dummy(2),dummy(3))

					theta_start(1) = theta_start(2) - pi/5.d0
					theta_start(3) = theta_start(2) + pi/5.d0

					dist = brent(theta_start(1), theta_start(2), theta_start(3), dist_fun, tolf, thetamin)

					call radius_theta(thetamin,r,xvec(1),xvec(2))

					if((bc_type==7).or.(bc_type==8)) then

						xb = xvec(1)
						zb = xvec(2)
						th = atan2(zb,xb-rmajor)
						if(th<0.d0) th = th + 2.d0*pi

						coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
										ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )

					endif

					dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

					if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

					else

						! proceed as usual

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

						ind(m_bound,1) = i
						ind(m_bound,2) = j

						inner = .false.

		!				coord(m_bound,1) = (1.d0+1.d0) * xvec(1) - 1.d0 * R_P
		!				coord(m_bound,2) = (1.d0+1.d0) * xvec(2) - 1.d0 * z_P

						coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
						coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

						call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
									ind(m_bound,3:4) )

						call check_points(ind(m_bound,3),ind(m_bound,4),n,dx,dz,inner)

						if(inner) then

							coord(m_bound,3) = 1.d0

						else

							coord(m_bound,3) = -1.d0

						endif

						coord(m_bound,4) = 0.d0 !distance from point to surface
						coord(m_bound,5) = dist !distance from surface to interpolated point


						if(bc_type==5) then

							call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
											dummy(3),zone,dummy(4),dummy(5),dummy(6))
							! NOTE: this works only because the inner boundary is a circle,
							! in the general case one should use the coordinates of the 
							! point on the boundary

						endif

					endif

				else
				! the point REALLY is on the boundary (numerically),
				! this is extremely unlikely to occur and dealt with easily

					zero_dist = .true.

					continue

				endif

			elseif(bound) then
				! (external) boundary point

				m_bound = m_bound+1

				R_P = x_coord(i)
				z_P = z_coord(j)

				call radius_1_3(R_P,z_P,ex,ez,theta_start(2),r,zone,dummy(1),dummy(2),dummy(3))

				theta_start(1) = theta_start(2) - pi/5.d0
				theta_start(3) = theta_start(2) + pi/5.d0

				dist = brent(theta_start(1), theta_start(2), theta_start(3), dist_fun, tolf, thetamin)

				call radius_theta(thetamin,r,xvec(1),xvec(2))

				if((bc_type==7).or.(bc_type==8)) then

					xb = xvec(1)
					zb = xvec(2)
					th = atan2(zb,xb-rmajor)
					if(th<0.d0) th = th + 2.d0*pi

					coord(m_bound,6) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
									ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )

				endif

				dist0 = ( (xvec(1) - R_P)**2 + (xvec(2) - z_P)**2 )**0.5d0

				if(dist0==0.d0) then

						! somehow we would not want this to happen,
						! but we'll treat the point as "real boundary"

						zero_dist = .true.

				else

					if(dist0<0.2d0*sqrt(dx_a(i)**2+dz_a(j)**2)) then

						dist = 0.5d0*sqrt(dx_a(i)**2+dz_a(j)**2)

					else

						dist = dist0

					endif

					ind(m_bound,1) = i
					ind(m_bound,2) = j

					inner = .false.

					coord(m_bound,1) = (1.d0 + dist/dist0) * (xvec(1) - R_P) + R_P
					coord(m_bound,2) = (1.d0 + dist/dist0) * (xvec(2) - z_P) + z_P

					call point1(i,j,coord(m_bound,1),coord(m_bound,2),  &
								ind(m_bound,3:4) )

					call check_points(ind(m_bound,3),ind(m_bound,4),n,dx,dz,inner)

					if(inner) then

						coord(m_bound,3) = 1.d0

					else

						coord(m_bound,3) = -1.d0

					endif

					coord(m_bound,4) = dist0 !distance from point to surface
					coord(m_bound,5) = dist !distance from surface to interpolated point

					if(bc_type==5) then

						call radius_1_3(R_P,z_P,dummy(1),dummy(2),coord(m_bound,6),  &
										dummy(3),zone,dummy(4),dummy(5),dummy(6))
						! NOTE: this works only because the inner boundary is a circle,
						! in the general case one should use the coordinates of the 
						! point on the boundary

					endif

					continue

				endif

				continue

			endif

			if(zero_dist) then
			! take care of the boundary point
			! (what follows in coord(m,1:2) does not matter, it's just to avoid "funny" numbers)

				ind(m_bound,1) = i
				ind(m_bound,2) = j

				coord(m_bound,1) = rmajor
				coord(m_bound,2) = 0.d0
				coord(m_bound,3) = 2.d0

			endif

       end do
    end do

	if(tri_type==11) then
		allocate(coord_bound(m_bound,6))
		allocate(ind_bound(m_bound,5))
	else
		allocate(ind_bound(m_bound,4))
		if((bc_type==7).or.(bc_type==8)) then
			allocate(coord_bound(m_bound,6))
		else
			allocate(coord_bound(m_bound,5))
		endif
	endif

	do i=1,m_bound

		do j=1,4
			ind_bound(i,j) = ind(i,j)
		enddo

		coord_bound(i,1) = 2.d0/dx_a(ind(i,3))*(coord(i,1) - x_coord(ind(i,3))) - 1.d0

		coord_bound(i,2) = 2.d0/dz_a(ind(i,4))*(coord(i,2) - z_coord(ind(i,4))) - 1.d0

		coord_bound(i,3) = coord(i,3) ! "in" or "out"

		coord_bound(i,4) = coord(i,4) ! distance from point to surface

		coord_bound(i,5) = coord(i,5) ! distance from surface to interpolated point

		if((bc_type==7).or.(bc_type==8)) 	coord_bound(i,6) = coord(i,6) ! psi on the boundary

		if(tri_type==11) then

			! what matters is the points outside!

			ind_bound(i,5) = sort_grid(ind(i,1),ind(i,2),1)

			if(bc_type==5) coord_bound(i,6) =  gtheta_bound(coord(i,6))

		endif

		continue

	enddo

	continue

end subroutine initialize_bc_minimum

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dist_fun(theta) result(dist)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: theta, dist
	real(kind=dkind) :: r, Rloc, Zloc

	call radius_theta(theta,r,Rloc,Zloc)

	dist = sqrt((Rloc-R_P)**2+(Zloc-z_P)**2)

end function dist_fun

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho0(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho

	if (bc_type==1) then 

		call bc_psi_rho1(psi,rho,nx,nz)

	elseif (bc_type==2) then

		call bc_psi_rho2(psi,rho,nx,nz)

	elseif (bc_type==3) then

		call bc_psi_rho3(psi,rho,nx,nz)

	elseif (bc_type==4) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho3(psi,rho,nx,nz)

	elseif (bc_type==5) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho4(psi,rho,nx,nz)

	elseif (bc_type==6) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho5(psi,rho,nx,nz)

	elseif ((bc_type==7).or.(bc_type==8)) then

		call bc_psi_rho7(psi,rho,nx,nz)

	elseif (bc_type==13) then

		call bc_psi_rho13(psi,rho,nx,nz)

	elseif (bc_type==14) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho13(psi,rho,nx,nz)

	elseif (bc_type==23) then

		call bc_psi_rho23(psi,rho,nx,nz)

	elseif (bc_type==24) then

		if(nx<bc_switch) call bc_psi_rho1(psi,rho,nx,nz)
		if(nx>=bc_switch) call bc_psi_rho23(psi,rho,nx,nz)

	else

		print*, 'unknown option for bc_type:'
		print*, 'bc_type =     ',bc_type

	endif

	if(tri_type==10) call bc_psi_rho1_bis(psi,nx,nz)

	continue

end subroutine bc_psi_rho0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho1(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord

	integer :: nx,nz
!	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k
	real(kind=dkind) :: x,z,dummy(1:3),alpha
	integer :: zone


    den = dofpsi(0.0d0)

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

    do j=1, nz
       do i=1, nx

		  if(tri_type==11) then
		  ! LDX case

				x = x_coord(i)
				z = z_coord(j)

				call radius_1_3(x,z,ex,ez,dummy(1),rminor,zone,alpha,dummy(2),dummy(3))

				if(((ex*alpha)**2 + (ez*alpha)**2) >= rminor**2) then

					if(zone==4) then
					! inner point

						rho(i,j) = den_in
						psi(i,j) = psi_in

					elseif(zone==1) then
					! outer point

						rho(i,j) = den_out
						psi(i,j) = psi_out

					endif

				endif

			else

				  if(sort_grid(i,j,0)<=0) then
					 rho(i,j) = den
					 psi(i,j) = 0.0d0
				  end if

			endif

       end do
    end do

	continue

end subroutine bc_psi_rho1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho2(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! WARNING: this is obsolete, it will not work
  ! this still to fix for general grid case

	use constant, only : dx,dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: den
	real (kind=dkind) :: ex2,ez2
	logical :: boundary
    integer :: i,j,k,ii,jj
	real (kind=dkind) :: rminor2,f_in,f_out


    den = dofpsi(0.0d0)

    do j=2, nz-1
       do i=2, nx-1

		  if(sort_grid(i,j,0)<=0) then

             rho(i,j) = den
             psi(i,j) = 0.0d0
			 boundary = .false.

			 ! dx,dz distance section

			 if(dx>=dz) then

				 ii=i
				 jj=j+1
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i
				 jj=j-1
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i-1
				 jj=j
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i+1
				 jj=j
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

			 else

				 ii=i+1
				 jj=j
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i-1
				 jj=j
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i
				 jj=j-1
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

				 ii=i
				 jj=j+1
				 if(sort_grid(ii,jj,0)==1) then
				 boundary = .true.
				 goto 666
				 endif

			endif

			 ! sqrt(dx**2+dz**2) distance section
			 ii=i+1
			 jj=j+1
			 if(sort_grid(ii,jj,0)==1) then
			 boundary = .true.
			 goto 666
			 endif

			 ii=i+1
			 jj=j-1
			 if(sort_grid(ii,jj,0)==1) then
			 boundary = .true.
			 goto 666
			 endif

			 ii=i-1
			 jj=j-1
			 if(sort_grid(ii,jj,0)==1) then
			 boundary = .true.
			 goto 666
			 endif

			 ii=i-1
			 jj=j+1
			 if(sort_grid(ii,jj,0)==1) then
			 boundary = .true.
			 goto 666
			 endif

666			 continue

			 if(boundary) then
			 ! this part erased, do not use bc_type=2 anymore

!				f_in = rminor2-dsqrt(ex2**2+ez2**2)	
!				f_out = -rminor+dsqrt(ex**2+ez**2)	
!				if(dabs(f_out/f_in)<.9d0) then
!					psi(i,j) = -psi(ii,jj)*f_out/f_in 
!				else
!					psi(i,j) = -(1.d0-1.d-1)*psi(ii,jj)
!				endif

				continue

			 endif

          end if
       end do
    end do


end subroutine bc_psi_rho2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho3(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif


	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

		rho(iQ,jQ) = den

		if(tri_type==11) then
			if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/50.d0 ) then
				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iQ,jQ))>min(psic/50.d0,1.d-4)) then
				psi(iQ,jQ) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iQ,jQ) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/50.d0 ) then
					psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iQ,jQ))>min(psic/50.d0,1.d-4)) then
!					psi(iQ,jQ) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho3

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho4(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! note: the code will only get here for the LDX case (tri_type=11)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	integer :: zone
	real(kind=dkind) :: max_dev = 2.d-2
	integer :: itemax = 10000

	den_in = dofpsi(psi_in/psic)
	den_out = dofpsi(psi_out/psic)


	! first pass

	do i=1,m_bound

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		zone = ind_bound(i,5)

		if( (zone==1).or.(zone==2) ) then
		! assign psi on the outer boundary

			psi_val = psi_out
			den = den_out

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

  			if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic*max_dev ) then
				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic*max_dev)
			endif

		elseif( (zone==3).or.(zone==4) ) then
		! assign grad psi on the inner boundary

			den = den_in

			psiP =  psi(iloc,jloc)*fi(1) + psi(iloc+1,jloc)*fi(2) +  &
						psi(iloc+1,jloc+1)*fi(3) + psi(iloc,jloc+1)*fi(4) 

			psi(iQ,jQ) = psiP - ( coord_bound(i,4) + coord_bound(i,5) ) *  &
										coord_bound(i,6)

		else
		!this should never happen
			print*, 'error in bc_psi_rho4'
			pause
			stop
		endif

		rho(iQ,jQ) = den

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if(coord_bound(i,3)==1.d0) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
			! assign psi on the outer boundary

				psi_val = psi_out

						!--------------!
				if( (iloc==0).or.(jloc==0) ) then
					psiloc(1) = psi_val
				else
					psiloc(1) = psi(iloc,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==0) ) then
					psiloc(2) = psi_val
				else
					psiloc(2) = psi(iloc+1,jloc)
				endif
						!--------------!
				if( (iloc==nx).or.(jloc==nz) ) then
					psiloc(3) = psi_val
				else
					psiloc(3) = psi(iloc+1,jloc+1)
				endif
						!--------------!
				if( (iloc==0).or.(jloc==nz) ) then
					psiloc(4) = psi_val
				else
					psiloc(4) = psi(iloc,jloc+1)
				endif
						!--------------!

				psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
							psiloc(3)*fi(3) + psiloc(4)*fi(4) 

				psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

   				if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic*max_dev ) then
					psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic*max_dev)
				endif

			elseif( (zone==3).or.(zone==4) ) then
			! assign grad psi on the inner boundary

				psiP =  psi(iloc,jloc)*fi(1) + psi(iloc+1,jloc)*fi(2) +  &
							psi(iloc+1,jloc+1)*fi(3) + psi(iloc,jloc+1)*fi(4) 

				psi(iQ,jQ) = psiP - ( coord_bound(i,4) + coord_bound(i,5) ) *  &
											coord_bound(i,6)

				if(dabs((psi(iQ,jQ)-psiold)/(psi(iq,jq)+psiold))>max_dev) then

					psi(iQ,jQ) = sign(1.d0,psi(iQ,jQ)) *  &
							min(abs(psiold),abs(psi(iQ,jQ)))*(1.d0+max_dev)/(1.d0-max_dev)

				endif

			else
			!this should never happen
				print*, 'error in bc_psi_rho4'
				pause
				stop
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			continue

		enddo

		if(error<2.d-12) exit
!		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho4

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho5(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! note: the code will only get here for the LDX case (tri_type=11)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	integer :: zone
	real(kind=dkind) :: max_dev = 2.d-2
	integer :: itemax = 10000

	den_in = dofpsi(psi_in/psic)
	den_out = dofpsi(psi_out/psic)

	! first pass

	do i=1,m_bound

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		zone = ind_bound(i,5)

		if( (zone==1).or.(zone==2) ) then
		! assign psi on the outer boundary

			psi_val = psi_out
			den = den_out

		elseif( (zone==3).or.(zone==4) ) then
		! assign psi on the inner boundary

			psi_val = psi_in
			den = den_in

		else
		!this should never happen
			print*, 'error in bc_psi_rho5'
			pause
			stop
		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!  		if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic*max_dev ) then
  		if( abs(psi(iQ,jQ)-psi_val)>psic*max_dev ) then
!			psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic*max_dev)
			psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val+psic*max_dev)
		endif

		rho(iQ,jQ) = den

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if(coord_bound(i,3)==1.d0) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			if((iQ==61).and.(jQ==39)) then
				continue
			endif

			psiold = psi(iQ,jQ)

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
			! assign psi on the outer boundary

				psi_val = psi_out

			elseif( (zone==3).or.(zone==4) ) then
			! assign psi on the inner boundary

				psi_val = psi_in

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!  			if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic*max_dev ) then
  			if( abs(psi(iQ,jQ)-psi_val)>psic*max_dev ) then
!				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic*max_dev)
				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val+psic*max_dev)
			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			continue

		enddo

		if(error<2.d-12) exit
!		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho5

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho7(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP
	real(kind=dkind) :: th
	integer :: itemax = 10000
	integer :: imin, imax, istep

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	! For this case there should not be any need to interpolate
!	if(nx<bc_switch) then
	! no interpolation

		do j=1, nz
		   do i=1, nx

			  if(sort_grid(i,j,0)<=0) then

				th = atan2(z_coord(j),x_coord(i)-rmajor)
				if(th<0.d0) th = th + 2.d0*pi

				 rho(i,j) = den
				 psi(i,j) = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )
				 !note that this is not accurate, but it should do for the coarse grids

				 continue

			  end if

		   end do
		end do

		continue
		return

!	endif

! 	this would be "else"
! 	first pass
! 
! 	do i=1,m_bound
! 
! 		psi_val = coord_bound(i,6)
! 
! 		xQ = coord_bound(i,1)
! 		zQ = coord_bound(i,2)
! 
! 		iQ = ind_bound(i,1)
! 		jQ = ind_bound(i,2)
! 		iloc = ind_bound(i,3)
! 		jloc = ind_bound(i,4)
! 
! 		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
! 		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
! 		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
! 		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)
! 
! 				--------------!
! 		if( (iloc==0).or.(jloc==0) ) then
! 			psiloc(1) = psi_val
! 		else
! 			psiloc(1) = psi(iloc,jloc)
! 		endif
! 				--------------!
! 		if( (iloc==nx).or.(jloc==0) ) then
! 			psiloc(2) = psi_val
! 		else
! 			psiloc(2) = psi(iloc+1,jloc)
! 		endif
! 				--------------!
! 		if( (iloc==nx).or.(jloc==nz) ) then
! 			psiloc(3) = psi_val
! 		else
! 			psiloc(3) = psi(iloc+1,jloc+1)
! 		endif
! 				--------------!
! 		if( (iloc==0).or.(jloc==nz) ) then
! 			psiloc(4) = psi_val
! 		else
! 			psiloc(4) = psi(iloc,jloc+1)
! 		endif
! 				--------------!
! 
! 		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
! 					psiloc(3)*fi(3) + psiloc(4)*fi(4) 
! 
! 		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
! 
! 		rho(iQ,jQ) = den
! 
! 		if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/2.5d0 ) then
! 			psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/2.5d0)	! 0.d0 !
! 		endif
! 
! 		continue
! 
! 	enddo
! 
! 	second pass
! 	error = 0.d0
! 
! 	k = 0
! 
! 	imin = 1
! 	imax = m_bound
! 	istep = 1
! 
! 	do
! 
! 		k = k+1
! 
! 		do i=imin, imax, istep
! 
! 			if(coord_bound(i,3)==1.d0) cycle !(this is an "internal" point)
! 
! 			psi_val = coord_bound(i,6)
! 
! 			xQ = coord_bound(i,1)
! 			zQ = coord_bound(i,2)
! 
! 			iQ = ind_bound(i,1)
! 			jQ = ind_bound(i,2)
! 			iloc = ind_bound(i,3)
! 			jloc = ind_bound(i,4)
! 
! 			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
! 			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
! 			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
! 			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)
! 
! 			psiold = psi(iQ,jQ)
! 
! 					--------------!
! 			if( (iloc==0).or.(jloc==0) ) then
! 				psiloc(1) = psi_val
! 			else
! 				psiloc(1) = psi(iloc,jloc)
! 			endif
! 					--------------!
! 			if( (iloc==nx).or.(jloc==0) ) then
! 				psiloc(2) = psi_val
! 			else
! 				psiloc(2) = psi(iloc+1,jloc)
! 			endif
! 					--------------!
! 			if( (iloc==nx).or.(jloc==nz) ) then
! 				psiloc(3) = psi_val
! 			else
! 				psiloc(3) = psi(iloc+1,jloc+1)
! 			endif
! 					--------------!
! 			if( (iloc==0).or.(jloc==nz) ) then
! 				psiloc(4) = psi_val
! 			else
! 				psiloc(4) = psi(iloc,jloc+1)
! 			endif
! 					--------------!
! 
! 			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
! 						psiloc(3)*fi(3) + psiloc(4)*fi(4) 
! 
! 			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))
! 
! 			if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/2.5d0 ) then
! 					psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/2.5d0)	! 0.d0 !
! 			endif
! 
! 			if(tri_type==11) then
! 				if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/50.d0 ) then
! 					psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
! 				endif
! 			endif
! 
! 			if(Broot==3) then
! 
! 				if(abs(psi(iQ,jQ))>min(psic/50.d0,1.d-4)) then
! 					psi(iQ,jQ) = -min(psic/50.d0,1.d-4)	! 0.d0 !
! 				endif
! 
! 			endif
! 
! 			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )
! 
! 			continue
! 
! 		enddo
! 
! 	if(error<2.d-12) exit
! 		if(error<1.d-13) exit
! 		if( (k>itemax).and.(itemax<1000001) ) then
! 			print*, 'bc not converged, error:', error
! 			exit
! 		endif
! 
! 		if (imin == 1) then
! 
! 			imin = m_bound
! 			imax = 1
! 			istep = -1
! 
! 		else
! 
! 			imin = 1
! 			imax = m_bound
! 			istep = 1
! 
! 		endif
! 
! 		error = 0.d0
! 
! 	enddo

	continue

	return

end subroutine bc_psi_rho7

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho1_bis(psi,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi
    real (kind=dkind) :: ex,ez
    integer :: i,j

	psi_ext = psic*(1.d0-0.5d0*x_size/a_elps)
	psi_ext = -.1d0

    do j=1, nz
       do i=1, nx

		  call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)

		  if((ex*ex + ez*ez) >= radius_ext**2) then

             psi(i,j) = psi_ext

          end if

       end do
    end do

end subroutine bc_psi_rho1_bis


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho13(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! I messed up with P and Q, they are inverted here with respect to the initializaion
! let's agree that "Q" is the interpolated point and "P" is the point on the boundary (to be verified in the inizialization)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc, fi_x, fi_z
    real (kind=dkind) :: xQ, zQ, xP, zP, dist, dist0, coordQ(1:2)
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iP, jP
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiQ, rhoQ, sum_RHS, cpoint
	integer :: P_ind, i_diff, j_diff
	real(kind=dkind) :: cx, cz, d
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(clean_edge_rho) then

		do j = 1, nz
		do i = 1,nx

			if(sort_grid(i,j,0)<0) then
				rho(i,j) = den
			endif

		enddo
		enddo

		clean_edge_rho = .false.

	endif

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)
			psi(iP,jP) = psi_val
			rho(iP,jP) = den

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
		zQ = coord_bound(i,2)

		coordQ(1) = bound_13(i,5)
		coordQ(2) = bound_13(i,6)

		iP = ind_bound(i,1)
		jP = ind_bound(i,2)

		dist0 = coord_bound(i,4) !distance from point to surface (P)
		dist = coord_bound(i,5) !distance from surface to interpolated point (Q)
		d = dist+dist0

		iloc = ind_bound(i,3)	! this is the corner point
		jloc = ind_bound(i,4)

		! normal unit vector and "P" point coordinate and index
		cx = bound_13(i,3)
		cz = bound_13(i,4)

		P_ind = nint(bound_13(i,0))	! use nint just to be safe
		xP = bound_13(i,1)
		zP = bound_13(i,2)
		! end of "P" point stuff

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		fi_x(1) = -0.25d0*(1.d0-zQ)/(dx/2.d0)
		fi_x(2) = 0.25d0*(1.d0-zQ)/(dx/2.d0)
		fi_x(3) = 0.25d0*(1.d0+zQ)/(dx/2.d0)
		fi_x(4) = -0.25d0*(1.d0+zQ)/(dx/2.d0)

		fi_z(1) = -0.25d0*(1.d0-xQ)/(dz/2.d0)
		fi_z(2) = -0.25d0*(1.d0+xQ)/(dz/2.d0)
		fi_z(3) = 0.25d0*(1.d0+xQ)/(dz/2.d0)
		fi_z(4) = 0.25d0*(1.d0-xQ)/(dz/2.d0)

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!		rho(iP,jP) = den

		! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = den
		else
			psiloc(1) = rho(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = den
		else
			psiloc(2) = rho(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = den
		else
			psiloc(3) = rho(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = den
		else
			psiloc(4) = rho(iloc,jloc+1)
		endif
				!--------------!

		rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4)

		sum_RHS = 0.d0

		do j = 1, 4

			if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

		enddo

		if(P_ind<0) then
		! easy case, internal Q point

			rho(iP,jP) = rhoQ - abs(sum_RHS) * d
			rho(iP,jP) = rhoQ + sum_RHS * d
			continue
			if(rho(iP,jP)>rhoQ) then
				continue
			endif

		else
		! need to take into account the value in P

			rho(iP,jP) = (rhoQ + sum_RHS*d) / (1.d0 - psiloc(P_ind) * (cx*fi_x(P_ind) + cz*fi_z(P_ind))*d)
			if(rho(iP,jP)>rhoQ) then
				continue
			endif

		endif

		! avoid negative densities
		if(rho(iP,jP)<0.d0) then
			rho(iP,jP) = rhoQ/2.d0
		endif

!!$		! we set grad n_den = 0, so n_denQ = n_denP
!!$		rho(iP,jP) = rhoQ
!!$		! NEED TO CHECK IF THIS IS RIGHT!!!!

		if(tri_type==11) then
			if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
				psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
				psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iP,jP) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)

			dist0 = coord_bound(i,4) !distance from point to surface (P)
			dist = coord_bound(i,5) !distance from surface to interpolated point (Q)
			d = dist+dist0

			iloc = ind_bound(i,3)	! this is the corner point
			jloc = ind_bound(i,4)

			! normal unit vector and "P" point coordinate and index
			cx = bound_13(i,3)
			cz = bound_13(i,4)

			P_ind = nint(bound_13(i,0))	! use nint just to be safe
			xP = bound_13(i,1)
			zP = bound_13(i,2)
			! end of "P" point stuff

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			fi_x(1) = -0.25d0*(1.d0-zQ)/dx
			fi_x(2) = 0.25d0*(1.d0-zQ)/dx
			fi_x(3) = 0.25d0*(1.d0+zQ)/dx
			fi_x(4) = -0.25d0*(1.d0+zQ)/dx

			fi_z(1) = -0.25d0*(1.d0-xQ)/dz
			fi_z(2) = -0.25d0*(1.d0+xQ)/dz
			fi_z(3) = 0.25d0*(1.d0+xQ)/dz
			fi_z(4) = 0.25d0*(1.d0-xQ)/dz

			psiold = psi(iP,jP)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
					psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
!					psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iP,jP)-psiold)/psiold ) )

			! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

			psiold = rho(iP,jP)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = rho(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = rho(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = rho(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = rho(iloc,jloc+1)
			endif
					!--------------!

			rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			sum_RHS = 0.d0

			do j = 1, 4

				if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

			enddo

			if(P_ind<0) then
			! easy case, internal Q point

				rho(iP,jP) = rhoQ - abs(sum_RHS) * d
				if(rho(iP,jP)>rhoQ) then
					continue
				endif

			else
			! need to take into account the value in P

				rho(iP,jP) = (rhoQ + sum_RHS*d) / (1.d0 - psiloc(P_ind) * (cx*fi_x(P_ind) + cz*fi_z(P_ind))*d)
				if(rho(iP,jP)>rhoQ) then
					continue
				endif

			endif

			! avoid negative densities
			if(rho(iP,jP)<0.d0) then
				rho(iP,jP) = rhoQ/2.d0
			endif

			error = max( error,abs( (rho(iP,jP)-psiold)/psiold ) )

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho13


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho23(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! I messed up with P and Q, they are inverted here with respect to the initializaion
! let's agree that "Q" is the interpolated point and "P" is the point on the boundary (to be verified in the inizialization)
! this version uses triangular elements for density extrapolation (only one pass needed)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc
	real (kind=dkind), dimension(1:3) :: fi3, fi3_x, fi3_z
    real (kind=dkind) :: xQ, zQ, xP, zP, dist, dist0, coordQ(1:2), coordQ2(1:2)
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iP, jP
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiQ, rhoQ, sum_RHS, cpoint
	integer :: P_ind, i_diff, j_diff
	integer :: ind_tri(1:3,1:2)
	real(kind=dkind) :: cx, cz, d
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(clean_edge_rho_psi) then

		do j = 1, nz
		do i = 1,nx

			if(sort_grid(i,j,0)<0) then
				rho(i,j) = den
				psi(i,j) = 0.d0
			endif

		enddo
		enddo

		clean_edge_rho_psi = .false.

	endif

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)
			psi(iP,jP) = psi_val
			rho(iP,jP) = den

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
		zQ = coord_bound(i,2)

		coordQ(1) = bound_13(i,5)
		coordQ(2) = bound_13(i,6)

		iP = ind_bound(i,1)
		jP = ind_bound(i,2)

		dist0 = coord_bound(i,4) !distance from point to surface (P)
		dist = coord_bound(i,5) !distance from surface to interpolated point (Q)
		d = dist+dist0

		iloc = ind_bound(i,3)	! this is the corner point
		jloc = ind_bound(i,4)

		! normal unit vector and "P" point coordinate and index
		cx = bound_13(i,3)
		cz = bound_13(i,4)

		P_ind = nint(bound_13(i,0))	! use nint just to be safe
		xP = bound_13(i,1)
		zP = bound_13(i,2)
		! end of "P" point stuff

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!		rho(iP,jP) = den

		! update rho extrapolating rho at the boundary
		! us point Q2, on the P-boundary normal, but in a triangular
		! element only containing internal points

		do j = 1, 3
		do k = 1, 2
			ind_tri(j,k) = bound_tri(i,j,k)
		enddo
		enddo

		xT1 = x_coord(ind_tri(1,1))
		xT2 = x_coord(ind_tri(2,1))
		xT3 = x_coord(ind_tri(3,1))

		zT1 = z_coord(ind_tri(1,2))
		zT2 = z_coord(ind_tri(2,2))
		zT3 = z_coord(ind_tri(3,2))

		two_A_tri = 2.d0 * coord_tri(i,4)

		do j = 1, 3

			psiloc(j) = rho(ind_tri(j,1),ind_tri(j,2))

		enddo

		coordQ2(1) = coord_tri(i,1)
		coordQ2(2) = coord_tri(i,2)

		fi3(1) = N1_tri(coordQ2(1),coordQ2(2))
		fi3(2) = N2_tri(coordQ2(1),coordQ2(2))
		fi3(3) = N3_tri(coordQ2(1),coordQ2(2))

		fi3_x(1) = N1_tri_x(coordQ2(1),coordQ2(2))
		fi3_x(2) = N2_tri_x(coordQ2(1),coordQ2(2))
		fi3_x(3) = N3_tri_x(coordQ2(1),coordQ2(2))

		fi3_z(1) = N1_tri_z(coordQ2(1),coordQ2(2))
		fi3_z(2) = N2_tri_z(coordQ2(1),coordQ2(2))
		fi3_z(3) = N3_tri_z(coordQ2(1),coordQ2(2))

		rhoQ =  psiloc(1)*fi3(1) + psiloc(2)*fi3(2) +  &
					psiloc(3)*fi3(3)

		sum_RHS = 0.d0

		do j = 1, 3
		! since all points in the triangleare internal, P cannot be one of them

			sum_RHS = sum_RHS + psiloc(j) * (cx*fi3_x(j) + cz*fi3_z(j))

		enddo

		! easy case by construction, internal Q2 point

		rho(iP,jP) = rhoQ - abs(sum_RHS) * coord_tri(i,3)
		rho(iP,jP) = rhoQ + sum_RHS * coord_tri(i,3)
		continue

		if(rho(iP,jP)>rhoQ) then
			continue
		endif

		! avoid negative densities (this could happen at the start of a new grid)
		if(rho(iP,jP)<0.d0) then
			rho(iP,jP) = rhoQ/2.d0
		endif

		if(tri_type==11) then
			if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
				psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
				psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iP,jP) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass (only psi needs to be updated)
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)

			dist0 = coord_bound(i,4) !distance from point to surface (P)
			dist = coord_bound(i,5) !distance from surface to interpolated point (Q)
			d = dist+dist0

			iloc = ind_bound(i,3)	! this is the corner point
			jloc = ind_bound(i,4)

			! normal unit vector and "P" point coordinate and index
			cx = bound_13(i,3)
			cz = bound_13(i,4)

			P_ind = nint(bound_13(i,0))	! use nint just to be safe
			xP = bound_13(i,1)
			zP = bound_13(i,2)
			! end of "P" point stuff

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iP,jP)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
					psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
!					psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iP,jP)-psiold)/psiold ) )

			! no more need to update rho

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho23

!--------------------------------------FE ruotins for interpolations--------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N1_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT2*zT3-xT3*zT2) + (zT2-zT3)*x + (xT3-xT2)*z
	answer = answer / two_A_tri

end function N1_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N2_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT3*zT1-xT1*zT3) + (zT3-zT1)*x + (xT1-xT3)*z
	answer = answer / two_A_tri

end function N2_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N3_tri(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT1*zT2-xT2*zT1) + (zT1-zT2)*x + (xT2-xT1)*z
	answer = answer / two_A_tri

end function N3_tri

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N1_tri_x(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (zT2-zT3)
	answer = answer / two_A_tri

end function N1_tri_x

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N2_tri_x(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (zT3-zT1)
	answer = answer / two_A_tri

end function N2_tri_x

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N3_tri_x(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (zT1-zT2)
	answer = answer / two_A_tri

end function N3_tri_x

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N1_tri_z(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT3-xT2)
	answer = answer / two_A_tri

end function N1_tri_z

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N2_tri_z(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT1-xT3)
	answer = answer / two_A_tri

end function N2_tri_z

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function N3_tri_z(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind) :: x, z, answer

	answer = (xT2-xT1)
	answer = answer / two_A_tri

end function N3_tri_z

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho13_old3(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! I messed up with P and Q, they are inverted here with respect to the initializaion
! let's agree that "Q" is the interpolated point and "P" is the point on the boundary (to be verified in the inizialization)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc, fi_x, fi_z, fi_x_P, fi_z_P
    real (kind=dkind) :: xQ, zQ, xP, zP, dist, dist0, coordQ(1:2)
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iP, jP
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiQ, rhoQ, sum_RHS, cpoint
	integer :: P_ind, i_diff, j_diff
	real(kind=dkind) :: cx, cz
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(clean_edge_rho) then

		do j = 1, nz
		do i = 1,nx

			if(sort_grid(i,j,0)<0) then
				rho(i,j) = den
			endif

		enddo
		enddo

		clean_edge_rho = .false.

	endif

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)
			psi(iP,jP) = psi_val
			rho(iP,jP) = den

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
		zQ = coord_bound(i,2)

		coordQ(1) = bound_13(i,5)
		coordQ(2) = bound_13(i,6)

		iP = ind_bound(i,1)
		jP = ind_bound(i,2)

		dist0 = coord_bound(i,4) !distance from point to surface (P)
		dist = coord_bound(i,5) !distance from surface to interpolated point (Q)

		iloc = ind_bound(i,3)	! this is the corner point
		jloc = ind_bound(i,4)

		! normal unit vector and "P" point coordinate and index
		cx = bound_13(i,3)
		cz = bound_13(i,4)

		P_ind = nint(bound_13(i,0))	! use nint just to be safe
		xP = bound_13(i,1)
		zP = bound_13(i,2)
		! end of "P" point stuff

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		fi_x(1) = -0.25d0*(1.d0-zQ)
		fi_x(2) = 0.25d0*(1.d0-zQ)
		fi_x(3) = 0.25d0*(1.d0+zQ)
		fi_x(4) = -0.25d0*(1.d0+zQ)

		fi_z(1) = -0.25d0*(1.d0-xQ)
		fi_z(2) = -0.25d0*(1.d0+xQ)
		fi_z(3) = 0.25d0*(1.d0+xQ)
		fi_z(4) = 0.25d0*(1.d0-xQ)

		if(P_ind>0) then
		! also need derivative in P

			fi_x_P(1) = -0.25d0*(1.d0-zP)
			fi_x_P(2) = 0.25d0*(1.d0-zP)
			fi_x_P(3) = 0.25d0*(1.d0+zP)
			fi_x_P(4) = -0.25d0*(1.d0+zP)

			fi_z_P(1) = -0.25d0*(1.d0-xP)
			fi_z_P(2) = -0.25d0*(1.d0+xP)
			fi_z_P(3) = 0.25d0*(1.d0+xP)
			fi_z_P(4) = 0.25d0*(1.d0-xP)

		endif

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!		rho(iP,jP) = den

		! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = den
		else
			psiloc(1) = rho(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = den
		else
			psiloc(2) = rho(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = den
		else
			psiloc(3) = rho(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = den
		else
			psiloc(4) = rho(iloc,jloc+1)
		endif
				!--------------!

		rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4)

		sum_RHS = 0.d0

		do j = 1, 4

			if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

		enddo

		if(P_ind<0) then
		! easy case, internal point

			rho(iP,jP) = rhoQ - abs(sum_RHS) * (dist+dist0)

		else
		! need to take into account the value in P

			do j = 1, 4

				if(j/=P_ind) sum_RHS = sum_RHS - psiloc(j) * (cx*fi_x_P(j) + cz*fi_z_P(j))

			enddo

			rho(iP,jP) = sum_RHS / (cx*(fi_x_P(P_ind)-fi_x(P_ind)) + cz*(fi_z_P(P_ind)-fi_z(P_ind)))

		endif

!!$		! we set grad n_den = 0, so n_denQ = n_denP
!!$		rho(iP,jP) = rhoQ
!!$		! NEED TO CHECK IF THIS IS RIGHT!!!!

		if(tri_type==11) then
			if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
				psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
				psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iP,jP) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)

			dist0 = coord_bound(i,4) !distance from point to surface (P)
			dist = coord_bound(i,5) !distance from surface to interpolated point (Q)

			iloc = ind_bound(i,3)	! this is the corner point
			jloc = ind_bound(i,4)

			! normal unit vector and "P" point coordinate and index
			cx = bound_13(i,3)
			cz = bound_13(i,4)

			P_ind = nint(bound_13(i,0))	! use nint just to be safe
			xP = bound_13(i,1)
			zP = bound_13(i,2)
			! end of "P" point stuff

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			fi_x(1) = -0.25d0*(1.d0-zQ)
			fi_x(2) = 0.25d0*(1.d0-zQ)
			fi_x(3) = 0.25d0*(1.d0+zQ)
			fi_x(4) = -0.25d0*(1.d0+zQ)

			fi_z(1) = -0.25d0*(1.d0-xQ)
			fi_z(2) = -0.25d0*(1.d0+xQ)
			fi_z(3) = 0.25d0*(1.d0+xQ)
			fi_z(4) = 0.25d0*(1.d0-xQ)

			if(P_ind>0) then
			! also need derivative in P

				fi_x_P(1) = -0.25d0*(1.d0-zP)
				fi_x_P(2) = 0.25d0*(1.d0-zP)
				fi_x_P(3) = 0.25d0*(1.d0+zP)
				fi_x_P(4) = -0.25d0*(1.d0+zP)

				fi_z_P(1) = -0.25d0*(1.d0-xP)
				fi_z_P(2) = -0.25d0*(1.d0+xP)
				fi_z_P(3) = 0.25d0*(1.d0+xP)
				fi_z_P(4) = 0.25d0*(1.d0-xP)

			endif

			psiold = psi(iP,jP)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
					psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
!					psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iP,jP)-psiold)/psiold ) )

			! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

			psiold = rho(iP,jP)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = rho(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = rho(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = rho(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = rho(iloc,jloc+1)
			endif
					!--------------!

			rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			sum_RHS = 0.d0

			do j = 1, 4

				if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

			enddo

			if(P_ind<0) then
			! easy case, internal point

				rho(iP,jP) = rhoQ - abs(sum_RHS) * (dist+dist0)

			else
			! need to take into account the value in P

				do j = 1, 4

					if(j/=P_ind) sum_RHS = sum_RHS - psiloc(j) * (cx*fi_x_P(j) + cz*fi_z_P(j))

				enddo

				rho(iP,jP) = sum_RHS / (cx*(fi_x_P(P_ind)-fi_x(P_ind)) + cz*(fi_z_P(P_ind)-fi_z(P_ind)))

			endif

			error = max( error,abs( (rho(iP,jP)-psiold)/psiold ) )

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho13_old3

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho13_old2(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! I messed up with P and Q, they are inverted here with respect to the initializaion
! let's agree that "Q" is the interpolated point and "P" is the point on the boundary (to be verified in the inizialization)

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi, psiloc, fi_x, fi_z, fi_x_P, fi_z_P
    real (kind=dkind) :: xQ, zQ, xP, zP
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iP, jP
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiQ, rhoQ, sum_RHS, cpoint
	integer :: P_ind, i_diff, j_diff
	real(kind=dkind) :: cx, cz, RR_Q, zz_Q, RR_P, zz_P, dist, dist0
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)
			psi(iP,jP) = psi_val
			rho(iP,jP) = den

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)	! this is officially "Q" (the interpolated point)
		zQ = coord_bound(i,2)

		iP = ind_bound(i,1)
		jP = ind_bound(i,2)

		dist0 = coord_bound(i,4) !distance from point to surface (P)
		dist = coord_bound(i,5) !distance from surface to interpolated point (Q)

		RR_P = x_coord(iP)
		zz_P = z_coord(jP)

!		RR_Q = RR_P + (xQ-RR_P)/(1.d0 + dist/dist0)
!		zz_Q = zz_P + (zQ-zz_P)/(1.d0 + dist/dist0)
		! this apparently should work, but it doesn't: I'd like to know why...

		RR_Q = RR_P+(1.d0+xQ)*dx/2.d0
		zz_Q = zz_P+(1.d0+zQ)*dz/2.d0

		cx = sqrt((RR_Q-RR_P)**2)/(dist+dist0)
		cz = sqrt((zz_Q-zz_P)**2)/(dist+dist0)

		iloc = ind_bound(i,3)	! this is the corner point
		jloc = ind_bound(i,4)

		i_diff = iP - iloc
		j_diff = jP - jloc

		! this determines which point of the cell corresponds to P
		! if the point is internal obviously it's not a cell point
		if((i_diff==0).and.(j_diff==0)) then
			P_ind = 1
			xP = -1.d0
			zP = -1.d0
		elseif((i_diff==1).and.(j_diff==0)) then
			P_ind = 2
			xP = 1.d0
			zP = -1.d0
		elseif((i_diff==1).and.(j_diff==1)) then
			P_ind = 3
			xP = 1.d0
			zP = 1.d0
		elseif((i_diff==0).and.(j_diff==1)) then
			P_ind = 4
			xP = -1.d0
			zP = 1.d0
		else
			P_ind = -1
		endif

		cx = (RR_P-RR_q)/(dist+dist0)
		cz = (zz_P-zz_q)/(dist+dist0)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		fi_x(1) = -0.25d0*(1.d0-zQ)
		fi_x(2) = 0.25d0*(1.d0-zQ)
		fi_x(3) = 0.25d0*(1.d0+zQ)
		fi_x(4) = -0.25d0*(1.d0+zQ)

		fi_z(1) = -0.25d0*(1.d0-xQ)
		fi_z(2) = -0.25d0*(1.d0+xQ)
		fi_z(3) = 0.25d0*(1.d0+xQ)
		fi_z(4) = 0.25d0*(1.d0-xQ)

		if(P_ind>0) then
		! also need derivative in P

			fi_x_P(1) = -0.25d0*(1.d0-zP)
			fi_x_P(2) = 0.25d0*(1.d0-zP)
			fi_x_P(3) = 0.25d0*(1.d0+zP)
			fi_x_P(4) = -0.25d0*(1.d0+zP)

			fi_z_P(1) = -0.25d0*(1.d0-xP)
			fi_z_P(2) = -0.25d0*(1.d0+xP)
			fi_z_P(3) = 0.25d0*(1.d0+xP)
			fi_z_P(4) = 0.25d0*(1.d0-xP)

		endif

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!		rho(iP,jP) = den

		! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = den
		else
			psiloc(1) = rho(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = den
		else
			psiloc(2) = rho(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = den
		else
			psiloc(3) = rho(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = den
		else
			psiloc(4) = rho(iloc,jloc+1)
		endif
				!--------------!

		rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4)

		sum_RHS = 0.d0

		do j = 1, 4

			if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

		enddo

		if(P_ind<0) then
		! easy case, internal point

			rho(iP,jP) = rhoQ - abs(sum_RHS) * (dist+dist0)

		else
		! need to take into account the value in P

			do j = 1, 4

				if(j/=P_ind) sum_RHS = sum_RHS - psiloc(j) * (cx*fi_x_P(j) + cz*fi_z_P(j))

			enddo

			rho(iP,jP) = sum_RHS / (cx*(fi_x_P(j)-fi_x(j)) + cz*(fi_z_P(j)-fi_z(j)))

		endif

!!$		! we set grad n_den = 0, so n_denQ = n_denP
!!$		rho(iP,jP) = rhoQ
!!$		! NEED TO CHECK IF THIS IS RIGHT!!!!

		if(tri_type==11) then
			if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
				psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
				psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iP,jP) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iP = ind_bound(i,1)
			jP = ind_bound(i,2)

			dist0 = coord_bound(i,4) !distance from point to surface (P)
			dist = coord_bound(i,5) !distance from surface to interpolated point (Q)

			RR_P = x_coord(iP)
			zz_P = z_coord(jP)

	!		RR_Q = RR_P + (xQ-RR_P)/(1.d0 + dist/dist0)
	!		zz_Q = zz_P + (zQ-zz_P)/(1.d0 + dist/dist0)
			! this apparently should work, but it doesn't: I'd like to know why...

			RR_Q = RR_P+(1.d0+xQ)*dx/2.d0
			zz_Q = zz_P+(1.d0+zQ)*dz/2.d0

			cx = sqrt((RR_Q-RR_P)**2)/(dist+dist0)
			cz = sqrt((zz_Q-zz_P)**2)/(dist+dist0)

			iloc = ind_bound(i,3)	! this is the corner point
			jloc = ind_bound(i,4)

			i_diff = iP - iloc
			j_diff = jP - jloc

			! this determines which point of the cell corresponds to P
			! if the point is internal obviously it's not a cell point
			if((i_diff==0).and.(j_diff==0)) then
				P_ind = 1
				xP = -1.d0
				zP = -1.d0
			elseif((i_diff==1).and.(j_diff==0)) then
				P_ind = 2
				xP = 1.d0
				zP = -1.d0
			elseif((i_diff==1).and.(j_diff==1)) then
				P_ind = 3
				xP = 1.d0
				zP = 1.d0
			elseif((i_diff==0).and.(j_diff==1)) then
				P_ind = 4
				xP = -1.d0
				zP = 1.d0
			else
				P_ind = -1
			endif

			cx = (RR_P-RR_q)/(dist+dist0)
			cz = (zz_P-zz_q)/(dist+dist0)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			fi_x(1) = -0.25d0*(1.d0-zQ)
			fi_x(2) = 0.25d0*(1.d0-zQ)
			fi_x(3) = 0.25d0*(1.d0+zQ)
			fi_x(4) = -0.25d0*(1.d0+zQ)

			fi_z(1) = -0.25d0*(1.d0-xQ)
			fi_z(2) = -0.25d0*(1.d0+xQ)
			fi_z(3) = 0.25d0*(1.d0+xQ)
			fi_z(4) = 0.25d0*(1.d0-xQ)

			if(P_ind>0) then
			! also need derivative in P

				fi_x_P(1) = -0.25d0*(1.d0-zP)
				fi_x_P(2) = 0.25d0*(1.d0-zP)
				fi_x_P(3) = 0.25d0*(1.d0+zP)
				fi_x_P(4) = -0.25d0*(1.d0+zP)

				fi_z_P(1) = -0.25d0*(1.d0-xP)
				fi_z_P(2) = -0.25d0*(1.d0+xP)
				fi_z_P(3) = 0.25d0*(1.d0+xP)
				fi_z_P(4) = 0.25d0*(1.d0-xP)

			endif

			psiold = psi(iP,jP)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iP,jP) = psiQ + (psi_val - psiQ) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iP,jP))-abs(psi_val))>psic/50.d0 ) then
					psi(iP,jP) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iP,jP))>min(psic/50.d0,1.d-4)) then
!					psi(iP,jP) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iP,jP)-psiold)/psiold ) )

			! update rho extrapolating rho at the boundary (set grad rho in P = grad rho on Q)

			psiold = rho(iP,jP)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = rho(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = rho(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = rho(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = rho(iloc,jloc+1)
			endif
					!--------------!

			rhoQ =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			sum_RHS = 0.d0

			do j = 1, 4

				if(j/=P_ind) sum_RHS = sum_RHS + psiloc(j) * (cx*fi_x(j) + cz*fi_z(j))

			enddo

			if(P_ind<0) then
			! easy case, internal point

				rho(iP,jP) = rhoQ - abs(sum_RHS) * (dist+dist0)

			else
			! need to take into account the value in P

				do j = 1, 4

					if(j/=P_ind) sum_RHS = sum_RHS - psiloc(j) * (cx*fi_x_P(j) + cz*fi_z_P(j))

				enddo

				rho(iP,jP) = sum_RHS / (cx*(fi_x_P(j)-fi_x(j)) + cz*(fi_z_P(j)-fi_z(j)))

			endif

			error = max( error,abs( (rho(iP,jP)-psiold)/psiold ) )

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho13_old2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bc_psi_rho13_old(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho
    real (kind=dkind), dimension(1:4) :: fi,psiloc
    real (kind=dkind) :: xQ, zQ
    real (kind=dkind) :: den, den_in, den_out
    integer :: i,j,k, iloc, jloc, iQ, jQ
	real(kind=dkind) :: error,psiold
	real(kind=dkind) :: x,ex,ez,psi_val, psiP, rhoP
	integer :: zone
	integer :: itemax = 10000

    den = dofpsi(0.0d0)
	psi_val = 0.d0

	if(tri_type==11) then
	! LDX case

		den_in = dofpsi(psi_in/psic)
		den_out = dofpsi(psi_out/psic)

	endif

	! first pass

	do i=1,m_bound

		if(coord_bound(i,3)==2.d0) then
		! the point is really on the boundary

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			psi(iQ,jQ) = psi_val
			rho(iQ,jQ) = den

			! then of course we can skip the rest
			cycle

		endif

		xQ = coord_bound(i,1)
		zQ = coord_bound(i,2)

		iQ = ind_bound(i,1)
		jQ = ind_bound(i,2)
		iloc = ind_bound(i,3)
		jloc = ind_bound(i,4)

		fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
		fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
		fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
		fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

		if(tri_type==11) then
		! LDX case

			zone = ind_bound(i,5)

			if( (zone==1).or.(zone==2) ) then
				psi_val = psi_out
				den = den_out
			elseif( (zone==3).or.(zone==4) ) then
				psi_val = psi_in
				den = den_in
			else
			!this should never happen
				psi_val = 0.d0
			endif

		endif

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = psi_val
		else
			psiloc(1) = psi(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = psi_val
		else
			psiloc(2) = psi(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = psi_val
		else
			psiloc(3) = psi(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = psi_val
		else
			psiloc(4) = psi(iloc,jloc+1)
		endif
				!--------------!

		psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4) 

		psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

!		rho(iQ,jQ) = den

		! update rho setting grad rho = 0 at the boundary

				!--------------!
		if( (iloc==0).or.(jloc==0) ) then
			psiloc(1) = den
		else
			psiloc(1) = rho(iloc,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==0) ) then
			psiloc(2) = den
		else
			psiloc(2) = rho(iloc+1,jloc)
		endif
				!--------------!
		if( (iloc==nx).or.(jloc==nz) ) then
			psiloc(3) = den
		else
			psiloc(3) = rho(iloc+1,jloc+1)
		endif
				!--------------!
		if( (iloc==0).or.(jloc==nz) ) then
			psiloc(4) = den
		else
			psiloc(4) = rho(iloc,jloc+1)
		endif
				!--------------!

		rhoP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
					psiloc(3)*fi(3) + psiloc(4)*fi(4)

		! we set grad n_den = 0, so n_denQ = n_denP
		rho(iQ,jQ) = rhoP
		! NEED TO CHECK IF THIS IS RIGHT!!!!

		if(tri_type==11) then
			if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/50.d0 ) then
				psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/50.d0)
			endif
		endif

		if(Broot==3) then

			if(abs(psi(iQ,jQ))>min(psic/50.d0,1.d-4)) then
				psi(iQ,jQ) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				psi(iQ,jQ) = 0.d0 !	-psic/50.d0	! 
			endif

		endif

		continue

	enddo

	! second pass
	error = 0.d0
	k = 0

	do

		k = k+1

		do i=1,m_bound

			if((coord_bound(i,3)==1.d0).or.(coord_bound(i,3)==2.d0)) cycle !(this is an "internal" point)

			xQ = coord_bound(i,1)
			zQ = coord_bound(i,2)

			iQ = ind_bound(i,1)
			jQ = ind_bound(i,2)
			iloc = ind_bound(i,3)
			jloc = ind_bound(i,4)

			fi(1) = 0.25d0*(1.d0-xQ)*(1.d0-zQ)
			fi(2) = 0.25d0*(1.d0+xQ)*(1.d0-zQ)
			fi(3) = 0.25d0*(1.d0+xQ)*(1.d0+zQ)
			fi(4) = 0.25d0*(1.d0-xQ)*(1.d0+zQ)

			psiold = psi(iQ,jQ)

			if(tri_type==11) then
			! LDX case

				zone = ind_bound(i,5)

				if( (zone==1).or.(zone==2) ) then
					psi_val = psi_out
				elseif( (zone==3).or.(zone==4) ) then
					psi_val = psi_in
				else
				!this should never happen
					psi_val = 0.d0
				endif

			endif

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = psi_val
			else
				psiloc(1) = psi(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = psi_val
			else
				psiloc(2) = psi(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = psi_val
			else
				psiloc(3) = psi(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = psi_val
			else
				psiloc(4) = psi(iloc,jloc+1)
			endif
					!--------------!

			psiP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4) 

			psi(iQ,jQ) = psiP + (psi_val - psiP) * (1.d0 + coord_bound(i,4)/coord_bound(i,5))

			if(tri_type==11) then
				if( abs(abs(psi(iQ,jQ))-abs(psi_val))>psic/50.d0 ) then
					psi(iQ,jQ) = 2.d0*psi_val-abs(psi_val-psic/50.d0)	! 0.d0 !
				endif
			endif

			if(Broot==3) then

				if(abs(psi(iQ,jQ))>min(psic/50.d0,1.d-4)) then
!					psi(iQ,jQ) = -min(psic/50.d0,1.d-4)	! 0.d0 !
				endif

			endif

			error = max( error,abs( (psi(iQ,jQ)-psiold)/psiold ) )

			! update rho setting grad rho = 0 at the boundary

			psiold = rho(iQ,jQ)

					!--------------!
			if( (iloc==0).or.(jloc==0) ) then
				psiloc(1) = den
			else
				psiloc(1) = rho(iloc,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==0) ) then
				psiloc(2) = den
			else
				psiloc(2) = rho(iloc+1,jloc)
			endif
					!--------------!
			if( (iloc==nx).or.(jloc==nz) ) then
				psiloc(3) = den
			else
				psiloc(3) = rho(iloc+1,jloc+1)
			endif
					!--------------!
			if( (iloc==0).or.(jloc==nz) ) then
				psiloc(4) = den
			else
				psiloc(4) = rho(iloc,jloc+1)
			endif
					!--------------!

			rhoP =  psiloc(1)*fi(1) + psiloc(2)*fi(2) +  &
						psiloc(3)*fi(3) + psiloc(4)*fi(4)

			! we set grad n_den = 0, so n_denQ = n_denP
			rho(iQ,jQ) = rhoP
			! NEED TO CHECK IF THIS IS RIGHT!!!!

			error = max( error,abs( (rho(iQ,jQ)-psiold)/psiold ) )

			continue

		enddo

!		if(error<2.d-12) exit
		if(error<1.d-13) exit
		if( (k>itemax).and.(itemax<1000001) ) then
			print*, 'bc not converged, error:', error
			exit
		endif

		error = 0.d0

	enddo

	continue

	return

end subroutine bc_psi_rho13_old

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function gtheta_bound(th) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this function computes grad psi on the inner boundary
  ! for the LDX case

	real(kind=dkind), intent(in) :: th
	real(kind=dkind) :: answer
	real(kind=dkind) :: smallr, bigR, z, K_ell, E_ell, k2, k, R02, R0
	real(kind=dkind) :: elle, ellf

!!$!	answer = -1.d-2 * (1.d0 + 0.1d0*cos(th))
!!$
!!$	R02 = rmajor**2
!!$
!!$	smallr = dbsval(th, r_ord, r_data(1:theta_points2+r_ord,6),  &
!!$	theta_points2, r_cscoef(2,1:theta_points2) )
!!$	!inner radius
!!$
!!$	bigR = rmajor + smallr*cos(th)
!!$	z = smallr*sin(th)
!!$
!!$	k2 = 4.d0*rmajor*bigR/((rmajor+bigR)**2+z**2)
!!$	k = 2.d0*sqrt(rmajor*bigR/((rmajor+bigR)**2+z**2))
!!$
!!$	K_ell = ellf(0.5d0*pi,k)
!!$	E_ell = elle(0.5d0*pi,k)
!!$
!!$!	answer = mu_mag/(pi*k)*I_coil*sqrt(rmajor*bigR)*  &
!!$!			((1.d0-0.5*k2)*K_ell-E_ell)
!!$
!!$	answer = mu_mag*I_coil/(2.d0*pi)/sqrt((bigR+rmajor)**2+z**2) *  &
!!$			( (bigR*cos(th)+z*sin(th))*K_ell + (bigR*cos(th)-z*sin(th))*  &
!!$				(R02+bigR**2+z**2)/((bigR-rmajor)**2+z**2)*E_ell )
!!$
!!$	R0 = rmajor + delta_coil
!!$	R02 = R0**2
!!$
!!$	bigR = rmajor + smallr*cos(th)
!!$	z = smallr*sin(th)
!!$
!!$	k2 = 4.d0*R0*bigR/((R0+bigR)**2+z**2)
!!$	k = 2.d0*sqrt(R0*bigR/((R0+bigR)**2+z**2))
!!$
!!$	K_ell = ellf(0.5d0*pi,k)
!!$	E_ell = elle(0.5d0*pi,k)
!!$
!!$!	answer = mu_mag/(pi*k)*I_coil*sqrt(rmajor*bigR)*  &
!!$!			((1.d0-0.5*k2)*K_ell-E_ell)
!!$
!!$	answer = mu_mag*I_coil/(2.d0*pi)/sqrt((bigR+R0)**2+z**2) *  &
!!$			( (bigR*cos(th)+z*sin(th))*K_ell + (bigR*cos(th)-z*sin(th))*  &
!!$				(R02+bigR**2+z**2)/((bigR-R0)**2+z**2)*E_ell )
!!$
!!$	answer = -1.d0

	answer = dbsval(th, psiprim_ord, psiprim_data(1:ibreak_pp,3),  &
		psiprim_dim, psiprim_cscoef(1,1:psiprim_dim) )


!	answer = 1.d1*answer

	continue

end function gtheta_bound

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer,apsi
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag) then
		answer = d_loc
		return
	elseif(psi==psic_flag) then
		answer = dc_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_n) then

			apsi=-dabs(psi)

			if(apsi/psic<d_data(1,1)) then
				answer = d_data(1,2)
			else
				answer = dbsval(apsi/psic, d_ord, d_data(1:ibreak_d,3),  &
					ibreak_d-d_ord, d_cscoef(1,1:ibreak_d-d_ord) )
			endif

			if(answer>0.d0) then
				continue
			else
				continue
			endif

			return

		else
			apsi = 0.d0
		endif

    endif

	if (numerical_n) then

		if(apsi/psic<d_data(1,1)) then
			answer = d_data(1,2)
		elseif (apsi>=psic) then
			answer = d_data(ibreak_d-d_ord,2)
		else

			answer = dbsval(apsi/psic, d_ord, d_data(1:ibreak_d,3),  &
					ibreak_d-d_ord, d_cscoef(1,1:ibreak_d-d_ord) )

		endif

	else

!	    if(dabs(psi) <= 0.0d0) then
!	    if(dabs(psi) <= dabs(psic*fraction) ) then
	    if(psi <= psic*fraction ) then
		   answer = dedge
		else if(dabs(psi) > psic) then
	       answer = dcenter
		else
		   answer = dedge + (dcenter - dedge)*  &
				dabs( psi/psic-fraction )**alpha_rho / ( 1.d0-fraction )**alpha_rho
!			answer = dedge + (dcenter - dedge)* dabs( psi/psic )**alpha_rho
		end if

	endif


	continue

end function dofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dddpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
 	integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer,apsi

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag) then
		answer = dp_loc
		return
	elseif(psi==psic_flag) then
		answer = dcp_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_n) then

			apsi=-dabs(psi)

			if(apsi/psic<d_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psic, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(1,:) )
			endif

			answer=answer/psic
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif
    
    if (numerical_n) then

			if ((apsi/psic)<d_data(1,1)) then

				answer = dbsder(1,d_data(1,1), d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(1,:) )

			elseif (apsi>=psic) then

				answer = dbsder(1,1.d0-1.d-6, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(1,:) )

			else

				answer = dbsder(1,apsi/psic, d_ord, d_data(1:ibreak_d,3),  &
								ibreak_d-d_ord, d_cscoef(1,:) )

			endif

			answer=answer/psic

!!$		if(apsi/psic>0.99d0) then
!!$
!!$			call splinderiv(d_data(1:data_dim,1),d_data(1:data_dim,2),  &
!!$							d_data(1:data_dim,3),0.99d0,answer,data_dim)
!!$			answer=answer/psic
!!$
!!$		else
!!$
!!$			call splinderiv(d_data(1:data_dim,1),d_data(1:data_dim,2),  &
!!$							d_data(1:data_dim,3),apsi/psic,answer,data_dim)
!!$			answer=answer/psic
!!$
!!$		endif

	else

!		if(dabs(psi) <= 0.0d0) then
!	    if(dabs(psi) <= dabs(psic*fraction) ) then
	    if(psi <= psic*fraction ) then
		   answer = 0.0d0
	    else if(dabs(psi) > psic) then
		   answer = 0.0d0
	    else
		   answer = alpha_rho*(dcenter - dedge)/psic  &
					*dabs( psi/psic-fraction )**(alpha_rho-1.d0)/( 1.d0-fraction )**alpha_rho
!			answer = alpha_rho*(dcenter - dedge)/psic  &
!					*dabs( psi/psic )**(alpha_rho-1.d0)
		end if

	endif

	if(abs(psi)>0.d0) answer = answer*psi/abs(psi)

	continue

end function dddpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function pofpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer,apsi
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag) then
		answer = p_loc
		return
	elseif(psi==psic_flag) then
		answer = pc_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

    if (zone_loc==-1) then

		if(numerical_p_iso) then

			apsi=-dabs(psi)

			if(apsi/psic<p_iso_data(1,1)) then
				answer = p_iso_data(1,2)
			else
				answer = dbsval(apsi/psic, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
						ibreak_p-p_iso_ord, p_iso_cscoef(1,1:ibreak_p-p_iso_ord) )
			endif

			return

		else
			apsi = 0.d0
		endif

    endif

!!$	answer = pedge + 0.5d0* ( b_phi_zero**2 - bzero(psi)**2 )
!!$	return

	if (tri_type==9) then

		answer = 2.d0/(rmajor**3*q0) * psi   + pedge !

    elseif (numerical_p_iso) then

			if(apsi/psic<p_iso_data(1,1)) then
				answer = p_iso_data(1,2)
			elseif (apsi>=psic) then
				answer = p_iso_data(ibreak_p-p_iso_ord,2)
			else

				answer = dbsval(apsi/psic, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
						ibreak_p-p_iso_ord, p_iso_cscoef(1,1:ibreak_p-p_iso_ord) )

			endif

	elseif(Broot==3) then

		if(dabs(apsi/psic) > 1.0d0) then
			answer =  pcenter
			print*, apsi/psic
!		elseif(apsi/psic > fraction) then
		elseif(psi/psic > fraction) then
			answer = pedge + (pcenter-pedge)* (   apsi/ ((1.d0-fraction)*psic) - &
												fraction/ (1.d0-fraction)   )**alpha
		else
			answer = pedge
		end if

	else

		if (p_opt == 4) then
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

			answer = pcenter * (volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0)

		elseif (p_opt == 5) then
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

			answer = pcenter * (volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0) *  &
								pauxofpsi(psi)

		elseif (p_opt == 6) then
		! p = K(1-psi)/(1+psi)**(gammahut+1)

			answer = pcenter * (1.d0 - apsi/psic) / (1.d0 + apsi/psic)**gammahut *  &
								pauxofpsi(psi)

		elseif (p_opt == 7) then

			if(apsi>psi_pres) then

				answer = 0.d0

			else

				answer = pcenter * (1.d0 - apsi/psi_pres) * (apsi/psi_pres)*gammahut

			endif

		elseif (p_opt == 8) then

			if((apsi>psi_pres).or.(apsi<1.d-6*psic)) then

				answer = 0.d0

			else

				answer = pcenter * (1.d0 - apsi/psi_pres)**aux_fact * (apsi/psi_pres)**gammahut

			endif

		elseif (p_opt == 9) then

			if(apsi>psi_pres) then

				answer = 0.d0

			elseif(apsi<psi_out-1.d-12*psic) then

				answer = pedge

			else

				answer = pcenter * (1.d0 - apsi/psi_pres)**aux_fact * (apsi/psi_pres)**gammahut

			endif

		elseif (p_opt == 10) then
		! RFP profile

			answer = pedge + (pcenter-pedge) * (apsi/psic_13)**2 *  &
				 (6.d0 - 8.d0*apsi/psic_13 + 3.d0*(apsi/psic_13)**2)


		else

			if(apsi/psic_13 > 1.0d0) then
				answer = pcenter
	!			print*, 'apsi/psic =   ',apsi/psic
	!		elseif(apsi/psic > 0.0d0) then
!			elseif(dabs(psi) >= dabs(psic*fraction) ) then
			elseif(psi >= (psic*fraction) ) then
				answer = pedge + (pcenter-pedge)* (   (apsi/psic_13 - &
												fraction)/ (1.d0-fraction)   )**alpha
			else
				answer = pedge
			end if

		endif



	endif

!	answer = 4.d0*(pcenter-pedge)*(psi/psic)*(1.d0-psi/psic) + pedge

	if(answer<=0.d0) answer = pedge*1.d-10


	continue

end function pofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dpdpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: p_l,p_r,p_c,apsi
	
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag) then
		answer = pp_loc
		return
	elseif(psi==psic_flag) then
		answer = pcp_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_p_iso) then

			apsi=-dabs(psi)

			if(apsi/psic<p_iso_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psic, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(1,:) )
			endif

			answer=answer/psic
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif


	if(tri_type==13) then
		apsi = apsi - psi_e
	endif

	if(tri_type==9) then

		answer = 2.d0/(rmajor**3*q0)

    elseif (numerical_p_iso) then

			if ((apsi/psic)<p_iso_data(1,1)) then

				answer = dbsder(1,p_iso_data(1,1), p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(1,:) )

			elseif (apsi>=psic) then

				answer = dbsder(1,1.d0-1.d-6, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(1,:) )

			else

				answer = dbsder(1,apsi/psic, p_iso_ord, p_iso_data(1:ibreak_p,3),  &
								ibreak_p-p_iso_ord, p_iso_cscoef(1,:) )

			endif

		answer=answer/psic

	elseif(Broot==3) then

!		if((apsi/psic > fraction).AND.(apsi/psic<=1.d0)) then
		if((psi/psic > fraction).AND.(apsi/psic<=1.d0)) then
!		if(apsi/psic > fraction) then
		   answer = (pcenter-pedge)/psic/(1.d0-fraction) * alpha *   &
											(   apsi/ ((1.d0-fraction)*psic) - &
												fraction/ (1.d0-fraction)   )**(alpha-1.d0)
		else
			answer = 0.0d0
		end if

	else

		if (p_opt == 4) then
		! LDX profile
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = pcenter *  &
						((2.d0+gammahut)*vol1-gammahut*volofpsi(psi)) /  &
						(volofpsi(psi)+vol1)**(gammahut+2.d0) *  &
						dvdpsi(psi)

		elseif (p_opt == 5) then
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = pcenter *  &
						(  &
						((2.d0+gammahut)*vol1-gammahut*volofpsi(psi)) /  &
						(volofpsi(psi)+vol1)**(gammahut+2.d0) *  &
						dvdpsi(psi) * pauxofpsi(psi) +  &
						(volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0) *  &
								dpauxdpsi(psi)  &
													)

		elseif (p_opt == 6) then
		! p = K(1-psi)/(1+psi)**(gammahut+1)

!				answer = pcenter * (1.d0 + apsi/psic)**(gammahut-2.d0) *  &
!						(  &
!							( 1.d0 - (psi/psic)**2 ) * dpauxdpsi(psi)  &
!							- (2.d0 + gammahut*(apsi/psic-1.d0)) *  &
!							pauxofpsi(psi) / psic  &
!							)

				answer = pcenter / (1.d0 + apsi/psic)**gammahut *  &
						(  &
							(1.d0-apsi/psic) * dpauxdpsi(psi) -  &
							( 1.d0 + gammahut*(1.d0-apsi/psic)/(1.d0+apsi/psic) )  &
							*pauxofpsi(psi) / psic  &
							)

		elseif (p_opt==7) then

			if(apsi>psi_pres) then

				answer = 0.d0

			else

				answer = pcenter/psi_pres * (apsi/psi_pres)**(gammahut-1.d0) *  &
						(gammahut-(gammahut+1.d0)*(apsi/psi_pres))

			endif

		elseif (p_opt==8) then


			if((apsi>=psi_pres).or.(apsi<1.d-6*psic)) then

				answer = 0.d0

			else

				answer = pcenter/psi_pres * &
							(apsi/psi_pres)**(gammahut-1.d0) *  &
						(1.d0-apsi/psi_pres)**(aux_fact-1.d0)*  &
						(gammahut-(gammahut+aux_fact)*(apsi/psi_pres))

			endif

		elseif (p_opt==9) then


			if((apsi>=psi_pres).or.(apsi<psi_out-1.d-12*psic)) then

				answer = 0.d0

			else

				answer = pcenter/psi_pres * &
							(apsi/psi_pres)**(gammahut-1.d0) *  &
						(1.d0-apsi/psi_pres)**(aux_fact-1.d0)*  &
						(gammahut-(gammahut+aux_fact)*(apsi/psi_pres))

			endif

		elseif (p_opt == 10) then
		! RFP profile

			answer = 12.d0 * (pcenter-pedge) * (1.d0-apsi/psic_13)**2 * apsi/psic_13 /psic_13


		else

			if(dabs(psi) > dabs(psic_13).or.(dabs(psi/psic)==0.d0) ) then
				answer = 0.d0
	!		if(apsi/psic > 0.0d0) then
!			elseif(dabs(psi) >= dabs(psic*fraction) ) then
			elseif(psi >= psic*fraction ) then
			   answer = alpha*(pcenter-pedge)/psic_13*  &
						( apsi/psic_13-fraction )**(alpha-1.0d0)/ (1.d0-fraction)**alpha
	!		   answer = alpha*(pcenter-pedge)/psic*(apsi/psic)**(alpha-1.0d0)
			else
				answer = 0.0d0
			end if

		endif

	endif

!	answer = (pcenter-pedge)*4.0d0*(1.d0-2.d0*apsi/psic) / psic


	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dpdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function bzero(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi,x,apsim
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep) then
		answer = b0_loc
		return
	elseif(psi==psic_flag) then
		answer = b0c_loc
		return
	endif

	apsi = dabs(psi)
	apsim=min(apsi,psic)

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_F) then

			apsi=-dabs(psi)

			if(apsi/psic<F_data(1,1)) then
				answer = F_data(1,2)
			else
				answer = dbsval(apsi/psic, F_ord, F_data(1:ibreak_B0,3),  &
					ibreak_B0-F_ord, F_cscoef(1,1:ibreak_B0-F_ord) )
			endif

			return

		else
			apsi = 0.d0
		endif

    endif
    
	if(tri_type==13) then
		apsi = apsi - psi_e
	endif


	if(tri_type==9) then

		answer = 1.d0/rmajor

	elseif(numerical_F) then

		if ((apsi/psic)<F_data(1,1)) then
			answer = F_data(1,2)
		elseif((apsi/psic)>=1.d0) then 
			answer = F_data(ibreak_B0-F_ord,2)
		else

			answer = dbsval(apsi/psic, F_ord, F_data(1:ibreak_B0,3),  &
					ibreak_B0-F_ord, F_cscoef(1,1:ibreak_B0-F_ord) )



		endif

	else

		if (F_opt==0) then
		! no toroidal field, for levitated dipole

			answer = 0.d0

		elseif (F_opt==1) then

			if (apsi/psic<fraction) then

				answer = fvaccuum / rmajor

			elseif (apsi>psic) then

				answer = fcenter / rmajor

			else

				answer = (fvaccuum + (fcenter - fvaccuum)*  &
							( apsi/psic-fraction )**kappa/ (1.d0-fraction)**kappa)/rmajor 
!				answer = (fvaccuum + (fcenter - fvaccuum)*(apsi/psic)**kappa)/rmajor	

			endif

		elseif(F_opt==2) then

			answer = dsqrt(fvaccuum**2-eta_P*2.d0*mu_mag*rmajor**2*pofpsi(psi))/rmajor

		elseif(F_opt==3) then

			x = pmult*(pfact-1.d0)* (p_eps+(psic-apsim)/psic)**pexp *  &
						(apsim/psic)**(pexp+1.d0) +pmult !replaces pmult

			answer = ((1.d0-ffrac)*fvaccuum + (dsqrt(ffrac*fvaccuum**2 - x*2.d0*mu_mag*rmajor**2*pofpsi(psi))  &
					+ amult* dsqrt(2.d0*(mach_alf_0 **2-1.d0)*b_phi_zero*apsim) ))/rmajor

!			answer = ((1.d0-ffrac)*fvaccuum + (dsqrt(ffrac*fvaccuum**2 - x*2.d0*mu_mag*rmajor**2*pofpsi(psi))  &
!					+ amult* dsqrt(2.d0*(mach_alf_0 **2-1.d0)*b_phi_zero*apsim/psic) ))/rmajor

		elseif(F_opt==4) then

!			answer = b_phi_zero - pmult*dsqrt( mu_mag*pofpsi(psi) )
			answer = b_phi_zero - 0.5d0*pmult*(   1.d0+tanh(  pfact*( -.75d0+  &
								dsqrt(pofpsi(psi)/pofpsi(psic)) )  )   )

			answer = b_phi_zero-0.5d0*pmult *beta_center*( 1.d0 + tanh(pfact*  &
								(dsqrt(pofpsi(psi)/pofpsi(psic))-0.5d0) ) )

		elseif(F_opt==5) then
		! RFP profile

			if(abs(psi/psic)>1.d0) apsi = psic
			! to avoid problems with the powers

			answer = b_phi_zero * ( 1.d0 + mu_RFP * ( apsi/psic_13 - 1.d0 +  &
						(1.d0 - apsi/psic_13)**(kappa+1.d0)/(kappa+1.d0) ) )


		else

			print*, 'wrong option for B_zero:   ',F_opt
			pause
			stop

		endif

	    if(psi < psic*fraction ) answer = b_phi_zero
!		if(dabs(psi) <= dabs(psic*fraction) ) answer = 0.d0


	endif

	continue

end function bzero

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dbzerodpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi,x,y,apsim
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep) then
		answer = b0p_loc
		return
	elseif(psi==psic_flag) then
		answer = b0cp_loc
		return
	endif

	apsi = dabs(psi)
	apsim=min(apsi,psic)

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_F) then

			apsi=-dabs(psi)

			if(apsi/psic<F_data(1,1)) then
				answer = 0.d0
			else
				answer = dbsder(1,apsi/psic, F_ord, F_data(1:ibreak_B0,3),  &
								ibreak_B0-F_ord, F_cscoef(1,:) )
			endif

			answer=answer/psic
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif
    
	if(tri_type==13) then
		apsi = apsi - psi_e
	endif


	if(tri_type==9) then

		answer = 0.d0

	elseif(numerical_F) then


!			if (apsi/psic<=1.d-6) then
			if (apsi/psic<F_data(1,1)) then

				answer = 0.d0

!				answer = dbsder(1,1.d-6, F_ord, F_data(1:ibreak_B0,3),  &
!							ibreak_B0-F_ord, F_cscoef(1,:) )

			elseif (apsi>=psic) then

				answer = dbsder(1,1.d0-1.d-6, F_ord, F_data(1:ibreak_B0,3),  &
							ibreak_B0-F_ord, F_cscoef(1,:) )

			else

				answer = dbsder(1,apsi/psic, F_ord, F_data(1:ibreak_B0,3),  &
								ibreak_B0-F_ord, F_cscoef(1,:) )

			endif

		answer=answer/psic

	else


		if (F_opt==0) then
		! no toroidal field, for levitated dipole

			answer = 0.d0

		elseif (F_opt==1) then

		   answer = kappa*(fcenter - fvaccuum)/psic* &
				 ( apsi/psic-fraction )**(kappa-1.0d0)/rmajor/ (1.d0-fraction)**kappa
!		   answer = kappa*(fcenter - fvaccuum)/psic* &
!				 (apsi/psic)**(kappa-1.0d0)/rmajor

		elseif(F_opt==2) then

			answer = -mu_mag*rmajor**2*eta_P*dpdpsi(psi)/  &
				dsqrt(fvaccuum**2-eta_P*2.d0*mu_mag*rmajor**2*pofpsi(psi))/rmajor

			if(apsi>0.d0) answer = answer*psi/apsi

		elseif(F_opt==3) then

			x = pmult*(pfact-1.d0)* (p_eps+(psic-apsim)/psic)**pexp *  &
						(apsim/psic)**(pexp+1.d0) +pmult !replaces pmult
			y = pmult*(pfact-1.d0)/psic * ( &
					(1.d0+pexp)*(apsim/psic)**pexp*(p_eps+(psic-apsim)/psic)**pexp -  &
					pexp*(apsim/psic)**(pexp+1.d0)*(p_eps+(psic-apsim)/psic)**(pexp-1.d0)  &
												)

				answer =( amult * dsqrt(b_phi_zero*(mach_alf_0**2-1.d0) )  &
						/(2.d0*apsi**0.5)  &
						- (  mu_mag * rmajor * ( x*dpdpsi(psi) + y*pofpsi(psi) )   )  &
						/ dsqrt( ffrac*fvaccuum**2-2.d0*mu_mag*pmult*rmajor*pofpsi(psi) )  )  &
							/rmajor

!				answer =( amult * dsqrt(b_phi_zero*(mach_alf_0**2-1.d0) * apsi/psic)  &
!						/(apsi*2.d0**0.5)  &
!						- (  mu_mag * rmajor * ( x*dpdpsi(psi) + y*pofpsi(psi) )   )  &
!						/ dsqrt( ffrac*fvaccuum**2-2.d0*mu_mag*pmult*rmajor*pofpsi(psi) )  )  &
!							/rmajor

		elseif(F_opt==4) then

!			answer = -pmult*mu_mag*dpdpsi(psi)/(2.d0*dsqrt(mu_mag*pofpsi(psi)))

			answer = -pfact*beta_center**pmult*dpdpsi(psi)/( 4.d0*dsqrt(pofpsi(psi)*pofpsi(psic)) *  &
											cosh(pfact*(dsqrt(pofpsi(psi)/pofpsi(psic))-0.5d0)) )

			if(apsi>0.d0) answer = answer*psi/apsi

		elseif(F_opt==5) then
		! RFP profile

			if(abs(psi/psic)>1.d0) apsi = psic
			! to avoid problems with the powers

			answer = - b_phi_zero*mu_RFP * ( - 1.d0 + (1.d0 - apsi/psic_13)**kappa ) / psic_13


		endif

		if(psi < psic*fraction ) answer = 0.d0
!		if(dabs(psi) <= dabs(psic*fraction) ) answer = 0.d0


	endif

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dbzerodpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function bpolzero(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi

!!$	if(bpol_av==0.d0) then
!!$
!!$		answer = dabs( 2.d0*(psi_in-psi_out)/rmajor/x_size )
!!$
!!$	else
!!$
!!$		answer = bpol_av
!!$
!!$	endif

!!$	answer = bpol0_temp * (1.d0 + ana_fac*(psi/psic-psi_B_min)**alpha_th)

!!$	if(ana_fac>0.d0) then
!!$		answer = ana_fac*dabs( 2.d0*(psi_in-psi_out)/rmajor/x_size )
!!$	else
!!$		answer = dabs( 2.d0*(psi_in-psi_out)/rmajor/x_size ) * 1.d-6
!!$	endif

	if((eq3_opt==9).or.(p_opt==9)) then

		answer = B_pmax

	else

		if(ana_fac>0.d0) then
			answer = ana_fac
		else
			answer = 1.d-6
		endif

	endif

	if(answer<0.d0) then
		continue
	endif

	continue

end function bpolzero

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dbpolzerodpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    real (kind=dkind) :: answer
	real(kind=dkind) :: apsi

	apsi = dabs(psi)

	answer = 0.d0

!!$	answer = alpha_th * bpol0_temp * ana_fac / psic *  &
!!$				(psi/psic-psi_B_min)**(alpha_th-1.d0)

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dbpolzerodpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_theta(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi !,M0
   
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_dep2) then
		answer = mth_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psic<mtheta_data(1,1)) then
				answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
			else
				answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1,1:ibreak_th-mtheta_ord) )
			endif

			return

		else
			apsi = 0.d0
		endif

    endif

	if((Broot==4).or.(Broot==5)) then
	! this case treted first to avoid messy "if" groups

		answer = mach_thetahat(psi)*bzero(psi)/b_phi_zero
		return

	endif

	if (numerical_mtheta) then


		if ((apsi/psic)<=mtheta_data(1,1)) then
			answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
		elseif((apsi/psic)>=1.d0) then 
			answer =  mach_theta_max/mach_theta_num *  &
					mtheta_data(ibreak_th-mtheta_ord,2)
		else
			answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1,1:ibreak_th-mtheta_ord) )
		endif

	else

	   if (Broot==3) then


			if(apsi==0.d0) then

				answer = M0

			else

	!			 M0 = mach_alf_0*dsqrt(2.d0/(beta_center*gamma))

	!				answer = M0
				 answer = M0*dsqrt( pofpsi(psic)*dofpsi(psic)/pofpsi(psi)/dofpsi(psi) )  &
									* bzero(psi)/bzero(psic)

	!			 answer = M0*dsqrt( pofpsi(psic)*dofpsi(psic)/pofpsi(psi)/dofpsi(psi) )  &
	!								* abs(bzero(psi)/bzero(psic))

			endif

		else

				! -------------------- t-2t shape --------------------

			if(apsi<=0.d0) then
				answer = mach_theta_edge
			elseif(apsi/psic>=2.d0*t_mth) then
				answer = 0.d0
			elseif(apsi/psic<=t_mth) then
				answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
						(2.d0/t_mth*apsi/psic - (apsi/psic/t_mth)**2)
			elseif(apsi/psic>t_mth) then
				answer = mach_theta_max *  &
						(2.d0*t_mth-apsi/psic)**2 * (2.d0*apsi/psic-t_mth) / t_mth**3
			endif

!!$				! -------------------- t shape --------------------
!!$
!!$			if(apsi<=0.d0) then
!!$				answer = mach_theta_edge
!!$			elseif(apsi>=psic) then
!!$				answer = mach_theta_edge
!!$			elseif(apsi/psic<=t_mth) then
!!$				answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
!!$						(2.d0/t_mth*apsi/psic - (apsi/psic/t_mth)**2)
!!$			elseif(apsi/psic>t_mth) then
!!$				answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
!!$						( 1.d0/(1.d0-t_mth)**2 * (1.d0-apsi/psic) *  &
!!$							(apsi/psic + 1.d0 - 2.d0*t_mth) )
!!$			endif

!!$							!------------------ normal shape!----------------------
!!$							if(apsi <= fraction*dabs(psic) .or. psi >= psic) then
!!$							   answer = mach_theta_edge !0.0d-2
!!$							else
!!$							   answer = mach_theta_max*4.0d0*( apsi/psic-fraction)  &
!!$										*( 1.0d0 - apsi/psic )
!!$				!!$			   answer = mach_theta_max*4.0d0*(apsi/psic)*(1.0d0 - apsi/psic)  
!!$							   answer = mach_theta_edge + (mach_theta_max-mach_theta_edge)*  &
!!$										4.0d0*(apsi/psic)*(1.0d0 - apsi/psic)  
!!$							end if

			!--------------------for stability----------------------
!!$			if(apsi <= 0.1d0*dabs(psic) .or. psi >= psic) then
!!$!			   answer = mach_theta_edge !0.0d-2
!!$				answer = 0.d0
!!$			else
!!$!			   answer = mach_theta_max*4.0d0*( apsi/psic-0.1d0)  &
!!$!						*( 1.0d0 - apsi/psic )
!!$!!!				answer = mach_theta_max*(-1.d3/108.d0)*(apsi/psic-0.1d0)**2*(apsi/psic-1.d0)
!!$				answer = mach_theta_max*(apsi/psic-0.1d0)**6*(1.d0-apsi/psic)
!!$!!				answer = mach_theta_max*(apsi/psic-0.1d0)**2*(1.d0-apsi/psic)**4
!!$
!!$			end if

		endif

	endif

	if (eq_type==3) answer = 0.d0

	continue

end function mach_theta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_thetadpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi
	real (kind=dkind) :: b0,d,p,b0c,dc,pc
	
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep2) then
		answer = mthp_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psic<mtheta_data(1,1)) then
				answer = 0.d0
			else
				answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )
			endif

			answer=answer/psic
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif


	if((Broot==4).or.(Broot==5)) then
	! this case treted first to avoid messy "if" groups

		answer = ( mach_thetahat(psi)*dbzerodpsi(psi)  + &
							dmach_thetahatdpsi(psi)*bzero(psi) )  &
					/b_phi_zero

		return

	endif

	if (numerical_mtheta) then

		if (apsi/psic<=mtheta_data(1,1)) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )


		elseif (apsi>=psic) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d0-1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )

		else

			answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )

		endif

		answer=answer/psic

	 else

		if(Broot==3) then

	!		 M0 = mach_alf_0*dsqrt(2.d0/(beta_center*gamma))
			 b0 = bzero(psi)
			 d = dofpsi(psi)
			 p = pofpsi(psi)
			 b0c = bzero(psic)
			 dc = dofpsi(psic)
			 pc = pofpsi(psic)

			 answer = -0.5d0*M0 * (dc*pc)**0.5/b0c/(d*p)**1.5*  &
				 (b0*p*dddpsi(psi) - 2.d0*d*p*dbzerodpsi(psi) + b0*d*dpdpsi(psi))

	!!		answer = -0.5d0*M0 * (dc*pc)**0.5/b0c**2/(d*p)**1.5*  &
	!!				(  2.d0*d*p*sign(1.d0,b0/b0c)*dbzerodpsi(apsi) -  &
	!!					abs(b0c*b0)* ( p*dddpsi(apsi) + d*dpdpsi(apsi) )  )


	!		 answer = -0.5d0*M0 * (dc*pc)**0.5/(d*p)**1.5*  &
	!			 (b0*p*dddpsi(psi) + b0*d*dpdpsi(psi))

	!		answer = 0.d0 


			if(apsi>0.d0) answer = answer*psi/apsi

		else

				! -------------------- t-2t shape --------------------

			if(apsi<=0.d0) then
				answer = 0.d0
			elseif(apsi/psic>=2.d0*t_mth) then
				answer = 0.d0
			elseif(apsi/psic<=t_mth) then
				answer = (mach_theta_max-mach_theta_edge) *  &
						2.d0/t_mth/psic*(1.d0 - (apsi/psic/t_mth))
			elseif(apsi/psic>t_mth) then
				answer = 6.d0*mach_theta_max *  &
						(apsi/psic-2.d0*t_mth)*(apsi/psic-t_mth)/t_mth**3/psic
			endif

!!$				! -------------------- t shape --------------------
!!$
!!$			if(apsi<=0.d0) then
!!$				answer = 0.d0
!!$			elseif(apsi>=psic) then
!!$				answer = 0.d0
!!$			elseif(apsi/psic<=t_mth) then
!!$				answer = (mach_theta_max-mach_theta_edge) *  &
!!$						2.d0/t_mth/psic*(1.d0 - (apsi/psic/t_mth))
!!$			elseif(apsi/psic>t_mth) then
!!$				answer = -2.d0*(mach_theta_max-mach_theta_edge) /  &
!!$						psic/(1.d0-t_mth)**2 * (apsi/psic-t_mth)
!!$			endif


							!------------------ normal shape!----------------------
!!$							if(apsi <= fraction*dabs(psic) .or. apsi > psic) then
!!$							   answer = 0.0d0
!!$							else
!!$							   answer = mach_theta_max*4.0d0*  &
!!$							   (fraction - 2.d0*apsi/psic + 1.d0) / psic
!!$				!			   answer = mach_theta_max*4.0d0*(1.0d0 - 2.0d0*apsi/psic)/psic
!!$					!!$		   answer = (mach_theta_max-mach_theta_edge)  &
!!$					!!$					*4.0d0*(1.0d0 - 2.0d0*apsi/psic)/psic
!!$							end if

			!--------------------for stability----------------------
!!$			if(apsi <= 0.1d0*dabs(psic) .or. apsi > psic) then
!!$			   answer = 0.0d0
!!$			else
!!$!			   answer = mach_theta_max*4.0d0*  &
!!$!			   (0.1d0 - 2.d0*apsi/psic + 1.d0) / psic
!!$
!!$!!				answer =  mach_theta_max*(-1.d3/36.d0)*(apsi/psic-0.1d0)*(apsi/psic-.7d0)/psic
!!$				answer =  mach_theta_max*(6.1d0-7.d0*apsi/psic)*(apsi/psic-0.1d0)**5/psic
!!$!!!				answer =  mach_theta_max*(  &
!!$!!!					2.d0*(1.d0-apsi/psic)**4*(apsi/psic-0.1d0)  &
!!$!!!					- 4.d0*(1.d0-apsi/psic)**3*(apsi/psic-0.1d0)**2 )/psic
!!$			end if



		endif

	endif

	if (eq_type==3) answer = 0.d0

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dmach_thetadpsi



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_thetahat(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi, xx

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep2) then
		answer = mth_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psic<mtheta_data(1,1)) then
				answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
			else
				answer = mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1,1:ibreak_th-mtheta_ord) )
			endif

			return

		else
			apsi = 0.d0
		endif

    endif

	if (numerical_mtheta) then

		if ((apsi/psic)<=mtheta_data(1,1)) then
			answer = mach_theta_max/mach_theta_num *  &
						mtheta_data(1,2)
		elseif((apsi/psic)>=1.d0) then
			answer =  mach_theta_max/mach_theta_num *  &
					mtheta_data(ibreak_th-mtheta_ord,2)
		else
			answer =  mach_theta_max/mach_theta_num *  &
					dbsval(apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
					ibreak_th-mtheta_ord, mtheta_cscoef(1,1:ibreak_th-mtheta_ord) )
		endif

	else

			! -------------------- t shape --------------------

		xx = apsi/psic

		if(xx<=t_mth) then
			answer = mach_theta_edge
		elseif(xx >= t_mth + 2.d0*w_mth) then
			answer = 0.d0
		elseif(xx > t_mth + w_mth) then !f2(x)
			answer = (mach_theta_max) *  &
					 (2.d0*xx-2.d0*t_mth-w_mth)*(t_mth+2.d0*w_mth-xx)**2/w_mth**3
		elseif(xx>t_mth) then !f1(x)
			answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
							(2.d0*t_mth+3.d0*w_mth-2.d0*xx)*(t_mth-xx)**2/w_mth**3
		endif

!!$		if(xx<=t_mth) then
!!$			answer = mach_theta_edge
!!$		elseif(xx>=1.d0-t_mth) then
!!$			answer = 0.d0
!!$		elseif(xx>2.d0*t_mth) then !f2(x)
!!$			answer = (mach_theta_max) *  &
!!$					 (1.d0-7.d0*t_mth+2.d0*xx)*(xx+t_mth-1.d0)**2/(1.d0-3.d0*t_mth)**3
!!$		elseif(apsi/psic>t_mth) then !f1(x)
!!$			answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
!!$							(5.d0*t_mth-2.d0*xx)*(t_mth-xx)**2/t_mth**3
!!$		endif

!!$		if(xx<=t_mth) then
!!$			answer = mach_theta_edge
!!$		elseif(xx>=1.d0-t_mth) then
!!$			answer = mach_theta_edge
!!$		elseif(xx>2.d0*t_mth) then !f2(x)
!!$			answer = mach_theta_edge +(mach_theta_max-mach_theta_edge) *  &
!!$					 (1.d0-7.d0*t_mth+2.d0*xx)*(xx+t_mth-1.d0)**2/(1.d0-3.d0*t_mth)**3
!!$		elseif(apsi/psic>t_mth) then !f1(x)
!!$			answer = mach_theta_edge + (mach_theta_max-mach_theta_edge) *  &
!!$							(5.d0*t_mth-2.d0*xx)*(t_mth-xx)**2/t_mth**3
!!$		endif

	endif

	continue

end function mach_thetahat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_thetahatdpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: apsi, xx
	real (kind=dkind) :: b0,d,p,b0c,dc,pc

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep2) then
		answer = mthp_loc
		return
	endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (zone_loc==-1) then

		if(numerical_mtheta) then

			apsi=-dabs(psi)

			if(apsi/psic<mtheta_data(1,1)) then
				answer = 0.d0
			else
				answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )
			endif

			answer=answer/psic
			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif
    
	if (numerical_mtheta) then

		if (apsi/psic<=mtheta_data(1,1)) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )


		elseif (apsi>=psic) then

			answer = mach_theta_max/mach_theta_num *  &
						dbsder(1,1.d0-1.d-6, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
						ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )

		else

			answer = mach_theta_max/mach_theta_num *  &
							dbsder(1,apsi/psic, mtheta_ord, mtheta_data(1:ibreak_th,3),  &
							ibreak_th-mtheta_ord, mtheta_cscoef(1,:) )

		endif

		answer=answer/psic

	 else


			! -------------------- t shape --------------------
		xx = apsi/psic

		if(xx<=t_mth) then
			answer = 0.d0
		elseif(xx >= t_mth + 2.d0*w_mth) then
			answer = 0.d0
		elseif(xx > t_mth + w_mth) then !f2(x)
			answer = - (mach_theta_max) *  &
					 6.d0*(2.d0*w_mth-xx+t_mth)*(xx-w_mth-t_mth)/w_mth**3
		elseif(xx>t_mth) then !f1(x)
			answer = (mach_theta_max-mach_theta_edge) *  &
							6.d0*(xx-t_mth)*(t_mth+w_mth-xx)/w_mth**3
		endif

!!$		if(xx<=t_mth) then
!!$			answer = 0.d0
!!$		elseif(xx>=1.d0-t_mth) then
!!$			answer = 0.d0
!!$		elseif(xx>2.d0*t_mth) then !f2(x)
!!$			answer = - (mach_theta_max) *  &
!!$					 6.d0*(2.d0*t_mth-xx)*(xx+t_mth-1.d0)/(1.d0-3.d0*t_mth)**3
!!$		elseif(apsi/psic>t_mth) then !f1(x)
!!$			answer = - (mach_theta_max-mach_theta_edge) *  &
!!$							6.d0*(2.d0*t_mth**2-3.d0*t_mth*xx+xx**2)/t_mth**3
!!$		endif

!!$		if(xx<=t_mth) then
!!$			answer = 0.d0
!!$		elseif(xx>=1.d0-t_mth) then
!!$			answer = 0.d0
!!$		elseif(xx>2.d0*t_mth) then !f2(x)
!!$			answer = - (mach_theta_max-mach_theta_edge) *  &
!!$					 6.d0*(2.d0*t_mth-xx)*(xx+t_mth-1.d0)/(1.d0-3.d0*t_mth)**3
!!$		elseif(apsi/psic>t_mth) then !f1(x)
!!$			answer = - (mach_theta_max-mach_theta_edge) *  &
!!$							6.d0*(2.d0*t_mth**2-3.d0*t_mth*xx+xx**2)/t_mth**3
!!$		endif

		answer = answer/psic

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

  continue

end function dmach_thetahatdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function mach_phi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: cs,omeg,apsi, mth

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep2) then
		answer = mph_loc
		return
	endif

     if (zone_loc==-1) then

		if((numerical_omega).or.(numerical_mphi)) then
			apsi=-dabs(psi)
			
			if(numerical_omega) then

				if(apsi/psic<mphi_data(1,1)) then
					omeg = mphi_data(1,2)
				else
					omeg = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )
				endif

				mth = mach_theta(psi,zone_loc)

				cs = dsqrt(gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc))

				answer = rmajor*omeg/cs + mth

			elseif(numerical_mphi) then

				if(apsi/psic<mphi_data(1,1)) then
					answer = mphi_data(1,2)
				else

					answer = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )

				endif

			endif

			return
			
		else
			apsi = 0.d0
		endif
		
    endif

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if ((numerical_omega).or.(Kepler_Omega)) then

		if(numerical_omega) then

			if((apsi/psic) >= 1.d0) then
				omeg = mphi_data(ibreak_fi-mphi_ord,2)
			elseif((apsi/psic) < mphi_data(1,1)) then
				omeg = mphi_data(1,2)
			else
				omeg = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
						ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )

			endif

		elseif(Kepler_Omega) then

			omeg = sqrt(G_gravity*M_gravity/rmajor**3)

		endif

		mth = mach_theta(psi)

		if(eq_type==1) then
			cs = dsqrt(gamma*pofpsi(psi)/dofpsi(psi))
		elseif(eq_type==3) then
			cs = dsqrt(tparofpsi(psi))
		endif

		answer = rmajor*omeg/cs + mth

    elseif (numerical_mphi) then

		if((apsi/psic) > 1.d0) then
			answer = mphi_data(ibreak_fi-mphi_ord,2)
		elseif((apsi/psic) < mphi_data(1,1)) then
			answer = mphi_data(1,2)
		else

			answer = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
					ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )

		endif

	else


		if(psi/psic <= fraction) then
!		if(psi/psic <= fraction .or. apsi/psic > 1.d0) then
!!		if(apsi/psic <= fraction .or. apsi/psic > 1.d0) then

		   answer = 0.d0

		else

			if(Broot==3) then

!					answer = mach_phi_max*(1.d0/(1.d0-fraction)*apsi/psic - fraction/(1.d0-fraction)) + 0.001d0
!					answer = 0.d0
				answer = mach_phi_max*(apsi/psic+mphi_min)**alpha_mphi

			else

					answer = mach_phi_max*  &
							( apsi/psic-fraction +mphi_min)**alpha_mphi/ (1.d0-fraction)**alpha_mphi
				!				answer = mach_phi_max*(apsi/psic+mphi_min)**alpha_mphi

				! --------------------------for stability-------------------------- !
!!$				if(apsi>0.1d0*psic) then 
!!$					answer = mach_phi_max*( (apsi/psic-0.1d0)/0.9d0)**alpha_mphi
!!$				else
!!$					answer = 0.d0
!!$				endif

!				answer = 0.5d0*  mach_phi_max * (1.d0+tanh( (apsi/psic-1.d-1)/1.d-2 ) )/  &
!								dsqrt(gamma*pofpsi(psi)/dofpsi(psi))  &
!								*dsqrt(gamma*pofpsi(psic)/dofpsi(psic)) 

!				answer = mach_phi_max*sin(1.d1*pi*apsi/psic)

!				answer = 0.5d0*  mach_phi_max * (1.d0+tanh( (apsi/psic-1.d-1)/1.d-2 ) )

				! --------------------------for stability-------------------------- !

			endif

		end if

!		answer = mach_phi_max*(1.d0-apsi/psic)

	endif

	continue

end function mach_phi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dmach_phidpsi(psi, zone) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: omeg,d_omeg
	real (kind=dkind) :: cs,dcs,p,d,apsi, dmth

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_dep2) then
		answer = mphp_loc
		return
	endif

    if (zone_loc==-1) then

		if((numerical_omega).or.(numerical_mphi)) then

			apsi=-dabs(psi)

			if (numerical_omega) then

				if(apsi/psic<mphi_data(1,1)) then
					answer = 0.d0
				else
					omeg = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
					ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )

					d_omeg = dbsder(1,apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
						ibreak_fi-mphi_ord, mphi_cscoef(1,:) )
				endif

				d_omeg=d_omeg/psic

				p = pofpsi(psi,zone_loc)
				d = dofpsi(psi,zone_loc)
				cs = dsqrt(gamma*p/d)
				dcs = gamma*(-p*dddpsi(psi,zone_loc)+d*dpdpsi(psi,zone_loc))/  &
					(2.d0*d**2*cs)


				if (apsi>0.d0) dcs = dcs*psi/apsi

				dmth = dmach_thetadpsi(psi,zone_loc)

				answer = rmajor*(d_omeg/cs-omeg/cs**2*dcs) + dmth

			elseif (numerical_mphi) then

				if(apsi/psic<mphi_data(1,1)) then

					answer = 0.d0

				else

					answer = dbsder(1,apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
									ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

				endif

				answer=answer/psic

			endif

			if(apsi/=0.d0) answer = answer*psi/apsi
			return

		else
			apsi = 0.d0
		endif

    endif

    
	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if ((numerical_omega).or.(Kepler_Omega)) then


		if (numerical_omega) then

			if((apsi/psic) >= 1.d0) then

				omeg = mphi_data(ibreak_fi-mphi_ord,2)

				d_omeg = dbsder(1,1.d0, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

			elseif((apsi/psic) < mphi_data(1,1)) then

				omeg = omega_data(1,2)

				d_omeg = 0.d0

			else

				omeg = dbsval(apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(1,1:ibreak_fi-mphi_ord) )

				d_omeg = dbsder(1,apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
				ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

			endif

		elseif(Kepler_Omega) then

			omeg = sqrt(G_gravity*M_gravity/rmajor**3)
			d_omeg = 0.d0

		endif

		d_omeg=d_omeg/psic

		if(eq_type==1) then
			p = pofpsi(psi)
			d = dofpsi(psi)
			cs = dsqrt(gamma*p/d)
			dcs = gamma*(-p*dddpsi(psi)+d*dpdpsi(psi))/  &
				(2.d0*d**2*cs)
		elseif(eq_type==3) then
			cs=dsqrt(tparofpsi(psi))
			dcs=dtpardpsi(psi)/(2.d0*cs)
		endif
		if (apsi>0.d0) dcs = dcs*psi/apsi

		dmth = dmach_thetadpsi(psi)

		answer = rmajor*(d_omeg/cs-omeg/cs**2*dcs) + dmth
			
	elseif (numerical_mphi) then

!!$		if(apsi/psic>0.99d0) then
!!$
!!$			call splinderiv(mphi_data(1:data_dim,1),mphi_data(1:data_dim,2),  &
!!$							mphi_data(1:data_dim,3),0.99d0,answer,data_dim)
!!$
!!$		else	
!!$
!!$			call splinderiv(mtheta_data(1:data_dim,1),mtheta_data(1:data_dim,2),  &
!!$							mtheta_data(1:data_dim,3),apsi/psic,answer,data_dim)
!!$
!!$		endif


!		answer = dcsder(1,apsi/psic, ibreak_fi-1, mphi_data(1:ibreak_fi,3),  &
!							mphi_cscoef(:,1:ibreak_fi) )

		if(apsi/psic>1.d0) then

			answer = dbsder(1,1.d0, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

		elseif(apsi/psic<mphi_data(1,1)) then

			answer = dbsder(1,0.d0, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

		else

			answer = dbsder(1,apsi/psic, mphi_ord, mphi_data(1:ibreak_fi,3),  &
							ibreak_fi-mphi_ord, mphi_cscoef(1,:) )

		endif

		answer=answer/psic

	 else


		if(psi <= fraction*psic .or. apsi > psic) then
!		if(apsi <= fraction*dabs(psic) .or. apsi > psic) then

			answer = 0.0d0

		else

			if(Broot==3) then

!					answer = mach_phi_max/(1.d0-fraction)/psic
!!!					answer = 0.d0
				answer = alpha_mphi * mach_phi_max * (apsi/psic+mphi_min)**(alpha_mphi-1.d0) /psic

			else

									answer = alpha_mphi * mach_phi_max *  &
										( apsi/psic-fraction + mphi_min)**(alpha_mphi-1.d0) / psic/ (1.d0-fraction)**alpha_mphi
					!				answer = alpha_mphi * mach_phi_max * (apsi/psic+mphi_min)**(alpha_mphi-1.d0) /psic

				! --------------------------for stability-------------------------- !
!!$				if(apsi>0.1d0*psic) then 
!!$					answer = mach_phi_max*alpha_mphi*( (apsi/psic-0.1d0)/0.9d0)**(alpha_mphi-1.d0)/0.9d0
!!$				else
!!$					answer = 0.d0
!!$				endif
!				answer = 5.d1*mach_phi_max / (cosh( (apsi/psic-1.d-1)/1.d-2 ))**2/psic

!				answer = (mach_phi_max*(1.d0/psic*(50.d0*gamma*dofpsi(psi)*pofpsi(psi)/(cosh(1.d1-1.d2*psi/psic))**2)  &
!					-0.25d0*(1.d0+tanh(1.d2*(-1.d-1+psi/psic)))*  &
!					gamma*(-pofpsi(psi)*dddpsi(psi)+dofpsi(psi)*dpdpsi(psi)))  /  &
!					(dofpsi(psi)**2*(gamma/dofpsi(psi)*pofpsi(psi))**1.5d0)        )&
!								*dsqrt(gamma*pofpsi(psic)/dofpsi(psic))  

!					answer = mach_phi_max*cos(1.d1*pi*apsi/psic)*1.d1/psic

!				answer = mach_phi_max*50.d0/(cosh(1.d1-1.d2*psi/psic))**2/psic

				! --------------------------for stability-------------------------- !

			endif

		end if

!		answer = -mach_phi_max/psic

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dmach_phidpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function pparofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer,apsi

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (eq_type == 1) then

	    if(dabs(psi) > 0.0d0) then
		   answer = qparedge + qparcenter*(dabs(psi)/psic)**alpha 
		else
			answer = qparedge
		end if

	elseif (eq_type == 3) then

	    if (numerical_p_par) then

			if ((apsi/psic)<=1.d-6) then
				answer = p_par_data(1,2)
			elseif (apsi>psic) then
				answer = p_par_data(ibreak_ppar-ppar_ord,2)
			else

				answer = dbsval(apsi/psic, ppar_ord, p_par_data(1:ibreak_ppar,3),  &
						ibreak_ppar-ppar_ord, ppar_cscoef(1,1:ibreak_ppar-ppar_ord) )

			endif

!!$			if ((apsi/psic)<=1.d-3) then
!!$				answer = p_par_data(1,2)
!!$			elseif ((apsi/psic)>=1.d0) then
!!$				answer = p_par_data(data_dim,2)
!!$			else
!!$				call splint(p_par_data(1:data_dim,1),p_par_data(1:data_dim,2),  &
!!$						p_par_data(1:data_dim,3),data_dim,apsi/psic,answer)
!!$			endif

		else

			if (eq3_opt == 1) then

				answer = dofpsi(psi) * tparofpsi(psi)

			elseif (eq3_opt == 2) then

				if(apsi > psic) then
					answer = qparcenter
				elseif(dabs(psi) > 0.0d0) then
					answer = qparedge + (qparcenter-qparedge)*(apsi/psic)**alpha 
				else
					answer = qparedge
				end if

			elseif (eq3_opt == 3) then
			! LDX profile

!				answer = qparedge + (qparcenter-qparedge)*  &
!					(psi_in-apsi)*(apsi+abs(psi_in-psi_out)/2.d0)**alpha /psic**(alpha+1.d0)

				answer = qparedge + (qparcenter-qparedge)*  &
					(psi_in-apsi)*(apsi+abs(psi_in-psi_out)*psifrac)**alpha /psic**(alpha+1.d0)

				if(answer<=qparedge) answer = qparedge

			elseif (eq3_opt == 4) then
			! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = qparcenter * (volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0)

			elseif (eq3_opt == 5) then
			! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = qparcenter * (volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0) *  &
									pauxofpsi(psi)

			elseif (eq3_opt == 6) then
			! p = K(1-psi)/(1+psi)**(gammahut+1)

				answer = qparcenter * (1.d0 - apsi/psic) / (1.d0 + apsi/psic)**gammahut *  &
									pauxofpsi(psi)


			elseif (eq3_opt == 9) then

				if(apsi>psi_pres) then

					answer = 0.d0

				elseif(apsi<psi_out-1.d-12*psic) then

					answer = qparedge

				else

					answer = qparcenter * (1.d0 - apsi/psi_pres)**aux_fact * (apsi/psi_pres)**gammahut

				endif

			endif

		endif

	endif

	if(answer<=0.d0) answer = qparedge*1.d-6

	continue

end function pparofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function pperpofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer,apsi

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (eq_type == 1) then

		if(dabs(psi) > 0.0d0) then
			answer = qperpedge + qperpcenter*(dabs(psi)/psic)**alpha
		else
			answer = qperpedge
		end if

	elseif (eq_type == 3) then

	    if (numerical_p_perp) then

			if ((apsi/psic)<=1.d-6) then
				answer = p_perp_data(1,2)
			elseif (apsi>psic) then
				answer = p_perp_data(ibreak_pperp-pperp_ord,2)
			else

				answer = dbsval(apsi/psic, pperp_ord, p_perp_data(1:ibreak_pperp,3),  &
						ibreak_pperp-pperp_ord, pperp_cscoef(1,1:ibreak_pperp-pperp_ord) )

			endif

!!$			if ((apsi/psic)<=1.d-3) then
!!$				answer = p_perp_data(1,2)
!!$			elseif ((apsi/psic)>=1.d0) then
!!$				answer = p_perp_data(data_dim,2)
!!$			else
!!$				call splint(p_perp_data(1:data_dim,1),p_perp_data(1:data_dim,2),  &
!!$						p_perp_data(1:data_dim,3),data_dim,apsi/psic,answer)
!!$			endif

		else

			if(eq3_opt==1) then

				answer = dofpsi(psi) * tperpofpsi(psi)

			elseif (eq3_opt==2) then

				if(dabs(psi) > psic) then
					answer = qperpcenter
				elseif(dabs(psi) > 0.0d0) then
					answer = qperpedge + (qperpcenter-qperpedge)*(apsi/psic)**alpha	
				else
					answer = qperpedge
				end if

			elseif (eq3_opt >= 3) then
			! LDX profile, this is never used

				answer = 0.d0

			endif

		endif

	endif

	continue

end function pperpofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dppardpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (numerical_p_par) then

		if ((apsi/psic)<=1.d-3) then

			answer = dbsder(1,0.d0, ppar_ord, p_par_data(1:ibreak_ppar,3),  &
							ibreak_ppar-ppar_ord, ppar_cscoef(1,:) )

		elseif (apsi>=psic) then

			answer = dbsder(1,1.d0, ppar_ord, p_par_data(1:ibreak_ppar,3),  &
							ibreak_ppar-ppar_ord, ppar_cscoef(1,:) )

		else

			answer = dbsder(1,apsi/psic, ppar_ord, p_par_data(1:ibreak_ppar,3),  &
							ibreak_ppar-ppar_ord, ppar_cscoef(1,:) )

		endif

!!$		if(apsi/psic>0.99d0) then
!!$
!!$			call splinderiv(p_par_data(1:data_dim,1),p_par_data(1:data_dim,2),  &
!!$							p_par_data(1:data_dim,3),0.99d0,answer,data_dim)
!!$
!!$		else
!!$
!!$			call splinderiv(p_par_data(1:data_dim,1),p_par_data(1:data_dim,2),  &
!!$							p_par_data(1:data_dim,3),apsi/psic,answer,data_dim)
!!$
!!$		endif

		answer=answer/psic

	else
		if((dabs(psi) > 0.0d0).AND.(eq_type==3).AND.(eq3_opt==1)) then
			answer = dddpsi(psi)*tparofpsi(psi)+dofpsi(psi)*dtpardpsi(psi)
			answer = answer*apsi/psi
		elseif((dabs(psi) > 0.0d0).AND.(eq_type==3).AND.(eq3_opt==2)) then
			answer = alpha*(qparcenter-qparedge)/psic**alpha*apsi**(alpha-1)

		elseif (eq3_opt == 3) then
		! LDX profile

!			answer =  -(qparcenter-qparedge)*  &
!				((1.d0+alpha)*apsi-alpha*psi_in+abs(psi_in-psi_out)/2.d0)*  &
!				(apsi+abs(psi_in-psi_out)/2.d0)**(alpha-1.d0) /psic**(alpha+1.d0)

			answer =  -(qparcenter-qparedge)*  &
				((1.d0+alpha)*apsi-alpha*psi_in+abs(psi_in-psi_out)*psifrac)*  &
				(apsi+abs(psi_in-psi_out)*psifrac)**(alpha-1.d0) /psic**(alpha+1.d0)

		elseif (eq3_opt == 4) then
		! LDX profile
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = qparcenter *  &
						((2.d0+gammahut)*vol1-gammahut*volofpsi(psi)) /  &
						(volofpsi(psi)+vol1)**(gammahut+2.d0) *  &
						dvdpsi(psi)

		elseif (eq3_opt == 5) then
		! p = K(vol-vol1)/(vol+vol1)**(gamma+1)

				answer = qparcenter *  &
						(  &
						((2.d0+gammahut)*vol1-gammahut*volofpsi(psi)) /  &
						(volofpsi(psi)+vol1)**(gammahut+2.d0) *  &
						dvdpsi(psi) * pauxofpsi(psi) +  &
						(volofpsi(psi)-vol1)/(volofpsi(psi)+vol1)**(gammahut+1.d0) *  &
								dpauxdpsi(psi)  &
													)
		elseif (eq3_opt == 6) then
		! p = K(1-psi)/(1+psi)**(gammahut+1)


				answer = qparcenter / (1.d0 + apsi/psic)**gammahut *  &
						(  &
							(1.d0-apsi/psic) * dpauxdpsi(psi) -  &
							( 1.d0 + gammahut*(1.d0-apsi/psic)/(1.d0+apsi/psic) )  &
							*pauxofpsi(psi) / psic  &
							)

		elseif (eq3_opt==9) then


			if((apsi>=psi_pres).or.(apsi<psi_out-1.d-12*psic)) then

				answer = 0.d0

			else

				answer = qparcenter/psi_pres * &
							(apsi/psi_pres)**(gammahut-1.d0) *  &
						(1.d0-apsi/psi_pres)**(aux_fact-1.d0)*  &
						(gammahut-(gammahut+aux_fact)*(apsi/psi_pres))

			endif

		else
		   answer = 0.d0
		end if

		if(dabs(psi) > psic) answer = 0.d0
  
	endif

	if (apsi>0.d0) answer = answer*apsi/psi

  continue


end function dppardpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dpperpdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

    if (numerical_p_perp) then

		if ((apsi/psic)<=1.d-3) then

			answer = dbsder(1,0.d0, ppar_ord, p_perp_data(1:ibreak_pperp,3),  &
							ibreak_pperp-pperp_ord, pperp_cscoef(1,:) )

		elseif (apsi>=psic) then

			answer = dbsder(1,1.d0, pperp_ord, p_perp_data(1:ibreak_pperp,3),  &
							ibreak_pperp-pperp_ord, pperp_cscoef(1,:) )

		else

			answer = dbsder(1,apsi/psic, pperp_ord, p_perp_data(1:ibreak_pperp,3),  &
							ibreak_pperp-pperp_ord, pperp_cscoef(1,:) )

		endif

!!$		if(apsi/psic>0.99d0) then
!!$
!!$			call splinderiv(p_perp_data(1:data_dim,1),p_perp_data(1:data_dim,2),  &
!!$						p_perp_data(1:data_dim,3),0.99d0,answer,data_dim)
!!$
!!$		else
!!$
!!$			call splinderiv(p_perp_data(1:data_dim,1),p_perp_data(1:data_dim,2),  &
!!$						p_perp_data(1:data_dim,3),apsi/psic,answer,data_dim)
!!$
!!$		endif

		answer=answer/psic

	else

		if((dabs(psi) > 0.0d0).AND.(eq_type==3).AND.(eq3_opt==1)) then
			answer = dddpsi(psi)*tperpofpsi(psi)+dofpsi(psi)*dtperpdpsi(psi) 
			answer = answer*apsi/psi
		elseif((dabs(psi) > 0.0d0).AND.(eq_type==3).AND.(eq3_opt==2)) then
			answer = alpha*(qperpcenter-qperpedge)/psic**alpha*apsi**(alpha-1)
		elseif(eq3_opt==3) then
			answer = 0.d0
		else
			answer = 0.d0
		end if

		if(dabs(psi) > psic) answer = 0.d0

	endif

	if (apsi>0.d0) answer = answer*apsi/psi

  continue

end function dpperpdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function tparofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

	if (numerical_p_par) then

		answer = pparofpsi(psi)/dofpsi(psi)

	else

		if (eq3_opt == 1) then

			if(dabs(psi) > 0.0d0) then
				answer = tparedge + tparcenter*(dabs(psi)/psic)**alpha 
			else
				answer = tparedge
			end if

		elseif (eq3_opt >= 2) then

			answer = pparofpsi(psi)/dofpsi(psi)

		endif

	endif

	continue

end function tparofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function tperpofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: b0,theta,tpar



	if (numerical_p_perp) then

		answer = pperpofpsi(psi)/dofpsi(psi)

	else

		if((eq_type==3).and.(eq3_opt==1)) then

			b0 = bzero(psi)
			theta = thetaofpsi(psi)
			tpar = tparofpsi(psi)
			answer = tpar*b0/(b0-theta*tpar)

		else

			answer = pperpofpsi(psi)/dofpsi(psi)

		endif

	endif

	continue

  end function tperpofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dtpardpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

 	if (numerical_p_par) then

		answer = dppardpsi(psi)/dofpsi(psi) -  &
				pparofpsi(psi)*dddpsi(psi)/dofpsi(psi)**2

	else

	   if (eq3_opt == 1 ) then

			if(dabs(psi) > 0.0d0) then
				answer = alpha*tparcenter/psic**alpha*dabs(psi)**(alpha-1) 
			else
				answer = 0.d0
			end if

			if (dabs(psi)>0.d0) answer = answer*dabs(psi)/psi

		elseif (eq3_opt >= 2) then

			answer = dppardpsi(psi)/dofpsi(psi) -  &
					pparofpsi(psi)*dddpsi(psi)/dofpsi(psi)**2

		endif

	endif

	continue

  end function dtpardpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dtperpdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: tpar,b0,dtpar,db0,thet,dthet


	if (numerical_p_perp) then

		answer = dpperpdpsi(psi)/dofpsi(psi) -  &
				pperpofpsi(psi)*dddpsi(psi)/dofpsi(psi)**2

	else

		if((eq_type==3).and.(eq3_opt==1)) then

			if(dabs(psi) > 0.0d0) then

				tpar = tparofpsi(psi)
				dtpar = dtpardpsi(psi)
				b0 = bzero(psi)
				db0 = dbzerodpsi(psi)
				thet = thetaofpsi(psi)
				dthet = dthetadpsi(psi)
				answer = dtpar*b0/(b0-thet*tpar) &
						+ db0*tpar/(b0-thet*tpar) &
						- tpar*b0/(b0-thet*tpar)**2*(db0-thet*dtpar-dthet*tpar)
			else
				answer = 0.d0
			end if

		else

			answer = dpperpdpsi(psi)/dofpsi(psi) -  &
					pperpofpsi(psi)*dddpsi(psi)/dofpsi(psi)**2

		endif

	endif

	continue

  end function dtperpdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function thetaofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: psibar

	if(tri_type==11) then

		if (numerical_p_perp) then

			answer = (1.d0-tparofpsi(psi)/tperpofpsi(psi))  &
						*bpolzero(psi)/tparofpsi(psi)
!			answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
!						*bpol_min/tparofpsi(psi)

		else

		! toroidal field is 0

!		answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
!						*bpolzero(psi)/tparofpsi(psi)

!!		answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
!!						*bpolzero(psi)/tparofpsi(psi)  &
!!						*0.5d0*(1.d0+tanh(20.d0*(psi/psic-0.2d0)))

!		answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
!						*bpolzero(psi)/tparofpsi(psi)  &
!						*abs((psi_in-psi)*(psi-psi_out))* &
!						(4.d0/(psi_in-psi_out)**2-psi_in*psi_out)/psic**2

!		answer = 2.d0*bpol_min/tparofpsi(psic)*(betaperp_center-betapar_center)/  &
!		(betaperp_center+betapar_center)

!		answer = 2.d0*bpol_min/tparofpsi(psic)*(betaperp_center-betapar_center)/  &
!					min(betaperp_center,betapar_center)*  &
!					4.d0/psic**2*(psi_in-psi)*(psi-psi_out)

!		answer = bpol_min/tparofpsi(psi)*  &
!						(betaperp_center-betapar_center)/betapar_center

!		answer = bpol_min/tparofpsi(psi)*(psi/psic)**alpha_th*  &
!						(betaperp_center-betapar_center)/betapar_center

			if(eq3_opt==9) then

				psibar = max(abs(psi_pmax-psi_in),abs(psi_pmax-psi_out))
!				answer = theteps * ( 1.d0 - ((psi-psi_pmax)/psibar)**alpha_th )
				answer = theteps * exp(-alpha_th*((psi-psi_pmax)/psibar)**2)

			else

				answer = theteps*bpol_min/tparofpsi(psi)*(psi_in-psi)/psic*(psi/psic)**alpha_th*  &
								(betaperp_center-betapar_center)/betapar_center

			endif

!			write(77,*) psi, psic, tparofpsi(psi), answer

		endif

	elseif (numerical) then

		answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
					*bzero(psi)/tparofpsi(psi)
		continue

	else

		if((eq_type==3).and.(eq3_opt==1)) then

			if(dabs(psi)>0.d0) then
				answer = theteps*bzero(psi)/tparofpsi(psi)
			else
				answer = 0.d0
			endif

		else

		answer = (tperpofpsi(psi)-tparofpsi(psi))/tperpofpsi(psi)  &
					*bzero(psi)/tparofpsi(psi)

		endif

	endif


	continue

  end function thetaofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dthetadpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: tpar,tperp,psim,psibar

	if (fcenter>fvaccuum) then
		psim = psi
	else
		psim = -psi
	endif

	tpar = tparofpsi(psi)

	if(tri_type==11) then

		! toroidal field is 0

!		tperp = tperpofpsi(psi)

!		answer = dbpolzerodpsi(psi)*(1.d0/tpar-1.d0/tperp) +  &
!				bpolzero(psi)*(dtperpdpsi(psi)/tperp**2-dtpardpsi(psi)/tpar**2)

!!		answer = dbpolzerodpsi(psi)*(1.d0/tpar-1.d0/tperp) +  &
!!				bpolzero(psi)*(dtperpdpsi(psi)/tperp**2-dtpardpsi(psi)/tpar**2)  &
!!				*0.5d0*(1.d0+tanh(20.d0*(psi/psic-0.2d0))) +  &
!!				(tperp-tpar)/tperp*bpolzero(psi)/tpar  &
!!				*10.d0/psic*(cosh(20.d0*(psi/psic-0.2d0)))**-2

!		answer = (  &
!				dbpolzerodpsi(psi)*(1.d0/tpar-1.d0/tperp) +  &
!				bpolzero(psi)*(dtperpdpsi(psi)/tperp**2-dtpardpsi(psi)/tpar**2)  &
!				*(psi_in-psi)*(psi-psi_out) +  &
!				(tperp-tpar)/tperp*bpolzero(psi)/tpar  &
!				*(-2.d0*psi+psi_in-psi_out)  &
!				)*(4.d0/(psi_in-psi_out)**2-psi_in*psi_out)/psic**2

!		answer = 0.d0

!		answer = 2.d0*bpol_min/tparofpsi(psic)*(betaperp_center-betapar_center)/  &
!					min(betaperp_center,betapar_center)*  &
!					4.d0/psic**2*(psi_in+psi_out - 2.d0*psi)

!		answer = -bpol_min*(betaperp_center-betapar_center)/betapar_center  &
!						*dtpardpsi(psi)/tparofpsi(psi)**2

!!		answer = bpol_min*(betaperp_center-betapar_center)/betapar_center/psic**alpha_th*  &
!!						( -dtpardpsi(psi)/tparofpsi(psi)**2 * psi**alpha_th +  &
!!						 alpha_th*psi**(alpha_th-1.d0)/tparofpsi(psi))

		if(eq3_opt==9) then

			psibar = max(abs(psi_pmax-psi_in),abs(psi_pmax-psi_out))
!			answer = -theteps * alpha_th * ((psi-psi_pmax)/psibar)**(alpha_th-1.d0) * whatever
			answer = -2.d0*theteps*alpha_th*exp(-alpha_th*((psi-psi_pmax)/psibar)**2) *  &
					(psi-psi_pmax)/psibar**2
		else

			answer = theteps*bpol_min*(betaperp_center-betapar_center)/betapar_center/psic**(alpha_th+1.d0)/  &
						tpar * ( (psi_in-psi)*(alpha_th - psi*dtpardpsi(psi)/tpar) - psi)

		endif

	elseif (numerical) then

		tperp = tperpofpsi(psi)
		answer = dbzerodpsi(psi)*(1.d0/tpar-1.d0/tperp) +  &
				bzero(psi)*(dtperpdpsi(psi)/tperp**2-dtpardpsi(psi)/tpar**2)

	else

		if((eq_type==3).and.(eq3_opt==1)) then

			if(psim>0.d0)	then
				answer = theteps*(dbzerodpsi(psi)/tpar-bzero(psi)*dtpardpsi(psi)/tpar**2)
			else
				answer = 0.d0
			endif

		else

			tperp = tperpofpsi(psi)
			answer = dbzerodpsi(psi)*(1.d0/tpar-1.d0/tperp) +  &
					bzero(psi)*(dtperpdpsi(psi)/tperp**2-dtpardpsi(psi)/tpar**2)

		endif

	endif

	continue

end function dthetadpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function deltaofpsi(psi,i,j,rho,bphi,scale,nxx,nzz) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent (in) :: i,j
	real (kind=dkind), intent (in) :: rho,bphi
	real (kind=dkind), intent (in) :: scale
	integer, intent (in) :: nxx,nzz
    real (kind=dkind), dimension(1:3,1:3), intent(in) :: psi	! only surrounding points
    real (kind=dkind) :: answer
	real (kind=dkind) :: bfield,spar,sperp,tpar,tperp,bthet

	bthet = bpol(psi,i,j,scale,nxx,nzz)
	bfield = dsqrt(bthet**2 + bphi**2)

	if (eq_type==3) then

		tpar = tparofpsi(psi(2,2))
		tperp = tpar*bfield/(bfield-thetaofpsi(psi(2,2))*tpar)

		answer = rho*(tpar-tperp)*mu_mag/bfield**2

		if(dabs(1.d0-answer)<0.d0) write(*,*) answer

	else

		answer = 0.d0

	endif

	continue

end function deltaofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function volofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	apsi = dabs(psi)


	if (dabs(psi/psic)<=1.d-6) then

		answer = dbsval(0.d0, vol_ord, vol_data(1:ibreak_v,3),  &
			vol_dim, vol_cscoef(1,1:vol_dim) )

	elseif (dabs(psi/psic)>=1.d0) then

		answer = dbsval(1.d0, vol_ord, vol_data(1:ibreak_v,3),  &
			vol_dim, vol_cscoef(1,1:vol_dim) )

	else

		answer = dbsval(apsi/psic, vol_ord, vol_data(1:ibreak_v,3),  &
			vol_dim, vol_cscoef(1,1:vol_dim) )

	endif

	continue

end function volofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dvdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	apsi = dabs(psi)

	if (dabs(psi/psic)<=1.d-6) then

		answer = dbsder(1,0.d0, vol_ord, vol_data(1:ibreak_v,3),  &
				vol_dim, vol_cscoef(1,:) )

	elseif (dabs(psi/psic)>=1.d0) then

		answer = dbsder(1,1.d0, vol_ord, vol_data(1:ibreak_v,3),  &
						vol_dim, vol_cscoef(1,:) )

	else

		answer = dbsder(1,apsi/psic, vol_ord, vol_data(1:ibreak_v,3),  &
						vol_dim, vol_cscoef(1,:) )

		answer=answer/psic

	endif

	if(apsi>0.d0) answer = answer*psi/apsi

	continue

end function dvdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function pauxofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	! this is for "step" functions

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	answer = 0.5d0 * (1.d0 - tanh( aux_fact*(apsi/psic - psifrac) ) )

!	answer = 1.d0


	continue

end function pauxofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dpauxdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer
	real (kind=dkind) :: apsi

	if((bc_type/=7).and.(bc_type/=8)) then
		apsi=dabs(psi)
	else
		apsi = psi
	endif

	answer = -0.5d0 * aux_fact/psic /  &
				( cosh( aux_fact*(apsi/psic - psifrac) ) )**2


	if(apsi>0.d0) answer = answer*psi/apsi

!	answer = 0.d0

	continue

end function dpauxdpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function bpol(psi,i_ext,j_ext,scale,nxx,nzz) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer, intent (in) :: i_ext,j_ext
	real (kind=dkind), intent (in) :: scale
	integer, intent (in) :: nxx,nzz
	real (kind=dkind) :: answer
    real (kind=dkind), dimension(1:3,1:3), intent(in) :: psi
	real (kind=dkind) :: rloc,bx,bz
	integer i,j

!	print*, 'warning, function bpol needs to be fixed'

	  i=2
	  j=2

	  rloc = x_coord(i_ext)


	  ! Calculate B_r
	  if(j_ext == 1) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_r as independent of z near the boundary.
		 bx = (psi(i,j) - psi(i,j+1))/dz_a(1)/rloc
	  else if(j_ext == nzz) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_r as independent of z near the boundary.
		 bx = (psi(i,j-1) - psi(i,j))/dz_a(nzz)/rloc
	  else
		 ! use a centered difference
	!             bx(j,k) = 0.5d0*(psi(j,k-1) - psi(j,k+1))/dz/x
		 bx = ( dz_a(j_ext-1)**2*psi(i,j+1) +  &
					 (dz_a(j_ext)**2-dz_a(j_ext-1)**2)*psi(i,j) -  &
					 dz_a(j_ext)**2*psi(i,j-1) ) /  &
					 (dz_a(j_ext)*dz_a(j_ext-1)*(dz_a(j_ext)+dz_a(j_ext-1)))/rloc
	  end if

	  ! Calculate Bz
	  if(i_ext == 1) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_z as independent of x near the boundary.
		 bz = (psi(i+1,j) - psi(i,j))/dx_a(1)/(rloc + 0.5d0*dx_a(1))
	  else if(i_ext == nxx) then
		 ! use a 1-sided difference, equivalent to extrapolating
		 ! B_z as independent of x near the boundary.
		 bz = (psi(i,j) - psi(i-1,j))/dx_a(nxx)/(rloc - 0.5d0*dx_a(nxx))
	  else
		 ! use a centered difference
	!             bz(j,k) = 0.5d0*(psi(j+1,k) - psi(j-1,k))/dx/x
		 bz = ( dx_a(i_ext-1)**2*psi(i+1,j) +  &
					 (dx_a(i_ext)**2-dx_a(i_ext-1)**2)*psi(i,j) -  &
					 dx_a(i_ext)**2*psi(i-1,j) ) /  &
					 (dx_a(i_ext)*dx_a(i_ext-1)*(dx_a(i_ext)+dx_a(i_ext-1)))/rloc
	  end if



!!$    i=2
!!$	j=2
!!$!	rloc = rmajor - 0.5d0*x_size + (i-1)*dxr
!!$	rloc = rcenter - 0.5d0*x_size + (i-1)*dxr
!!$
!!$    if(i_ext == 1) then
!!$      ! use a 1-sided difference
!!$      ! print *, "One sided x-dif left"
!!$        dpsidx = (psi(i+1,j) - psi(i,j))/dxx
!!$    else if(i_ext == nxx) then
!!$      ! use a 1-sided difference
!!$      ! print *, "One sided x-dif right"
!!$        dpsidx = (psi(i,j) - psi(i-1,j))/dxx
!!$    else
!!$      ! use a centered difference
!!$        dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dxx
!!$    end if
!!$
!!$    if(j_ext == 1) then
!!$      ! use a 1-sided difference
!!$      ! print *, "One sided z-dif left"
!!$        dpsidz = (psi(i,j+1) - psi(i,j))/dzz
!!$    else if(j_ext == nzz) then
!!$      ! use a 1-sided difference
!!$      ! print *, "One sided z-dif right"
!!$        dpsidz = (psi(i,j) - psi(i,j-1))/dzz
!!$    else
!!$      ! use a centered difference
!!$        dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dzz
!!$    end if

!!$	answer = dsqrt((dpsidx**2+dpsidz**2))/rloc
	answer = (dsqrt(bx**2+bz**2))/scale
	!"scale" factors in the half/full grid size step

  continue

end function bpol

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function grav_potential(x,z) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only :   gravity_type, G_gravity, M_gravity

	real(kind=dkind) :: x, z, answer

	answer = 0.d0

	if(gravity_type==1) then
	! point mass

		answer = -G_gravity*M_gravity/sqrt(x**2+z**2)

	elseif(gravity_type==2) then
	! constant gravity (-R direction)

		answer = -G_gravity*M_gravity*x

	endif

end function grav_potential

  ! ------------------------------------------------------------------
  ! The following set of function are used to construct Hameiri's 
  ! functions of Psi from Betti & Freidberg's set.
  ! ------------------------------------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function sofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer, p, d
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = s_loc
		return
	endif

    p = pofpsi(psi,zone_loc)
    d = dofpsi(psi,zone_loc)

    answer = p/d**gamma

    return
  end function sofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dsdpsi(psi, zone) result(answer)		
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: d
    real (kind=dkind) :: answer
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif


	if(psi==psi_flag_ham) then
		answer = sp_loc
		return
	endif

    d = dofpsi(psi,zone_loc)

    answer = (dpdpsi(psi,zone_loc) - gamma*pofpsi(psi,zone_loc)*dddpsi(psi,zone_loc)/d)/d**gamma

    return
  end function dsdpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function phiofpsi(psi, zone) result(answer)		
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: cs, b0

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_ham) then
		answer = phi_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

!			b0 = bzero(psi)
!			if(abs(b0)<b_phi_zero*1.d-2)	b0 = b_phi_zero*1.d-2
!			answer = dsqrt(mu_mag*gamma*pofpsi(psi)*dofpsi(psi))* &
!				 mach_thetahat(psi)/b0

			b0 = b_phi_zero
			answer = dsqrt(mu_mag*gamma*pofpsi(psi,zone_loc)*dofpsi(psi,zone_loc))* &
				 mach_thetahat(psi,zone_loc)/b0

		else

			b0 = bzero(psi,zone_loc)

			if (eq_type==1) then

!	print*, psi,pofpsi(psi),dofpsi(psi)

				answer = dsqrt(mu_mag*gamma*pofpsi(psi,zone_loc)*dofpsi(psi,zone_loc))* &
					 mach_theta(psi,zone_loc)/b0

			elseif (eq_type==3)	then
				! no poloidal flow allowed
				answer = 0.d0
			endif

		endif

	endif

    return

  end function phiofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dphidpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: d,p,b,m,cs,ppar,pperp, db
    real (kind=dkind) :: answer

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_ham) then
		answer = phip_loc
		return
	endif

	if(F_opt==0) then
	! bzero is 0, let's start with no poloidal flow

		answer = 0.d0

	else

		d = dofpsi(psi,zone_loc)
		p = pofpsi(psi,zone_loc)

		if((Broot==4).or.(Broot==5)) then
		! RFP case, bzero can be 0, better replace it with a constant

			b = b_phi_zero
!			b = bzero(psi)
!			db = dbzerodpsi(psi)

!			if(abs(b)<b_phi_zero*1.d-2) then
!				b = b_phi_zero*1.d-2
!				db = 0.d0
!			endif

			m = mach_thetahat(psi,zone_loc)

			answer = dsqrt(mu_mag*gamma*p*d)/b*( &
				dmach_thetahatdpsi(psi,zone_loc) &
				+ 0.5d0*m*(dddpsi(psi,zone_loc)/d + dpdpsi(psi,zone_loc)/p) &
				)

!					answer = dsqrt(mu_mag*gamma*p*d)/b*( &
!						dmach_thetadpsi(psi) &
!						+ 0.5d0*m*(dddpsi(psi)/d + dpdpsi(psi)/p) &
!						- db*m/b)

		else

			b = bzero(psi,zone_loc)

			if(d > 0.0d0 .and. dabs(b) > 0.0d0 .and. p > 0.0d0) then
			   m = mach_theta(psi,zone_loc)

			   if (eq_type==1) then

					answer = dsqrt(mu_mag*gamma*p*d)/b*( &
						dmach_thetadpsi(psi,zone_loc) &
						+ 0.5d0*m*(dddpsi(psi,zone_loc)/d + dpdpsi(psi,zone_loc)/p) &
						- dbzerodpsi(psi,zone_loc)*m/b)

				elseif (eq_type==3) then
					answer = 0.d0
				endif

			else
			   answer = 0.0d0
			end if

		endif

	endif

    return

end function dphidpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function omegaofpsi(psi, zone) result(answer)		
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
	real (kind=dkind) :: cs
	
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_ham) then
		answer = omega_loc
		return
	endif

	if (eq_type==1) then

	    answer = dsqrt(gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc))* &
		      (mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc))/rmajor

	elseif (eq_type==3) then

		cs = dsqrt(tparofpsi(psi))

		answer = cs*mach_phi(psi)/rmajor			

	endif

    return
  end function omegaofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function domegadpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: p,d,cs,dcs,ppa,ppe
    real (kind=dkind) :: answer

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_ham) then
		answer = omegap_loc
		return
	endif

    p = pofpsi(psi,zone_loc)
	ppa=pparofpsi(psi)
!	ppe=pperpofpsi(psi)

    if((p > 0.0d0).AND.(eq_type==1)) then

       d = dofpsi(psi,zone_loc)

       answer = (0.5d0*dsqrt(gamma/(p*d))*(dpdpsi(psi,zone_loc) -dddpsi(psi,zone_loc)*p/d)* &
            (mach_phi(psi,zone_loc) - mach_theta(psi,zone_loc)) + &
            dsqrt(gamma*p/d)*(dmach_phidpsi(psi,zone_loc) - dmach_thetadpsi(psi,zone_loc)))/ &
            rmajor

!	elseif(((ppa+ppe)>0.d0).AND.(eq_type==3)) then
	elseif((ppa>0.d0).AND.(eq_type==3)) then

		cs = dsqrt(tparofpsi(psi))
		dcs = 0.5d0*tparofpsi(psi)**(-0.5d0)*dtpardpsi(psi)


		answer = ( dcs*mach_phi(psi) +  &
					cs*dmach_phidpsi(psi)  )/rmajor

	else
       answer = 0.0d0
    end if

    return
  end function domegadpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function iofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer
    
    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif


	if(psi==psi_flag_ham) then
		answer = i_loc
		return
	endif

    answer = rmajor*bzero(psi,zone_loc)/dsqrt(mu_mag)

    return
  end function iofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function didpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: answer

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = ip_loc
		return
	endif

 		  answer = rmajor*dbzerodpsi(psi,zone_loc)/dsqrt(mu_mag)

    return
  end function didpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function hofpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: mt,mp,cs
    real (kind=dkind) :: answer

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif

	if(psi==psi_flag_ham) then
		answer = h_loc
		return
	endif

    mt = mach_theta(psi,zone_loc)
    mp = mach_phi(psi,zone_loc)

	if (eq_type==1) then

		answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &
			gamma*pofpsi(psi,zone_loc)/dofpsi(psi,zone_loc)

	elseif (eq_type==3) then

		cs = dsqrt(tparofpsi(psi))
		answer = -(cs**2)* 0.5d0*mp**2  +  &
					wofpsi(psi)

	endif

    if((answer <= 0.0d0).AND.(eq_type.ne.3)) then
		if(Broot/=3) then
		   print *, "Error H(Psi) = ",answer
		   pause
		   stop
	    endif
    end if

    return
  end function hofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dhdpsi(psi, zone) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent (in) :: psi
    integer, optional, intent (in) :: zone
    integer :: zone_loc
    real (kind=dkind) :: d,p,mt,mp
	real (kind=dkind) :: ppa,ppe,cs,dcs
    real (kind=dkind) :: answer

    if(present(zone)) then
    	zone_loc = zone
    else
    	zone_loc = 0
    endif
    
	if(psi==psi_flag_ham) then
		answer = hp_loc
		return
	endif

    d = dofpsi(psi,zone_loc)
    p = pofpsi(psi,zone_loc)
    mt = mach_theta(psi,zone_loc)
    mp = mach_phi(psi,zone_loc)
	ppa=pparofpsi(psi)
!	ppe=pperpofpsi(psi)


    if((d>0.0d0).AND.(eq_type==1)) then

	   answer = (1.0d0/(gamma-1.0d0) + mt*mp - 0.5d0*mp**2)* &
            gamma*(dpdpsi(psi,zone_loc) - dddpsi(psi,zone_loc)*p/d)/d + &
            ((mt-mp)*dmach_phidpsi(psi,zone_loc) + mp*dmach_thetadpsi(psi,zone_loc))* &
            gamma*p/d

	elseif((d>0.d0).AND.(eq_type==3)) then

		cs = dsqrt(ppa/d)
		dcs = 0.5d0*tparofpsi(psi)**(-0.5d0)*dtpardpsi(psi)

		answer = -cs*( dcs*mp**2 + cs*mp*dmach_phidpsi(psi)) + &
					dwdpsi(psi)

	else
       answer = 0.0d0
    end if

    return
  end function dhdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function wofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: psi	!
	real (kind=dkind) :: answer, logarg

! IMPORTANT: THIS IS NOT W(psi), BUT SOMETHING SIMILAR TO IT USED TO DEFINE H(psi)


	if (eq_type==3) then

		if(tri_type==11) then

!			write(88,*) psi,psic,tparofpsi(psi),bpolzero(psi),thetaofpsi(psi),  &
!				dabs(bpolzero(psi)-thetaofpsi(psi)*tparofpsi(psi)) / bpolzero(psi)
!			write(88,*) '   '

			logarg = dabs(bpolzero(psi)-thetaofpsi(psi)*tparofpsi(psi)) / bpolzero(psi)

!!$			if(logarg>0.d0) then
!!$
!!$				answer = tparofpsi(psi) * log(logarg)
!!$
!!$			else
!!$
!!$				answer = 0.d0
!!$
!!$			endif

			if(logarg<=0.d0) then
				logarg = 1.d-22 - logarg
			endif

			answer = tparofpsi(psi) * log(logarg)


!			answer = tparofpsi(psi) * dlog( &
!					dabs(bpolzero(psi)-thetaofpsi(psi)*tparofpsi(psi)) / bpolzero(psi))

		else

			answer = tparofpsi(psi) * dlog( &
					dabs(bzero(psi)-thetaofpsi(psi)*tparofpsi(psi)) / bzero(psi))

		endif

	else

		answer = 0.d0

	endif

	continue

  end function wofpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dwdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: psi
	real (kind=dkind) :: answer
    real (kind=dkind) :: d,tpar,theta,b0,p,mp, b0p, logarg

    d = dofpsi(psi)
    p = pofpsi(psi)
    mp = mach_phi(psi)

	if (eq_type==3) then

		if (dabs(psi)>0.d0) then

			if(tri_type==11) then

				b0 = bpolzero(psi)
				b0p = dbpolzerodpsi(psi)

			else

				b0 = bzero(psi)
				b0p = dbzerodpsi(psi)

			endif

			theta = thetaofpsi(psi)
			d = dofpsi(psi)
			tpar = tparofpsi(psi)

			logarg = dabs(b0-theta*tpar)/b0

			if(logarg<=0.d0) logarg = 1.d-22 - logarg

				answer = ( dtpardpsi(psi) * dlog(logarg) + &
					tpar*b0 / ( dabs(b0-theta*tpar) ) * &
					( -b0p*dabs(b0-theta*tpar)/b0**2 &
					  +  sign(1.d0,(b0-theta*tpar))/b0 *  &
		!  				  +  dabs(sign(1.d0,(b0-theta*tpar))/b0) *  &
						( b0p-theta*dtpardpsi(psi)-dthetadpsi(psi)*tpar ) ) )  

!!$			if(logarg>0.d0) then
!!$
!!$				answer = ( dtpardpsi(psi) * dlog( dabs(b0-theta*tpar)/b0) + &
!!$					tpar*b0 / ( dabs(b0-theta*tpar) ) * &
!!$					( -b0p*dabs(b0-theta*tpar)/b0**2 &
!!$					  +  sign(1.d0,(b0-theta*tpar))/b0 *  &
!!$		!  				  +  dabs(sign(1.d0,(b0-theta*tpar))/b0) *  &
!!$						( b0p-theta*dtpardpsi(psi)-dthetadpsi(psi)*tpar ) ) )  
!!$			else
!!$
!!$				answer = 0.d0
!!$
!!$			endif


		else
			answer = 0.d0
		endif

	else

		answer = 0.d0

	endif

	continue

  end function dwdpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function sparofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

	answer=pparofpsi(psi)*bzero(psi)**2/dofpsi(psi)**3

	continue

  end function sparofpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function sperpofpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

	answer=pperpofpsi(psi)/dofpsi(psi)/bzero(psi)

	continue

  end function sperpofpsi


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dspardpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

	answer=dppardpsi(psi)*bzero(psi)**2/dofpsi(psi)**3  &
			-3.d0*pparofpsi(psi)*bzero(psi)**2/dofpsi(psi)**4*dddpsi(psi)  &
			+2.d0*pparofpsi(psi)*bzero(psi)/dofpsi(psi)**3*dbzerodpsi(psi)

  continue

  end function dspardpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dsperpdpsi(psi) result (answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real (kind=dkind), intent (in) :: psi
	real (kind=dkind) :: answer

	answer=dpperpdpsi(psi)/(dofpsi(psi)*bzero(psi))	&
			- pperpofpsi(psi)/(dofpsi(psi)**2*bzero(psi))*dddpsi(psi)  &
			- pperpofpsi(psi)/(dofpsi(psi)*bzero(psi)**2)*dbzerodpsi(psi)  

	continue

  end function dsperpdpsi



  ! ------------------------------------------------------------------


  ! This function calculates the residual of the Bernoulli equation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual,residual1,residual2,residual3
	real (kind=dkind) :: term1,term2,term3
	real (kind=dkind) ::term4,term5,termperp,termpar
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli]: rho = ",rho
       print *, "[bernoulli]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli]: rho must be positive"
       stop
    end if


	b_torc = g_b_torc
	b_polc = g_b_polc
	bfield = dsqrt(b_torc**2+b_polc**2)


	if (eq_type==1) then

		residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
					- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
					- 0.5d0/mu_mag*(g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2 &
					- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
					- g_Lambda

	elseif (eq_type==3) then
		residual = g_H + 0.5d0*(g_r*g_Omega)**2 &	! no poloidal flow allowed!
				  - g_tpar*dlog(rho/(g_D*bfield)*dabs(bfield-g_theta*g_tpar))
	endif

    return
  end function bernoulli

  ! ------------------------------------------------------------------

  ! This function calculates the residual of the Bernoulli equation
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_gauss(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual,residual1,residual2,residual3
	real (kind=dkind) :: term1,term2,term3
	real (kind=dkind) ::term4,term5,termperp,termpar
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli_gauss]: rho = ",rho
       print *, "[bernoulli_gauss]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli_gauss]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli_gauss]: rho must be positive"
       stop
    end if


	b_torc = g_b_torc
	b_polc = g_b_polc
	bfield = dsqrt(b_torc**2+b_polc**2)


		residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
					- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
					- 0.5d0/mu_mag*(g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2 &
					- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
					- g_Lambda  &
					- fBernmax * exp(-delta_Bern*((psi_Bern-psi_degen)/psic)**2)

    return

  end function bernoulli_gauss


  ! This function calculates the first derivative of the Bernoulli 
  ! equation with respect to rho.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_drho(rho) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: db
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[dbern_drho]: rho = ",rho
       print *, "[dbern_drho]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[dbern_drho]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[dbern_drho]: rho must be positive"
       stop
    end if

	b_torc = g_b_torc
	b_polc = g_b_polc
	bfield = dsqrt(b_polc**2+b_torc**2)

	if (eq_type==1) then
		db = - gamma*g_S*rho**(gamma-2.0d0) &
			  + (g_dpsidx**2 + g_dpsidz**2)*(g_Phi/(rho*g_r))**2/rho/mu_mag &
		      + (g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2/ &
			  (rho - g_Phi**2)
	elseif (eq_type==3) then
		db = -g_tpar/rho	! no poloidal flow allowed!
	endif

    return
  end function dbern_drho


  ! ------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_L(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual,residual1,residual2,residual3
	real (kind=dkind) :: term1,term2,term3
	real (kind=dkind) ::term4,term5,termperp,termpar
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli]: rho = ",rho
       print *, "[bernoulli]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli]: rho must be positive"
       stop
    end if


	b_torc = g_b_torc
	b_polc = sqrt(g_dpsidx_L**2 + g_dpsidz_L**2)
	bfield = dsqrt(b_torc**2+b_polc**2)

	residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
				- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
				- 0.5d0/mu_mag*(g_dpsidx_L**2 + g_dpsidz_L**2)*(g_Phi/(rho*g_r))**2 &
				- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
				- g_Lambda


    return
  end function bernoulli_L

  ! ------------------------------------------------------------------


  ! This function calculates the first derivative of the Bernoulli 
  ! equation with respect to rho.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_drho_L(rho) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: db
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[dbern_drho]: rho = ",rho
       print *, "[dbern_drho]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[dbern_drho]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[dbern_drho]: rho must be positive"
       stop
    end if

	b_torc = g_b_torc
	b_polc = sqrt(g_dpsidx_L**2 + g_dpsidz_L**2)
	bfield = dsqrt(b_polc**2+b_torc**2)

	db = - gamma*g_S*rho**(gamma-2.0d0) &
		  + (g_dpsidx_L**2 + g_dpsidz_L**2)*(g_Phi/(rho*g_r))**2/rho/mu_mag &
	      + (g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2/ &
		  (rho - g_Phi**2)

    return
  end function dbern_drho_L


  ! ------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function bernoulli_R(rho) result (residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: residual,residual1,residual2,residual3
	real (kind=dkind) :: term1,term2,term3
	real (kind=dkind) ::term4,term5,termperp,termpar
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[bernoulli]: rho = ",rho
       print *, "[bernoulli]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[bernoulli]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[bernoulli]: rho must be positive"
       stop
    end if


	b_torc = g_b_torc
	b_polc = sqrt(g_dpsidx_R**2 + g_dpsidz_R**2)
	bfield = dsqrt(b_torc**2+b_polc**2)

	residual = g_H + 0.5d0*(g_r*g_Omega)**2 &
				- gamma/(gamma-1.0d0)*g_S*rho**(gamma-1.0d0) &
				- 0.5d0/mu_mag*(g_dpsidx_R**2 + g_dpsidz_R**2)*(g_Phi/(rho*g_r))**2 &
				- 0.5d0*(g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2   &
				- g_Lambda

    return

  end function bernoulli_R

  ! ------------------------------------------------------------------


  ! This function calculates the first derivative of the Bernoulli 
  ! equation with respect to rho.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function dbern_drho_R(rho) result (db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind) :: db
	real (kind=dkind) :: b_torc,b_polc,bfield

    if(rho == g_Phi*g_Phi .or. rho == 0.0d0 ) then
       print *, "[dbern_drho]: rho = ",rho
       print *, "[dbern_drho]: rho(Alfvenic) = ",g_Phi*g_Phi
       print *, "[dbern_drho]: This is not yet handled"
       stop
    else if(rho <= 0.0d0) then
       print *, "[dbern_drho]: rho must be positive"
       stop
    end if

	b_torc = g_b_torc
	b_polc = sqrt(g_dpsidx_R**2 + g_dpsidz_R**2)
	bfield = dsqrt(b_polc**2+b_torc**2)

	db = - gamma*g_S*rho**(gamma-2.0d0) &
		  + (g_dpsidx_R**2 + g_dpsidz_R**2)*(g_Phi/(rho*g_r))**2/rho/mu_mag &
	      + (g_Phi*(g_I/g_r + g_r*g_Omega*g_Phi)/(rho - g_Phi**2))**2/ &
		  (rho - g_Phi**2)

    return

  end function dbern_drho_R


  ! ------------------------------------------------------------------


  ! This function calculates the residual of the Bernoulli equation
  ! as well as the first derivative of the Bernoulli equation with 
  ! respect to rho for use in a Newton-Raphson solver.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap(rho,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to rho.

    b = bernoulli(rho)

    db = dbern_drho(rho)

  end subroutine newt_rap


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine newt_rap_gauss(rho,b,db)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: rho
    real (kind=dkind), intent(inout) :: b, db
    ! b is the residual of the bernoulli function
    ! and db is the derivative of this function with respect to rho.

    b = bernoulli_gauss(rho)

    db = dbern_drho(rho)

!	print*, rho, b, db
!	pause

  end subroutine newt_rap_gauss


  ! ------------------------------------------------------------------


  ! Initialize u(1..nx , 1..nz) with an initial guess to the 
  ! solution which respects the boundary conditions u = 0.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine guess_soln(u,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: u
    integer :: i,j,k
    real (kind=dkind) :: ex,ez,th,r_i,r_o,psie
	integer :: dummy_int
	real(kind=dkind), dimension(1:3) :: dummy



	if(tri_type==11) psic = max(psi_in,psi_out)

    ! Put some peaked initial distribution in the center
    do i=1, nx
       do j=1, nz
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
              u(i,j) = 0.0d0
          else

		  if((tri_type==11).or.(tri_type==13)) then

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,th,rminor,  &
								dummy_int,dummy(1),dummy(2),dummy(3))
			!this is required to get rminor

			r_o = dbsval(th, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )

			r_i = dbsval(th, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

		  else

			call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
			!this is required to get rminor

		  endif



			if(tri_type==11) then

!				u(i,j)=(psi_in+psi_out)/2.d0 * (1.d0 +  &
!						1.d-0*(rminor**2 - ex*ex - ez*ez)/rminor**2)

				u(i,j) = ( psi_in - (psi_in-psi_out)/(r_o-r_i) *  &
								(sqrt(ex**2+ez**2)-r_i) )

			elseif((bc_type==7).or.(bc_type==8)) then

				 psie = dbsval(th, psib_ord, psib_data(1:ibreak_pb,3),  &
							ibreak_pb-psib_ord, psib_cscoef(1,1:ibreak_pb-psib_ord) )

				if(tri_type==13) then

					u(i,j) = (psic + (psie-psic)*(ex**2+ez**2)/rminor**2) *  &
							abs(r_i-sqrt(ex**2+ez**2))/max(abs(r_i-sqrt(ex**2+ez**2)),1.d-2)

				else

					u(i,j) = psic + (psie-psic)*(ex**2+ez**2)/rminor**2

				endif

			else
				u(i,j) = psic *(rminor**2 - ex*ex - ez*ez)/rminor**2
!             u(i,j) = psic*( (rminor**2 - ex*ex - ez*ez) )**2/rminor**4
!             u(i,j) = psic*dsqrt( (rminor**2 - ex*ex - ez*ez) )/rminor

!				u(i,j) = (x_coord(i)/rmajor)**8*psic *(rminor**2 - ex*ex - ez*ez)/rminor**2


			endif

          end if
       end do
    end do

	psic = maxval(u)

!!$	u = 0.d0
!!$
!!$	do i = 2*nx/5,3*nx/5
!!$	do j = 2*nz/5,3*nz/5
!!$
!!$		u(i,j) = psic
!!$
!!$	enddo
!!$	enddo

	continue

  end subroutine guess_soln


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine initialize_density(psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: rho
    integer :: i,j

    do i=1, nx
       do j=1, nz
          rho(i,j) = dofpsi(psi(i,j))
       end do
    end do

	call bc_psi_rho0(psi,rho,nx,nz)

	continue

  end subroutine initialize_density

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine initialize_b(b,psi,rho,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    use constant, only: x_coord, dx, dz
	implicit none
	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,rho
    real (kind=dkind), dimension(1:nx,1:nz), intent(out) :: b
	real (kind=dkind) :: r,ex,ez
    integer :: i,j,k

	if(eq_type==1) then

		b(1:nx,1:nz) = 0.d0
		return

	endif



	do i=1,nx

		r = x_coord(i)

		do j=1,nz


			if(sort_grid(i,j,0)==1) then

				b(i,j) = dsqrt(mu_mag)*(iofpsi(psi(i,j))/r + r*Omegaofpsi(psi(i,j))  &
								* phiofpsi(psi(i,j)))  &
						/ (1.d0-(phiofpsi(psi(i,j))**2)/rho(i,j))


			else

				b(i,j) = 0.d0

			endif

		enddo
	enddo

	continue
	return

  end subroutine initialize_b

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine update_b(psi,rho,b,nx,nz,orp,anorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    use constant
	use triangularity
	implicit none
	integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,rho
	real (kind=dkind), intent(in) :: orp
	real (kind=dkind), intent(inout) :: anorm
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: b
	real (kind=dkind) :: r,ex,ez,res
	real (kind=dkind), dimension(1:3,1:3) :: psi3x3
    integer :: i,j,iii,jjj,k

    if(eq_type==1) return

	if (F_opt==0) then
	! no toroidal field, for levitated dipole

		anorm = 0.d0
		return

	endif


	do i=1,nx

		r = x_coord(i)

		do j=1,nz

		    if(sort_grid(i,j,0)==-1) then
             			cycle
			end if


			if (i==1) then

				do iii=1,3
				psi3x3(1,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo
					do iii=1,2
					do jjj=1,2
					psi3x3(iii+1,jjj+1)=psi(iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo
					do iii=1,2
					do jjj=0,1
					psi3x3(iii+1,jjj+1)=psi(iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii+1,jjj+2)=psi(iii,j+jjj)
					enddo
					enddo
				endif

				elseif (i==nx) then

				do iii=1,3
				psi3x3(3,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo
					do iii=1,2
					do jjj=1,2
					psi3x3(iii,jjj+1)=psi(i-2+iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo
					do iii=1,2
					do jjj=0,1
					psi3x3(iii,jjj+1)=psi(i-2+iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii,jjj+2)=psi(i-2+iii,j+jjj)
					enddo
					enddo
				endif

				else

				if (j==1) then
					do iii=-1,1
						psi3x3(iii+2,1)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj+1)=psi(i+iii,jjj)
						enddo
					enddo
				elseif (j==nz) then
					do iii=-1,1
						psi3x3(iii+2,3)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj)=psi(i+iii,j-2+jjj)
						enddo
					enddo
				else
					do iii=-1,1
					do jjj=-1,1
						psi3x3(iii+2,jjj+2)=psi(i+iii,j+jjj)
					enddo
					enddo

				endif

			endif

			if (eq_type==3) then

				res= b(i,j) - (dsqrt(mu_mag)*iofpsi(psi(i,j))/r)  &
							/ ( 1.d0 - deltaofpsi(psi3x3,i,j,  &
									rho(i,j),b(i,j),1.d0,nx,nz) )
				b(i,j) = b(i,j) - orp*res
				anorm = anorm + dabs(res)

			endif

			continue
		enddo
	enddo

	continue
	return

  end subroutine update_b

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine get_psi_pres(psi,nn)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !for LDX equilibria with eq3_opt=7

	integer :: nn
	real(kind=dkind), dimension(nn) :: psi
	integer :: i
	real(kind=dkind) :: temp_l, temp_r

	i = 0

	if(bc_type==6) then

		psi_pres = psi_in
		return

	endif

	do

		i = i+1
		if(sort_grid(i,(nn+1)/2,1)==4) then
			i = i-1
			temp_l = psi(i)
			exit
		endif

	enddo

	do

		i = i+1
		if(sort_grid(i,(nn+1)/2,1)<4) then
			temp_r = psi(i)
			exit
		endif

	enddo

	psi_pres = min(temp_l,temp_r)

	continue

  end subroutine get_psi_pres


  ! Test for symmetry in z
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine symmetric(q,nx,nz,name)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: q 
    character (len=*), intent(in) :: name
    integer :: riz,liz
    integer :: ir
    real (kind=dkind) :: delta, mx

    mx = 0.0d0 ! Initialize the max variable

    do liz=1, nz
       riz = nz + 1 - liz
       if(riz <= liz) exit

       do ir=1, nx
          delta = dabs(q(ir,liz) - q(ir,riz))
          mx = dmax1(mx,delta)
       end do
    end do

    if(mx > 0.0d0) then
       print *, name," is NOT symmetric, Max | Difference | = ",mx
    else
       print *, name," is symmetric"
    end if

  end subroutine symmetric


  ! van Leer Slope Limiter
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function van_Leer_slope(l,c,r,dx) result (slope)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: l,c,r,dx
    real (kind=dkind) :: slope, up, down

    up = (r-c)/dx
    down = (c-l)/dx
    if((up==0.0d0) .or. (down==0.0)) then
       slope = 0.0d0
    else
       slope = (up*dabs(down) + down*dabs(up))/(dabs(up)+dabs(down))
    endif

    return
  end function van_Leer_slope

  ! ------------------------------------------------------------------

  ! van Leer Slope Limiter
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function van_Leer_slope_new(l,c,r,dxl,dxr) result (slope)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real (kind=dkind), intent(in) :: l,c,r,dxl,dxr
    real (kind=dkind) :: slope, up, down

    up = (r-c)/dxr
    down = (c-l)/dxl
    if((up==0.0d0) .or. (down==0.0)) then
       slope = 0.0d0
    else
       slope = (up*dabs(down) + down*dabs(up))/(dabs(up)+dabs(down))
    endif

    return
  end function van_Leer_slope_new



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic
	real (kind=dkind), dimension (1:3,1:3) :: psi3x3
	integer :: i_zone = 0
	integer iii,jjj,k

    ! -----------------------------------------------------


	if(eq_type==1) then
		if(Broot<4) then
			mtm_limit = 1.d0
		else
			mtm_limit = 1.d2
		endif
	endif

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

!		  rho(i,j) = dofpsi(psi(i,j))
!		  cycle

!!$		  if( ( (tri_type==13).and.(sort_grid(i,j,0)==1) )  &
!!$				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
!!$				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
			! Edited March 18 2021, to allow for the regular solution in the outer region when bc_type==7
			! Edited again Jannuary 21 2022, to allow for the regular solution in the outer region when bc_type==7 and tri_type==13
		  if ( (tri_type==13).and.(sort_grid(i,j,0)<1) )  then
!		  if ( (tri_type==13).and.(sort_grid(i,j,0)==1) )  then
!!!		  if( ( (tri_type==13).and.(sort_grid(i,j,0)==1) )  &
!!!				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  )then
				if(gravity_type==0) then
					rho(i,j) = dofpsi(0.d0)
				else
					rho(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(i,j,0)-2)
				i_zone = max(-1,i_zone)
			endif		

          !  "Setup Global Parameters"
          phic = phiofpsi(psi(i,j),i_zone)
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j),i_zone)
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j),i_zone)
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j),i_zone)
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j),i_zone)

			g_D = dofpsi(psi(i,j),i_zone)		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  if (eq_type==3) then


		  if (i==1) then

				do iii=1,3
				psi3x3(1,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo
					do iii=1,2
					do jjj=1,2
					psi3x3(iii+1,jjj+1)=psi(iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo			
					do iii=1,2
					do jjj=0,1
					psi3x3(iii+1,jjj+1)=psi(iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii+1,jjj+2)=psi(iii,j+jjj)
					enddo
					enddo
				endif

			elseif (i==nx) then

				do iii=1,3
				psi3x3(3,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo			
					do iii=1,2
					do jjj=1,2
					psi3x3(iii,jjj+1)=psi(i-2+iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo			
					do iii=1,2
					do jjj=0,1
					psi3x3(iii,jjj+1)=psi(i-2+iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii,jjj+2)=psi(i-2+iii,j+jjj)
					enddo
					enddo
				endif

			else

				if (j==1) then
					do iii=-1,1
						psi3x3(iii+2,1)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj+1)=psi(i+iii,jjj)
						enddo
					enddo
				elseif (j==nz) then
					do iii=-1,1
						psi3x3(iii+2,3)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj)=psi(i+iii,j-2+jjj)
						enddo
					enddo
				else
					do iii=-1,1
					do jjj=-1,1
						psi3x3(iii+2,jjj+2)=psi(i+iii,j+jjj)
					enddo
					enddo

				endif

			endif

!		  	g_delta = deltaofpsi(psi3x3,i,j,  &
!							rho(i,j),b_phi(i,j),dx,dx,dz,nx,nz)
		  	g_delta = deltaofpsi(psi3x3,i,j,  &
							rho(i,j),b_phi(i,j),1.d0,nx,nz)
			g_wofpsi = wofpsi(psi(i,j)) ! should be	g_wofpsi = wofpsi(psi(i,j),i_zone) if that was implemented
			term1=g_wofpsi
!			g_spar = sparofpsi(psi(i,j))
!			g_sperp = sperpofpsi(psi(i,j))
			g_mtheta=mach_theta(psi(i,j),i_zone)
			g_theta = thetaofpsi(psi(i,j)) ! should be g_theta = thetaofpsi(psi(i,j),i_zone) if that was implemented
			g_tpar = tparofpsi(psi(i,j)) ! should be g_tpar = tparofpsi(psi(i,j),i_zone) if that was implemented
			g_D = dofpsi(psi(i,j),i_zone)

		  endif

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

		  if (eq_type==3) then
			g_b_torc = b_phi(i,j)
		  elseif (eq_type == 1) then
			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))
		  endif

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho

		  if (eq_type==1) then

				rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  elseif (eq_type==3) then

 		  		rhomax = g_D*g_bfield/dabs(g_bfield-g_theta*g_tpar) &
							*dexp((2.d0*g_H + (g_r*g_omega)**2)/(2.d0*g_tpar))

		  endif

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(rhomin,rhomax,1000)
!!!					pause
!!!					rho(i,j) = dofpsi(psi(i,j))
!!!					cycle

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j),i_zone)
                   g_Omega = omegaofpsi(psi(i,j),i_zone)
                   g_H = hofpsi(psi(i,j),i_zone)

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)


!!$					   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz),i_zone)		
	   g_I = iofpsi(psi(min_ix,min_iz),i_zone)
	   g_D = dofpsi(psi(min_ix,min_iz),i_zone)
       ! Calculate the derivatives of psi
!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )

       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz),i_zone)
          g_Omega = omegaofpsi(psi(min_ix,min_iz),i_zone)
          g_H = hofpsi(psi(min_ix,min_iz),i_zone)


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

!!$			   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz),i_zone)
          g_Omega = omegaofpsi(psi(min_ix,min_iz),i_zone)
          g_H = hofpsi(psi(min_ix,min_iz),i_zone)
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

  end subroutine update_rho

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_frozen_March_17_2021(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! We are updating the Bernoulli solver ro use sort_grid instead of local calculations. Freeze the previous version
! for safety and possibly allow the use of either.
	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic

	real (kind=dkind), dimension (1:3,1:3) :: psi3x3
	integer iii,jjj,k

    ! -----------------------------------------------------


	if(eq_type==1) then
		if(Broot<4) then
			mtm_limit = 1.d0
		else
			mtm_limit = 1.d2
		endif
	endif

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

!		  rho(i,j) = dofpsi(psi(i,j))
!		  cycle

		  if( ( (tri_type==13).and.(sort_grid(i,j,0)==1) )  &
				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
				.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then
!		  if( ( (tri_type==13).and.(sort_grid(i,j,0)==1) )  &
!				.or. (((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  )then
				if(gravity_type==0) then
					rho(i,j) = dofpsi(0.d0)
				else
					rho(i,j) = ( (gamma-1.d0)/gamma *  &
								(hofpsi(0.d0) - grav_potential(x_coord(i),z_coord(j)) + 0.5d0*x_coord(i)**2*omegaofpsi(0.d0)**2)  &
								/sofpsi(0.d0)  &
								)**(1.d0/(gamma-1.d0))
				endif
				 cycle
          end if


          !  "Setup Global Parameters"
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))

		  if (eq_type==3) then


		  if (i==1) then

				do iii=1,3
				psi3x3(1,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo
					do iii=1,2
					do jjj=1,2
					psi3x3(iii+1,jjj+1)=psi(iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo			
					do iii=1,2
					do jjj=0,1
					psi3x3(iii+1,jjj+1)=psi(iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii+1,jjj+2)=psi(iii,j+jjj)
					enddo
					enddo
				endif

			elseif (i==nx) then

				do iii=1,3
				psi3x3(3,iii)=0.d0
				enddo
				if (j==1) then
					do jjj=1,3
					psi3x3(jjj,1)=0.d0
					enddo			
					do iii=1,2
					do jjj=1,2
					psi3x3(iii,jjj+1)=psi(i-2+iii,jjj)
					enddo
					enddo
				elseif (j==nz) then
					do jjj=1,3
					psi3x3(jjj,3)=0.d0
					enddo			
					do iii=1,2
					do jjj=0,1
					psi3x3(iii,jjj+1)=psi(i-2+iii,j-1+jjj)
					enddo
					enddo
				else
					do iii=1,2
					do jjj=-1,1
					psi3x3(iii,jjj+2)=psi(i-2+iii,j+jjj)
					enddo
					enddo
				endif

			else

				if (j==1) then
					do iii=-1,1
						psi3x3(iii+2,1)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj+1)=psi(i+iii,jjj)
						enddo
					enddo
				elseif (j==nz) then
					do iii=-1,1
						psi3x3(iii+2,3)=0.d0
						do jjj=1,2
							psi3x3(iii+2,jjj)=psi(i+iii,j-2+jjj)
						enddo
					enddo
				else
					do iii=-1,1
					do jjj=-1,1
						psi3x3(iii+2,jjj+2)=psi(i+iii,j+jjj)
					enddo
					enddo

				endif

			endif

!		  	g_delta = deltaofpsi(psi3x3,i,j,  &
!							rho(i,j),b_phi(i,j),dx,dx,dz,nx,nz)
		  	g_delta = deltaofpsi(psi3x3,i,j,  &
							rho(i,j),b_phi(i,j),1.d0,nx,nz)
			g_wofpsi = wofpsi(psi(i,j))
			term1=g_wofpsi
!			g_spar = sparofpsi(psi(i,j))
!			g_sperp = sperpofpsi(psi(i,j))
			g_mtheta=mach_theta(psi(i,j))
			g_theta = thetaofpsi(psi(i,j))
			g_tpar = tparofpsi(psi(i,j))
			g_D = dofpsi(psi(i,j))

		  endif

		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

		  if (eq_type==3) then
			g_b_torc = b_phi(i,j)
		  elseif (eq_type == 1) then
			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))
		  endif

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho

		  if (eq_type==1) then

				rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  elseif (eq_type==3) then

 		  		rhomax = g_D*g_bfield/dabs(g_bfield-g_theta*g_tpar) &
							*dexp((2.d0*g_H + (g_r*g_omega)**2)/(2.d0*g_tpar))

		  endif

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(rhomin,rhomax,1000)
!!!					pause
!!!					rho(i,j) = dofpsi(psi(i,j))
!!!					cycle

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)


!!$					   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))		
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi
!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )

       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

!!$			   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

  end subroutine update_rho_frozen_March_17_2021




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_one_sided(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic

	real (kind=dkind), dimension (1:3,1:3) :: psi3x3
	integer iii,jjj,k

    ! -----------------------------------------------------


	mtm_limit = 1.d0

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

          !  "Setup Global Parameters"
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))


		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x

          if((Broot==0).and.(abs(psi(i,j)-psi_degen)<psic/10.).and.(psi(i,j)/=psi_degen)  &
			.and.(abs(min_ix-i)>3).and.(abs(min_iz-j)>3)) then
          ! check for one-sided derivatives

			if((psi(i+1,j)-psi_degen)*(psi(i-1,j)-psi_degen) > 0.d0) then

				! no crossing in x, use centered differences
				g_dpsidx = (psi(i+1,j)-psi(i-1,j))/dx/2.d0
				g_dpsidx_R = g_dpsidx
				g_dpsidx_L = g_dpsidx

			elseif((psi(i+1,j)-psi_degen)*(psi(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				g_dpsidx = (psi(i+1,j)-psi(i,j))/dx

			elseif((psi(i,j)-psi_degen)*(psi(i-1,j)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				g_dpsidx = (psi(i,j)-psi(i-1,j))/dx

			else

				print*, 'error in update_rho_one_sided, x derivative,', i, j

			endif

			if((psi(i,j+1)-psi_degen)*(psi(i,j-1)-psi_degen) > 0.d0) then

				! no crossing in z, use centered differences
				g_dpsidz = (psi(i,j+1)-psi(i,j-1))/dz/2.d0

			elseif((psi(i,j+1)-psi_degen)*(psi(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				g_dpsidz = (psi(i,j+1)-psi(i,j))/dz

			elseif((psi(i,j)-psi_degen)*(psi(i,j-1)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				g_dpsidz = (psi(i,j)-psi(i,j-1))/dz

			else

				print*, 'error in update_rho_one_sided, z derivative,', i, j

			endif

		else


          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

         endif

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho


			rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
				     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

!!!			       print *, "Failed to find root to Bernoulli eqn."
!!!                print *, "Coordinate (",i,",",j,"), repeat =",repeat
!!!					call b_plot(rhomin,rhomax,1000)
!!!					pause
!!!					rho(i,j) = dofpsi(psi(i,j))
!!!					cycle

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)


!!$					   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
					   ! not efficient, but there is a problem with the reversal...

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))		
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi


!!!!$          if((Broot==0).and.(abs(psi(min_ix,min_iz)-psi_degen)<psic/10.).and.(psi(min_ix,min_iz)/=psi_degen)) then
!!!!$          if((Broot==0).and.(abs(psi(min_ix,j)-psi_degen)<psic/10.).and.(psi(min_ix,j)/=psi_degen)  &
!!!!$			.and.(abs(min_ix-i)>3).and.(abs(min_iz-j)>3)) then
!!!!$          ! check for one-sided derivatives
!!!!$
!!!!$			if((psi(min_ix+1,min_iz)-psi_degen)*(psi(min_ix-1,min_iz)-psi_degen) > 0.d0) then
!!!!$
!!!!$				! no crossing in x, use centered differences
!!!!$				g_dpsidx = (psi(min_ix+1,min_iz)-psi(min_ix-1,min_iz))/dx/2.d0
!!!!$
!!!!$			elseif((psi(min_ix+1,min_iz)-psi_degen)*(psi(min_ix,min_iz)-psi_degen) >= 0.d0) then
!!!!$
!!!!$				! use one-sided right derivative
!!!!$				g_dpsidx = (psi(min_ix+1,min_iz)-psi(min_ix,min_iz))/dx
!!!!$
!!!!$			elseif((psi(min_ix,j)-psi_degen)*(psi(min_ix-1,min_iz)-psi_degen) >= 0.d0) then
!!!!$
!!!!$				! use one-sided left derivative
!!!!$				g_dpsidx = (psi(min_ix,min_iz)-psi(min_ix-1,min_iz))/dx
!!!!$
!!!!$			else
!!!!$
!!!!$				print*, 'error in update_rho_one_sided, minimum search, x derivative,', min_ix, min_iz
!!!!$
!!!!$			endif
!!!!$
!!!!$			if((psi(min_ix,min_iz+1)-psi_degen)*(psi(min_ix,min_iz-1)-psi_degen) > 0.d0) then
!!!!$
!!!!$				! no crossing in z, use centered differences
!!!!$				g_dpsidz = (psi(min_ix,min_iz+1)-psi(min_ix,min_iz-1))/dz/2.d0
!!!!$
!!!!$			elseif((psi(min_ix,min_iz+1)-psi_degen)*(psi(min_ix,min_iz)-psi_degen) >= 0.d0) then
!!!!$
!!!!$				! use one-sided right derivative
!!!!$				g_dpsidz = (psi(min_ix,min_iz+1)-psi(min_ix,min_iz))/dz
!!!!$
!!!!$			elseif((psi(min_ix,min_iz)-psi_degen)*(psi(min_ix,min_iz-1)-psi_degen) >= 0.d0) then
!!!!$
!!!!$				! use one-sided left derivative
!!!!$				g_dpsidz = (psi(min_ix,min_iz)-psi(min_ix,min_iz-1))/dz
!!!!$
!!!!$			else
!!!!$
!!!!$				print*, 'error in update_rho_one_sided, minimum search, z derivative,', min_ix, min_iz
!!!!$
!!!!$			endif
!!!!$
!!!!$		else

!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )

!!!!$        endif


       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

!!$			   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
			   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


!!$		   if(Broot==5) rhomin = 1.d-10 * dofpsi(0.d0)
		   ! not efficient, but there is a problem with the reversal...

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

	continue

    ! -----------------------------------------------------

  end subroutine update_rho_one_sided


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_gauss(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic
	real (kind=skind), dimension(nx,nz) :: out

	integer iii,jjj,k

    ! -----------------------------------------------------

!	print*, 'in Bernoulli Gauss'

	mtm_limit = 1.d0

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j,0)>0)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo



!  if(allocated(fmax_2D)) then
!  	if(size(fmax_2D,1)==nx) then
!  		continue
!  	else
!	  	deallocate(fmax_2D)
!  	endif
!  endif

  if(allocated(fmax_2D)) then
	continue
  else
  	allocate(fmax_2D(1:nx,1:nz))
  	fmax_2D = 0.d0
  endif



!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

          !  "Setup Global Parameters"
		  psi_Bern = psi(i,j)
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))


		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho


			rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
				     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif

!			print*, i,j,rhomin, rhomax
!			print*, '      '

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

		! find Bernmax (maximum of Bernoulli function between the roots)
        Bernmax = rtbis(dbern_drho,light,heavy,1.0d-10)
		fBernmax = bernoulli(Bernmax)
		fmax_2D(i,j) = fBernmax

		! now solve Bernoulli again, using the modified Bernoulli
!          light_gauss = rtsafe(newt_rap_gauss,xb1(1),xb2(1),1.0d-14,10000)
!          heavy_gauss = rtsafe(newt_rap_gauss,xb1(2),xb2(2),1.0d-14,10000)
          light = rtsafe(newt_rap_gauss,rhomin,Bernmax,1.0d-14,10000)
          light_gauss = light
          heavy = rtsafe(newt_rap_gauss,Bernmax,rhomax,1.0d-14,10000)
          heavy_gauss = heavy


          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy_gauss
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light_gauss
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy_gauss
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light_gauss
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

!			print*, i,j, rho(i,j), heavy_gauss, light_gauss

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))		
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi

!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )


       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if


!	out(1:nx,1:nz) = fmax_2D(1:nx,1:nz)
!	call radial_plot(out,psi,nx,nz,"fmax",nz/2)

!	pause


	continue

    ! -----------------------------------------------------

  end subroutine update_rho_gauss



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_rho_gauss_first(psi_in,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho,  &
						min_ix,min_iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi_in,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
	real (kind=dkind), dimension(1:nx,1:nz) :: psi
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    real (kind=dkind) :: light_gauss, heavy_gauss ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho
    integer, intent(inout) :: min_ix, min_iz
    ! We use rho_ext to store rho for which Bernoulli function is a maximum
    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
!	external rtsafe, rtbis
    ! VERY IMPORTANT
    ! psi_degen records the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: dmtm, mtm_soln
    ! Maximum iteration loop for Bisection search
    integer, parameter :: mtm_jmax = 40
    real (kind=dkind) :: mtm_limit 
    real (kind=dkind) :: phic

	integer iii,jjj,k

    ! -----------------------------------------------------

!	print*, 'in Bernoulli'

	mtm_limit = 1.d0

	delta_Bern = 1.d0*(nx-min_ix)

	do j=1,nz
		do i=1,nx
				if((bc_type/=7).and.(bc_type/=8)) then
					psi(i,j) = dabs(psi_in(i,j))
				else
					psi(i,j) = psi_in(i,j)
				endif
		enddo
	enddo

!!$   "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue

    min_drho = -1.0d0 ! Initialize to non-sensical values
    min_ix = 1
    min_iz = 1

    ! Now we want to update the density
    do j=1, nz
       do i=1, nx
          ! Only solve the inner region problem

		  if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

          !  "Setup Global Parameters"
		  psi_Bern = psi(i,j)
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))

			g_D = dofpsi(psi(i,j))		!

			g_Lambda = grav_potential(x_coord(i),z_coord(j))


		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             !  "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
          end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             !  "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             !  "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

          end if

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho


			rhomax = ( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
				     (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  if (rhomax<=0.d0) then
				write(*,*) 'error in rhomax = ',rhomax
				write(*,*) 'i,j = ', i,j
				pause
				stop
		  endif

          if((phic == 0.0d0).AND.(eq_type==1)) then
             ! There is only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  elseif(eq_type==3) then
             ! No poloidal flow allowed:
			 ! there is ALWAYS only 1 root given by rhomax
			 rho(i,j) = rhomax
			 continue
             cycle
		  end if

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
          rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
		   if(eq_type==1) then
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r	
				if(rhomin > 1.1d0*g_Phi**2) then
				! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				if(rhomin>rhomax) then
					rhomin = 1.1d0*g_phi**2
				endif
				rhomin = dmax1(rhomin,1.0d-31)
		   else
				rhomin = 1.d-32
		   endif

			if(rhomax<rhomin) then
				continue
!				rhomin = rhomin/10.
!				rhomax = rhomax*10.
			endif


!			print*, i,j,rhomin, rhomax
!			print*, '      '


          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,12,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots! in ",i,j
!             call b_plot(rhomin,rhomax,1000)
!			 pause
!!$			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                min_ix = i
                min_iz = j
                cycle
             else if(tmp > 0.0d0) then
                ! "We found the two separate roots"
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

				if(Broot==2) then

					! The light root
					xb1(1) = rhomin
					xb2(1) = rho_ext
					! The heavy root
					xb1(2) = rho_ext
					xb2(2) = rhomax

					if(nx>=65) then

						print*, 'problem in Bernoulli solver: '
						print*, 'no solution found with prescribed Mach theta!'
						print*, i, j
						print*, '--------------------------'

					endif

					cycle

				endif

                ! We need to reduce Mach Theta Max
                ! Store the current value of mach_theta_max
                mtm_soln = mach_theta_max
                ! Choose an initial increment
                dmtm = mtm_acc
	!			write(*,*) mach_theta_max
                do
                   ! Reduce Mach Theta Max
                   mach_theta_max = mtm_soln - dmtm
				   if (mach_theta_max<0.) then
						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
						pause
						stop
				   endif
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))

				   ! Calculate bounds for rho

						rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
							    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))


                   ! Calculate a minimum density
						rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r 
						if(rhomin > 1.1d0*g_Phi**2) then
						! "Only a sub-Alfvenic Root!"
						else
							rhomin = 1.1d0*g_Phi**2
						end if
						rhomin = dmax1(rhomin,1.0d-31)

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
				   ! Increase the separation
				   if(dmtm>mtm_soln/20.d0) then
						dmtm = dmtm*1.05d0
				   else
						dmtm = 10.0d0*dmtm 
				   endif
                 end do
                goto 101 ! repeat the density calculation
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
          ! Choose the last bracket (highest density)
          heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)

		! find Bernmax (maximum of Bernoulli function between the roots)
        Bernmax = rtbis(dbern_drho,light,heavy,1.0d-10)
		fBernmax = bernoulli(Bernmax)
!		print*, i,j , Bernmax, fBernmax, delta_Bern

		! now solve Bernoulli again, using the modified Bernoulli
!          light_gauss = rtsafe(newt_rap_gauss,xb1(1),xb2(1),1.0d-14,10000)
!          heavy_gauss = rtsafe(newt_rap_gauss,xb1(2),xb2(2),1.0d-14,10000)
          light = rtsafe(newt_rap_gauss,rhomin,Bernmax,1.0d-14,10000)
          light_gauss = light
          heavy = rtsafe(newt_rap_gauss,Bernmax,rhomax,1.0d-14,10000)
          heavy_gauss = heavy


          if((Broot == 0).or.(Broot == 5)) then
             ! Put a low density (rapidly moving) region 
             ! around the high density (slowly moving) core 
             if(psi(i,j) >= psi_degen) then
                ! Choose the last bracket (highest density)
                rho(i,j) = heavy_gauss
             else
                ! Choose the first bracket (lowest density)
                rho(i,j) = light_gauss
             end if
          else if((Broot == 1).or.(Broot == 4)) then
             ! Choose the heavy root
             rho(i,j) = heavy_gauss
          else if(Broot == 2) then
             ! choose the light root
             rho(i,j) = light_gauss
          else
             print *, "I don't understand Broot =",Broot
			 pause
             stop
          endif

          ! Find the minimum separation between the Bernoulli roots
          if(min_drho < 0.0d0) then
             ! initialization
             mean_rho = 0.5d0*(heavy + light)
             min_drho = (heavy - light)/mean_rho
             min_ix = i
             min_iz = j
          else
             tmp = 0.5d0*(heavy + light)
             if((heavy - light)/tmp < min_drho) then
                mean_rho = tmp
                min_drho = (heavy - light)/mean_rho
                min_ix = i
                min_iz = j
             end if
          end if

!		print*, i, j, light_gauss, heavy_gauss
!		print*, i, j, rho(i,j)

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------

    ! Update psi degenerate
    psi_degen = psi(min_ix,min_iz)

    ! -----------------------------------------------------

    if (((Broot == 0).or.(Broot == 5)) .and. seek_mtm == 1 .and. repeat == 1) then
!!$          print *, "Seek Mach Theta Max: min drho/rho =",min_drho
!!$  Search for Mach Theta Max which causes the maximum of the
!!$  Bernoulli function to be > 0, but less than some epsilon,
!!$  say 1.0d-8 or so.  This search will be handled by bisection.

!!$  The first step is to bracket the root.  The current value of
!!$  mach_theta_max will serve as the lower bound.  We now need to
!!$  find an upper bound for which the maximum of the Bernoulli
!!$  function is less than zero.

       ! Store the current value of mach_theta_max as our 
       ! best guess for the solution.  This is used below!!
       mtm_soln = mach_theta_max

!!$  Set the globals for the Bernouli eqn. which do not depend on
!!$  the value of Mach Theta.
	   g_r = x_coord(min_ix)
       ! g_Phi = phiofpsi(psi(min_ix,min_iz))
	   g_S = sofpsi(psi(min_ix,min_iz))		
	   g_I = iofpsi(psi(min_ix,min_iz))
	   g_D = dofpsi(psi(min_ix,min_iz))
       ! Calculate the derivatives of psi

!       g_dpsidx = 0.5d0*(psi(min_ix+1,min_iz) - psi(min_ix-1,min_iz))/dx_a(min_ix)
		 g_dpsidx = ( dx_a(min_ix-1)**2*psi(min_ix+1,min_iz) +  &
						(dx_a(min_ix)**2-dx_a(min_ix-1)**2)*psi(min_ix,min_iz) -  &
						dx_a(min_ix)**2*psi(min_ix-1,min_iz) ) /  &
						( dx_a(min_ix)*dx_a(min_ix-1)*(dx_a(min_ix)+dx_a(min_ix-1)) )

!       g_dpsidz = 0.5d0*(psi(min_ix,min_iz+1) - psi(min_ix,min_iz-1))/dz_a(min_iz)
		g_dpsidz = ( dz_a(min_iz-1)**2*psi(min_ix,min_iz+1) +  &
						(dz_a(min_iz)**2-dz_a(min_iz-1)**2)*psi(min_ix,min_iz) -  &
						dz_a(min_iz)**2*psi(min_ix,min_iz-1) ) /  &
						( dz_a(min_iz)*dz_a(min_iz-1)*(dz_a(min_iz)+dz_a(min_iz-1)) )


       g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)

	   g_Lambda = grav_potential(x_coord(i),z_coord(j))

       ! Adjust mach_theta_max to push Bernoulli max below zero.
       dmtm = mtm_acc
       do
          ! Increase mach_theta_max
          mach_theta_max = mtm_soln + dmtm
          ! Limit mach theta max
          mach_theta_max = dmin1(mach_theta_max,mtm_limit)

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))


		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

		  ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
					 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)


          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp < 0.0d0) exit ! upper bound on mach_theta_max found

          if(mach_theta_max == mtm_limit) then
             ! we can't go any higher, just repeat the density calc.
             print *, "Hit the mach_theta_max limit"
             repeat = 0
             goto 101
          end if

          ! Otherwise we try again
		  ! Increase the separation
		  if(dmtm > 0.1d0*mtm_soln) then
				dmtm = dmtm*1.05d0
		  else
				dmtm = 10.0d0*dmtm 
		  endif 
       end do

       ! Okay, we found a good upper bound for mach_theta_max
       ! Set dmtm, the range of mach theta max under consideration
       dmtm = mach_theta_max - mtm_soln

       ! Now we enter the bisection phase of this search
       do j = 1, mtm_jmax
          dmtm = 0.5d0*dmtm
          mach_theta_max = mtm_soln + dmtm

          ! Set the globals which depend on Mach Theta
          g_Phi = phiofpsi(psi(min_ix,min_iz))
          g_Omega = omegaofpsi(psi(min_ix,min_iz))
          g_H = hofpsi(psi(min_ix,min_iz))
		  g_Lambda = grav_potential(x_coord(min_ix),z_coord(min_iz))

		  ! Calculate bounds for rho

				rhomax = 1.01d0*( (g_H + 0.5d0*(g_r*g_Omega)**2 - g_Lambda)* &
					    (gamma-1.0d0)/(gamma*g_S) )**(1.0d0/(gamma-1.0d0))

          ! Calculate a minimum density
				rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
							 (2.0d0*g_H + (g_r*g_Omega)**2) - 2.d0*g_Lambda)/g_r
				if(rhomin > 1.1d0*g_Phi**2) then
					! "Only a sub-Alfvenic Root!"
				else
					rhomin = 1.1d0*g_Phi**2
				end if
				rhomin = dmax1(rhomin,1.0d-31)

          ! Now locate the maximum of the Bernoulli function.
          rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
          ! Evaluate the Bernoulli function 
          tmp = bernoulli(rho_ext)

          if(tmp .ge. 0.0d0) mtm_soln = mach_theta_max
          ! Compare dmtm to the desired accuracy and the solution
          if(dabs(dmtm) .lt. mtm_acc .or. tmp .eq. 0.0d0) exit
       end do
       if (j > mtm_jmax + 1) then
          print *, '[ngs_solve]: too many bisections for Mach Theta Max'
       end if
       ! Set mach_theta_max and repeat the density calculation once
       mach_theta_max = mtm_soln
       ! print *, "NB: Mach Theta Max =",mach_theta_max
       repeat = 0
       goto 101
    end if

!	print*, 'done with Bernoulli'

	continue

    ! -----------------------------------------------------

  end subroutine update_rho_gauss_first



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine update_rho_super(psi,rho,b_phi,nx,nz,seek_mtm,mtm_acc,min_drho)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    integer, intent(in) :: nx,nz,seek_mtm
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi,b_phi
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: rho
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), intent(in) :: mtm_acc
    integer :: i,j
    real (kind=dkind) :: light, heavy ! Temp. for Bernoulli roots
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: mean_rho, tmp
    real (kind=dkind), intent(inout) :: min_drho

    real (kind=dkind) :: rho_ext ! rho extremum
    integer :: repeat ! Logical for repeating the density calc.
    real (kind=dkind) :: ex,ez
    integer :: nb
    real (kind=dkind), dimension(1:10) :: xb1,xb2
    real (kind=dkind) :: rhomax,rhomin
	real (kind=dkind) :: rhomax1,rhomax2,rhomax3
!    real (kind=dkind) :: rtsafe, rtbis
	real (kind=dkind) term1
	real(kind=dkind) :: plocg, Alf, disc
!	external rtsafe, rtbis
!!$    ! VERY IMPORTANT
!!$    ! psi_degen records the value of psi where the solution to the
!!$    ! Bernoulli equation is approximately degenerate.  It is 
!!$    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
!!$    ! The next set of variables are used for the bisection search for
!!$    ! locating the degenerate root for Mach Theta Max
!!$    real (kind=dkind) :: dmtm, mtm_soln
!!$    ! Maximum iteration loop for Bisection search
!!$    integer, parameter :: mtm_jmax = 40
!!$    real (kind=dkind), parameter :: mtm_limit = 1.0d0
    real (kind=dkind) :: phic

	real (kind=dkind), dimension (1:3,1:3) :: psi3x3
	integer iii,jjj,k
	real (kind=dkind) :: bigR,x

    ! -----------------------------------------------------


!!$    print *, "Calculate The Density"

    ! -----------------------------------------------------
    repeat = 1

    ! This next line allows us to repeat the density calculation
    ! using a goto in case of failure of the Bernoulli solver
101 continue


    ! Now we want to update the density
    do j=1, nz
		outside = .false.
		inside = .false.
       do i=1, nx
          ! Only solve the inner elipse problem
!          ex = ((i-1)*dx - 0.5d0*x_size)/a_elps
!          ez = ((j-1)*dz - 0.5d0*z_size)/b_elps
!          if((ex*ex + ez*ez) > 1.0d0) then
!		  if((i==9).and.(j==9)) pause

                x = x_coord(i)

			if(sort_grid(i,j,0)<=0) then
				 cycle
          end if

		  if((psi(i,j)/psic)<fraction) then
!		  if(dabs(psi(i,j)/psic)<fraction) then
			rho(i,j) = 1.d-10	!	dofpsi(0.0d0)
			if(inside) outside=.true.
			cycle
		  endif

!		if(x>rmajor+a_elps/3.d0) then
!		  	rho(i,j) = dofpsi(0.0d0)
!			cycle
!		endif !trucco

		  inside = .true.

!		  if((inside).AND.(outside)) then
!				rho(i,j) = dofpsi(0.0d0)	
!				cycle
!		  endif

!		  bigR = rmajor - 0.5d0*x_size + (i-1)*dx

!		  if (bigR>5.8d0) cycle
          ! print *, "Passed Rho Elipse test"



          ! print *, "Setup Global Parameters"
          phic = phiofpsi(psi(i,j))
          ! Set the globals for solving Bernouli eqn.
!          g_r = rmajor - 0.5d0*x_size + (i-1)*dx
		  g_r = x_coord(i)
          g_Phi = phic
          ! print *, "g_S"
          g_S = sofpsi(psi(i,j))
          ! print *, "g_I"
          g_I = iofpsi(psi(i,j))
          ! print *, "g_Omega"
          g_Omega = omegaofpsi(psi(i,j))
          ! print *, "g_H"
          g_H = hofpsi(psi(i,j))


		  g_dx=dx_a(i)
		  g_dz=dz_a(j)
		  g_indi=i
		  g_indj=j
		  g_nx=nx
		  g_nz=nz


          ! Calculate the derivatives of psi
          ! g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
          ! g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
          ! d Psi / d x
          if(i == 1) then
             ! use a 1-sided difference
             ! print *, "One sided x-dif left"
             g_dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(i)
          else if(i == nx) then
             ! use a 1-sided difference
             ! print *, "One sided x-dif right"
             g_dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(i)
          else
             ! use a centered difference
!             g_dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx_a(i)
 			 g_dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
							(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
							dx_a(i)**2*psi(i-1,j) ) /  &
							( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
         end if
          ! d Psi / d z
          if(j == 1) then
             ! use a 1-sided difference
             ! print *, "One sided z-dif left"
             g_dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
          else if(j == nz) then
             ! use a 1-sided difference
             ! print *, "One sided z-dif right"
             g_dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
          else
             ! use a centered difference
!             g_dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz_a(j)
			g_dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
							(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
							dz_a(j)**2*psi(i,j-1) ) /  &
							( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )
          end if

			g_b_torc = dsqrt(mu_mag)*(g_i/g_r + g_r*g_phi*g_Omega)/ &
                     (1.0d0 - g_phi*g_phi/rho(i,j))

		  g_b_polc =  (dsqrt(g_dpsidx**2+g_dpsidz**2)/g_r)
		  g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)


          ! Calculate bounds for rho
          ! print *, "Calculate density bounds"


			rhomax = g_phi**2 * (1.d0-1.d-9)

          ! Otherwise, increase rhomax a little mach_theta_max == mtm_limit
 !         rhomax = 1.01d0*rhomax

          ! Calculate a minimum density
!!$          rhomin = dsqrt((g_dpsidx**2 + g_dpsidz**2)*g_phi**2/ &
!!$               (2.0d0*g_H + (g_r*g_Omega)**2))/g_r
!!$          if(rhomin > 1.1d0*g_Phi**2) then
!!$             ! print *, "Only a sub-Alfvenic Root!"
!!$          else
!!$             rhomin = 1.1d0*g_Phi**2
!!$          end if
!!$          rhomin = dmax1(rhomin,1.0d-31)
			rhomin = 1.0d-31		! bookmark
!		  call b_plot(rhomin,rhomax,1000)
 !         call p_plot(rhomin,rhomax,1000)

          ! Now bracket the roots in the simple case
          nb = 10 ! indicate the space available
          call zbrak(bernoulli,rhomin,rhomax,120,xb1,xb2,nb)

          ! Check for a crazy error...
          if(nb == 1 .or. nb > 2) then
             print *, "Found ",nb,"roots!"
!             call b_plot(rhomin,rhomax,5000)
!			 pause
!			 stop
          end if

          ! If that root bracketing failed try a different approach.
          ! The Bernoulli function is sufficiently smooth that we can
          ! locate a maximum (the Bernoulli function is concave down)
          ! by looking for zero or root of the first derivative of the
          ! Bernoulli function with respect to rho.  We use bisection
          ! to locate the zero of d(Bernoulli)/d(rho).
       !   if(nb .eq. 0) then
		  if(nb .eq. 0) then
!		  if((nb .eq. 0).or.(nb==1)) then
             rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
             ! Next we evaluate the Bernoulli function 
             tmp = bernoulli(rho_ext)
             ! Now if we have a degenerate root (tmp = 0.0d0)
             if(tmp .eq. 0.0d0) then		! what about tmp<tiny?
			 ! "=" sign doesn't make numerical sense...
                ! print *, "We found 2 degenerate roots!"
                ! we have 2 degenerate roots
                rho = rho_ext
                mean_rho = rho_ext
                min_drho = 0.0d0
                cycle
             else if(tmp > 0.0d0) then
                ! print *, "We found the two separate roots"
                ! we have 2 roots
                ! Set xb1(1,2), xb2(1,2) and nb = 2
                nb = 2
                ! The light root
                xb1(1) = rhomin
                xb2(1) = rho_ext
                ! The heavy root
                xb1(2) = rho_ext
                xb2(2) = rhomax
             else
                ! mach theta max is too large, there are no roots.
                if(repeat == 0) then
                   ! This should never happen!
!!a                   print *, "Failed to find root to Bernoulli eqn."
!!a                   print *, "Coordinate (",i,",",j,"), repeat =",repeat
                endif

	!			for Alfvenic equilibrium
				if (rho_ext>(dofpsi(psic)+dofpsi(0.d0))/2) then
					rho(i,j) = rho_ext
				else
					rho(i,j) = 55.d-2	!dofpsi(psi(i,j))*rmajor/x	!		
				endif
				print *, 'warning: Bernoulii solver has failed in',i,j
!!$	            call b_plot(rhomin,rhomax,100)
!!$				pause
!!$				g_phi = phic*10.d0
!!$	            call b_plot(rhomin,rhomax,100)
!!$				pause
!!$				g_phi = phic/10.d0
!!$	            call b_plot(rhomin,rhomax,100)
!!$				pause

				cycle


                do
                   ! Reduce Mach Theta Max
 !                  mach_theta_max = mtm_soln - dmtm
!				   if (mach_theta_max<0.) then
!						write(*,*) 'error: Mach_theta_max < 0'
!						call b_plot(rhomin,rhomax,1000)
!						pause
	!					stop
!				   endif
!!$                      print *, "FB: Mach Theta Max =",mach_theta_max
!					if((i==12).and.(j==6)) call b_plot(rhomin,rhomax,1000)
                   ! Set the globals which depend on Mach Theta
                   g_Phi = phiofpsi(psi(i,j))
                   g_Omega = omegaofpsi(psi(i,j))
                   g_H = hofpsi(psi(i,j))
!				   g_D = dofpsi(psi(i,j))
		   		   if (eq_type==2) then
						if ((g_phi==0.d0).OR.(rho(i,j)==0.d0)) then
							g_b_torc=g_i/g_r
						else
 							g_b_torc = (g_i/g_r+g_r*g_Omega*g_phi)/ &
										(1.d0 - (g_phi**2)/rho(i,j) -g_delta)
						endif
			!			write(*,*) g_b_torc,g_b_polc
			!			write(*,*) g_phi,rho(i,j)
						g_bfield = dsqrt(g_b_torc**2+g_b_polc**2)
				   endif

				   ! Calculate bounds for rho

					rhomax=	(1.d0-1d-9)*g_Phi**2


                   ! Calculate a minimum density

					rhomin = 1.0d-31		! bookmark

                   ! Now bracket the roots:
                   ! Again we accomplish this by first locating the
                   ! maximum of the Bernoulli function.
                   rho_ext = rtbis(dbern_drho,rhomin,rhomax,1.0d-10)
                   ! Next we evaluate the Bernoulli function 
                   tmp = bernoulli(rho_ext)
                   if(tmp > 0.0d0) then
                      ! print *, "We found the two separate roots"
                      ! we have 2 roots
                      ! Set xb1(1,2), xb2(1,2) and nb = 2
                      nb = 2
                      ! The light root
                      xb1(1) = rhomin
                      xb2(1) = rho_ext
                      ! The heavy root
                      xb1(2) = rho_ext
                      xb2(2) = rhomax
                      ! Exit the do loop and repeat the density calc.
                      exit
                   end if
                   ! Otherwise we try again
!				   call b_plot(rhomin,rhomax,1000)
!                   dmtm = 10.0d0*dmtm ! Increase the separation
                end do
                goto 101 ! repeat the density calculation
!!$
!!$                if(nx > 33) then
!!$                   call b_plot(rhomin,rhomax,1000)
!!$                else 
!!$                   goto 101 ! repeat the density calculation
!!$                end if
             end if
          end if

          ! Solve for the two roots
          ! Choose the first bracket (lowest density)
!		  print *, rhomin, rhomax, i,j

		  if(nb==1) then
	          light = 1.d-20
		      ! Choose the last bracket (highest density)
			  heavy = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
			else
	          light = rtsafe(newt_rap,xb1(1),xb2(1),1.0d-14,10000)
		      ! Choose the last bracket (highest density)
			  heavy = rtsafe(newt_rap,xb1(2),xb2(2),1.0d-14,10000)
		  endif

             ! Choose the heavy root
             rho(i,j) = heavy	!	light	! 

!		  print *, rho(i,j), i,j

		  Alf = g_phi**2/rho(i,j)	! Alfvén square
		  plocg = gamma*g_S*rho(i,j)**gamma
!!$		  disc = -(1.d0-Alf)**2* ( Alf*(g_bfield**2+plocg)-plocg)  &
		  disc = - ( Alf*(g_bfield**2+plocg)-plocg)  &
				/(Alf**2*g_b_polc**2 - Alf*(g_bfield**2+plocg)+plocg)

		  if(disc<0.d0) then

			print*, 'system is hyperbolic in',i,j,disc
			continue
!			pause

		  else

			continue

		  endif

       end do
    end do
    ! print *, "Finished: Rho Update"
    ! -----------------------------------------------------


	continue

    ! -----------------------------------------------------

  end subroutine update_rho_super
  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve_wrapper(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! depending on input, calls a different routine to solve the GS-Bernoulli system

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi,rho,residual,b_phi
    real (kind=dkind) :: eps
    real (kind=dkind) :: orp

	if(((jump_option==1).or.(jump_option==2).or.(jump_option==3).or.(jump_option==4)).and.(Broot==0)) then                                                                                                                                                                                                                                                                                                                                                                        
	     call ngs_solve_jump(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif((jump_option==0).or.(Broot>0)) then
	     call ngs_solve(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif(jump_option==-1) then
	     call ngs_solve_no_limiter(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif(jump_option==-2) then
	     call ngs_solve_grad_psi(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif(jump_option==-3) then
	     call ngs_solve_grad_psi_consistent(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif(jump_option==-4) then
		if(nx>inter_switch) then
		     call ngs_solve_grad_psi_consistent_smooth_transition(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
		 else
			 call ngs_solve(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
		 endif
	elseif(jump_option==-5) then
	     call ngs_solve_grad_psi_gauss(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	elseif(jump_option==-7) then
	     call ngs_solve_all_gauss(psi,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,orp)
	endif

	continue

end subroutine ngs_solve_wrapper


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
    real (kind=dkind), intent(in) :: eps
    ! Input Over Relaxation Parameter
    ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
    ! else we set orp = in_orp
    real (kind=dkind), intent(in) :: in_orp
    integer ::ipass,i,isw,j,jsw,k,p,h,alloc_stat
    real (kind=dkind) :: res, den, x
    ! res -> residual
    ! den -> denominator of the Newton-Gauss-Seidel update
    ! x -> Radial position of grid point (i,j)
    real (kind=dkind) :: dx2,dz2,mx
    real (kind=dkind) :: anorm, eps_anorm
	real (kind=dkind) :: anorm2, eps_anorm2	! these 2 are fo b_phi
    real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
    real (kind=dkind) :: orp, std_orp ! over relaxation parameter
    ! Phi @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: phirx,philx,phirz,philz
    ! Density @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    ! The square of the Alfvenic Mach Number
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
    real (kind=dkind) :: rhoc,phic,omegac,deltac,thetac,tparc,dc
    ! by is the phi component of the magnetic field
    ! b_dot_v is the inner product of B and V
    real (kind=dkind) :: by, b_dot_v
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: last_mtm ! last mach theta max
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: min_drho, tmp
    integer :: min_ix, min_iz
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: mtm_soln
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), parameter :: mtm_acc = 1.0d-12
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    real (kind=dkind) :: drho_ds ! Used with a slope limiter
	! The following is used to allow the possibility of anisotropic pressure
	real (kind=dkind) :: deltarx,deltalx,deltarz,deltalz
	real (kind=dkind) :: bphirx,bphilx,bphirz,bphilz
	real (kind=dkind), dimension (1:3,1:3) :: psi_delta,psi3x3
	real (kind=dkind) :: b_field_l,b_pol_l,dpsidx,dpsidz,psinew
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0
	real(kind=dkind) :: bpol_min_temp,bpol_max_temp, psi_pmax_temp, B_pmax_temp
	real (kind=dkind), dimension (-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: last_anorm
	real(kind=dkind) :: inorm = 0.d0

    real (kind=dkind), dimension(:,:,:), allocatable :: psi123
	integer :: h_steps=1
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: tic,toc
	integer :: pippa
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: orp_3step = 5.d-2
	integer :: i_zone = 1 !This will ony be changed if bc_type==7
!	logical, save :: initialize_zones = .true.

	integer iii,jjj

	if((tri_type==13).and.(nx>=inter_switch).and.(nz>=inter_switch)) then
!	if((tri_type==13).and.(nx>inter_switch).and.(nz>inter_switch)) then

		h_steps = 3

	endif

	if(allocated(psi123)) deallocate(psi123)
	allocate(psi123(1:h_steps,1:nx,1:nz),stat = alloc_stat)
	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi123"
		 pause
		 stop
	endif

	psi123 = 0.d0

	do h=1,h_steps
	do j=1,nz
	do i=1,nx

		psi123(h,i,j) = psi_in(i,j)

	enddo
	enddo
	enddo

	if((eq_type==3).and.(tri_type==11)) then
		bpol0_temp = bpol0_fix
		if(nx>n_min) then
			if(ana_fac==0.d0) then
				ana_fac = ana_fac0*1.d-3
			else
				ana_fac = ana_fac0
			endif
		endif
	endif

    if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
    last_mtm = mach_theta_max
    eps_anorm = 0.0d0
	eps_anorm2 = 0.0d0

	  if((tri_type==11).and.(((eq3_opt==7).or.(eq3_opt==8)).or.  &
			((p_opt==7).or.(p_opt==8)))) call get_psi_pres(psi123(1,:,(nz+1)/2),nx)

	fix_orp0 = fix_orp

!!$	 	call step_output(nx,psi,rho,residual)

    ! The under/over relaxation parameter
    ! Standard Case:
    if(in_orp <= 0.0d0) then
       ! Note that the parameters below are tuned for starting with 
       ! a 17 x 17 grid.
       orp = 1.0d0
       if(nx <= 5) then
          orp = 0.5d0
          rjac = 0.75d0
       else if(nx <= 9) then
          orp = 0.5d0
          rjac = 0.9d0
       else if(nx <= 17) then
          rjac = 0.98d0
       else if(nx <= 33) then
          rjac = 0.995d0
       else if(nx <= 65) then
          rjac = 0.9972d0
       else 
!!$ It would appear that 129 is about where the solution begins to converge
          rjac = 0.0d0
       end if
       std_orp = orp
    else
       print *, "Input Over Relaxation Parameter =",in_orp
       orp = in_orp
       std_orp = in_orp
    endif

   if (accelerate) then
		continue
   else
	    orp = fix_orp
   endif


    dx2 = dx*dx
    dz2 = dz*dz

!	if((initialize_zones).and.(bc_type==7)) then
	if(bc_type==7) then

		call update_sort_grid(psi123(1,:,:),nx,nz,inorm)
!		initialize_zones = .false.

	endif

    if(Broot/=3) then

	! Update rho before "relaxing psi" but do not seek mach theta max
	    call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
							min_ix,min_iz)

	elseif(Broot==3) then

	    call update_rho_super(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho)

	endif

! set up a few functions for Broot = 3

	if(Broot==3) then

		psic_flag = 1.d6

		dc_loc = dofpsi(psic)
		dcp_loc = dddpsi(psic)
		pc_loc = pofpsi(psic)
		pcp_loc = dpdpsi(psic)
		b0c_loc = bzero(psic)
		b0cp_loc = dbzerodpsi(psic)

		psic_flag = psic

	endif


    ! Iterate until convergence
    do k=1, max_it
!!$     "Update Psi"

!		call cpu_time(tic)

       ! Reset the norm and max of Psi
       mx = 0.0d0
       anorm = 0.0d0
	   anorm2 = 0.0d0

	   do j=2,nz-1
	   do i=2,nx-1
			if ( ((tri_type==13).and.(sort_grid(i,j,0)==2)).or.  &
				((tri_type/=13).and.(bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j,0)==1)) ) then
					psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
!!$				March 12 2021: This should NOT be necessary for bc_type==7
!!$ 			elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
!!$ 				psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			endif
	   enddo
	   enddo

	   bpol_max_temp = 0.d0
	   bpol_min_temp = 1.d22
	   psi_pmax_temp = 0.d0

		do h=1,h_steps
		! vertical stabilization cycle

       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
             do i=isw+1, nx-1, 2
                ! "Update Psi(",i,",",j,")"
                ! Only solve the inner region problem

				if(tri_type==10) then

!!$					if((ex*ex + ez*ez) >= radius_ext**2) then
!!$						cycle
!!$					endif

				else

					if(sort_grid(i,j,0)<0) then ! change back to <= March 17 2021
!					if((sort_grid(i,j,0)<0).or.((sort_grid(i,j,0)==0).and.(bc_type/=7))) then
					   cycle
					end if

				endif

				! set up local psi values
				! this should save considerable resources in looking up values

				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi123(h,i,j)
				if(i>1) psil = psi123(h,i-1,j)
				if(i<nx) psir = psi123(h,i+1,j)
				if(j>1) psid = psi123(h,i,j-1)
				if(j<nz) psiu = psi123(h,i,j+1)

				psi_around(-1,0) = psil
				psi_around(1,0) = psir
				psi_around(0,0) = psi0
				psi_around(0,-1) = psid
				psi_around(0,1) = psiu

				! Set up zone. May want to do the same for tri_type==13 here. March 12 2021

				if(bc_type==7) then
					! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
					i_zone = min(0,sort_grid(i,j,0)-2)
					i_zone = max(-1,i_zone)
				endif

				! set up functions of psi
				! NOTE: DIFFERENTIATE FOR EQ_TYPE=3 LATER ON

				psi_flag = 1.d9
				psi_flag_dep = 1.d9
				psi_flag_dep2 = 1.d9
				psi_flag_ham = 1.d9

				d_loc = dofpsi(psi0,i_zone)
				dp_loc = dddpsi(psi0,i_zone)

				p_loc = pofpsi(psi0,i_zone)
				pp_loc = dpdpsi(psi0,i_zone)

				psi_flag = psi0

				b0_loc = bzero(psi0,i_zone)
				b0p_loc = dbzerodpsi(psi0,i_zone)

				psi_flag_dep = psi0

				mth_loc = mach_theta(psi0,i_zone)
				mthp_loc = dmach_thetadpsi(psi0,i_zone)

				mph_loc = mach_phi(psi0,i_zone)
				mphp_loc = dmach_phidpsi(psi0,i_zone)

				psi_flag_dep2 = psi0

				s_loc = sofpsi(psi0,i_zone)
				sp_loc = dsdpsi(psi0,i_zone)

				phi_loc = phiofpsi(psi0,i_zone)
				phip_loc = dphidpsi(psi0,i_zone)

				omega_loc = omegaofpsi(psi0,i_zone)
				omegap_loc = domegadpsi(psi0,i_zone)

				i_loc = iofpsi(psi0,i_zone)
				ip_loc = didpsi(psi0,i_zone)

				h_loc = hofpsi(psi0,i_zone)
				hp_loc = dhdpsi(psi0,i_zone)

				psi_flag_ham = psi0

				! end of functions set up

                x = x_coord(i)

                ! Calculate rho, phi and omega at the current location
                rhoc = rho(i,j)
                phic = phiofpsi(psi0,i_zone)
                omegac = omegaofpsi(psi0,i_zone)

				if(eq_type==3) then 
					thetac = thetaofpsi(psi0)
					tparc = tparofpsi(psi0)
				endif

				dc = dofpsi(psi0,i_zone)
                ! Calculate B_phi = by
				if (eq_type == 1) then
					by = dsqrt(mu_mag)*(iofpsi(psi0,i_zone)/x + x*phic*omegac)/ &
						 (1.0d0 - phic*phic/rhoc)
				elseif (eq_type == 3) then
					deltac=deltaofpsi(psi123(h,i-1:i+1,j-1:j+1),i,j,rho(i,j),  &
									b_phi(i,j),1.d0,nx,nz)
					by = b_phi(i,j)
				endif

                ! -----------------------------------------------------

                ! Calculate Phi & Density, rho @ +/- 0.5*dx, +/- 0.5*dz
                ! and then the square of the Alfvenic Mach Number.
!!$  NOTE: The use of a slope limiter to interpolate rho is essential
!!$  for eliminating an oscillation in Bpol where the solution switches
!!$  from the sub-slow to super-slow root. -- 3/14/2002 -- T. Gardiner
                ! -----------------------------------------------------

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimx(0) = ( dx_a(i-1)**2*psir +  &
									(dx_a(i)**2-dx_a(i-1)**2)*psi0 -  &
									dx_a(i)**2*psil ) /  &
									( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )

					ma2c = phic*phic/rhoc

				endif

!                drho_ds = van_Leer_slope(rho(i-1,j),rho(i,j),rho(i+1,j),dx)
                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

				! March 12 2021
				! We are evaluating free functions in half grid points. Rather than testing for wich zone
				! each point is in, we will always assume that the are in the same zone as the central point. 

                ! Right x
                phirx = phiofpsi(0.5d0*(psir + psi0),i_zone)
!                rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                rhorx = rho(i,j) + 0.5d0*dx_a(i)*drho_ds
                ma2rx = phirx*phirx/rhorx
				if(eq_type==3) bphirx = (b_phi(i+1,j) + b_phi(i,j))*0.5d0
				psiprimx(1) = (psir - psi0)/dx_a(i)

                ! Left x
                philx = phiofpsi(0.5d0*(psil + psi0),i_zone)
!                rholx = rho(i,j) - 0.5d0*dx*drho_ds
                rholx = rho(i,j) - 0.5d0*dx_a(i-1)*drho_ds
                ma2lx = philx*philx/rholx
				if(eq_type==3) bphilx = (b_phi(i-1,j) + b_phi(i,j))*0.5d0
				psiprimx(-1) = (psi0 - psil)/dx_a(i-1)

                ! -----------------------------------------------------

!                drho_ds = van_Leer_slope(rho(i,j-1),rho(i,j),rho(i,j+1),dz)
                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimz(0) = ( dz_a(j-1)**2*psiu +  &
									(dz_a(j)**2-dz_a(j-1)**2)*psi0 -  &
									dz_a(j)**2*psid ) /  &
									( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

				endif

                ! Right z
                phirz = phiofpsi(0.5d0*(psi0 + psiu),i_zone)
!                rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                rhorz = rho(i,j) + 0.5d0*dz_a(j)*drho_ds
                ma2rz = phirz*phirz/rhorz
				if(eq_type==3) bphirz = (b_phi(i,j) + b_phi(i,j+1))*0.5d0
				psiprimz(1) = (psiu - psi0)/dz_a(j)


                ! Left z
                philz = phiofpsi(0.5d0*(psi0 + psid),i_zone)
!                rholz = rho(i,j) - 0.5d0*dz*drho_ds
                rholz = rho(i,j) - 0.5d0*dz_a(j)*drho_ds
                ma2lz = philz*philz/rholz
				if(eq_type==3) bphilz = (b_phi(i,j) + b_phi(i,j-1))*0.5d0
				psiprimz(-1) = (psi0 - psid)/dz_a(j-1)


				! -----------------------------------------------------


				! calculate the magnetic field

!				call cpu_time(tic)

!				do pippa=1,1000

!				if(eq_type==3) then

					call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
!					call psi_derivative(i,j,nx,nz,psi123(h,:,:),dpsidx,dpsidz)
					b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)
					b_field_l = dsqrt(by**2+b_pol_l**2)

!				endif
!				enddo

!				call cpu_time(toc)
!				print*, 'psi_derivative time', toc-tic
!				pause

				if(tri_type==11) then

					if (b_pol_l>bpol_max_temp) bpol_max_temp = b_pol_l
					if (b_pol_l<bpol_min_temp) then
						bpol_min_temp = b_pol_l
						i_bmin = i
						j_bmin = j
					endif

					if(((p_opt==9).or.(eq3_opt==9)).and.(j==1+nz/2)) then
						if(abs(psi_pmax-psi0)<abs(psi_pmax-psi_pmax_temp)) then
							B_pmax_temp = b_pol_l
							psi_pmax_temp = psi0
						endif
					endif

				endif

				! -----------------------------------------------------

				! Calculate Delta (function of p_par and p_perp)

				if(eq_type==3) then

				! Right x

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=5.d-1*(psi123(h,i-1,j+1)+psi0)
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=psi0
				psi_delta(2,2)=5.d-1*(psi0+psiu)
				psi_delta(2,3)=psiu
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=5.d-1*(psi0+psi123(h,i+1,j+1))
				psi_delta(3,3)=0.d0

				deltarx=deltaofpsi(psi_delta,i,j,rhorx,bphirx,0.5d0,nx,nz)

				! Left x

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=5.d-1*(psi123(h,i-1,j-1)+psi0)
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=psid
				psi_delta(2,2)=5.d-1*(psi0+psid)
				psi_delta(2,3)=psi0
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=5.d-1*(psi0+psi123(h,i+1,j-1))
				psi_delta(3,3)=0.d0

				deltalx=deltaofpsi(psi_delta,i,j,rholx,bphilx,0.5d0,nx,nz)

				! Right z

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=psi0
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=5.d-1*(psi0+psi123(h,i+1,j-1))
				psi_delta(2,2)=5.d-1*(psi0+psir)
				psi_delta(2,3)=5.d-1*(psi0+psi123(h,i+1,j+1))
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=psir
				psi_delta(3,3)=0.d0

				deltarz=deltaofpsi(psi_delta,i,j,rhorz,bphirz,0.5d0,nx,nz)

				! Left z

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=psil
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=5.d-1*(psi0+psi123(h,i-1,j-1))
				psi_delta(2,2)=5.d-1*(psi0+psil)
				psi_delta(2,3)=5.d-1*(psi0+psi123(h,i-1,j+1))
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=psi0
				psi_delta(3,3)=0.d0

				deltalz=deltaofpsi(psi_delta,i,j,rholz,bphilz,0.5d0,nx,nz)

				endif

                ! -----------------------------------------------------

                ! Calculate B dot v
!                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + &
!                     ( (psi(i+1,j) - psi(i-1,j))**2/dx2 &
!                     + (psi(i,j+1) - psi(i,j-1))**2/dz2 )*0.25d0/(x*x) )
                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)

                ! -----------------------------------------------------

				! March 12 2021: we don't need to distinguish the two zones this way for bc_type==7 anymore.
				! We are using the zone label "i_zone" instead.
				! January 21 2022: Adapted to the same option for tri_type==13
				if( ((tri_type/=13).and.(bc_type/=7).AND.((psi0/psic)<fraction))  &
								.OR.  &
!!!					((bc_type==7).AND.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))  &
!!!								.OR.  &
!					((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
					((tri_type==13).AND.(sort_grid(i,j,0)<1)) ) then
					! OUTER REGION

!!$					if(h_steps==1) then

						if(grid_type==0) then

							res = (1.d0/mu_mag)*(( (1.0d0 )*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 )*(psil-psi0)/(x - 0.5d0*dx) )/(x*dx2) &
								+ ( (1.0d0 )*(psiu-psi0) &
								  + (1.0d0 )*(psid-psi0) )/(x*x*dz2) ) 

						else

							res = (2.d0/mu_mag)*( (dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										- dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										(dz_a(j-1)**2*psiprimz(1)   &
										- dz_a(j)**2*psiprimz(-1) +  &
										(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif

!!$					else
!!$
!!$						res = 0.d0
!!$
!!$					endif






				else	!inner region

				    if (eq_type==1) then

!!$						if(h_steps==1) then

							if(grid_type==0) then

								term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
								  /(x*dx2) &
								+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
								  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

							else

								term0 = (2.d0/mu_mag)*( ((1.0d0 - ma2rx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
											-(1.0d0 - ma2lx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1))   &
											+(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x )   &
											/(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1))   &

											+((1.0d0 - ma2rz)*dz_a(j-1)**2*psiprimz(1)   &
											-(1.0d0 - ma2lz)*dz_a(j)**2*psiprimz(-1)   &
											+(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) )   &
											/(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

							endif

!!$						else
!!$
!!$							term0 = 0.d0
!!$
!!$						endif



						term1= b_dot_v*dphidpsi(psi0,i_zone)/dsqrt(mu_mag) 
						term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0,i_zone) 
						term3= by*didpsi(psi0,i_zone)/x/dsqrt(mu_mag) 
						term4= rhoc*dhdpsi(psi0,i_zone) 
						term5= rhoc**gamma*dsdpsi(psi0,i_zone)/(gamma-1.0d0)

						if(i_zone==-1) then
							continue
						endif

						res = term0+term1+term2+term3+term4-term5

					elseif (eq_type==3) then

!!$						if(h_steps==1) then

							term0 = (2.d0/mu_mag)*( ((1.0d0 - deltarx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))  &
										- (1.0d0 - deltalx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(1.0d0 - deltac)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										((1.0d0 - deltarz)*dz_a(j-1)**2*psiprimz(1)  &
										- (1.0d0 - deltalz)*dz_a(j)**2*psiprimz(-1) +  &
										(1.0d0 - deltac)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

!!$						else
!!$
!!$							term0 = 0.d0
!!$
!!$						endif

						! March 21 2021 Need to be careful with the anisotropic functions
						! (not sure they are defined with the zone selection, yet)
						res = term0  &

							+ rhoc*dtpardpsi(psi0)  &
							+ by*didpsi(psi0,i_zone)/x/dsqrt(mu_mag) &
							+ rhoc*x**2*omegac*domegadpsi(psi0,i_zone) &
							+ rhoc*dhdpsi(psi0) &
							- rhoc*(			&
									dtpardpsi(psi0)*dlog(rhoc/dc*dabs(b_field_l-thetac*tparc)/b_field_l) &
									- tparc*dddpsi(psi0)/dc  &
									- tparc  &
									  *(tparc*dthetadpsi(psi0)+thetac*dtpardpsi(psi0)) &
									  *sign(1.d0,(b_field_l-thetac*tparc)) )

						if(b_field_l-thetac*tparc<0.1d0*b_field_l) then

							print *, b_field_l-thetac*tparc, k, i, j
							reduce = .true.
							continue

						endif

					endif


				continue

				endif




				res = res*mu_mag


                ! Store the residual
                if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
                   residual(i,j) = res
                end if

				if(res>1.d0) then
					continue
				endif

                ! -----------------------------------------------------

                if (eq_type==1) then

!!$					den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
!!$				             (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
!!$					      -( (1.0d0 - ma2rz) + &
!!$						     (1.0d0 - ma2lz) )/(x*x*dz2)



					den = -2.d0*( ((1.0d0 - ma2rx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - ma2lx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - ma2rz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - ma2lz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

					continue


				elseif (eq_type==3) then

!!$					den = -( (1.0d0 - deltarx)/(x + 0.5d0*dx) + &
!!$							(1.0d0 - deltalx)/(x - 0.5d0*dx) )/(x*dx2) &
!!$							-( (1.0d0 - deltarz) + &
!!$							(1.0d0 - deltalz) )/(x*x*dz2)


					den = -2.d0*( ((1.0d0 - deltarx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - deltalx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - deltac)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - deltarz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - deltalz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - deltac)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

				endif

     !    				psinew=-psinew*mu_mag/den

			  !  psi(i,j) = psi(i,j)*(1.d0-orp) + orp*psinew

				if(h_steps==1) then
					psi123(1,i,j) = psi0 - orp*res/den
				elseif((h_steps==3).and.(h<3)) then

!!$					if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
!!$									.OR.  &
!!$						((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
!!$						! OUTER REGION
!!$
!!$						psi123(h+1,i,j) = psi0 - orp_3step*res/den
!!$
!!$					else
!!$						!INNER REGION

						psi123(h+1,i,j) = psi0 - res/den

!!$					endif


				elseif((h_steps==3).and.(h==3)) then

!!$					if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
!!$									.OR.  &
!!$						((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
!!$						! OUTER REGION
!!$
!!$						psi123(1,i,j) = psi123(2,i,j)
!!$
!!$					else
!!$						!INNER REGION

						psi123(1,i,j) = (1.d0-orp_3step)*psi123(1,i,j) + 2.d0*orp_3step*psi123(2,i,j)  &
										- orp_3step*psi123(3,i,j)

!!$					endif

				endif

				if(h==h_steps) then

					! Find the max absolute value of psi in the plasma
					if((tri_type==13).and.(sort_grid(i,j,0)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
!					elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
					elseif((bc_type==7).and.(sort_grid(i,j,0)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((tri_type/=13).and.(bc_type/=7)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					endif

				endif

                ! Calculate the norm of the residual error
!                if ( (tri_type==10).and.((ex*ex + ez*ez) > rminor**2)  ) then
!					continue
!				else
					anorm = anorm + dabs(res)
!				endif

             end do
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do

		if((h_steps==3).and.(h<3)) then

		   call bc_psi_rho0(psi123(h+1,:,:),rho,nx,nz)

		endif


	   enddo
	   ! end of the vertical stability cycle

!	   close(77)
!	   close(88)

       ! Set the new maximum value of psic */
       psic = mx
	   if(tri_type==13)  then
		psic_13 = psic - psi_e
	  else
		psic_13 = psic
	  endif
       ! "Finished: Psi Update"

	   	if((Broot==3).and.(psic/=psic_flag)) then

			psic_flag = 1.d6

			dc_loc = dofpsi(psic)
			dcp_loc = dddpsi(psic)
			pc_loc = pofpsi(psic)
			pcp_loc = dpdpsi(psic)
			b0c_loc = bzero(psic)
			b0cp_loc = dbzerodpsi(psic)

			psic_flag = psic

		endif

       ! -----------------------------------------------------

       ! Move the density calculation to a separate function

		if(Broot/=3) then

		! Update rho and seek mach theta max
			call update_rho(psi123(1,:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
								min_ix,min_iz)

		elseif(Broot==3) then

			call update_rho_super(psi123(1,:,:),rho,b_phi,nx,nz,0,mtm_acc,min_drho)

		endif


	   call update_b(psi123(1,:,:),rho,b_phi,nx,nz,orp,anorm2)

!   	   call bc_psi_rho0(psi123(1,:,:),rho,nx,nz)


		do j=1,nz
		do i=1,nx

			psi_in(i,j) = psi123(1,i,j)

		enddo
		enddo


   	   call bc_psi_rho0(psi_in,rho,nx,nz)

		do j=1,nz
		do i=1,nx

			psi123(1,i,j) = psi_in(i,j)

		enddo
		enddo


	   if(tri_type==11) then
	   ! update LDX stuff
			bpol_max = bpol_max_temp
			bpol_min = bpol_min_temp
			bpol_av = 0.5d0*(bpol_max+bpol_min)

			if((eq3_opt==9).or.(p_opt==9)) B_pmax = B_pmax_temp

!			qperpcenter =  betaperp_center*bpol_max*bpol_max/2.0d0/mu_mag 
!			qperpedge = qperpcenter * qpee_o_qpec
!			qparcenter = betapar_center*bpol_max*bpol_max/2.0d0/mu_mag
!			qparedge = qparcenter * qpae_o_qpac

			if(((eq3_opt==7).or.(eq3_opt==8)).or.  &
				((p_opt==7).or.(p_opt==8))) call get_psi_pres(psi123(1,:,(nz+1)/2),nx)

	   endif

       if(in_orp <= 0.0d0) then
          ! Chebyshev acceleration 
          if(nx >= 5) then
             std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)


          end if

       endif


	   if (accelerate) then
		    orp = std_orp
	        if (bc_type==2) orp = dmin1(orp,max_orp)
	   else
	        orp = fix_orp
	   endif


       ! -----------------------------------------------------

       if(k == 1) then
          eps_anorm = eps*anorm
		  eps_anorm2 = eps*anorm2
		  eps_anorm2 = dmax1(1.d-6,eps_anorm2)
		  if((eq_type==3).and.(tri_type==11)) last_anorm = anorm
       else
          if((in_orp <= 0.0d0).and.accelerate) then
             ! As we get closer to convergence, we want the solution
             ! to relax.  So as anorm approaches eps_anorm we want 
             ! the over relaxation parameter to go to some const. ~ 1
             ! Use x and mtm_soln as a temporary variable
             x = anorm/eps_anorm
             mtm_soln = 1.0d0 ! The limiting value for x ~ 1
             tmp = x/(x + orp/mtm_soln - 1.0d0)
             tmp = dmin1(tmp,1.0d0)
             orp = orp*tmp
             orp = dmax1(orp,mtm_soln)
          endif
       end if

	if( (k<=50).and.(k>=25) ) then 
!		fix_orp = (fix_orp1-fix_orp0)/25.d0*(k-25.d0) + fix_orp0
		continue
	endif

!!$ 			call step_output(nx,psi123(1,:,:),rho,residual)


       if(k == 1 .or. modulo(k,25) == 0) then
          print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		  print *, k,": anorm2 = ",real(anorm2,skind)," eps_anorm2 = ",real(eps_anorm2,skind)
          print *, "The Over-Relaxation Parameter =",orp,std_orp

 			call step_output(nx,psi123(1,:,:),rho,residual)
			if((tri_type==-2).and.(k==50).and.(.not.(tri_type_m2_ready))) then
				initialize_r_for_tri_type_m2 = .true.
			endif
			! Only update the plasma/vacuum interface every 25 iterations
			if ((k>25).and.(tri_type==13)) call update_interface(psi123(1,:,:),nx,inorm)
			if(bc_type==7) call update_sort_grid(psi123(1,:,:),nx,nz,inorm)
			continue
      end if



		write(111,*) k,anorm

	if( (tri_type==11).and.(eq_type==3)) then

		if((reduce).or.(anorm>last_anorm)) then

			bpol0_temp = bpol0_temp*0.96d0
			ana_fac = ana_fac*0.96d0

		else

			bpol0_temp = bpol0_temp/0.95d0
			if(bpol0_temp>bpol0_fix) bpol0_temp = bpol0_fix

			ana_fac = ana_fac/0.95d0
			if(ana_fac>ana_fac0) ana_fac = ana_fac0

		endif

		last_anorm = anorm
		reduce = .false.

		psi_B_min = psi123(1,i_bmin,j_bmin)/psic

	endif


!		call cpu_time(toc)
!
!		print*, 'iteration time = ', toc-tic
!		pause

!		call step_output(nx,psi,rho,residual)

       ! Check the convergence criteria
	   if(anorm < eps_anorm .and. k > min_it .and. anorm2 < eps_anorm2 .and. inorm<1.d-5) exit
    end do


	do j=1,nz
	do i=1,nx

		psi_in(i,j) = psi123(1,i,j)

	enddo
	enddo

	deallocate(psi123)


    print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
	print *, "anorm2 = ",anorm2," eps_anorm2 = ",eps_anorm2
    print *, "Average residual =",anorm/(nx*nz)

    if ((broot==0).or.(Broot==5)) print *, "Mach Theta Max =",mach_theta_max
    print *, "Final Solution has Psi Center = ",psic
    print *


	return

  end subroutine ngs_solve


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve_no_limiter(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
    real (kind=dkind), intent(in) :: eps
    ! Input Over Relaxation Parameter
    ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
    ! else we set orp = in_orp
    real (kind=dkind), intent(in) :: in_orp
    integer ::ipass,i,isw,j,jsw,k,p,h,alloc_stat
    real (kind=dkind) :: res, den, x
    ! res -> residual
    ! den -> denominator of the Newton-Gauss-Seidel update
    ! x -> Radial position of grid point (i,j)
    real (kind=dkind) :: dx2,dz2,mx
    real (kind=dkind) :: anorm, eps_anorm
	real (kind=dkind) :: anorm2, eps_anorm2	! these 2 are fo b_phi
    real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
    real (kind=dkind) :: orp, std_orp ! over relaxation parameter
    ! Phi @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: phirx,philx,phirz,philz
    ! Density @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    ! The square of the Alfvenic Mach Number
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
    real (kind=dkind) :: rhoc,phic,omegac,deltac,thetac,tparc,dc
    ! by is the phi component of the magnetic field
    ! b_dot_v is the inner product of B and V
    real (kind=dkind) :: by, b_dot_v
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: last_mtm ! last mach theta max
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: min_drho, tmp
    integer :: min_ix, min_iz
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: mtm_soln
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), parameter :: mtm_acc = 1.0d-12
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    real (kind=dkind) :: drho_ds ! Used with a slope limiter
	! The following is used to allow the possibility of anisotropic pressure
	real (kind=dkind) :: deltarx,deltalx,deltarz,deltalz
	real (kind=dkind) :: bphirx,bphilx,bphirz,bphilz
	real (kind=dkind), dimension (1:3,1:3) :: psi_delta,psi3x3
	real (kind=dkind) :: b_field_l,b_pol_l,dpsidx,dpsidz,psinew
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0
	real(kind=dkind) :: bpol_min_temp,bpol_max_temp, psi_pmax_temp, B_pmax_temp
	real (kind=dkind), dimension (-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: last_anorm
	real(kind=dkind) :: inorm = 0.d0

    real (kind=dkind), dimension(:,:,:), allocatable :: psi123
	integer :: h_steps=1
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: tic,toc
	integer :: pippa
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: orp_3step = 5.d-2

	integer iii,jjj

	if((tri_type==13).and.(nx>=inter_switch).and.(nz>=inter_switch)) then
!	if((tri_type==13).and.(nx>inter_switch).and.(nz>inter_switch)) then

		h_steps = 3

	endif

	if(allocated(psi123)) deallocate(psi123)
	allocate(psi123(1:h_steps,1:nx,1:nz),stat = alloc_stat)
	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi123"
		 pause
		 stop
	endif

	psi123 = 0.d0

	do h=1,h_steps
	do j=1,nz
	do i=1,nx

		psi123(h,i,j) = psi_in(i,j)

	enddo
	enddo
	enddo

	if((eq_type==3).and.(tri_type==11)) then
		bpol0_temp = bpol0_fix
		if(nx>n_min) then
			if(ana_fac==0.d0) then
				ana_fac = ana_fac0*1.d-3
			else
				ana_fac = ana_fac0
			endif
		endif
	endif

    if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
    last_mtm = mach_theta_max
    eps_anorm = 0.0d0
	eps_anorm2 = 0.0d0

	  if((tri_type==11).and.(((eq3_opt==7).or.(eq3_opt==8)).or.  &
			((p_opt==7).or.(p_opt==8)))) call get_psi_pres(psi123(1,:,(nz+1)/2),nx)

	fix_orp0 = fix_orp

!!$	 	call step_output(nx,psi,rho,residual)

    ! The under/over relaxation parameter
    ! Standard Case:
    if(in_orp <= 0.0d0) then
       ! Note that the parameters below are tuned for starting with 
       ! a 17 x 17 grid.
       orp = 1.0d0
       if(nx <= 5) then
          orp = 0.5d0
          rjac = 0.75d0
       else if(nx <= 9) then
          orp = 0.5d0
          rjac = 0.9d0
       else if(nx <= 17) then
          rjac = 0.98d0
       else if(nx <= 33) then
          rjac = 0.995d0
       else if(nx <= 65) then
          rjac = 0.9972d0
       else 
!!$ It would appear that 129 is about where the solution begins to converge
          rjac = 0.0d0
       end if
       std_orp = orp
    else
       print *, "Input Over Relaxation Parameter =",in_orp
       orp = in_orp
       std_orp = in_orp
    endif

   if (accelerate) then
		continue
   else
	    orp = fix_orp
   endif


    dx2 = dx*dx
    dz2 = dz*dz

    if(Broot/=3) then

	! Update rho before "relaxing psi" but do not seek mach theta max
	    call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
							min_ix,min_iz)

	elseif(Broot==3) then

	    call update_rho_super(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho)

	endif

! set up a few functions for Broot = 3

	if(Broot==3) then

		psic_flag = 1.d6

		dc_loc = dofpsi(psic)
		dcp_loc = dddpsi(psic)
		pc_loc = pofpsi(psic)
		pcp_loc = dpdpsi(psic)
		b0c_loc = bzero(psic)
		b0cp_loc = dbzerodpsi(psic)

		psic_flag = psic

	endif


    ! Iterate until convergence
    do k=1, max_it
!!$     "Update Psi"

!		call cpu_time(tic)

       ! Reset the norm and max of Psi
       mx = 0.0d0
       anorm = 0.0d0
	   anorm2 = 0.0d0

	   do j=2,nz-1
	   do i=2,nx-1
			if ( ((tri_type==13).and.(sort_grid(i,j,0)==2)).or.  &
				((tri_type/=13).and.(bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j,0)==1)) ) then
					psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
				psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			endif
	   enddo
	   enddo

	   bpol_max_temp = 0.d0
	   bpol_min_temp = 1.d22
	   psi_pmax_temp = 0.d0

!	   open(77,file='check.dat')
!	   open(88,file='checkw.dat')


		do h=1,h_steps
		! vertical stabilization cycle

       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
             do i=isw+1, nx-1, 2
                ! "Update Psi(",i,",",j,")"
                ! Only solve the inner region problem

				if(tri_type==10) then

!!$					if((ex*ex + ez*ez) >= radius_ext**2) then
!!$						cycle
!!$					endif

				else

					if(sort_grid(i,j,0)<=0) then
					   cycle
					end if

				endif

				! set up local psi values
				! this should save considerable resources in looking up values


				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi123(h,i,j)
				if(i>1) psil = psi123(h,i-1,j)
				if(i<nx) psir = psi123(h,i+1,j)
				if(j>1) psid = psi123(h,i,j-1)
				if(j<nz) psiu = psi123(h,i,j+1)

				psi_around(-1,0) = psil
				psi_around(1,0) = psir
				psi_around(0,0) = psi0
				psi_around(0,-1) = psid
				psi_around(0,1) = psiu

				! set up functions of psi
				! NOTE: DIFFERENTIATE FOR EQ_TYPE=3 LATER ON

				psi_flag = 1.d9
				psi_flag_dep = 1.d9
				psi_flag_dep2 = 1.d9
				psi_flag_ham = 1.d9

				d_loc = dofpsi(psi0)
				dp_loc = dddpsi(psi0)

				p_loc = pofpsi(psi0)
				pp_loc = dpdpsi(psi0)

				psi_flag = psi0

				b0_loc = bzero(psi0)
				b0p_loc = dbzerodpsi(psi0)

				psi_flag_dep = psi0

				mth_loc = mach_theta(psi0)
				mthp_loc = dmach_thetadpsi(psi0)

				mph_loc = mach_phi(psi0)
				mphp_loc = dmach_phidpsi(psi0)

				psi_flag_dep2 = psi0

				s_loc = sofpsi(psi0)
				sp_loc = dsdpsi(psi0)

				phi_loc = phiofpsi(psi0)
				phip_loc = dphidpsi(psi0)

				omega_loc = omegaofpsi(psi0)
				omegap_loc = domegadpsi(psi0)

				i_loc = iofpsi(psi0)
				ip_loc = didpsi(psi0)

				h_loc = hofpsi(psi0)
				hp_loc = dhdpsi(psi0)

				psi_flag_ham = psi0

				! end of functions set up

                x = x_coord(i)

                ! Calculate rho, phi and omega at the current location
                rhoc = rho(i,j)
                phic = phiofpsi(psi0)
                omegac = omegaofpsi(psi0)

				if(eq_type==3) then 
					thetac = thetaofpsi(psi0)
					tparc = tparofpsi(psi0)
				endif

				dc = dofpsi(psi0)
                ! Calculate B_phi = by
				if (eq_type == 1) then
					by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
						 (1.0d0 - phic*phic/rhoc)
				elseif (eq_type == 3) then
					deltac=deltaofpsi(psi123(h,i-1:i+1,j-1:j+1),i,j,rho(i,j),  &
									b_phi(i,j),1.d0,nx,nz)
					by = b_phi(i,j)
				endif

                ! -----------------------------------------------------

                ! Calculate Phi & Density, rho @ +/- 0.5*dx, +/- 0.5*dz
                ! and then the square of the Alfvenic Mach Number.
!!$  NOTE: The use of a slope limiter to interpolate rho is essential
!!$  for eliminating an oscillation in Bpol where the solution switches
!!$  from the sub-slow to super-slow root. -- 3/14/2002 -- T. Gardiner
                ! -----------------------------------------------------



				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimx(0) = ( dx_a(i-1)**2*psir +  &
									(dx_a(i)**2-dx_a(i-1)**2)*psi0 -  &
									dx_a(i)**2*psil ) /  &
									( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )

					ma2c = phic*phic/rhoc

				endif

!                drho_ds = van_Leer_slope(rho(i-1,j),rho(i,j),rho(i+1,j),dx)
                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                ! Right x
                phirx = phiofpsi(0.5d0*(psir + psi0))
!                rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                rhorx = (rho(i,j) + rho(i+1,j))/2.d0
                ma2rx = phirx*phirx/rhorx
				if(eq_type==3) bphirx = (b_phi(i+1,j) + b_phi(i,j))*0.5d0
				psiprimx(1) = (psir - psi0)/dx_a(i)

                ! Left x
                philx = phiofpsi(0.5d0*(psil + psi0))
!                rholx = rho(i,j) - 0.5d0*dx*drho_ds
                rholx = (rho(i,j)+rho(i-1,j))/2.d0
                ma2lx = philx*philx/rholx
				if(eq_type==3) bphilx = (b_phi(i-1,j) + b_phi(i,j))*0.5d0
				psiprimx(-1) = (psi0 - psil)/dx_a(i-1)

                ! -----------------------------------------------------

!                drho_ds = van_Leer_slope(rho(i,j-1),rho(i,j),rho(i,j+1),dz)
                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimz(0) = ( dz_a(j-1)**2*psiu +  &
									(dz_a(j)**2-dz_a(j-1)**2)*psi0 -  &
									dz_a(j)**2*psid ) /  &
									( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

				endif

                ! Right z
                phirz = phiofpsi(0.5d0*(psi0 + psiu))
!                rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                rhorz = (rho(i,j) + rho(i,j+1))/2.d0
                ma2rz = phirz*phirz/rhorz
				if(eq_type==3) bphirz = (b_phi(i,j) + b_phi(i,j+1))*0.5d0
				psiprimz(1) = (psiu - psi0)/dz_a(j)


                ! Left z
                philz = phiofpsi(0.5d0*(psi0 + psid))
!                rholz = rho(i,j) - 0.5d0*dz*drho_ds
                rholz = (rho(i,j) +rho(i,j-1))/2.d0
                ma2lz = philz*philz/rholz
				if(eq_type==3) bphilz = (b_phi(i,j) + b_phi(i,j-1))*0.5d0
				psiprimz(-1) = (psi0 - psid)/dz_a(j-1)


				! -----------------------------------------------------


				! calculate the magnetic field

!				call cpu_time(tic)

!				do pippa=1,1000

!				if(eq_type==3) then

					call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
!					call psi_derivative(i,j,nx,nz,psi123(h,:,:),dpsidx,dpsidz)
					b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)
					b_field_l = dsqrt(by**2+b_pol_l**2)

!				endif
!				enddo

!				call cpu_time(toc)
!				print*, 'psi_derivative time', toc-tic
!				pause

				if(tri_type==11) then

					if (b_pol_l>bpol_max_temp) bpol_max_temp = b_pol_l
					if (b_pol_l<bpol_min_temp) then
						bpol_min_temp = b_pol_l
						i_bmin = i
						j_bmin = j
					endif

					if(((p_opt==9).or.(eq3_opt==9)).and.(j==1+nz/2)) then
						if(abs(psi_pmax-psi0)<abs(psi_pmax-psi_pmax_temp)) then
							B_pmax_temp = b_pol_l
							psi_pmax_temp = psi0
						endif
					endif

				endif

				! -----------------------------------------------------

				! Calculate Delta (function of p_par and p_perp)

				if(eq_type==3) then

				! Right x

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=5.d-1*(psi123(h,i-1,j+1)+psi0)
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=psi0
				psi_delta(2,2)=5.d-1*(psi0+psiu)
				psi_delta(2,3)=psiu
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=5.d-1*(psi0+psi123(h,i+1,j+1))
				psi_delta(3,3)=0.d0

				deltarx=deltaofpsi(psi_delta,i,j,rhorx,bphirx,0.5d0,nx,nz)

				! Left x

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=5.d-1*(psi123(h,i-1,j-1)+psi0)
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=psid
				psi_delta(2,2)=5.d-1*(psi0+psid)
				psi_delta(2,3)=psi0
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=5.d-1*(psi0+psi123(h,i+1,j-1))
				psi_delta(3,3)=0.d0

				deltalx=deltaofpsi(psi_delta,i,j,rholx,bphilx,0.5d0,nx,nz)

				! Right z

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=psi0
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=5.d-1*(psi0+psi123(h,i+1,j-1))
				psi_delta(2,2)=5.d-1*(psi0+psir)
				psi_delta(2,3)=5.d-1*(psi0+psi123(h,i+1,j+1))
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=psir
				psi_delta(3,3)=0.d0

				deltarz=deltaofpsi(psi_delta,i,j,rhorz,bphirz,0.5d0,nx,nz)

				! Left z

				psi_delta(1,1)=0.d0
				psi_delta(1,2)=psil
				psi_delta(1,3)=0.d0
				psi_delta(2,1)=5.d-1*(psi0+psi123(h,i-1,j-1))
				psi_delta(2,2)=5.d-1*(psi0+psil)
				psi_delta(2,3)=5.d-1*(psi0+psi123(h,i-1,j+1))
				psi_delta(3,1)=0.d0
				psi_delta(3,2)=psi0
				psi_delta(3,3)=0.d0

				deltalz=deltaofpsi(psi_delta,i,j,rholz,bphilz,0.5d0,nx,nz)

				endif

                ! -----------------------------------------------------

                ! Calculate B dot v
!                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + &
!                     ( (psi(i+1,j) - psi(i-1,j))**2/dx2 &
!                     + (psi(i,j+1) - psi(i,j-1))**2/dz2 )*0.25d0/(x*x) )
                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)

                ! -----------------------------------------------------


				if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
								.OR.  &
					((bc_type==7).AND.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))  &
								.OR.  &
					((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
					! OUTER REGION

!!$					if(h_steps==1) then

						if(grid_type==0) then

							res = (1.d0/mu_mag)*(( (1.0d0 )*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 )*(psil-psi0)/(x - 0.5d0*dx) )/(x*dx2) &
								+ ( (1.0d0 )*(psiu-psi0) &
								  + (1.0d0 )*(psid-psi0) )/(x*x*dz2) ) 

						else

							res = (2.d0/mu_mag)*( (dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										- dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										(dz_a(j-1)**2*psiprimz(1)   &
										- dz_a(j)**2*psiprimz(-1) +  &
										(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif

!!$					else
!!$
!!$						res = 0.d0
!!$
!!$					endif






				else	!inner region

				    if (eq_type==1) then

!!$						if(h_steps==1) then

							if(grid_type==0) then

								term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
								  /(x*dx2) &
								+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
								  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

							else

								term0 = (2.d0/mu_mag)*( ((1.0d0 - ma2rx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
											-(1.0d0 - ma2lx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1))   &
											+(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x )   &
											/(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1))   &

											+((1.0d0 - ma2rz)*dz_a(j-1)**2*psiprimz(1)   &
											-(1.0d0 - ma2lz)*dz_a(j)**2*psiprimz(-1)   &
											+(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) )   &
											/(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

							endif

!!$						else
!!$
!!$							term0 = 0.d0
!!$
!!$						endif



						term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
						term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
						term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
						term4= rhoc*dhdpsi(psi0) 
						term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

						res = term0+term1+term2+term3+term4-term5

					elseif (eq_type==3) then

!!$						if(h_steps==1) then

							term0 = (2.d0/mu_mag)*( ((1.0d0 - deltarx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))  &
										- (1.0d0 - deltalx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(1.0d0 - deltac)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										((1.0d0 - deltarz)*dz_a(j-1)**2*psiprimz(1)  &
										- (1.0d0 - deltalz)*dz_a(j)**2*psiprimz(-1) +  &
										(1.0d0 - deltac)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

!!$						else
!!$
!!$							term0 = 0.d0
!!$
!!$						endif

						res = term0  &

							+ rhoc*dtpardpsi(psi0)  &
							+ by*didpsi(psi0)/x/dsqrt(mu_mag) &
							+ rhoc*x**2*omegac*domegadpsi(psi0) &
							+ rhoc*dhdpsi(psi0) &
							- rhoc*(			&
									dtpardpsi(psi0)*dlog(rhoc/dc*dabs(b_field_l-thetac*tparc)/b_field_l) &
									- tparc*dddpsi(psi0)/dc  &
									- tparc  &
									  *(tparc*dthetadpsi(psi0)+thetac*dtpardpsi(psi0)) &
									  *sign(1.d0,(b_field_l-thetac*tparc)) )

						if(b_field_l-thetac*tparc<0.1d0*b_field_l) then

							print *, b_field_l-thetac*tparc, k, i, j
							reduce = .true.
							continue

						endif

					endif


				continue

				endif




				res = res*mu_mag


                ! Store the residual
                if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
                   residual(i,j) = res
                end if

				if(res>1.d0) then
					continue
				endif

                ! -----------------------------------------------------

                if (eq_type==1) then

						if(grid_type==0) then

					den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
				             (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
					      -( (1.0d0 - ma2rz) + &
						     (1.0d0 - ma2lz) )/(x*x*dz2)

				else

					den = -2.d0*( ((1.0d0 - ma2rx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - ma2lx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - ma2rz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - ma2lz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

				endif

					continue


				elseif (eq_type==3) then

!!$					den = -( (1.0d0 - deltarx)/(x + 0.5d0*dx) + &
!!$							(1.0d0 - deltalx)/(x - 0.5d0*dx) )/(x*dx2) &
!!$							-( (1.0d0 - deltarz) + &
!!$							(1.0d0 - deltalz) )/(x*x*dz2)


					den = -2.d0*( ((1.0d0 - deltarx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
									(1.0d0 - deltalx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
									(1.0d0 - deltac)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
									(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

									((1.0d0 - deltarz)*dz_a(j-1)**2/dz_a(j) +  &
									(1.0d0 - deltalz)*dz_a(j)**2/dz_a(j-1) +  &
									(1.0d0 - deltac)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
									/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
									(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

				endif

     !    				psinew=-psinew*mu_mag/den

			  !  psi(i,j) = psi(i,j)*(1.d0-orp) + orp*psinew

				if(h_steps==1) then
					psi123(1,i,j) = psi0 - orp*res/den
				elseif((h_steps==3).and.(h<3)) then

!!$					if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
!!$									.OR.  &
!!$						((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
!!$						! OUTER REGION
!!$
!!$						psi123(h+1,i,j) = psi0 - orp_3step*res/den
!!$
!!$					else
!!$						!INNER REGION

						psi123(h+1,i,j) = psi0 - res/den

!!$					endif


				elseif((h_steps==3).and.(h==3)) then

!!$					if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
!!$									.OR.  &
!!$						((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
!!$						! OUTER REGION
!!$
!!$						psi123(1,i,j) = psi123(2,i,j)
!!$
!!$					else
!!$						!INNER REGION

						psi123(1,i,j) = (1.d0-orp_3step)*psi123(1,i,j) + 2.d0*orp_3step*psi123(2,i,j)  &
										- orp_3step*psi123(3,i,j)

!!$					endif

				endif

				if(h==h_steps) then

					! Find the max absolute value of psi in the plasma
					if((tri_type==13).and.(sort_grid(i,j,0)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((tri_type/=13).and.(bc_type/=7)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					endif

				endif

                ! Calculate the norm of the residual error
!                if ( (tri_type==10).and.((ex*ex + ez*ez) > rminor**2)  ) then
!					continue
!				else
					anorm = anorm + dabs(res)
!				endif

             end do
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do

		if((h_steps==3).and.(h<3)) then

		   call bc_psi_rho0(psi123(h+1,:,:),rho,nx,nz)

		endif


	   enddo
	   ! end of the vertical stability cycle

!	   close(77)
!	   close(88)

       ! Set the new maximum value of psic */
       psic = mx
	   if(tri_type==13)  then
		psic_13 = psic - psi_e
	  else
		psic_13 = psic
	  endif
       ! "Finished: Psi Update"

	   	if((Broot==3).and.(psic/=psic_flag)) then

			psic_flag = 1.d6

			dc_loc = dofpsi(psic)
			dcp_loc = dddpsi(psic)
			pc_loc = pofpsi(psic)
			pcp_loc = dpdpsi(psic)
			b0c_loc = bzero(psic)
			b0cp_loc = dbzerodpsi(psic)

			psic_flag = psic

		endif

       ! -----------------------------------------------------

       ! Move the density calculation to a separate function

		if(Broot/=3) then

		! Update rho and seek mach theta max
			call update_rho(psi123(1,:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
								min_ix,min_iz)

		elseif(Broot==3) then

			call update_rho_super(psi123(1,:,:),rho,b_phi,nx,nz,0,mtm_acc,min_drho)

		endif


	   call update_b(psi123(1,:,:),rho,b_phi,nx,nz,orp,anorm2)

!   	   call bc_psi_rho0(psi123(1,:,:),rho,nx,nz)


		do j=1,nz
		do i=1,nx

			psi_in(i,j) = psi123(1,i,j)

		enddo
		enddo


   	   call bc_psi_rho0(psi_in,rho,nx,nz)

		do j=1,nz
		do i=1,nx

			psi123(1,i,j) = psi_in(i,j)

		enddo
		enddo


	   if(tri_type==11) then
	   ! update LDX stuff
			bpol_max = bpol_max_temp
			bpol_min = bpol_min_temp
			bpol_av = 0.5d0*(bpol_max+bpol_min)

			if((eq3_opt==9).or.(p_opt==9)) B_pmax = B_pmax_temp

!			qperpcenter =  betaperp_center*bpol_max*bpol_max/2.0d0/mu_mag 
!			qperpedge = qperpcenter * qpee_o_qpec
!			qparcenter = betapar_center*bpol_max*bpol_max/2.0d0/mu_mag
!			qparedge = qparcenter * qpae_o_qpac

			if(((eq3_opt==7).or.(eq3_opt==8)).or.  &
				((p_opt==7).or.(p_opt==8))) call get_psi_pres(psi123(1,:,(nz+1)/2),nx)

	   endif

       if(in_orp <= 0.0d0) then
          ! Chebyshev acceleration 
          if(nx >= 5) then
             std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)


          end if

       endif


	   if (accelerate) then
		    orp = std_orp
	        if (bc_type==2) orp = dmin1(orp,max_orp)
	   else
	        orp = fix_orp
	   endif


       ! -----------------------------------------------------

       if(k == 1) then
          eps_anorm = eps*anorm
		  eps_anorm2 = eps*anorm2
		  eps_anorm2 = dmax1(1.d-6,eps_anorm2)
		  if((eq_type==3).and.(tri_type==11)) last_anorm = anorm
       else
          if((in_orp <= 0.0d0).and.accelerate) then
             ! As we get closer to convergence, we want the solution
             ! to relax.  So as anorm approaches eps_anorm we want 
             ! the over relaxation parameter to go to some const. ~ 1
             ! Use x and mtm_soln as a temporary variable
             x = anorm/eps_anorm
             mtm_soln = 1.0d0 ! The limiting value for x ~ 1
             tmp = x/(x + orp/mtm_soln - 1.0d0)
             tmp = dmin1(tmp,1.0d0)
             orp = orp*tmp
             orp = dmax1(orp,mtm_soln)
          endif
       end if

	if( (k<=50).and.(k>=25) ) then 
!		fix_orp = (fix_orp1-fix_orp0)/25.d0*(k-25.d0) + fix_orp0
		continue
	endif

!!$ 			call step_output(nx,psi123(1,:,:),rho,residual)


       if(k == 1 .or. modulo(k,25) == 0) then
          print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		  print *, k,": anorm2 = ",real(anorm2,skind)," eps_anorm2 = ",real(eps_anorm2,skind)
          print *, "The Over-Relaxation Parameter =",orp,std_orp

 			call step_output(nx,psi123(1,:,:),rho,residual)
			if ((k>25).and.(tri_type==13)) call update_interface(psi123(1,:,:),nx,inorm)
			continue
      end if



		write(111,*) k,anorm

	if( (tri_type==11).and.(eq_type==3)) then

		if((reduce).or.(anorm>last_anorm)) then

			bpol0_temp = bpol0_temp*0.96d0
			ana_fac = ana_fac*0.96d0

		else

			bpol0_temp = bpol0_temp/0.95d0
			if(bpol0_temp>bpol0_fix) bpol0_temp = bpol0_fix

			ana_fac = ana_fac/0.95d0
			if(ana_fac>ana_fac0) ana_fac = ana_fac0

		endif

		last_anorm = anorm
		reduce = .false.

		psi_B_min = psi123(1,i_bmin,j_bmin)/psic

	endif


!		call cpu_time(toc)
!
!		print*, 'iteration time = ', toc-tic
!		pause

!		call step_output(nx,psi,rho,residual)

       ! Check the convergence criteria
	   if(anorm < eps_anorm .and. k > min_it .and. anorm2 < eps_anorm2 .and. inorm<1.d-5) exit
    end do


	do j=1,nz
	do i=1,nx

		psi_in(i,j) = psi123(1,i,j)

	enddo
	enddo

	deallocate(psi123)


    print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
	print *, "anorm2 = ",anorm2," eps_anorm2 = ",eps_anorm2
    print *, "Average residual =",anorm/(nx*nz)

    if ((broot==0).or.(Broot==5)) print *, "Mach Theta Max =",mach_theta_max
    print *, "Final Solution has Psi Center = ",psic
    print *


	return

  end subroutine ngs_solve_no_limiter


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! a jump condition for the total pressure is enforced at the transonic interface
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine ngs_solve_jump(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a

    implicit none
	integer, intent(in) :: nx,nz,min_it,max_it,flag
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
	real(kind=dkind), dimension(1:nx,1:nz) :: Ptot, lambda, alpha_jump
	real(kind=dkind) :: lambda0, t0_ratio, den_ratio, dpsi_jump
	real(kind=dkind) :: lambda_relax = 1.d-1
	real(kind=dkind), dimension(1:nx,1:nz,2) :: gpsi_b, gPtor !grad psi for jump condition
    real (kind=dkind), intent(in) :: eps
    ! Input Over Relaxation Parameter
    ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
    ! else we set orp = in_orp
    real (kind=dkind), intent(in) :: in_orp
    integer ::ipass,i,isw,j,jsw,k,p,h,alloc_stat
    real (kind=dkind) :: res, den, x
    ! res -> residual
    ! den -> denominator of the Newton-Gauss-Seidel update
    ! x -> Radial position of grid point (i,j)
    real (kind=dkind) :: dx2,dz2,mx
    real (kind=dkind) :: anorm, eps_anorm
	real (kind=dkind) :: anorm2, eps_anorm2	! these 2 are fo b_phi
    real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
    real (kind=dkind) :: orp, std_orp ! over relaxation parameter
    ! Phi @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: phirx,philx,phirz,philz
    ! Density @ right 0.5*dx, left 0.5*dx, right 0.5*dz, left 0.5*dz 
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    ! The square of the Alfvenic Mach Number
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
    real (kind=dkind) :: rhoc,phic,omegac,dc
    ! by is the phi component of the magnetic field
    ! b_dot_v is the inner product of B and V
    real (kind=dkind) :: by, b_dot_v
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: last_mtm ! last mach theta max
    ! minimum delta rho for bernoulli roots, the mean val. and location
    real (kind=dkind) :: min_drho, tmp
    integer :: min_ix, min_iz
    ! The next set of variables are used for the bisection search for
    ! locating the degenerate root for Mach Theta Max
    real (kind=dkind) :: mtm_soln
    ! mtm_acc controls the bisection loop seeking mach theta max
    real (kind=dkind), parameter :: mtm_acc = 1.0d-12
    ! VERY IMPORTANT
    ! psi_degen record the value of psi where the solution to the
    ! Bernoulli equation is approximately degenerate.  It is 
    ! approximately equal to 0.5d0*psic but generally slightly less.
!!$    real (kind=dkind) :: psi_degen
    real (kind=dkind) :: drho_ds ! Used with a slope limiter
	! The following is used to allow the possibility of anisotropic pressure
	real (kind=dkind) :: b_field_l,b_pol_l,dpsidx,dpsidz,dpsidx_b,dpsidz_b,psinew
	real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term6, term00
	real(kind=dkind) :: bpol_min_temp,bpol_max_temp, psi_pmax_temp, B_pmax_temp
	real (kind=dkind), dimension (-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: last_anorm
	real(kind=dkind) :: inorm = 0.d0

    real (kind=dkind), dimension(:,:,:), allocatable :: psi123
	integer :: h_steps=1
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: tic,toc
	integer :: pippa
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: orp_3step = 5.d-2

	integer iii,jjj

	if((tri_type==13).and.(nx>=inter_switch).and.(nz>=inter_switch)) then
!	if((tri_type==13).and.(nx>inter_switch).and.(nz>inter_switch)) then

		h_steps = 3

	endif

	if(allocated(psi123)) deallocate(psi123)
	allocate(psi123(1:h_steps,1:nx,1:nz),stat = alloc_stat)
	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi123"
		 pause
		 stop
	endif

	psi123 = 0.d0

	do h=1,h_steps
	do j=1,nz
	do i=1,nx

		psi123(h,i,j) = psi_in(i,j)

	enddo
	enddo
	enddo

	if((eq_type==3).and.(tri_type==11)) then
		bpol0_temp = bpol0_fix
		if(nx>n_min) then
			if(ana_fac==0.d0) then
				ana_fac = ana_fac0*1.d-3
			else
				ana_fac = ana_fac0
			endif
		endif
	endif

    if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
    last_mtm = mach_theta_max
    eps_anorm = 0.0d0
	eps_anorm2 = 0.0d0

	  if((tri_type==11).and.(((eq3_opt==7).or.(eq3_opt==8)).or.  &
			((p_opt==7).or.(p_opt==8)))) call get_psi_pres(psi123(1,:,(nz+1)/2),nx)

	fix_orp0 = fix_orp

!!$	 	call step_output(nx,psi,rho,residual)

    ! The under/over relaxation parameter
    ! Standard Case:
    if(in_orp <= 0.0d0) then
       ! Note that the parameters below are tuned for starting with 
       ! a 17 x 17 grid.
       orp = 1.0d0
       if(nx <= 5) then
          orp = 0.5d0
          rjac = 0.75d0
       else if(nx <= 9) then
          orp = 0.5d0
          rjac = 0.9d0
       else if(nx <= 17) then
          rjac = 0.98d0
       else if(nx <= 33) then
          rjac = 0.995d0
       else if(nx <= 65) then
          rjac = 0.9972d0
       else 
!!$ It would appear that 129 is about where the solution begins to converge
          rjac = 0.0d0
       end if
       std_orp = orp
    else
       print *, "Input Over Relaxation Parameter =",in_orp
       orp = in_orp
       std_orp = in_orp
    endif

   if (accelerate) then
		continue
   else
	    orp = fix_orp
   endif


    dx2 = dx*dx
    dz2 = dz*dz

	lambda0 = 0.d0
	lambda = 0.d0

	! Update rho before "relaxing psi" but do not seek mach theta max
	call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
						min_ix,min_iz)


    ! Iterate until convergence
    do k=1, max_it
!!$     "Update Psi"

!		call cpu_time(tic)

       ! Reset the norm and max of Psi
       mx = 0.0d0
       anorm = 0.0d0
	   anorm2 = 0.0d0

	   do j=2,nz-1
	   do i=2,nx-1
			if ( ((tri_type==13).and.(sort_grid(i,j,0)==2)).or.  &
				((tri_type/=13).and.(bc_type/=7).and.(bc_type/=8).and.(sort_grid(i,j,0)==1)) ) then
					psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
				psi123(1,i,j) = dabs(psi123(1,i,j))	!	segno
			endif
	   enddo
	   enddo

	   bpol_max_temp = 0.d0
	   bpol_min_temp = 1.d22
	   psi_pmax_temp = 0.d0

!	   open(77,file='check.dat')
!	   open(88,file='checkw.dat')


		do h=1,h_steps
		! vertical stabilization cycle

		t0_ratio = 0.d0
		den_ratio = 0.d0

       jsw = 1
       do ipass = 1, 2
          isw=jsw;
          do j=2, nz-1
             do i=isw+1, nx-1, 2
                ! "Update Psi(",i,",",j,")"
                ! Only solve the inner region problem

				if(tri_type==10) then

!!$					if((ex*ex + ez*ez) >= radius_ext**2) then
!!$						cycle
!!$					endif

				else

					if(sort_grid(i,j,0)<=0) then
					   cycle
					end if

				endif

				! set up local psi values
				! this should save considerable resources in looking up values


				psil = 0.d0
				psir = 0.d0
				psiu = 0.d0
				psid = 0.d0

				psi0 = psi123(h,i,j)
				if(i>1) psil = psi123(h,i-1,j)
				if(i<nx) psir = psi123(h,i+1,j)
				if(j>1) psid = psi123(h,i,j-1)
				if(j<nz) psiu = psi123(h,i,j+1)

				psi_around(-1,0) = psil
				psi_around(1,0) = psir
				psi_around(0,0) = psi0
				psi_around(0,-1) = psid
				psi_around(0,1) = psiu

				! set up functions of psi
				! NOTE: DIFFERENTIATE FOR EQ_TYPE=3 LATER ON

				psi_flag = 1.d9
				psi_flag_dep = 1.d9
				psi_flag_dep2 = 1.d9
				psi_flag_ham = 1.d9

				d_loc = dofpsi(psi0)
				dp_loc = dddpsi(psi0)

				p_loc = pofpsi(psi0)
				pp_loc = dpdpsi(psi0)

				psi_flag = psi0

				b0_loc = bzero(psi0)
				b0p_loc = dbzerodpsi(psi0)

				psi_flag_dep = psi0

				mth_loc = mach_theta(psi0)
				mthp_loc = dmach_thetadpsi(psi0)

				mph_loc = mach_phi(psi0)
				mphp_loc = dmach_phidpsi(psi0)

				psi_flag_dep2 = psi0

				s_loc = sofpsi(psi0)
				sp_loc = dsdpsi(psi0)

				phi_loc = phiofpsi(psi0)
				phip_loc = dphidpsi(psi0)

				omega_loc = omegaofpsi(psi0)
				omegap_loc = domegadpsi(psi0)

				i_loc = iofpsi(psi0)
				ip_loc = didpsi(psi0)

				h_loc = hofpsi(psi0)
				hp_loc = dhdpsi(psi0)

				psi_flag_ham = psi0

				! end of functions set up

                x = x_coord(i)

                ! Calculate rho, phi and omega at the current location
                rhoc = rho(i,j)
                phic = phiofpsi(psi0)
                omegac = omegaofpsi(psi0)


				dc = dofpsi(psi0)
                ! Calculate B_phi = by
				by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
					 (1.0d0 - phic*phic/rhoc)

                ! -----------------------------------------------------

                ! Calculate Phi & Density, rho @ +/- 0.5*dx, +/- 0.5*dz
                ! and then the square of the Alfvenic Mach Number.
!!$  NOTE: The use of a slope limiter to interpolate rho is essential
!!$  for eliminating an oscillation in Bpol where the solution switches
!!$  from the sub-slow to super-slow root. -- 3/14/2002 -- T. Gardiner
                ! -----------------------------------------------------



				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimx(0) = ( dx_a(i-1)**2*psir +  &
									(dx_a(i)**2-dx_a(i-1)**2)*psi0 -  &
									dx_a(i)**2*psil ) /  &
									( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )

					ma2c = phic*phic/rhoc

				endif

!                drho_ds = van_Leer_slope(rho(i-1,j),rho(i,j),rho(i+1,j),dx)
                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                ! Right x
                phirx = phiofpsi(0.5d0*(psir + psi0))
!                rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                rhorx = rho(i,j) + 0.5d0*dx_a(i)*drho_ds
                ma2rx = phirx*phirx/rhorx
				psiprimx(1) = (psir - psi0)/dx_a(i)

                ! Left x
                philx = phiofpsi(0.5d0*(psil + psi0))
!                rholx = rho(i,j) - 0.5d0*dx*drho_ds
                rholx = rho(i,j) - 0.5d0*dx_a(i-1)*drho_ds
                ma2lx = philx*philx/rholx
				psiprimx(-1) = (psi0 - psil)/dx_a(i-1)

                ! -----------------------------------------------------

!                drho_ds = van_Leer_slope(rho(i,j-1),rho(i,j),rho(i,j+1),dz)
                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

				if((grid_type==1).or.(grid_type==2).or.(grid_type==3)) then

					psiprimz(0) = ( dz_a(j-1)**2*psiu +  &
									(dz_a(j)**2-dz_a(j-1)**2)*psi0 -  &
									dz_a(j)**2*psid ) /  &
									( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )

				endif

                ! Right z
                phirz = phiofpsi(0.5d0*(psi0 + psiu))
!                rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                rhorz = rho(i,j) + 0.5d0*dz_a(j)*drho_ds
                ma2rz = phirz*phirz/rhorz
				psiprimz(1) = (psiu - psi0)/dz_a(j)


                ! Left z
                philz = phiofpsi(0.5d0*(psi0 + psid))
!                rholz = rho(i,j) - 0.5d0*dz*drho_ds
                rholz = rho(i,j) - 0.5d0*dz_a(j)*drho_ds
                ma2lz = philz*philz/rholz
				psiprimz(-1) = (psi0 - psid)/dz_a(j-1)


				! -----------------------------------------------------


				! calculate the magnetic field

					call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
!					call psi_derivative(i,j,nx,nz,psi123(h,:,:),dpsidx,dpsidz)
					b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)
					b_field_l = dsqrt(by**2+b_pol_l**2)


				if(tri_type==11) then

					if (b_pol_l>bpol_max_temp) bpol_max_temp = b_pol_l
					if (b_pol_l<bpol_min_temp) then
						bpol_min_temp = b_pol_l
						i_bmin = i
						j_bmin = j
					endif

					if(((p_opt==9).or.(eq3_opt==9)).and.(j==1+nz/2)) then
						if(abs(psi_pmax-psi0)<abs(psi_pmax-psi_pmax_temp)) then
							B_pmax_temp = b_pol_l
							psi_pmax_temp = psi0
						endif
					endif

				endif

                ! -----------------------------------------------------

                ! Calculate B dot v
!                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + &
!                     ( (psi(i+1,j) - psi(i-1,j))**2/dx2 &
!                     + (psi(i,j+1) - psi(i,j-1))**2/dz2 )*0.25d0/(x*x) )
                b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)

                ! -----------------------------------------------------


				if( ((tri_type/=13).AND.((psi0/psic)<fraction))  &
								.OR.  &
					((bc_type==7).AND.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))  &
								.OR.  &
					((tri_type==13).AND.(sort_grid(i,j,0)==1)) ) then
					! OUTER REGION

!!$					if(h_steps==1) then

						if(grid_type==0) then

							res = (1.d0/mu_mag)*(( (1.0d0 )*(psir-psi0)/(x + 0.5d0*dx) &
								  + (1.0d0 )*(psil-psi0)/(x - 0.5d0*dx) )/(x*dx2) &
								+ ( (1.0d0 )*(psiu-psi0) &
								  + (1.0d0 )*(psid-psi0) )/(x*x*dz2) ) 

						else

							res = (2.d0/mu_mag)*( (dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
										- dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1)) +  &
										(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x ) /  &
										(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

										(dz_a(j-1)**2*psiprimz(1)   &
										- dz_a(j)**2*psiprimz(-1) +  &
										(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) ) /  &
										(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

						endif


				else	!inner region


					if(grid_type==0) then

						term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
						  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
						  /(x*dx2) &
						+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
						  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

					else

						term0 = (2.d0/mu_mag)*( ((1.0d0 - ma2rx)*dx_a(i-1)**2*psiprimx(1)/(x+0.5d0*dx_a(i))   &
									-(1.0d0 - ma2lx)*dx_a(i)**2*psiprimx(-1)/(x-0.5d0*dx_a(i-1))   &
									+(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)*psiprimx(0)/x )   &
									/(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1))   &

									+((1.0d0 - ma2rz)*dz_a(j-1)**2*psiprimz(1)   &
									-(1.0d0 - ma2lz)*dz_a(j)**2*psiprimz(-1)   &
									+(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)*psiprimz(0) )   &
									/(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

					endif


					term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
					term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
					term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
					term4= rhoc*dhdpsi(psi0) 
					term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

					if(lambda(i,j)==0.d0) then

						term6 = 0.d0
						term00 = 0.d0

					else

						dpsidx_b = gpsi_b(i,j,1)
						dpsidz_b = gpsi_b(i,j,2)
						term6 = lambda(i,j) * (  &
										dpsidx_b*gPtor(i,j,1) + dpsidz_b*gPtor(i,j,2)  &
										- dpsidx_b/(x**3*mu_mag)*(dpsidx_b**2+dpsidz_b**2) +  &
										1.d0/(mu_mag*x**2)* ( dpsidx_b**2*(psir+psil)/dx2 +  &
										2.d0*dpsidx_b*dpsidz_b*(psi123(h,i+1,j+1)+psi123(h,i-1,j-1)-psi123(h,i-1,j+1)-psi123(h,i+1,j-1))/(4.d0*dx*dz)  &
										+ dpsidz_b**2*(psiu+psid)/dz2 )  &
										- alpha_jump(i,j) )

						term00 = - 2.d0/(mu_mag*x**2)*lambda(i,j)  * (dpsidx_b**2/dx2 + dpsidz_b**2/dz2)*psi0
						term0 = term0 + term00

						t0_ratio = max(t0_ratio,abs(term00/term0))

					endif

					res = term0+term1+term2+term3+term4-term5+term6

				endif


				continue


				res = res*mu_mag


                ! Store the residual
                if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
                   residual(i,j) = res
                end if

				if(res>1.d0) then
					continue
				endif

                ! -----------------------------------------------------


				den = -2.d0*( ((1.0d0 - ma2rx)*dx_a(i-1)**2/(x+0.5d0*dx_a(i))/dx_a(i) +  &
								(1.0d0 - ma2lx)*dx_a(i)**2/(x-0.5d0*dx_a(i-1))/dx_a(i-1) +  &
								(1.0d0 - ma2c)*(dx_a(i)**2-dx_a(i-1)**2)**2/  &
								(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) ) /  &
								(x*(dx_a(i)+dx_a(i-1))*dx_a(i)*dx_a(i-1)) +  &

								((1.0d0 - ma2rz)*dz_a(j-1)**2/dz_a(j) +  &
								(1.0d0 - ma2lz)*dz_a(j)**2/dz_a(j-1) +  &
								(1.0d0 - ma2c)*(dz_a(j)**2-dz_a(j-1)**2)**2  &
								/((dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) ) /  &
								(x**2*(dz_a(j)+dz_a(j-1))*dz_a(j)*dz_a(j-1)) )

				if(lambda(i,j)/=0.d0) then

					den_ratio = max(den_ratio,abs(2.d0/(x**2)*lambda(i,j)  * (dpsidx_b**2/dx2 + dpsidz_b**2/dz2)/den))

					den = den - 2.d0/(x**2)*lambda(i,j)  * (dpsidx_b**2/dx2 + dpsidz_b**2/dz2)

				endif

				continue


     !    				psinew=-psinew*mu_mag/den

			  !  psi(i,j) = psi(i,j)*(1.d0-orp) + orp*psinew

				if(h_steps==1) then
					psi123(1,i,j) = psi0 - orp*res/den
				elseif((h_steps==3).and.(h<3)) then

					psi123(h+1,i,j) = psi0 - res/den

				elseif((h_steps==3).and.(h==3)) then

!!$						!INNER REGION

					psi123(1,i,j) = (1.d0-orp_3step)*psi123(1,i,j) + 2.d0*orp_3step*psi123(2,i,j)  &
									- orp_3step*psi123(3,i,j)

				endif

				if(h==h_steps) then

					! Find the max absolute value of psi in the plasma
					if((tri_type==13).and.(sort_grid(i,j,0)==2)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((bc_type==7).and.((i>nx/3).and.(i<2*nx/3).and.(j>nz/3).and.(j<2*nz/3))) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					elseif((tri_type/=13).and.(bc_type/=7)) then
						mx = dmax1(dabs(psi123(h,i,j)),mx)
					endif

				endif

				anorm = anorm + dabs(res)

             end do
             isw = 3 - isw
          end do
          jsw = 3 - jsw
       end do

		if((h_steps==3).and.(h<3)) then

		   call bc_psi_rho0(psi123(h+1,:,:),rho,nx,nz)

		endif


	   enddo
	   ! end of the vertical stability cycle

       ! Set the new maximum value of psic */
       psic = mx
	   if(tri_type==13)  then
		psic_13 = psic - psi_e
	  else
		psic_13 = psic
	  endif
       ! "Finished: Psi Update"

       ! -----------------------------------------------------

       ! Move the density calculation to a separate function

		! Update rho and seek mach theta max
		call update_rho(psi123(1,:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
							min_ix,min_iz)

		do j=1,nz
		do i=1,nx

			psi_in(i,j) = psi123(1,i,j)

		enddo
		enddo

   	   call bc_psi_rho0(psi_in,rho,nx,nz)

		do j=1,nz
		do i=1,nx

			psi123(1,i,j) = psi_in(i,j)

		enddo
		enddo

		! update lambda0
		if(nx>33) then
			if(jump_option==1) then
				if(lambda0==0.d0) lambda0 = -1.d-6
				if(t0_ratio==0.d0) then
					continue
				else
					if(t0_ratio>den_ratio) then
						lambda0 = lambda0*((1.d0-lambda_relax) + 1.d0/(t0_ratio)*lambda_relax)
					else
						lambda0 = lambda0*((1.d0-lambda_relax) + 1.d0/(den_ratio)*lambda_relax)
					endif
					if(max(t0_ratio,den_ratio)<.5d3) lambda0 = lambda0 *1.25
				endif
			elseif((jump_option==2).or.(jump_option==3).or.(jump_option==4)) then
				lambda0 = lambda0_fix
			endif
		else
			lambda0 = 0.d0
		endif

		call get_jump_arrays(nx, nz,psi_in,rho,lambda0,Ptot,lambda, alpha_jump,gpsi_b,gPtor)

       if(in_orp <= 0.0d0) then
          ! Chebyshev acceleration 
          if(nx >= 5) then
             std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)

          end if

       endif


	   if (accelerate) then
		    orp = std_orp
	        if (bc_type==2) orp = dmin1(orp,max_orp)
	   else
	        orp = fix_orp
	   endif


       ! -----------------------------------------------------

       if(k == 1) then
          eps_anorm = eps*anorm
		  eps_anorm2 = eps*anorm2
		  eps_anorm2 = dmax1(1.d-6,eps_anorm2)
       else
          if((in_orp <= 0.0d0).and.accelerate) then
             ! As we get closer to convergence, we want the solution
             ! to relax.  So as anorm approaches eps_anorm we want 
             ! the over relaxation parameter to go to some const. ~ 1
             ! Use x and mtm_soln as a temporary variable
             x = anorm/eps_anorm
             mtm_soln = 1.0d0 ! The limiting value for x ~ 1
             tmp = x/(x + orp/mtm_soln - 1.0d0)
             tmp = dmin1(tmp,1.0d0)
             orp = orp*tmp
             orp = dmax1(orp,mtm_soln)
          endif
       end if

	if( (k<=50).and.(k>=25) ) then 
!		fix_orp = (fix_orp1-fix_orp0)/25.d0*(k-25.d0) + fix_orp0
		continue
	endif

!!$ 			call step_output(nx,psi123(1,:,:),rho,residual)


       if(k == 1 .or. modulo(k,25) == 0) then
          print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
		  print *, k,": anorm2 = ",real(anorm2,skind)," eps_anorm2 = ",real(eps_anorm2,skind)
          print *, "The Over-Relaxation Parameter =",orp,std_orp

 			call step_output(nx,psi123(1,:,:),rho,residual)
			if ((k>25).and.(tri_type==13)) call update_interface(psi123(1,:,:),nx,inorm)
			continue
      end if

		write(111,*) k,anorm

!		call step_output(nx,psi,rho,residual)

       ! Check the convergence criteria
	   if(anorm < eps_anorm .and. k > min_it .and. anorm2 < eps_anorm2 .and. inorm<1.d-5) exit
    end do


	do j=1,nz
	do i=1,nx

		psi_in(i,j) = psi123(1,i,j)

	enddo
	enddo

	deallocate(psi123)


    print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
	print *, "anorm2 = ",anorm2," eps_anorm2 = ",eps_anorm2
    print *, "Average residual =",anorm/(nx*nz)

    print *, "Mach Theta Max =",mach_theta_max
    print *, "Final Solution has Psi Center = ",psic
    print *


	return

  end subroutine ngs_solve_jump












! ------------------------------------------------------------------
! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this version solves the GS equation with grad Ptot in the RHS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_grad_psi(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen record the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(1:nx,1:nz) :: vphi2_RHS, bphi2_RHS
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  real(kind=dkind), dimension(1:2,1:nx,1:nz) :: gradPtot_RHS
  integer, dimension(1:2,1:nx,1:nz) :: grad_dir
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real(kind=dkind) :: gpsi_switch, gpsi_frac, gpsi_frac_max
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

	gpsi_frac_max = lambda0_fix
 	gpsi_frac = 1.d-1
  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  call set_directions
	  call GS_RHS4


  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j,0)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo

     if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)

     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j,0)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)

              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then

                 ! get |grad psi| in surrounding points (using old iteration values)

                 if(grad_dir(1,i,j)==0) then
                    gpsirx = (gradpsi_RHS(0,i+1,j)+gradpsi_RHS(0,i,j))/2.d0
                    gpsilx= (gradpsi_RHS(0,i,j)+gradpsi_RHS(0,i-1,j))/2.d0
                 elseif(grad_dir(1,i,j)==1) then
                    gpsirx = (gradpsi_RHS(0,i+1,j)+gradpsi_RHS(0,i,j))/2.d0
                    gpsilx= (3.d0*gradpsi_RHS(0,i,j) - gradpsi_RHS(0,i+1,j))/2.d0
                 elseif(grad_dir(1,i,j)==-1) then
                    gpsirx = (3.d0*gradpsi_RHS(0,i,j) - gradpsi_RHS(0,i-1,j))/2.d0
                    gpsilx= (gradpsi_RHS(0,i,j)+gradpsi_RHS(0,i-1,j))/2.d0
                 endif

                 if(grad_dir(2,i,j)==0) then
                    gpsirz = (gradpsi_RHS(0,i,j+1)+gradpsi_RHS(0,i,j))/2.d0
                    gpsilz= (gradpsi_RHS(0,i,j)+gradpsi_RHS(0,i,j-1))/2.d0
                 elseif(grad_dir(2,i,j)==1) then
                    gpsirz = (gradpsi_RHS(0,i,j+1)+gradpsi_RHS(0,i,j))/2.d0
                    gpsilz= (3.d0*gradpsi_RHS(0,i,j) - gradpsi_RHS(0,i,j+1))/2.d0
                 elseif(grad_dir(2,i,j)==-1) then
                    gpsirz = (3.d0*gradpsi_RHS(0,i,j) - gradpsi_RHS(0,i,j-1))/2.d0
                    gpsilz= (gradpsi_RHS(0,i,j)+gradpsi_RHS(0,i,j-1))/2.d0
                 endif


                 term0 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res = term0+term1+term2+term3


                 den = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz -1.d0/gpsilz ) / dz2  &
                      )

				res = res/gradpsi_RHS(0,i,j)**2
				den = den/gradpsi_RHS(0,i,j)**2

				continue

              else

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = (rho(i,j) + rho(i+1,j))/2.d0
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = (rho(i,j)+rho(i-1,j))/2.d0
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = (rho(i,j) + rho(i,j+1))/2.d0
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = (rho(i,j) +rho(i,j-1))/2.d0
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 omegac = omegaofpsi(psi0)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res = term0+term1+term2+term3+term4-term5


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue

              endif


              !res = res*mu_mag


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,25) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					residual(i,j) = res
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
              anorm = anorm + dabs(res)

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
	 if(nx>=inter_switch) gpsi_frac = min(gpsi_frac*1.1,gpsi_frac_max)
     	if(nx>=inter_switch) call GS_RHS4

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,25) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return

contains

  !--------------------------------------------------

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
          gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

    ! write results for debugging

    open(11, file='Ptot_RHS.plt')
    open(12, file='gradPtot_RHS.plt')
    open(13, file='gradpsi_RHS.plt')
    open(14, file='bphi2_RHS.plt')
    open(15, file='vphi2_RHS.plt')

    write(11,*)'TITLE="solution of GS equation with flow"'
    write(11,*)'Variables =" X ","Y", "Ptot"'
    write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(12,*)'TITLE="solution of GS equation with flow"'
    write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
    write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(13,*)'TITLE="solution of GS equation with flow"'
    write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
    write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(14,*)'TITLE="solution of GS equation with flow"'
    write(14,*)'Variables =" X ","Y", "bphi2"'
    write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(15,*)'TITLE="solution of GS equation with flow"'
    write(15,*)'Variables =" X ","Y", "vphi2"'
    write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


    do j=1,nz

       z = z_coord(j)

       do i=1,nx

          x = x_coord(i)

          write(11,88) x,z,Ptot(i,j)
          write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
          write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
          write(14,88) x,z,bphi2_RHS(i,j)
          write(15,88) x,z,vphi2_RHS(i,j)

       end do
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS2
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

			if(abs(grad_dir(1,i,j))+abs(grad_dir(2,i,j))/=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
          gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

    ! write results for debugging

    open(11, file='Ptot_RHS.plt')
    open(12, file='gradPtot_RHS.plt')
    open(13, file='gradpsi_RHS.plt')
    open(14, file='bphi2_RHS.plt')
    open(15, file='vphi2_RHS.plt')

    write(11,*)'TITLE="solution of GS equation with flow"'
    write(11,*)'Variables =" X ","Y", "Ptot"'
    write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(12,*)'TITLE="solution of GS equation with flow"'
    write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
    write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(13,*)'TITLE="solution of GS equation with flow"'
    write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
    write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(14,*)'TITLE="solution of GS equation with flow"'
    write(14,*)'Variables =" X ","Y", "bphi2"'
    write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(15,*)'TITLE="solution of GS equation with flow"'
    write(15,*)'Variables =" X ","Y", "vphi2"'
    write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


    do j=1,nz

       z = z_coord(j)

       do i=1,nx

          x = x_coord(i)

          write(11,88) x,z,Ptot(i,j)
          write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
          write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
          write(14,88) x,z,bphi2_RHS(i,j)
          write(15,88) x,z,vphi2_RHS(i,j)

       end do
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS2

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS3
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,25)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS3




  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS4
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then

             if((abs((Ptot(i+1,j)-Ptot(i-1,j))/tdx)>abs((Ptot(i+1,j)-Ptot(i,j))/dx)).and.(abs((Ptot(i+1,j)-Ptot(i-1,j))/tdx)>abs((Ptot(i,j)-Ptot(i-1,j))/dx))) then
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
			elseif(abs((Ptot(i+1,j)-Ptot(i,j))/dx)>abs((Ptot(i,j)-Ptot(i-1,j))/dx)) then
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
			else
				gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
			endif

          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then

             if((abs((Ptot(i,j+1)-Ptot(i,j-1))/tdz)>abs((Ptot(i,j+1)-Ptot(i,j))/dz)).and.(abs((Ptot(i,j+1)-Ptot(i,j-1))/tdz)>abs((Ptot(i,j)-Ptot(i,j-1))/dz))) then
				gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
			elseif(abs((Ptot(i,j+1)-Ptot(i,j))/dz)>abs((Ptot(i,j)-Ptot(i,j-1))/dz)) then
				gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
			else
				gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
			endif
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,25)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS4


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine set_directions
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

    grad_dir = -10
	dir_filter = -10

    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R index
          if(sort_grid(i+1,j,0)>0) then

			dir_filter(i,j) = 0

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = 0
             else
                grad_dir(1,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

          ! Z index
          if(sort_grid(i,j+1,0)>0) then

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = 0
             else
                grad_dir(2,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

		if(dir_filter(i,j)==0) then

			if((i-dir_switch<1).or.(i+dir_switch>nx).or.(j-dir_switch<1).or.(j+dir_switch>nz)) then
				dir_filter(i,j) = -1
			else
				do ii = -dir_switch, dir_switch
					do jj = -dir_switch, dir_switch
						if(sort_grid(i+ii,j+jj,0)<=0) then
							dir_filter(i,j) = -1
						endif
					enddo
				enddo
			endif
		endif

       enddo
    enddo

	open(11,file='directions.plt')

	write(11,*)'TITLE="solution of GS equation with flow"'
	write(11,*)'Variables =" X ","Y", "grad_x", "grad_z", "filter"'
	write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		do j=1,nz
		   do i=1,nx

			  write(11,77) x_coord(i), z_coord(j), grad_dir(1,i,j), grad_dir(2,i,j), dir_filter(i,j)

			enddo
		enddo

	continue

77  format(2(e12.6, 3x), 3(I3, 3x))

  end subroutine set_directions


end subroutine ngs_solve_grad_psi
















! ------------------------------------------------------------------
! ------------------------------------------------------------------


! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this version solves the GS equation with grad Ptot in the RHS
! gradients in term0 are calculated in a consistent way
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_grad_psi_consistent(psi_in,rho,residual,b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen record the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(1:nx,1:nz) :: vphi2_RHS, bphi2_RHS, gradpsi_iph_RHS, gradpsi_jph_RHS
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  real(kind=dkind), dimension(1:2,1:nx,1:nz) :: gradPtot_RHS
  integer, dimension(1:2,1:nx,1:nz) :: grad_dir
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real(kind=dkind) :: gpsi_switch, gpsi_frac, gpsi_frac_max
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v, drho_ds


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

        call step_output(nx,psi_in(:,:),rho,residual)


	gpsi_frac_max = lambda0_fix
 	gpsi_frac = 1.d-1
  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  call set_directions
	  call GS_RHS6


  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j,0)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo

     if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)

     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j,0)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)
              omegac = omegaofpsi(psi0)

              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then

                 ! get |grad psi| in surrounding points (using old iteration values)

                gpsirx = gradpsi_iph_RHS(i,j)
				gpsilx = gradpsi_iph_RHS(i-1,j)

                gpsirz = gradpsi_jph_RHS(i,j)
				gpsilz = gradpsi_jph_RHS(i,j-1)


!                 term0 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
!                      x * ( (psir-psi0)/(gpsirx*(x+dx/2.d0)) - (psi0-psil)/(gpsilx*(x-dx/2.d0)) ) / dx2 +  &
!                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
!                      )


                 term0 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res = term0+term1+term2+term3


!                 den = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
!                      x*( -1.d0/(gpsirx*(x+dx/2.d0)) - 1.d0/(gpsilx*(x-dx/2.d0)) ) / dx2 +  &
!                      ( -1.d0/gpsirz -1.d0/gpsilz ) / dz2  &
!                      )

                 den = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz - 1.d0/gpsilz ) / dz2  &
                      )

				res = res/gradpsi_RHS(0,i,j)**2/x
				den = den/gradpsi_RHS(0,i,j)**2/x

				continue

              else

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
!                 rhorx = (rho(i,j) + rho(i+1,j))/2.d0
                rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
!                 rholx = (rho(i,j)+rho(i-1,j))/2.d0
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
!                 rhorz = (rho(i,j) + rho(i,j+1))/2.d0
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
!                 rholz = (rho(i,j) +rho(i,j-1))/2.d0
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res = term0+term1+term2+term3+term4-term5


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue

              endif


              !res = res*mu_mag


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,200) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					if(nx>=inter_switch) then
						residual(i,j) = res/den
					else
						residual(i,j) = res
					endif
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
				if(nx>=inter_switch) then
	              anorm = anorm + dabs(res/den)
				else
	              anorm = anorm + dabs(res)
				endif

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
	 if(nx>=inter_switch) gpsi_frac = min(gpsi_frac*1.1,gpsi_frac_max)
     	if(nx>=inter_switch) call GS_RHS6

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

        !!$ call step_output(nx,psi_in(:,:),rho,residual)

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,200) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return

contains

  !--------------------------------------------------

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
          gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

    ! write results for debugging

    open(11, file='Ptot_RHS.plt')
    open(12, file='gradPtot_RHS.plt')
    open(13, file='gradpsi_RHS.plt')
    open(14, file='bphi2_RHS.plt')
    open(15, file='vphi2_RHS.plt')

    write(11,*)'TITLE="solution of GS equation with flow"'
    write(11,*)'Variables =" X ","Y", "Ptot"'
    write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(12,*)'TITLE="solution of GS equation with flow"'
    write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
    write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(13,*)'TITLE="solution of GS equation with flow"'
    write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
    write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(14,*)'TITLE="solution of GS equation with flow"'
    write(14,*)'Variables =" X ","Y", "bphi2"'
    write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(15,*)'TITLE="solution of GS equation with flow"'
    write(15,*)'Variables =" X ","Y", "vphi2"'
    write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


    do j=1,nz

       z = z_coord(j)

       do i=1,nx

          x = x_coord(i)

          write(11,88) x,z,Ptot(i,j)
          write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
          write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
          write(14,88) x,z,bphi2_RHS(i,j)
          write(15,88) x,z,vphi2_RHS(i,j)

       end do
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS2
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

			if(abs(grad_dir(1,i,j))+abs(grad_dir(2,i,j))/=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
          gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

    ! write results for debugging

    open(11, file='Ptot_RHS.plt')
    open(12, file='gradPtot_RHS.plt')
    open(13, file='gradpsi_RHS.plt')
    open(14, file='bphi2_RHS.plt')
    open(15, file='vphi2_RHS.plt')

    write(11,*)'TITLE="solution of GS equation with flow"'
    write(11,*)'Variables =" X ","Y", "Ptot"'
    write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(12,*)'TITLE="solution of GS equation with flow"'
    write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
    write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(13,*)'TITLE="solution of GS equation with flow"'
    write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
    write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(14,*)'TITLE="solution of GS equation with flow"'
    write(14,*)'Variables =" X ","Y", "bphi2"'
    write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

    write(15,*)'TITLE="solution of GS equation with flow"'
    write(15,*)'Variables =" X ","Y", "vphi2"'
    write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


    do j=1,nz

       z = z_coord(j)

       do i=1,nx

          x = x_coord(i)

          write(11,88) x,z,Ptot(i,j)
          write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
          write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
          write(14,88) x,z,bphi2_RHS(i,j)
          write(15,88) x,z,vphi2_RHS(i,j)

       end do
    end do

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS2

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS3
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,25)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS3




  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS4
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then

             if((abs((Ptot(i+1,j)-Ptot(i-1,j))/tdx)>abs((Ptot(i+1,j)-Ptot(i,j))/dx)).and.(abs((Ptot(i+1,j)-Ptot(i-1,j))/tdx)>abs((Ptot(i,j)-Ptot(i-1,j))/dx))) then
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
			elseif(abs((Ptot(i+1,j)-Ptot(i,j))/dx)>abs((Ptot(i,j)-Ptot(i-1,j))/dx)) then
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
			else
				gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
			endif

          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then

             if((abs((Ptot(i,j+1)-Ptot(i,j-1))/tdz)>abs((Ptot(i,j+1)-Ptot(i,j))/dz)).and.(abs((Ptot(i,j+1)-Ptot(i,j-1))/tdz)>abs((Ptot(i,j)-Ptot(i,j-1))/dz))) then
				gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
			elseif(abs((Ptot(i,j+1)-Ptot(i,j))/dz)>abs((Ptot(i,j)-Ptot(i,j-1))/dz)) then
				gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
			else
				gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
			endif
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,25)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS4


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS5
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2

    real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)
          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt( ((psi_in(i+1,j)-psi_in(i,j))/dx)**2 +  &
											( ( (psi_in(i+1,j+1)+psi_in(i,j+1))/2.d0 - (psi_in(i+1,j-1)+psi_in(i,j-1))/2.d0)/tdz)**2)

		  gradpsi_jph_RHS(i,j) = sqrt(  ( ( (psi_in(i+1,j+1)+psi_in(i+1,j))/2.d0 - (psi_in(i-1,j+1)+psi_in(i-1,j))/2.d0)/tdx)**2) +  &
											((psi_in(i,j+1)-psi_in(i,j))/dz)**2

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + (gradpsi_RHS(0,i,j)/x)**2) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,25)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS5



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS6
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,200)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS6


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine set_directions
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: ii, jj

    grad_dir = -10
	dir_filter = -10

    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R index
          if(sort_grid(i+1,j,0)>0) then

			dir_filter(i,j) = 0

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = 0
             else
                grad_dir(1,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

          ! Z index
          if(sort_grid(i,j+1,0)>0) then

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = 0
             else
                grad_dir(2,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

		if(dir_filter(i,j)==0) then

			if((i-dir_switch<1).or.(i+dir_switch>nx).or.(j-dir_switch<1).or.(j+dir_switch>nz)) then
				dir_filter(i,j) = -1
			else
				do ii = -dir_switch, dir_switch
					do jj = -dir_switch, dir_switch
						if(sort_grid(i+ii,j+jj,0)<=0) then
							dir_filter(i,j) = -1
						endif
					enddo
				enddo
			endif
		endif

       enddo
    enddo

	open(11,file='directions.plt')

	write(11,*)'TITLE="solution of GS equation with flow"'
	write(11,*)'Variables =" X ","Y", "grad_x", "grad_z", "filter"'
	write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		do j=1,nz
		   do i=1,nx

			  write(11,77) x_coord(i), z_coord(j), grad_dir(1,i,j), grad_dir(2,i,j), dir_filter(i,j)

			enddo
		enddo

	continue

77  format(2(e12.6, 3x), 3(I3, 3x))

  end subroutine set_directions


end subroutine ngs_solve_grad_psi_consistent




! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this version solves the GS equation with grad Ptot in the RHS
! gradients in term0 are calculated in a consistent way.
! The two forms of the GS equation are mixed to get a smooth transition
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_grad_psi_consistent_smooth_transition(psi_in,rho,residual,  &
						b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  real (kind=dkind) :: res_1, den_1
  real (kind=dkind) :: res_2, den_2
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen record the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0
  real(kind=dkind) :: term1_1, term2_1, term3_1, term0_1
  real(kind=dkind) :: term1_2, term2_2, term3_2, term4_2, term5_2, term0_2
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(1:nx,1:nz) :: vphi2_RHS, bphi2_RHS, gradpsi_iph_RHS, gradpsi_jph_RHS
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  real(kind=dkind), dimension(1:2,1:nx,1:nz) :: gradPtot_RHS
  integer, dimension(1:2,1:nx,1:nz) :: grad_dir
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real(kind=dkind) :: gpsi_switch, gpsi_frac, gpsi_frac_max
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v, drho_ds
  real(kind=dkind) :: psi_distance ! how far to go from psi_degen in using the Ptot formulation
  real(kind=dkind) :: smooth_fact
  integer :: k_save = 100


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

        call step_output(nx,psi_in(:,:),rho,residual)


  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  psi_distance = 0.d0
	  call set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,k,k_save)
	  call GS_RHS8


  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j,0)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo

 !    if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)


     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j,0)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)
              omegac = omegaofpsi(psi0)

              if(dir_filter(i,j)==1) then
			  ! Ptot zone

                 ! get |grad psi| in surrounding points (using old iteration values)

                gpsirx = gradpsi_iph_RHS(i,j)
				gpsilx = gradpsi_iph_RHS(i-1,j)

                gpsirz = gradpsi_jph_RHS(i,j)
				gpsilz = gradpsi_jph_RHS(i,j-1)

                 term0 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res = term0+term1+term2+term3

                 den = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz - 1.d0/gpsilz ) / dz2  &
                      )

				res = res/gradpsi_RHS(0,i,j)**2/x
				den = den/gradpsi_RHS(0,i,j)**2/x

				continue

              elseif(dir_filter(i,j)==0) then
			  ! mixed zone

				smooth_fact = abs(psi_in(i,j)-psi_degen)/psi_distance - 1.d0

                 ! get |grad psi| in surrounding points (using old iteration values)

                gpsirx = gradpsi_iph_RHS(i,j)
				gpsilx = gradpsi_iph_RHS(i-1,j)

                gpsirz = gradpsi_jph_RHS(i,j)
				gpsilz = gradpsi_jph_RHS(i,j-1)


                 term0_1 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1_1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2_1 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3_1 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res_1 = term0_1+term1_1+term2_1+term3_1


                 den_1 = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz - 1.d0/gpsilz ) / dz2  &
                      )

				!--------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0_2 = (1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1_2 = b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2_2 = x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3_2 = by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4_2 = rhoc*dhdpsi(psi0) 
                 term5_2 = rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res_2 = term0_2 + term1_2 + term2_2 + term3_2 + term4_2 - term5_2


                 den_2 = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

				res = (1.d0-smooth_fact)*res_1/den_1 + smooth_fact*res_2/den_2
				den = 1.d0

				continue

              else
			  ! old formulation zone

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res = term0+term1+term2+term3+term4-term5


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue

              endif


              !res = res*mu_mag


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					if(nx>=inter_switch) then
						residual(i,j) = res/den
					else
						residual(i,j) = res
					endif
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
				if(nx>=inter_switch) then
	              anorm = anorm + dabs(res/den)
				else
	              anorm = anorm + dabs(res)
				endif

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
     	if(nx>=inter_switch) then
			psi_distance = min(psic/lambda0_fix/2.d0, psi_distance + 0.05d0*psic/lambda0_fix/2.d0)
			call set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,k,k_save)
			call GS_RHS8
		endif

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
	if((nx>=inter_switch).and.(Broot==0)) then
     call update_rho_one_sided(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
    else
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
	endif

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

        !!$ call step_output(nx,psi_in(:,:),rho,residual)

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,k_save) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return

contains

  !--------------------------------------------------


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS6_5
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS6_5



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS7
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          if((abs(psi_in(i,j)-psi_degen)<psic/10.).and.(psi_in(i,j)/=psi_degen)) then
          ! check for one-sided derivatives

			if((psi_in(i+1,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) > 0.d0) then

				! no crossing in x, use centered differences
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx

			elseif((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i,j))/dx

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(1,i,j) = (psi_in(i,j)-psi_in(i-1,j))/dx

			else

				print*, 'error in GS_RHS7, x derivative,', i, j

			endif

			if((psi_in(i,j+1)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! no crossing in z, use centered differences
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz

			elseif((psi_in(i,j+1)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j))/dz

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j)-psi_in(i,j-1))/dz

			else

				print*, 'error in GS_RHS7, z derivative,', i, j

			endif

		else

          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

		endif

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

!         if(abs(psi_in(i,j)-psi_degen)<psic/10.) then
!          ! check for one-sided derivatives
!
!			if((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) > 0.d0) then
!
!				! no crossing in x, use centered differences
!				gradpsi_iph_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
!
!			elseif((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) > 0.d0) then
!
!				! use one-sided right derivative
!				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i,j))/dx
!
!			elseif((psi_in(i,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) > 0.d0) then
!
!				! use one-sided left derivative
!				gradpsi_RHS(1,i,j) = (psi_in(i,j)-psi_in(i-1,j))/dx
!
!			else
!
!				print*, 'error in GS_RHS7, x derivative,', i, j
!
!			endif
!
!			if((psi_in(i,j+1)-psi_degen)*(psi_in(i,j-1)-psi_degen) > 0.d0) then
!
!				! no crossing in z, use centered differences
!				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
!
!			elseif((psi_in(i,j+1)-psi_degen)*(psi_in(i,j)-psi_degen) > 0.d0) then
!
!				! use one-sided right derivative
!				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j))/dz
!
!			elseif((psi_in(i,j)-psi_degen)*(psi_in(i,j-1)-psi_degen) > 0.d0) then
!
!				! use one-sided left derivative
!				gradpsi_RHS(2,i,j) = (psi_in(i,j)-psi_in(i,j-1))/dz
!
!			else
!
!				print*, 'error in GS_RHS7, z derivative,', i, j
!
!			endif
!
!		else

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

!		endif


       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS7


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS7_5
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          if((abs(psi_in(i,j)-psi_degen)<psic/10.).and.(psi_in(i,j)/=psi_degen)) then
          ! check for one-sided derivatives

			if((psi_in(i+1,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) > 0.d0) then

				! no crossing in x, use centered differences
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx

			elseif((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i,j))/dx

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(1,i,j) = (psi_in(i,j)-psi_in(i-1,j))/dx

			else

				print*, 'error in GS_RHS7_5, x derivative,', i, j

			endif

			if((psi_in(i,j+1)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! no crossing in z, use centered differences
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz

			elseif((psi_in(i,j+1)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j))/dz

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j)-psi_in(i,j-1))/dz

			else

				print*, 'error in GS_RHS7_5, z derivative,', i, j

			endif

		else

          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

		endif

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

!		endif


       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle
          if((abs(psi_in(i,j)-psi_degen)<psic/10.).and.(psi_in(i,j)/=psi_degen)) then
          ! check for one-sided derivatives

			if((psi_in(i+1,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) > 0.d0) then

				! no crossing in x, use centered differences
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx

			elseif((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx

			else

				print*, 'error in GS_RHS7, x derivative,', i, j

			endif

			if((psi_in(i,j+1)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! no crossing in z, use centered differences
	             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz

			elseif((psi_in(i,j+1)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
	             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz

			else

				print*, 'error in GS_RHS7, z derivative,', i, j

			endif

		else

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

		endif

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS7_5

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS8
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations
  ! grad Ptot is also calculated with interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS, bscoef_Ptot_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          if((abs(psi_in(i,j)-psi_degen)<psic/10.).and.(psi_in(i,j)/=psi_degen)) then
          ! check for one-sided derivatives

			if((psi_in(i+1,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) > 0.d0) then

				! no crossing in x, use centered differences
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx

			elseif((psi_in(i+1,j)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i,j))/dx

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i-1,j)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(1,i,j) = (psi_in(i,j)-psi_in(i-1,j))/dx

			else

				print*, 'error in GS_RHS7, x derivative,', i, j

			endif

			if((psi_in(i,j+1)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! no crossing in z, use centered differences
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz

			elseif((psi_in(i,j+1)-psi_degen)*(psi_in(i,j)-psi_degen) >= 0.d0) then

				! use one-sided right derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j))/dz

			elseif((psi_in(i,j)-psi_degen)*(psi_in(i,j-1)-psi_degen) >= 0.d0) then

				! use one-sided left derivative
				gradpsi_RHS(2,i,j) = (psi_in(i,j)-psi_in(i,j-1))/dz

			else

				print*, 'error in GS_RHS7, z derivative,', i, j

			endif

		else

          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

		endif

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

!		endif


       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

	! now interpolate Ptot

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),Ptot,nx,  &
				ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_Ptot_RHS)


    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          gradPtot_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_Ptot_RHS)
          gradPtot_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_Ptot_RHS)

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS8


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,  &
!					k, k_save)
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!    integer :: nx, nz
!	real (kind=dkind), dimension(1:nx,1:nz) :: psi_in
!	integer, dimension(1:2,1:nx,1:nz) :: grad_dir
!    integer :: dir_filter(1:nx,1:nz)
!	real(kind=dkind) :: psi_degen,psi_distance
!	integer :: k, k_save

!	integer :: i, j
	integer :: ii, jj

    grad_dir = -10
	dir_filter = -10

	if(nx<inter_switch) return

    do j = 1, nz
       loop_grid_i: do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R index
          if(sort_grid(i+1,j,0)>0) then

			dir_filter(i,j) = 1

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = 0
             else
                grad_dir(1,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

          ! Z index
          if(sort_grid(i,j+1,0)>0) then

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = 0
             else
                grad_dir(2,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

		if((dir_filter(i,j)==1).or.(dir_filter(i,j)==0)) then

			if((i-dir_switch<1).or.(i+dir_switch>nx).or.(j-dir_switch<1).or.(j+dir_switch>nz)) then
				dir_filter(i,j) = -1
				cycle loop_grid_i
			else
				do ii = -dir_switch, dir_switch
					do jj = -dir_switch, dir_switch
						if(sort_grid(i+ii,j+jj,0)<=0) then
							dir_filter(i,j) = -1
							cycle loop_grid_i
						endif
					enddo
				enddo
			endif

			if(abs(psi_in(i,j)-psi_degen)<psi_distance) then
				dir_filter(i,j) = 1
			elseif(abs(psi_in(i,j)-psi_degen)<2.d0*psi_distance) then
				dir_filter(i,j) = 0
			else
				dir_filter(i,j) = -1
			endif

		endif


       enddo loop_grid_i
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

		open(11,file='directions.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "grad_x", "grad_z", "filter"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		do j=1,nz
		   do i=1,nx

			  write(11,77) x_coord(i), z_coord(j), grad_dir(1,i,j), grad_dir(2,i,j), dir_filter(i,j)

			enddo
		enddo

	endif

	continue

77  format(2(e12.6, 3x), 3(I3, 3x))

  end subroutine set_directions_4


end subroutine ngs_solve_grad_psi_consistent_smooth_transition





! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! this version solves the GS equation with grad Ptot in the RHS
! the Bernoulli equation is modified to have a continuous solution
! a smooth transition between the two zones is still used
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_grad_psi_gauss(psi_in,rho,residual,  &
						b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  real (kind=dkind) :: res_1, den_1
  real (kind=dkind) :: res_2, den_2
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen record the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0
  real(kind=dkind) :: term1_1, term2_1, term3_1, term0_1
  real(kind=dkind) :: term1_2, term2_2, term3_2, term4_2, term5_2, term0_2
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(1:nx,1:nz) :: vphi2_RHS, bphi2_RHS, gradpsi_iph_RHS, gradpsi_jph_RHS
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  real(kind=dkind), dimension(1:2,1:nx,1:nz) :: gradPtot_RHS
  integer, dimension(1:2,1:nx,1:nz) :: grad_dir
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real(kind=dkind) :: gpsi_switch, gpsi_frac, gpsi_frac_max
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v, drho_ds
  real(kind=dkind) :: psi_distance ! how far to go from psi_degen in using the Ptot formulation
  real(kind=dkind) :: smooth_fact
  integer :: k_save = 100


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

        call step_output(nx,psi_in(:,:),rho,residual)


  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  psi_distance = 0.d0
	  call set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,k,k_save)
	  call GS_RHS_G1


  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j,0)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo

 !    if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)


     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j,0)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)
              omegac = omegaofpsi(psi0)

              if(dir_filter(i,j)==1) then
			  ! Ptot zone

                 ! get |grad psi| in surrounding points (using old iteration values)

                gpsirx = gradpsi_iph_RHS(i,j)
				gpsilx = gradpsi_iph_RHS(i-1,j)

                gpsirz = gradpsi_jph_RHS(i,j)
				gpsilz = gradpsi_jph_RHS(i,j-1)

                 term0 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res = term0+term1+term2+term3

                 den = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz - 1.d0/gpsilz ) / dz2  &
                      )

				res = res/gradpsi_RHS(0,i,j)**2/x
				den = den/gradpsi_RHS(0,i,j)**2/x

				continue

              elseif(dir_filter(i,j)==0) then
			  ! mixed zone

				smooth_fact = abs(psi_in(i,j)-psi_degen)/psi_distance - 1.d0

                 ! get |grad psi| in surrounding points (using old iteration values)

                gpsirx = gradpsi_iph_RHS(i,j)
				gpsilx = gradpsi_iph_RHS(i-1,j)

                gpsirz = gradpsi_jph_RHS(i,j)
				gpsilz = gradpsi_jph_RHS(i,j-1)


                 term0_1 = (1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( (psir-psi0)/gpsirx - (psi0-psil)/gpsilx ) / dx2 +  &
                      ( (psiu-psi0)/gpsirz - (psi0-psid)/gpsilz ) / dz2  &
                      )

                 term1_1 = - rhoc * x * vphi2_RHS(i,j) * gradpsi_RHS(1,i,j)

                 term2_1 = x**2 * (  &
                      gradpsi_RHS(1,i,j)*gradPtot_RHS(1,i,j) +  &
                      gradpsi_RHS(2,i,j)*gradPtot_RHS(2,i,j) )

                 term3_1 =  bphi2_RHS(i,j) * x * gradpsi_RHS(1,i,j) 


                 res_1 = term0_1+term1_1+term2_1+term3_1


                 den_1 = 		(1.d0-phic**2/rhoc)*gradpsi_RHS(0,i,j)**3 * (  &
                      ( -1.d0/gpsirx - 1.d0/gpsilx ) / dx2 +  &
                      ( -1.d0/gpsirz - 1.d0/gpsilz ) / dz2  &
                      )

				!--------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0_2 = (1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1_2 = b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2_2 = x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3_2 = by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4_2 = rhoc*dhdpsi(psi0) 
                 term5_2 = rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res_2 = term0_2 + term1_2 + term2_2 + term3_2 + term4_2 - term5_2


                 den_2 = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

				res = (1.d0-smooth_fact)*res_1/den_1 + smooth_fact*res_2/den_2
				den = 1.d0

				continue

              else
			  ! old formulation zone

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)


                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

                 res = term0+term1+term2+term3+term4-term5


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue

              endif


              !res = res*mu_mag


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					if(nx>=inter_switch) then
						residual(i,j) = res/den
					else
						residual(i,j) = res
					endif
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
				if(nx>=inter_switch) then
	              anorm = anorm + dabs(res/den)
				else
	              anorm = anorm + dabs(res)
				endif

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
     	if(nx>=inter_switch) then
			psi_distance = min(psic/lambda0_fix/2.d0, psi_distance + 0.05d0*psic/lambda0_fix/2.d0)
			call set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,k,k_save)
			call GS_RHS_G1
		endif

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
	if((nx>=inter_switch).and.(Broot==0)) then
     call update_rho_gauss(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
    else
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
	endif

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

        !!$ call step_output(nx,psi_in(:,:),rho,residual)

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,k_save) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return

contains

  !--------------------------------------------------


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_G0
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle


          gradpsi_RHS(1,i,j) = 	dbs2dr(1,0,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)
          gradpsi_RHS(2,i,j) = 	dbs2dr(0,1,x_coord(i),z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)

          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_G0


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_G1
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_psi_RHS
	integer, parameter :: ord_RHS = 2
	real(kind=dkind) :: xknot_RHS(1:nx+ord_RHS), zknot_RHS(1:nz+ord_RHS)
	real(kind=dkind) :: Ptot(1:nx,1:nz)
    real(kind=dkind) :: tdx, tdz
    real(kind=dkind) :: gpsi_max
    real(kind=dkind) :: z

	print*, ord_RHS

    tdx = dx*2.d0
    tdz = dz*2.d0

    Ptot = 0.d0
    gradPtot_RHS = 0.d0
    gradpsi_RHS = 0.d0
	gradpsi_iph_RHS = 0.d0
	gradpsi_jph_RHS = 0.d0
    bphi2_RHS = 0.d0
    vphi2_RHS = 0.d0

    gpsi_max = 0.d0

	call DBSNAK(nx,x_coord(1:nx),ord_RHS,xknot_RHS)
	call DBSNAK(nz,z_coord(1:nz),ord_RHS,zknot_RHS)

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord(1:nz),psi_in,nx,  &
					ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,bscoef_psi_RHS)

    ! first fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)
		  if(dir_filter(i,j)==0) gpsi_max = max(gradpsi_RHS(0,i,j),gpsi_max )

		  gradpsi_iph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i)+dx/2.d0,z_coord(j),ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

		  gradpsi_jph_RHS(i,j) = sqrt(dbs2dr(1,0,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 +  &
													dbs2dr(0,1,x_coord(i),z_coord(j)+dz/2.d0,ord_RHS,ord_RHS,xknot_RHS,zknot_RHS,nx,nz,bscoef_psi_RHS)**2 )

       enddo
    enddo

    gpsi_switch = gpsi_max / gpsi_frac

    ! next, everything else except grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          psi0 = psi_in(i,j)
          x = x_coord(i)

          ! Calculate rho, phi and omega at the current location
          rhoc = rho(i,j)
          phic = phiofpsi(psi0)
          omegac = omegaofpsi(psi0)

          by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
               (1.0d0 - phic*phic/rhoc)
          bphi2_RHS(i,j) = by**2

          vphi2_RHS(i,j) = (x*omegac + (phic/rhoc)*by/sqrt(mu_mag))**2

          Ptot(i,j) = Sofpsi(psi0)*rhoc**gamma + (bphi2_RHS(i,j)/2.d0 + ((gradpsi_RHS(0,i,j)/x)**2)/2.d0) / mu_mag

       enddo
    enddo

    ! last, get grad Ptot
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R gradient
          if(grad_dir(1,i,j) == 0) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i-1,j))/tdx
          elseif(grad_dir(1,i,j) == 1) then
             gradPtot_RHS(1,i,j) = (Ptot(i+1,j)-Ptot(i,j))/dx
          elseif(grad_dir(1,i,j) == -1) then
             gradPtot_RHS(1,i,j) = (Ptot(i,j)-Ptot(i-1,j))/dx
          endif

          ! Z gradient
          if(grad_dir(2,i,j) == 0) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j-1))/tdz
          elseif(grad_dir(2,i,j) == 1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j+1)-Ptot(i,j))/dz
          elseif(grad_dir(2,i,j) == -1) then
             gradPtot_RHS(2,i,j) = (Ptot(i,j)-Ptot(i,j-1))/dz
          endif

       enddo
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

    ! write results for debugging

		open(11, file='Ptot_RHS.plt')
		open(12, file='gradPtot_RHS.plt')
		open(13, file='gradpsi_RHS.plt')
		open(14, file='bphi2_RHS.plt')
		open(15, file='vphi2_RHS.plt')
		open(16, file='gradpsi_half_grid_RHS.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "Ptot"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(12,*)'TITLE="solution of GS equation with flow"'
		write(12,*)'Variables =" X ","Y", "Ptot_x", "Ptot_z"'
		write(12,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(13,*)'TITLE="solution of GS equation with flow"'
		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(14,*)'TITLE="solution of GS equation with flow"'
		write(14,*)'Variables =" X ","Y", "bphi2"'
		write(14,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(15,*)'TITLE="solution of GS equation with flow"'
		write(15,*)'Variables =" X ","Y", "vphi2"'
		write(15,*)'ZONE I=',nx,',J=',nz,',F=Point'

		write(16,*)'TITLE="solution of GS equation with flow"'
		write(16,*)'Variables =" X ","Y", "gpsi_iph", "gpsi_jph"'
		write(16,*)'ZONE I=',nx,',J=',nz,',F=Point'


		do j=1,nz

		   z = z_coord(j)

		   do i=1,nx

			  x = x_coord(i)

			  write(11,88) x,z,Ptot(i,j)
			  write(12,89) x,z,gradPtot_RHS(1,i,j),gradPtot_RHS(2,i,j)
			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
			  write(14,88) x,z,bphi2_RHS(i,j)
			  write(15,88) x,z,vphi2_RHS(i,j)
			  write(16,89) x,z,gradpsi_iph_RHS(i,j),gradpsi_jph_RHS(i,j)

		   end do
		end do

		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)

	endif

88  format(3(e12.6,3x))
89  format(4(e12.6,3x))
90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_G1





  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine set_directions_4 !(nx,nz,psi_in,grad_dir,dir_filter,psi_degen,psi_distance,  &
!					k, k_save)
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!    integer :: nx, nz
!	real (kind=dkind), dimension(1:nx,1:nz) :: psi_in
!	integer, dimension(1:2,1:nx,1:nz) :: grad_dir
!    integer :: dir_filter(1:nx,1:nz)
!	real(kind=dkind) :: psi_degen,psi_distance
!	integer :: k, k_save

!	integer :: i, j
	integer :: ii, jj

    grad_dir = -10
	dir_filter = -10

	if(nx<inter_switch) return

    do j = 1, nz
       loop_grid_i: do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          ! R index
          if(sort_grid(i+1,j,0)>0) then

			dir_filter(i,j) = 1

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = 0
             else
                grad_dir(1,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i-1,j,0)>0) then
                grad_dir(1,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

          ! Z index
          if(sort_grid(i,j+1,0)>0) then

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = 0
             else
                grad_dir(2,i,j) = 1
				dir_filter(i,j) = -1
             endif

          else

             if(sort_grid(i,j-1,0)>0) then
                grad_dir(2,i,j) = -1
				dir_filter(i,j) = -1
             else
                ! this of course should never happen...
                print*, 'error in set_directions'
             endif

          endif

		if((dir_filter(i,j)==1).or.(dir_filter(i,j)==0)) then

			if((i-dir_switch<1).or.(i+dir_switch>nx).or.(j-dir_switch<1).or.(j+dir_switch>nz)) then
				dir_filter(i,j) = -1
				cycle loop_grid_i
			else
				do ii = -dir_switch, dir_switch
					do jj = -dir_switch, dir_switch
						if(sort_grid(i+ii,j+jj,0)<=0) then
							dir_filter(i,j) = -1
							cycle loop_grid_i
						endif
					enddo
				enddo
			endif

			if(abs(psi_in(i,j)-psi_degen)<psi_distance) then
				dir_filter(i,j) = 1
			elseif(abs(psi_in(i,j)-psi_degen)<2.d0*psi_distance) then
				dir_filter(i,j) = 0
			else
				dir_filter(i,j) = -1
			endif

		endif


       enddo loop_grid_i
    enddo

	if((k==1).or.(modulo(k,k_save)==0)) then

		open(11,file='directions.plt')

		write(11,*)'TITLE="solution of GS equation with flow"'
		write(11,*)'Variables =" X ","Y", "grad_x", "grad_z", "filter"'
		write(11,*)'ZONE I=',nx,',J=',nz,',F=Point'

		do j=1,nz
		   do i=1,nx

			  write(11,77) x_coord(i), z_coord(j), grad_dir(1,i,j), grad_dir(2,i,j), dir_filter(i,j)

			enddo
		enddo

	endif

	continue

77  format(2(e12.6, 3x), 3(I3, 3x))

  end subroutine set_directions_4


end subroutine ngs_solve_grad_psi_gauss










! Red-Black Newton-Gauss-Seidel relaxation.  The current
! value of the solution psi[1..nx][1..nz] is updated.
! max_it -> the maximum number of iterations
! eps -> Fraction by which to reduce the residual
! the Bernoulli equation is modified to have a continuous solution
! Hameiri's GS equation (modified with the new Bernoulli) is used everywhere
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine ngs_solve_all_gauss(psi_in,rho,residual,  &
						b_phi,nx,nz,min_it,max_it,eps,flag,in_orp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dx, dz, dx_a, dz_a

  implicit none
  integer, intent(in) :: nx,nz,min_it,max_it,flag
  real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: psi_in,rho,residual,b_phi
  real (kind=dkind), intent(in) :: eps
  ! Input Over Relaxation Parameter
  ! if in_orp <= 0.0d0 then we use Chebyshev Acceleration
  ! else we set orp = in_orp
  real (kind=dkind), intent(in) :: in_orp
  integer ::ipass,i,isw,j,jsw,k,alloc_stat
  real (kind=dkind) :: res, den, x !, xl, xr
  real (kind=dkind) :: res_1, den_1
  real (kind=dkind) :: res_2, den_2
  ! res -> residual
  ! den -> denominator of the Newton-Gauss-Seidel update
  ! x -> Radial position of grid point (i,j)
  real (kind=dkind) :: dx2,dz2,mx
  real (kind=dkind) :: anorm, eps_anorm
  real (kind=dkind) :: rjac ! Spectral Radius of the Jacobi Iteration
  real (kind=dkind) :: orp, std_orp ! over relaxation parameter
  real (kind=dkind) :: rhoc,phic,omegac
  ! by is the phi component of the magnetic field
  real (kind=dkind) :: by
  real (kind=dkind) :: last_mtm ! last mach theta max
  ! minimum delta rho for bernoulli roots, the mean val. and location
  real (kind=dkind) :: min_drho, tmp
  integer :: min_ix, min_iz
  ! The next set of variables are used for the bisection search for
  ! locating the degenerate root for Mach Theta Max
  real (kind=dkind) :: mtm_soln
  ! mtm_acc controls the bisection loop seeking mach theta max
  real (kind=dkind), parameter :: mtm_acc = 1.0d-12
  ! VERY IMPORTANT
  ! psi_degen record the value of psi where the solution to the
  ! Bernoulli equation is approximately degenerate.  It is 
  ! approximately equal to 0.5d0*psic but generally slightly less.
  real (kind=dkind) :: gpsirx, gpsilx, gpsirz, gpsilz
  real(kind=dkind) :: term1, term2, term3, term4, term5, term0, term6
  real(kind=dkind) :: last_anorm
  real(kind=dkind), dimension(0:2,1:nx,1:nz) :: gradpsi_RHS
  integer :: dir_filter(1:nx,1:nz)
  real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
  real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
  real (kind=dkind) :: rhorx,rholx,rhorz,rholz
  real (kind=dkind) :: phirx,philx,phirz,philz
  real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
  real(kind=dkind) :: dpsidx, dpsidz, b_pol_l, b_dot_v, drho_ds
  real(kind=dkind) :: psi_distance ! how far to go from psi_degen in using the Ptot formulation
  real(kind=dkind) :: smooth_fact
  integer :: k_save = 100


  if(psi_degen==0.d0) psi_degen = 0.5d0*psic ! Initialization
  last_mtm = mach_theta_max
  eps_anorm = 0.0d0

  fix_orp0 = fix_orp

  ! The under/over relaxation parameter
  ! Standard Case:
  if(in_orp <= 0.0d0) then
     ! Note that the parameters below are tuned for starting with 
     ! a 17 x 17 grid.
     orp = 1.0d0
     if(nx <= 5) then
        orp = 0.5d0
        rjac = 0.75d0
     else if(nx <= 9) then
        orp = 0.5d0
        rjac = 0.9d0
     else if(nx <= 17) then
        rjac = 0.98d0
     else if(nx <= 33) then
        rjac = 0.995d0
     else if(nx <= 65) then
        rjac = 0.9972d0
     else
        ! It would appear that 129 is about where the solution begins to converge
        rjac = 0.0d0
     end if
     std_orp = orp
  else
     print *, "Input Over Relaxation Parameter =",in_orp
     orp = in_orp
     std_orp = in_orp
  endif

  if (accelerate) then
     continue
  else
     orp = fix_orp
  endif


  dx2 = dx*dx
  dz2 = dz*dz
  
  if(allocated(fmax_2D)) then
!  	if(size(fmax_2D,1)==nx) then
  		continue
 ! 	else
	  	deallocate(fmax_2D)
  !	endif
  endif

!  if(allocated(fmax_2D)) then
!	continue
 ! else
  	allocate(fmax_2D(1:nx,1:nz))
  	fmax_2D = 0.d0
  !endif


  ! Update rho before "relaxing psi" but do not seek mach theta max
  call update_rho(psi_in,rho,b_phi,nx,nz,0,mtm_acc,min_drho,  &
       min_ix,min_iz)

        call step_output(nx,psi_in(:,:),rho,residual)


  ! set up gradient directions (fixed on the grid) and initial RHS (this is not going to be red-black)
	  call GS_RHS_all_Gauss



  ! Iterate until convergence
  do k=1, max_it
     !     "Update Psi"

     !		call cpu_time(tic)

     ! Reset the norm and max of Psi
     mx = 0.0d0
     anorm = 0.0d0

     do j=2,nz-1
        do i=2,nx-1
           if(sort_grid(i,j,0)>0) psi_in(i,j) = dabs(psi_in(i,j))	!	segno
        enddo
     enddo


 !    if(nx>=inter_switch) call step_output(nx,psi_in(:,:),rho,residual)

	if(delta_Bern_fact==0.d0) delta_Bern_fact = 1.d0

	delta_Bern = delta_Bern_fact*(nx-min_ix)



     jsw = 1
     do ipass = 1, 2
        isw=jsw;
        do j=2, nz-1
           do i=isw+1, nx-1, 2
              ! "Update Psi(",i,",",j,")"
              ! Only solve the inner region problem

              if(sort_grid(i,j,0)<=0) cycle

              ! set up local psi values

              psil = 0.d0
              psir = 0.d0
              psiu = 0.d0
              psid = 0.d0

              psi0 = psi_in(i,j)
              if(i>1) psil = psi_in(i-1,j)
              if(i<nx) psir = psi_in(i+1,j)
              if(j>1) psid = psi_in(i,j-1)
              if(j<nz) psiu = psi_in(i,j+1)


              x = x_coord(i)
!              xl = (x_coord(i)+x_coord(i-1))/2.d0
 !             xr = (x_coord(i+1)+x_coord(i))/2.d0

              ! Calculate rho, phi and omega at the current location
              rhoc = rho(i,j)
              phic = phiofpsi(psi0)
              omegac = omegaofpsi(psi0)

			  ! old formulation with new Bernoulli

                drho_ds = van_Leer_slope_new(rho(i-1,j),rho(i,j),rho(i+1,j),dx_a(i-1),dx_a(i))

                 ! Right x
                 phirx = phiofpsi(0.5d0*(psir + psi0))
                 rhorx = rho(i,j) + 0.5d0*dx*drho_ds
                 ma2rx = phirx*phirx/rhorx

                 ! Left x
                 philx = phiofpsi(0.5d0*(psil + psi0))
                 rholx = rho(i,j) - 0.5d0*dx*drho_ds
                 ma2lx = philx*philx/rholx

                 ! -----------------------------------------------------

                drho_ds = van_Leer_slope_new(rho(i,j-1),rho(i,j),rho(i,j+1),dz_a(j-1),dz_a(j))

                 ! Right z
                 phirz = phiofpsi(0.5d0*(psi0 + psiu))
                 rhorz = rho(i,j) + 0.5d0*dz*drho_ds
                 ma2rz = phirz*phirz/rhorz

                 ! Left z
                 philz = phiofpsi(0.5d0*(psi0 + psid))
                 rholz = rho(i,j) - 0.5d0*dz*drho_ds
                 ma2lz = philz*philz/rholz

                 ! -----------------------------------------------------

                 psi_around(-1,0) = psil
                 psi_around(1,0) = psir
                 psi_around(0,0) = psi0
                 psi_around(0,-1) = psid
                 psi_around(0,1) = psiu

                 by = dsqrt(mu_mag)*(iofpsi(psi0)/x + x*phic*omegac)/ &
                      (1.0d0 - phic*phic/rhoc)

                 call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)
                 b_pol_l =  (dsqrt(dpsidx**2+dpsidz**2) / x)

                 b_dot_v = x*omegac*by + (phic/rhoc)*( by*by + b_pol_l*b_pol_l )/sqrt(mu_mag)



                 term0 =(1.d0/mu_mag)*(( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &
                      + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
                      /(x*dx2) &
                      + ( (1.0d0 - ma2rz)*(psiu-psi0) &
                      + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )


                 term1= b_dot_v*dphidpsi(psi0)/dsqrt(mu_mag) 
                 term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi0) 
                 term3= by*didpsi(psi0)/x/dsqrt(mu_mag) 
                 term4= rhoc*dhdpsi(psi0) 
                 term5= rhoc**gamma*dsdpsi(psi0)/(gamma-1.0d0)

				if((nx>=inter_switch).and.(Broot==0)) then

					 term6 = -(gradpsi_RHS(1,i,j)*(fmax_2D(i+1,j)-fmax_2D(i-1,j))/dx/2.d0 +  &
                 			gradpsi_RHS(2,i,j)*(fmax_2D(i,j+1)-fmax_2D(i,j-1))/dz/2.d0)/gradpsi_RHS(0,i,j)**2 +  &
                 			2.d0*delta_Bern*(psi_in(i,j)-psi_degen)/psic**2*fmax_2D(i,j)

					 term6 = rhoc*term6 * exp(-delta_Bern*((psi_in(i,j)-psi_degen)/psic)**2)

				 else

					term6 = 0.d0

				endif

!				if((nx>=inter_switch).and.(term6/=0.d0)) then
!					print*, i,j, term0, term1, term2, term3, term4, term5, term6
!					print*, '   '
!					pause
!				endif

                 res = term0+term1+term2+term3+term4-term5+term6


                 den = -( (1.0d0 - ma2rx)/(x + 0.5d0*dx) + &
                      (1.0d0 - ma2lx)/(x - 0.5d0*dx) )/(x*dx2) &
                      -( (1.0d0 - ma2rz) + &
                      (1.0d0 - ma2lz) )/(x*x*dz2)

		              continue



              res = res*mu_mag
			  ! Not sure what happened here...


              ! Store the residual
              if( (flag == 1).or.(k == max_it).or. (modulo(k,k_save) == 0).or.(k==1) ) then
!              if((gradpsi_RHS(0,i,j)>gpsi_switch).and.(dir_filter(i,j)==0)) then
!					residual(i,j) = res/gradpsi_RHS(0,i,j)**2
!				else
					if(nx>=inter_switch) then
						residual(i,j) = res/den
					else
						residual(i,j) = res
					endif
!				endif
              end if

              if(res>1.d0) then
                 continue
              endif

              ! -----------------------------------------------------

              continue

              psi_in(i,j) = psi0 - orp*res/den

              ! Find the max absolute value of psi in the plasma
              mx = dmax1(dabs(psi_in(i,j)),mx)

              ! Calculate the norm of the residual error
				if(nx>=inter_switch) then
	              anorm = anorm + dabs(res/den)
				else
	              anorm = anorm + dabs(res)
				endif

           end do
           isw = 3 - isw
        end do
        jsw = 3 - jsw
     end do

     ! update RHS
     	if(nx>=inter_switch) then
			call GS_RHS_all_Gauss
		endif

     ! Set the new maximum value of psic */
     psic = mx
     psic_13 = psic
     ! "Finished: Psi Update"

     ! -----------------------------------------------------

     ! Move the density calculation to a separate function
     ! Update rho and seek mach theta max
	if((nx>=inter_switch).and.(Broot==0)) then
     call update_rho_gauss(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
    else
     call update_rho(psi_in(:,:),rho,b_phi,nx,nz,1,mtm_acc,min_drho,  &
          min_ix,min_iz)
	endif

     call bc_psi_rho0(psi_in,rho,nx,nz)

     if(in_orp <= 0.0d0) then
        ! Chebyshev acceleration 
        if(nx >= 5) then
           std_orp = 1.0d0/(1.0d0 - 0.25d0*rjac*rjac*std_orp)
        end if
     endif

     if ((accelerate).and.(nx<inter_switch)) then
        orp = std_orp
     else
		orp = fix_orp
     endif

        !!$ call step_output(nx,psi_in(:,:),rho,residual)

     ! -----------------------------------------------------

     if(k == 1) then
        eps_anorm = eps*anorm
     else
        if((in_orp <= 0.0d0).and.accelerate.and.(nx<inter_switch)) then
           ! As we get closer to convergence, we want the solution
           ! to relax.  So as anorm approaches eps_anorm we want 
           ! the over relaxation parameter to go to some const. ~ 1
           ! Use x and mtm_soln as a temporary variable
           x = anorm/eps_anorm
           mtm_soln = 1.0d0 ! The limiting value for x ~ 1
           tmp = x/(x + orp/mtm_soln - 1.0d0)
           tmp = dmin1(tmp,1.0d0)
           orp = orp*tmp
           orp = dmax1(orp,mtm_soln)
        endif
     end if

     if(k == 1 .or. modulo(k,k_save) == 0) then
        print *, k,": anorm = ",real(anorm,skind)," eps_anorm = ",real(eps_anorm,skind)
        print *, "The Over-Relaxation Parameter =",orp,std_orp
        call step_output(nx,psi_in(:,:),rho,residual)
        continue
     end if

     write(111,*) k,anorm

     ! Check the convergence criteria
     if(anorm < eps_anorm .and. k > min_it) exit

  end do

	deallocate(fmax_2D)

  print *, k,": anorm = ",anorm," eps_anorm = ",eps_anorm
  print *, "Average residual =",anorm/(nx*nz)

  print *, "Mach Theta Max =",mach_theta_max
  print *, "Final Solution has Psi Center = ",psic
  print *


  return


contains



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine GS_RHS_all_Gauss
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ! this version calculates |grad psi| in i+1/2, j+1/2 using interpolations

	use pseudo_IMSl, only : DBSNAK, DBS2IN, dbs2dr

    real(kind=dkind) :: tdx, tdz


    tdx = dx*2.d0
    tdz = dz*2.d0

    gradpsi_RHS = 1.d-10
    

    ! fill grad psi
    do j = 1, nz
       do i = 1, nx

          if(sort_grid(i,j,0)<=0) cycle

          gradpsi_RHS(1,i,j) = (psi_in(i+1,j)-psi_in(i-1,j))/tdx
          gradpsi_RHS(2,i,j) = (psi_in(i,j+1)-psi_in(i,j-1))/tdz
          gradpsi_RHS(0,i,j) = sqrt(gradpsi_RHS(1,i,j)**2+gradpsi_RHS(2,i,j)**2)

          if(gradpsi_RHS(0,i,j)<1.d-6) gradpsi_RHS(0,i,j)=1.d-6 !(better patch with interpolation later...)

       enddo
    enddo


!	if((k==1).or.(modulo(k,k_save)==0)) then
!
!    ! write results for debugging
!
!		open(13, file='gradpsi_RHS.plt')
!
!		write(13,*)'TITLE="solution of GS equation with flow"'
!		write(13,*)'Variables =" X ","Y", "gpsi", "gpsi_x", "gpsi_z"'
!		write(13,*)'ZONE I=',nx,',J=',nz,',F=Point'
!
!		do j=1,nz
!
!		   z = z_coord(j)
!
!		   do i=1,nx
!
!			  x = x_coord(i)
!
!			  write(13,90) x,z,gradpsi_RHS(0,i,j),gradpsi_RHS(1,i,j),gradpsi_RHS(2,i,j)
!
!		   end do
!		end do
!
!		close(13)
!
!	endif

90  format(5(e12.6,3x))

    continue

  end subroutine GS_RHS_all_Gauss


end subroutine ngs_solve_all_gauss




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_jump_arrays(nx,nz,psi_last,rho,lambda0,Ptot,lambda, alpha_jump,gpsi_b,gPtor)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	! calculates the arrays needed by jump condition

	integer :: nx, nz
    real (kind=dkind), dimension(1:nx,1:nz) :: psi_last, rho
	real(kind=dkind), dimension(1:nx,1:nz) :: Ptot, lambda, alpha_jump
	real(kind=dkind) :: lambda0, dpsi_jump
	real(kind=dkind), dimension(1:nx,1:nz,2) :: gpsi_b, gPtor !grad psi for jump condition

	real(kind=dkind), dimension(1:nx,1:nz,2) :: gpsi, gPtot
	real(kind=dkind), dimension(1:nx,1:nz) :: Ptor
	real(kind=dkind) :: bphi_l
	real(kind=dkind), dimension(-1:1,-1:1) :: psi_around
	real(kind=dkind) :: dpsidx, dpsidz

	integer :: i, j

	if(lambda0==0.d0) then

		lambda = 0.d0
		alpha = 0.d0
		gpsi_b = 0.d0
		gPtor = 0.d0

		return

	endif

	! first get Ptot (and grad psi as a bonus)
	do j = 1, nz
	do i = 1, nx

		if(sort_grid(i,j,0)<0) then

			Ptot(i,j) = 0.d0

		else

			bphi_l = dsqrt(mu_mag)*(iofpsi(psi_last(i,j))/x_coord(i) + x_coord(i)*phiofpsi(psi_last(i,j))*omegaofpsi(psi_last(i,j)))/ &
				   (1.0d0 - phiofpsi(psi_last(i,j))**2/rho(i,j))

			psi_around(-1,0) = psi_last(i-1,j)
			psi_around(1,0) = psi_last(i+1,j)
			psi_around(0,0) = psi_last(i,j)
			psi_around(0,-1) = psi_last(i,j-1)
			psi_around(0,1) = psi_last(i,j+1)

			call psi_derivative_new(i,j,nx,nz,psi_around,dpsidx,dpsidz)

			gpsi(i,j,1) = dpsidx
			gpsi(i,j,2) = dpsidz

			Ptot(i,j) = Sofpsi(psi_last(i,j))*rho(i,j)**gamma +  &
							( ( (dpsidx**2+dpsidz**2) / x_coord(i)**2) + bphi_l**2) / (2.d0*mu_mag)
			Ptor(i,j) = Sofpsi(psi_last(i,j))*rho(i,j)**gamma +  &
							(bphi_l**2) / (2.d0*mu_mag)

		endif

	enddo
	enddo


	! then get grad Ptot
	do j = 1, nz
	do i = 1, nx

		if(sort_grid(i,j,0)<0) then

			gPtot(i,j,:) = 0.d0
			gPtor(i,j,:) = 0.d0

		else

			gPtot(i,j,1) = (Ptot(i+1,j)-Ptot(i-1,j))/(2.d0*dx)
			gPtot(i,j,2) = (Ptot(i,j+1)-Ptot(i,j-1))/(2.d0*dz)

			gPtor(i,j,1) = (Ptor(i+1,j)-Ptor(i-1,j))/(2.d0*dx)
			gPtor(i,j,2) = (Ptor(i,j+1)-Ptor(i,j-1))/(2.d0*dz)

		endif

	enddo
	enddo

	! next get alpha_jump (in two steps)
	do j = 1, nz
	do i = 1, nx

		if(sort_grid(i,j,0)<0) then

			alpha_jump(i,j) = 0.d0

		else

			alpha_jump(i,j) = gpsi(i,j,1)*gPtot(i,j,1)+gpsi(i,j,2)*gPtot(i,j,2)

		endif

	enddo
	enddo

	do j = 2, nz-1
	do i = 2, nx-1

		alpha_jump(i,j) = (alpha_jump(i+1,j+1)+alpha_jump(i+1,j-1)+alpha_jump(i-1,j+1)+alpha_jump(i-1,j-1))/4.d0

	enddo
	enddo

	! next get lambda

	dpsi_jump = abs(psic-psi_degen)*50./nx

	do j = 1, nz
	do i = 1, nx

		if((abs(psi_last(i,j)-psi_degen)>dpsi_jump).or.(sort_grid(i,j,0)<=0)) then

			lambda(i,j) = 0.d0

		else

			if(jump_option==3) then
				lambda(i,j) = lambda0
			elseif(jump_option==4) then
				lambda(i,j) = lambda0 * sign(1.d0,psi_last(i,j)-psi_degen)
			else
				lambda(i,j) = lambda0 * (1.d0 - ((psi_last(i,j)-psi_degen)/dpsi_jump)**2 )
			endif

		endif

	enddo
	enddo

	! finally, get gpsi_b
	do j = 1, nz
	do i = 1, nx

		if(lambda(i,j) == 0) then

			cycle

		else

			if(abs(gPtot(i,j,1))>abs(gPtot(i,j,2))) then

				gpsi_b(i,j,1) = (alpha_jump(i,j) - gPtot(i,j,2)*gpsi(i,j,2))/gPtot(i,j,1)
				gpsi_b(i,j,2) = gpsi(i,j,2)

			else

				gpsi_b(i,j,1) = gpsi(i,j,1)
				gpsi_b(i,j,2) = (alpha_jump(i,j) - gPtot(i,j,1)*gpsi(i,j,1))/gPtot(i,j,2)

			endif

		endif

	enddo
	enddo

	continue

	end subroutine get_jump_arrays


  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine step_output(n,psi,rho,residual)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer :: n
	real (kind=skind), dimension(n,n) :: out
	real (kind=dkind), dimension(n,n) :: rho, psi,residual
	integer :: i,j
	real(kind=dkind) :: X0loc, theta, mth, mphi,ex,ez,dx,dz, x

	out(1:n,1:n) = psi(1:n,1:n)
	call radial_plot(out,psi,n,n,"psi_step",n/2)

	out(1:n,1:n) = rho(1:n,1:n)
	call radial_plot(out,psi,n,n,"rho_step",n/2)

	out(1:n,1:n) = residual(1:n,1:n)
	call radial_plot(out,psi,n,n,"residual_step",n/2)


	continue

	return

  end subroutine step_output
  ! ------------------------------------------------------------------


  ! ------------------------------------------------------------------
  ! These routines are responsible for outputing the various physical 
  ! quantities such as the density, pressure, velocity, etc.
  ! Recall that (x, y, z) -> (r, phi, z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine physical_var(psi,rho,nx,nz,vx,vy,vz,bx,by,bz,p,ppar,pperp,  &
							malfven2,mslow,mcusp,mpol,beta,betapar,betaperp, &
							j_phi,temp,tpar,tperp,mtor,j_par, j_x, j_z,cs,csp,hyper, &
							br_2,bz_2,br_3,bz_3,gg, br_gg, bz_gg, j_phi_new)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx, dz, dx_a, dz_a
	use pseudo_IMSL, only : DBSNAK, DBS2VL, DBS2IN

    integer, intent(in) :: nx,nz
    real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi, rho
    real (kind=dkind), dimension(1:nx,1:nz), intent(inout) :: &
         vx,vy,vz,bx,by,bz,p,ppar,pperp,malfven2,mslow,mcusp,mpol,beta,  &
		 betapar,betaperp,temp,tpar,tperp,mtor,j_phi,j_par, j_x, j_z,cs,csp,  &
		 hyper,br_2,bz_2,br_3,bz_3, br_gg, bz_gg, j_phi_new
	real (kind=dkind), dimension(1:nx,1:nz,1:3), intent(inout) :: gg ! (gg_R,gg_Z,div gg)
    ! malfven_2 is the square of the alfvenic mach number
!	real (kind=dkind), dimension(1:nx,1:nz) :: j_x, j_z
!    real (kind=dkind) :: rtsafe
!    external rtsafe
    real (kind=dkind) :: phic ! Temp Var. for phiofpsi(psi(j,k))
    real (kind=dkind) :: omegac ! Temp Var. for omegaofpsi(psi(j,k))
    real (kind=dkind) :: vp2 ! Square of the poloidal velocity
    real (kind=dkind) :: ex,ez
    real (kind=dkind) :: dx2,dz2
    integer :: i,j,k
    real (kind=dkind) :: x,rloc,bfield,delstar
    real (kind=dkind) :: cf2,cs2,a2,ca2,ct2,cusp2
	real(kind=dkind), dimension(-1:1) :: psiprimx, psiprimz
	real(kind=dkind) :: term1,term2,term3,term4,term5,b_dot_v, term6
	real(kind=dkind) :: Alf, plocg, bpol
	real (kind=dkind), dimension(2) :: Gloc ! (1-M_A^2) grad psi /R^2
    real (kind=dkind) :: phirx,philx,phirz,philz
    real (kind=dkind) :: rhorx,rholx,rhorz,rholz
    real (kind=dkind) :: ma2rx, ma2lx, ma2rz, ma2lz, ma2c
	real(kind=dkind) :: psi0,psir,psil,psiu,psid ! local psi, right, left, up and down
	real(kind=dkind) :: drho_ds
	integer, parameter :: ord_temp = 2
	real(kind=dkind), dimension(1:nx,1:nz) :: bscoef_gg_r, bscoef_gg_z
	real(kind=dkind), dimension(1:nx+ord_temp) :: xknot_loc
	real(kind=dkind), dimension(1:nz+ord_temp) :: zknot_loc
	real(kind=dkind), dimension(1:nx+ord_temp) :: xknot_loc_m05
	real(kind=dkind), dimension(1:nz+ord_temp) :: zknot_loc_m05
	real(kind=dkind), dimension(1:nx) :: x_coord_m05
	real(kind=dkind), dimension(1:nz) :: z_coord_m05
	integer :: i_zone

    dx2 = dx*dx
    dz2 = dz*dz

	i_zone = 0 ! This will only be changed for bc_type==7

    ! Calculate the poloidal magnetic field
    do k = 1, nz
       do j = 1, nx

		 if(sort_grid(j,k,0)>=1) then

          ! Calculate the position of the (j,k) grid point
          x = x_coord(j)


			psiprimx(0) = ( dx_a(j-1)**2*psi(j+1,k) +  &
							(dx_a(j)**2-dx_a(j-1)**2)*psi(j,k) -  &
							dx_a(j)**2*psi(j-1,k) ) /  &
							( dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)) )

			psiprimx(1) = (psi(j+1,k) - psi(j,k))/dx_a(j)

			psiprimx(-1) = (psi(j,k) - psi(j-1,k))/dx_a(j-1)

			psiprimz(0) = ( dz_a(k-1)**2*psi(j,k+1) +  &
							(dz_a(k)**2-dz_a(k-1)**2)*psi(j,k) -  &
							dz_a(k)**2*psi(j,k-1) ) /  &
							( dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)) )

			psiprimz(1) = (psi(j,k+1) - psi(j,k))/dz_a(k)


			psiprimz(-1) = (psi(j,k) - psi(j,k-1))/dz_a(k-1)


          ! Calculate B_r
          if(k == 1) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_r as independent of z near the boundary.
!             bx(j,k) = (psi(j,k) - psi(j,k+1))/dz_a(1)/x !older version had the correct sign
             bx(j,k) = (psi(j,k+1) - psi(j,k))/dz_a(1)/x
          else if(k == nz) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_r as independent of z near the boundary.
!             bx(j,k) = (psi(j,k-1) - psi(j,k))/dz_a(nz)/x !older version had the correct sign
             bx(j,k) = (psi(j,k) - psi(j,k-1))/dz_a(nz)/x
          else
             ! use a centered difference
!             bx(j,k) = 0.5d0*(psi(j,k-1) - psi(j,k+1))/dz/x
			 bx(j,k) = ( dz_a(k-1)**2*psi(j,k+1) +  &
						 (dz_a(k)**2-dz_a(k-1)**2)*psi(j,k) -  &
						 dz_a(k)**2*psi(j,k-1) ) /  &
						 (dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)))/x
          end if

		  bx(j,k) = - bx(j,k) ! WARNING: the derivatives are taken with the PLUS sign for convenience!

          ! Calculate Bz
          if(j == 1) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_z as independent of x near the boundary.
             bz(j,k) = (psi(j+1,k) - psi(j,k))/dx_a(1)/(x + 0.5d0*dx_a(1))
          else if(j == nx) then
             ! use a 1-sided difference, equivalent to extrapolating
             ! B_z as independent of x near the boundary.
             bz(j,k) = (psi(j,k) - psi(j-1,k))/dx_a(nx)/(x - 0.5d0*dx_a(nx))
          else
             ! use a centered difference
!             bz(j,k) = 0.5d0*(psi(j+1,k) - psi(j-1,k))/dx/x
			 bz(j,k) = ( dx_a(j-1)**2*psi(j+1,k) +  &
						 (dx_a(j)**2-dx_a(j-1)**2)*psi(j,k) -  &
						 dx_a(j)**2*psi(j-1,k) ) /  &
						 (dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)))/x
          end if

		  ! -------------------------------------differentiate here for external points-------------------------------------
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma in the
		  ! "inner region" of the calculation

		if( (((tri_type/=13).AND.(bc_type/=7)).AND.((psi(j,k)/psic)<fraction))  &
						.OR.  &
!			((bc_type==7).AND.(sqrt((x_coord(j)-rmajor)**2+z_coord(k)**2)>0.4d0*z_size))  &
!						.OR.  &
			((tri_type==13).AND.(sort_grid(j,k,0)==1)) ) then

			! external point

			phic = 0.d0

			vx(j,k) = 0.d0
			vz(j,k) = 0.d0
			by(j,k) = dsqrt(mu_mag)*iofpsi(0.d0)/x_coord(j)

			vy(j,k) = 0.d0

			! Finally we calculate the pressure
			if (eq_type==1) then
				p(j,k) = sofpsi(0.d0)*rho(j,k)**gamma		
			elseif (eq_type==3) then
				bfield = dsqrt(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
				ppar(j,k) = tparofpsi(0.d0)*rho(j,k)
				if ((bfield-thetaofpsi(0.d0)*tparofpsi(0.d0))>0.d0) then
					pperp(j,k) = tparofpsi(0.d0)*rho(j,k)*bfield  &
						/(bfield-thetaofpsi(0.d0))*tparofpsi(0.d0)
				else
					pperp(j,k) = 0.d0
				endif
			endif	
			! Calculate the plasma beta
			if (eq_type==1) then
				beta(j,k) = 2.0d0*mu_mag*p(j,k)/(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
			elseif (eq_type==3) then

				if ((bx(j,k)**2  + bz(j,k)**2)>0.d0) then
					betapar(j,k) = 2.0d0*mu_mag*ppar(j,k)/  &
						(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
					betaperp(j,k) = 2.0d0*mu_mag*pperp(j,k)/  &
						(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
					beta(j,k) = (2.d0*betaperp(j,k)+betapar(j,k))/3.d0

				continue
				else
					betapar(j,k) = 0.d0
					betaperp(j,k) = 0.d0
				endif
			endif
			! And the plasma temperature
			if (eq_type==1) then
				temp(j,k) = mass*p(j,k)/rho(j,k)/eV
			elseif (eq_type==3) then
				tpar(j,k) = mass*ppar(j,k)/rho(j,k)/eV
				tperp(j,k) = mass*pperp(j,k)/rho(j,k)/eV
				temp(j,k) = (2.d0*tperp(j,k)+tpar(j,k))/3.d0
			endif

			! Now the toroidal current density

			if ((j>1).AND.(k>1).AND.(j<nx).AND.(k<nz)) then

				delstar = 0.D0

				j_phi(j,k) = x*delstar ! "+" is the correct sign
				j_phi_new(j,k) = 0.d0

			endif


			! Now let's calculate some Mach numbers
			! First some intermediate quantities

			vp2 = 0.0d0

			! Sound Speed
			if (eq_type==1) then
			a2 = gamma*p(j,k)/rho(j,k)
			elseif (eq_type==3) then
			a2 = ppar(j,k)/rho(j,k)
			endif			
			! Poloidal Alfven Speed Squared
			if(phic .eq. 0.0d0) then
			ca2 = (bx(j,k)*bx(j,k) + bz(j,k)*bz(j,k))/mu_mag/rho(j,k)
			endif
			! Toroidal Alfven Speed Squared
			ct2 = by(j,k)*by(j,k)/mu_mag/rho(j,k)
			! Fast Mode Speed Squared
			cf2 = 0.5d0*(a2 + ca2 + ct2 &
			+ dsqrt((a2 + ct2 - ca2)**2 + 4.0d0*ca2*ct2) )
			! Slow Mode Speed Squared
			cs2 = a2*ca2/cf2

			! Calculate the square of the Alfvenic Mach Number
			malfven2(j,k) = phic**2/rho(j,k)
			!          malfven2(j,k) = mu_mag*phic**2/rho(j,k)

			! Now the slow mode Mach number
			mslow(j,k) = dsqrt(malfven2(j,k)*cf2/a2)

			mslow(j,k) = sqrt(vp2/a2)*sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

			! Now the cusp Mach number
			mcusp(j,k) = mslow(j,k)*dsqrt(1.0d0 + cs2/cf2)

			mpol(j,k) = sqrt((vx(j,k)**2+vz(j,k)**2)/a2) *sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

			! Now the toroidal Mach number
			mtor(j,k) = dsqrt(vy(j,k)**2/a2)

			! sound speeds
			cs(j,k) = sqrt(a2)
			csp(j,k) = sqrt(a2*(bx(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+by(j,k)**2+bz(j,k)**2))

			hyper(j,k) = 1.d0

			cycle

		endif


		  ! -------------------------------------proceed from here for inner points-------------------------------------
		  ! October 8 2021 Added further differentiation for bc_type==7 and external region with plasma

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(j,k,0)-2)
				i_zone = max(-1,i_zone)
			endif

          ! Initialize phic
          phic = phiofpsi(psi(j,k),i_zone)
          omegac = omegaofpsi(psi(j,k),i_zone)

          ! Now let's calculate the poloidal velocity components
          vx(j,k) = phic*bx(j,k)/rho(j,k)/dsqrt(mu_mag)
          vz(j,k) = phic*bz(j,k)/rho(j,k)/dsqrt(mu_mag)

          ! Next we calculate B_phi
		  if (eq_type==1) then
			by(j,k) = dsqrt(mu_mag)*(iofpsi(psi(j,k),i_zone)/x + x*phic*omegac)/ &
				   (1.0d0 - phic**2/rho(j,k))
!			by(j,k) = (rmajor*bzero(psi(j,k))+x**2*sqrt(mu_mag)*phic*omegac)/(1.d0-phic**2/rho(j,k))/x
		  endif


          ! Next we calculate V_phi
          vy(j,k) = (phic*iofpsi(psi(j,k),i_zone)/rho(j,k)/x + x*omegac)/ &
               (1.0d0 - phic**2/rho(j,k))

          ! Finally we calculate the pressure
          if (eq_type==1) then
				p(j,k) = sofpsi(psi(j,k),i_zone)*rho(j,k)**gamma		
		  elseif (eq_type==3) then ! We have not implemented the free-boundary version of this yet
				bfield = dsqrt(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
				ppar(j,k) = tparofpsi(psi(j,k))*rho(j,k)
					if ((bfield-thetaofpsi(psi(j,k))*tparofpsi(psi(j,k)))>0.d0) then
						pperp(j,k) = tparofpsi(psi(j,k))*rho(j,k)*bfield  &
							/(bfield-thetaofpsi(psi(j,k))*tparofpsi(psi(j,k)))
					else
						pperp(j,k) = 0.d0
					endif
		  endif	
          ! Calculate the plasma beta
		  if (eq_type==1) then
			  beta(j,k) = 2.0d0*mu_mag*p(j,k)/(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
		  continue
		  elseif (eq_type==3) then

				if ((bx(j,k)**2  + bz(j,k)**2)>0.d0) then
					betapar(j,k) = 2.0d0*mu_mag*ppar(j,k)/  &
								(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
					betaperp(j,k) = 2.0d0*mu_mag*pperp(j,k)/  &
								(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)
					beta(j,k) = (2.d0*betaperp(j,k)+betapar(j,k))/3.d0

					continue
				else
					betapar(j,k) = 0.d0
					betaperp(j,k) = 0.d0
				endif
		  endif
		  ! And the plasma temperature
		  if (eq_type==1) then
			temp(j,k) = mass*p(j,k)/rho(j,k)/eV
		  elseif (eq_type==3) then
			tpar(j,k) = mass*ppar(j,k)/rho(j,k)/eV
			tperp(j,k) = mass*pperp(j,k)/rho(j,k)/eV
			temp(j,k) = (2.d0*tperp(j,k)+tpar(j,k))/3.d0
		  endif

          ! Now the toroidal current density

		 if ((j>1).AND.(k>1).AND.(j<nx).AND.(k<nz)) then


			if(eq_type==1) then

				b_dot_v = x*omegaofpsi(psi(j,k),i_zone)*by(j,k) + (phiofpsi(psi(j,k),i_zone)/rho(j,k))*  &
						( bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2 )/sqrt(mu_mag)

				term1= b_dot_v*dphidpsi(psi(j,k),i_zone)/dsqrt(mu_mag) 
				term2= x*(phiofpsi(psi(j,k),i_zone)*by(j,k)/dsqrt(mu_mag) +  &
						rho(j,k)*x*omegaofpsi(psi(j,k),i_zone))*domegadpsi(psi(j,k),i_zone)
				term3= by(j,k)*didpsi(psi(j,k),i_zone)/x/dsqrt(mu_mag)
				term4= rho(j,k)*dhdpsi(psi(j,k),i_zone)
				term5= rho(j,k)**gamma*dsdpsi(psi(j,k),i_zone)/(gamma-1.0d0)

				delstar = term1+term2+term3+term4-term5

				term6 = (psi(j+1,k)-psi(j-1,k)) *  &
								(phiofpsi(psi(j+1,k),i_zone)**2/rho(j+1,k)-phiofpsi(psi(j-1,k),i_zone)**2/rho(j-1,k)) / (4.d0*dx2) +  &
							(psi(j,k+1)-psi(j,k-1)) *  &
								(phiofpsi(psi(j,k+1))**2/rho(j,k+1)-phiofpsi(psi(j,k-1),i_zone)**2/rho(j,k-1)) / (4.d0*dz2)
				term6 = term6/x**2

			elseif(eq_type==3) then ! No free-boundary option yet

				bfield = dsqrt(bx(j,k)**2 + by(j,k)**2 + bz(j,k)**2)

				delstar =  &

                    + rho(j,k)*dtpardpsi(psi(j,k))  &
					+ by(j,k)*didpsi(psi(j,k))/x/dsqrt(mu_mag) &
					+ rho(j,k)*x**2*omegaofpsi(psi(j,k))*domegadpsi(psi(j,k)) &
					+ rho(j,k)*dhdpsi(psi(j,k)) &
					- rho(j,k)*(			&
							dtpardpsi(psi(j,k))*dlog(rho(j,k)/dofpsi(psi(j,k))*  &
									dabs(bfield-thetaofpsi(psi(j,k))*tparofpsi(psi(j,k)))/bfield) &
							- tparofpsi(psi(j,k))*dddpsi(psi(j,k))/dofpsi(psi(j,k))  &
							- tparofpsi(psi(j,k))  &
							  *(tparofpsi(psi(j,k))*dthetadpsi(psi(j,k))+thetaofpsi(psi(j,k))*dtpardpsi(psi(j,k))) &
							  *sign(1.d0,(bfield-thetaofpsi(psi(j,k))*tparofpsi(psi(j,k)))) )

				term6 = 0.d0

			endif

!!$
!!$	                term1= b_dot_v*dphidpsi(psi(i,j))/dsqrt(mu_mag) 
!!$                    term2= x*(phic*by/dsqrt(mu_mag) + rhoc*x*omegac)*domegadpsi(psi(i,j)) 
!!$                    term3= by*didpsi(psi(i,j))/x/dsqrt(mu_mag) 
!!$                    term4= rhoc*dhdpsi(psi(i,j)) 
!!$                    term5= rhoc**gamma*dsdpsi(psi(i,j))/(gamma-1.0d0)


			j_phi(j,k) = x * ( delstar - term6) / (1.d0-phic**2/rho(j,k))
			! "+" is the correct sign

			j_phi_new(j,k) = ( -( (psi(j+1,k)-psi(j,k))/(x+dx/2.d0) - (psi(j,k)-psi(j-1,k))/(x-dx/2.d0) )/dx2  &
										- (psi(j,k+1)-2.d0*psi(j,k)+psi(j,k-1))/dz2/x ) / mu_mag

		 endif


		  ! Now let's calculate some Mach numbers
          ! First some intermediate quantities

          ! Rather than using the poloidal velocities just calculated
          ! which involve differentiation, let's use the Bernoulli 
          ! equation to determine the square of the poloidal velocity
          ! which is really all we need for the Mach numbers
          if(phic .eq. 0.0d0) then
             vp2 = 0.0d0
          else
			if (eq_type==1) then
				vp2 = 2.0d0*hofpsi(psi(j,k),i_zone) + (x*omegac)**2 &
					  - (phic*by(j,k)/rho(j,k))**2 &
					  - 2.0d0*gamma/(gamma-1.0d0)*p(j,k)/rho(j,k)
			elseif (eq_type==3) then
				 vp2 = 0.d0
			endif

             ! Because of roundoff error this can be negative.
             vp2 = dmax1(0.0d0,vp2)
          endif

          ! Sound Speed
		  if (eq_type==1) then
		    a2 = gamma*p(j,k)/rho(j,k)
		  elseif (eq_type==3) then
			a2 = ppar(j,k)/rho(j,k)
		  endif			
          ! Poloidal Alfven Speed Squared
          if(phic .eq. 0.0d0) then
             ca2 = (bx(j,k)*bx(j,k) + bz(j,k)*bz(j,k))/mu_mag/rho(j,k)
          else 
             ca2 = rho(j,k)*vp2/(mu_mag*phic**2)
          endif
          ! Toroidal Alfven Speed Squared
          ct2 = by(j,k)*by(j,k)/mu_mag/rho(j,k)
          ! Fast Mode Speed Squared
          cf2 = 0.5d0*(a2 + ca2 + ct2 &
               + dsqrt((a2 + ct2 - ca2)**2 + 4.0d0*ca2*ct2) )
          ! Slow Mode Speed Squared
          cs2 = a2*ca2/cf2

          ! Calculate the square of the Alfvenic Mach Number
          malfven2(j,k) = phic**2/rho(j,k)
!          malfven2(j,k) = mu_mag*phic**2/rho(j,k)

          ! Now the slow mode Mach number
          mslow(j,k) = dsqrt(malfven2(j,k)*cf2/a2)

		  mslow(j,k) = sqrt(vp2/a2)*sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

          ! Now the cusp Mach number
          mcusp(j,k) = mslow(j,k)*dsqrt(1.0d0 + cs2/cf2)

			mpol(j,k) = sqrt((vx(j,k)**2+vz(j,k)**2)/a2) *sqrt((bx(j,k)**2+by(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+bz(j,k)**2))

		  ! Now the toroidal Mach number
          mtor(j,k) = dsqrt(vy(j,k)**2/a2)

		  ! sound speeds
		  cs(j,k) = sqrt(a2)
		  csp(j,k) = sqrt(a2*(bx(j,k)**2+bz(j,k)**2)/(bx(j,k)**2+by(j,k)**2+bz(j,k)**2))

!!$          ! Note: cs2 is only zero when the poloidal field goes to zero.
!!$          ! Since rho*vp = Phi * Bp and Phi = Phi(psi) is finite, 
!!$          ! this implies that vp is also zero. The correct limiting 
!!$          ! value should be treated better later.

			! check for hyperbolic regions

		  Alf = malfven2(j,k)	! Alfvén square
		  plocg = gamma*Sofpsi(psi(j,k),i_zone)*rho(j,k)**gamma
		  bpol = sqrt(bx(j,k)**2+bz(j,k)**2)
		  bfield=sqrt(bx(j,k)**2+by(j,k)**2+bz(j,k)**2)
		  hyper(j,k) = - ( Alf*(bfield**2+plocg*mu_mag)-plocg*mu_mag)  &
				/(Alf**2*bpol**2 - Alf*(bfield**2+plocg*mu_mag)+plocg*mu_mag)

			! recalculate B_pol using the LHS of the GS equation

			Gloc(1) = (1.d0-Alf**2)/x**2*(psi(j+1,k)-psi(j-1,k))/(2.d0*dx)
			Gloc(2) = (1.d0-Alf**2)/x**2*(psi(j,k+1)-psi(j,k-1))/(2.d0*dz)

			br_2(j,k) = -(psi(j,k+1) - psi(j,k))/dz_a(k)/x
			br_3(j,k) = -(psi(j,k) - psi(j,k-1))/dz_a(k)/x

			bz_2(j,k) = (psi(j+1,k) - psi(j,k))/dx_a(j)/(x + 0.5d0*dx_a(j))
			bz_3(j,k) = (psi(j,k) - psi(j-1,k))/dx_a(j)/(x - 0.5d0*dx_a(j))

			continue

       endif

	   end do
    end do

	! need to calculate the current at the end, because it requires the fields (in particular, b_phi)
    do k = 1, nz
       do j = 1, nx

		 if(sort_grid(j,k,0)>=1) then

          ! Calculate the position of the (j,k) grid point
          x = x_coord(j)

		 ! now the other components of the current
		 ! first J_R

          if( (k == 1).or.(sort_grid(j,k-1,0)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(j,k) = (by(j,k+1) - by(j,k))/dz_a(k)
          else if( (k == nz).or.(sort_grid(j,k+1,0)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_x(j,k) = (by(j,k) - by(j,k-1))/dz_a(k)
          else
             ! use a centered difference
			 j_x(j,k) = ( dz_a(k-1)**2*by(j,k+1) +  &
						 (dx_a(k)**2-dx_a(k-1)**2)*by(j,k) -  &
						 dz_a(k)**2*by(j,k-1) ) /  &
						 (dz_a(k)*dz_a(k-1)*(dz_a(k)+dz_a(k-1)))
          end if

		 j_x(j,k) =  - j_x(j,k) / mu_mag ! WARNING: the derivatives are taken with the PLUS sign for convenience!

		 ! then J_Z

          if( (j == 1).or.(sort_grid(j-1,k,0)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(j,k) = (by(j+1,k) - by(j,k))/dx_a(j)
          else if( (j == nx).or.(sort_grid(j+1,k,0)<=0) ) then
             ! use a 1-sided difference, equivalent to extrapolating
             j_z(j,k) = (by(j,k) - by(j-1,k))/dx_a(j)
          else
             ! use a centered difference
			 j_z(j,k) = ( dx_a(j-1)**2*by(j+1,k) +  &
						 (dx_a(j)**2-dx_a(j-1)**2)*by(j,k) -  &
						 dx_a(j)**2*by(j-1,k) ) /  &
						 (dx_a(j)*dx_a(j-1)*(dx_a(j)+dx_a(j-1)))
          end if

		 j_z(j,k) = (j_z(j,k) + by(j,k)/x) / mu_mag

		 ! then use the various components of the current to get J_par

		 j_par(j,k) = ( j_x(j,k)*bx(j,k) + j_z(j,k)*bz(j,k) + j_phi(j,k)*by(j,k) ) /  &
						dsqrt(bx(j,k)**2+by(j,k)**2+bz(j,k)**2)

		endif

	  enddo
	enddo

	! calculate staggered-grid stuff after the rest (i,j) -> (i-1/2,j-1/2)

    do k = 2, nz-1
       do j = 2, nx-1

			x = x_coord(j)

		if( (((tri_type/=13).AND.(bc_type/=7)).AND.((psi(j,k)/psic)<fraction))  &
							.OR.  &
!				((bc_type==7).AND.(sqrt((x_coord(j)-rmajor)**2+z_coord(k)**2)>0.4d0*z_size))  &
!							.OR.  &
				((tri_type==13).AND.(sort_grid(j,k,0)==1)) ) then

				cycle

			endif

			if(bc_type==7) then
				! for now we want the index to be 0 for the plasma region and -1 for the vacuum region
				i_zone = min(0,sort_grid(j,k,0)-2)
				i_zone = max(-1,i_zone)
			endif

			psi0 = psi(j,k)
			if(j>1) psil = psi(j-1,k)
			if(j<nx) psir = psi(j+1,k)
			if(k>1) psid = psi(j,k-1)
			if(k<nz) psiu = psi(j,k+1)

			drho_ds = van_Leer_slope_new(rho(j-1,k),rho(j,k),rho(j+1,k),dx_a(j-1),dx_a(j))

			phirx = phiofpsi(0.5d0*(psir + psi0),i_zone)
			if(jump_option==-1) then
				rhorx = (rho(j,k) + rho(j+1,k))/2.d0
			else
				rhorx = rho(j,k) + 0.5d0*dx_a(j)*drho_ds
			endif
			ma2rx = phirx*phirx/rhorx

			philx = phiofpsi(0.5d0*(psil + psi0),i_zone)
			if(jump_option==-1) then
				rholx = (rho(j,k) + rho(j-1,k))/2.d0
			else
				rholx = rho(j,k) - 0.5d0*dx_a(j-1)*drho_ds
			endif
			ma2lx = philx*philx/rholx

			! -----------------------------------------------------

			drho_ds = van_Leer_slope_new(rho(j,k-1),rho(j,k),rho(j,k+1),dz_a(k-1),dz_a(k))

			phirz = phiofpsi(0.5d0*(psi0 + psiu),i_zone)
			if(jump_option==-1) then
				rhorz = (rho(j,k)+rho(j,k+1))/2.d0
			else
				rhorz = rho(j,k) + 0.5d0*dz_a(k)*drho_ds
			endif
			ma2rz = phirz*phirz/rhorz

			philz = phiofpsi(0.5d0*(psi0 + psid),i_zone)
			if(jump_option==-1) then
				rholz = (rho(j,k)+rho(j,k-1))/2.d0
			else
				rholz = rho(j,k) - 0.5d0*dz_a(k)*drho_ds
			endif
			ma2lz = philz*philz/rholz

			gg(j,k,1) = (-(1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx))/x/dx	! (j-1/2,k)

			gg(j,k,2) = (-(1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz)	! (j,k-1/2)

			gg(j,k,3) = (( (1.0d0 - ma2rx)*(psir-psi0)/(x + 0.5d0*dx) &	! (j,k)
						  + (1.0d0 - ma2lx)*(psil-psi0)/(x - 0.5d0*dx) )  &
						  /(x*dx2) &
						+ ( (1.0d0 - ma2rz)*(psiu-psi0) &
						  + (1.0d0 - ma2lz)*(psid-psi0) )/(x*x*dz2) )

			br_gg(j,k) = -gg(j,k,2)/(1.0d0 - ma2lz)*x	! (j,k-1/2)

			bz_gg(j,k) = gg(j,k,1)/(1.0d0 - ma2lx)*x	! (j-1/2,k)


!			br_3(j,k) = -(psi(j,k) - psi(j,k-1))/dz_a(k)/x

!			bz_2(j,k) = (psi(j+1,k) - psi(j,k))/dx_a(j)/(x + 0.5d0*dx_a(j))


	  enddo
	enddo

	! interpolate gg

	! define half-grid step coordinates

	do i = 1, nx
		x_coord_m05(i) = x_coord(i) - dx/2.d0
	enddo

	do i = 1, nz
		z_coord_m05(i) = z_coord(i) - dz/2.d0
	enddo


	call DBSNAK(nx,x_coord(1:nx),ord_temp,xknot_loc)
!	print*, 'one'
	call DBSNAK(nz,z_coord(1:nz),ord_temp,zknot_loc)
!	print*, 'two'
	call DBSNAK(nx,x_coord_m05(1:nx),ord_temp,xknot_loc_m05)
!	print*, 'three'
	call DBSNAK(nz,z_coord_m05(1:nz),ord_temp,zknot_loc_m05)
!	print*, 'four'

	call DBS2IN(nx,x_coord_m05(1:nx),nz,z_coord(1:nz),gg(:,:,1),nx,  &
					ord_temp,ord_temp,xknot_loc_m05,zknot_loc,bscoef_gg_r)
!	print*, 'five'

	call DBS2IN(nx,x_coord(1:nx),nz,z_coord_m05(1:nz),gg(:,:,2),nx,  &
					ord_temp,ord_temp,xknot_loc,zknot_loc_m05,bscoef_gg_z)
!	print*, 'six'

	do j = 1, nz
	do i = 1, nx

		if((sort_grid(i,j,0)==0).or.(i==1).or.(j==1).or.(i==nx).or.(j==nz)) then

			br_gg(i,j) = 0.d0
			bz_gg(i,j) = 0.d0
			cycle

		endif

		br_gg(i,j) = DBS2VL(x_coord(i),z_coord(j),ord_temp,ord_temp,xknot_loc, &
						zknot_loc_m05,n,n,bscoef_gg_z)

		br_gg(i,j) = -br_gg(i,j) * x_coord(i) /(1.d0-malfven2(i,j))

		bz_gg(i,j) = DBS2VL(x_coord(i),z_coord(j),ord_temp,ord_temp,xknot_loc_m05, &
						zknot_loc,n,n,bscoef_gg_r)

		bz_gg(i,j) = bz_gg(i,j) * x_coord(i) /(1.d0-malfven2(i,j))

	enddo
	enddo

	continue

  end subroutine physical_var

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_trapped(nx, nz, psi, bpol, bphi, trapped)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer, intent(in) :: nx, nz

	real(kind=dkind), dimension(1:nx,1:nz) :: psi, bpol, bphi, trapped
	real(kind=dkind) :: psiloc, Bloc, thcrit, Bmax

	integer :: i,j

	do j = 1, nz
	do i = 1, nx

		if(sort_grid(i,j,0)<0) then

			trapped(i,j) = -1.d0
			cycle

		endif

		psiloc = psi(i,j)

		if(psiloc<modB_coef(1,1)) then

			trapped(i,j) = -1.d0
			cycle

		endif

		if(psiloc>modB_coef(1,enq)) then

			trapped(i,j) = 0.d0
			cycle

		endif

		Bmax = dbsval(psiloc,modB_ord, modB_coef(3,:),  &
					enq, modB_coef(4,1:enq) )

		Bloc = sqrt( bpol(i,j)**2 + bphi(i,j)**2 )

		if(Bloc>Bmax) then
		!this should not happen...

			thcrit = 0.d0
			continue

		else

			thcrit = sqrt(1.d0-Bloc/Bmax)

		endif

		trapped(i,j) = thcrit

	enddo
	enddo

	deallocate(modB_coef)

	continue

	return

end subroutine get_trapped


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine check_flow
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!This subroutine verifies whether there is any flow in the input

	real(kind=dkind) :: vmax	!maximum Mach number (only matters whether it's =0 or /=0)
	real(kind=dkind) :: x,y1,y2

	integer :: i

	vmax = 0.d0
	static_equi = .true.

	do i=0,201

		x=i/201.d0*psic
		y1 = mach_theta(x)
		y2 = mach_phi(x)
		vmax = max(vmax,abs(y1),abs(y2))

		if(vmax>0.d0) then

			static_equi = .false.
			exit

		endif

	enddo

	continue

end subroutine check_flow


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psiout
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!This subroutine writes the most important functions of Psi
	!and their derivatives in a text format
	!this is mainly for debugging purposes, and will not 
	! be used in the normal operation of the code

	integer i
	real (kind=dkind) x, y, xmin
	real(kind=dkind) :: gtheta_bound_dump
	integer :: zone_temp
	integer :: nsave = 201

!!$	return
	xmin = 0.d0

	!-----------------------------------------------------
	open(44,file='d.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	if(numerical_n) xmin = d_data(1,1)

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dofpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='dddpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dddpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	xmin = 0.d0

	!-----------------------------------------------------
	open(44,file='pofpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	if(numerical_p_iso) xmin = p_iso_data(1,1)

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=pofpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='dpdpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dpdpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	xmin = 0.d0

	if (eq_type==3) then

		!-----------------------------------------------------
		open(44,file='pparofpsi.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=pparofpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

		!-----------------------------------------------------
		open(44,file='pperpofpsi.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=pperpofpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

		!-----------------------------------------------------
		open(44,file='dppardpsi.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dppardpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

		!-----------------------------------------------------
		open(44,file='dpperpdpsi.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dpperpdpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

	endif

	!-----------------------------------------------------
	open(44,file='mach_phi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	if((numerical_omega).or.(numerical_mphi)) xmin = mphi_data(1,1)

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=mach_phi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

!	-----------------------------------------------------
	open(44,file='dmach_phidpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dmach_phidpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	xmin = 0.d0

	!-----------------------------------------------------
	open(44,file='mach_theta.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	if(numerical_mtheta) xmin = mtheta_data(1,1)

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=mach_theta(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

!	-----------------------------------------------------
	open(44,file='dmach_thetadpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dmach_thetadpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	xmin = 0.d0

	!-----------------------------------------------------
	open(44,file='bzero.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	if(numerical_F) xmin = F_data(1,1)

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=bzero(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='dbzerodpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dbzerodpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	xmin = 0.d0
!

	! Set xmin to a negative constant for Hameiri's functions
	!-----------------------------------------------------
	open(44,file='omega.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	xmin = -0.5d0

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=omegaofpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='domegadpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=domegadpsi(x,zone_temp)
		write(44,*) x/psic,y
	enddo

	close(44)

!	!-----------------------------------------------------
	open(44,file='h.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=hofpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='dhdpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'	
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dhdpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

!
!	!-----------------------------------------------------
!	open(44,file='i.plt',status='unknown',action='write')
!
!	write(44,*)'TITLE="solution of GS equation with flow"'
!    write(44,*)'Variables =" psi " "var"'		!, 	fname
!    write(44,*)'ZONE I=',201,',F=Point'
!
!	do i=0, nsave
!		x=i*psic/nsave
!		y=iofpsi(x)
!		write(44,*) x,y
!	enddo
!	close(44)
!
!	!-----------------------------------------------------
!	open(44,file='didpsi.plt',status='unknown',action='write')
!
!	write(44,*)'TITLE="solution of GS equation with flow"'
!    write(44,*)'Variables =" psi " "var"'		!, 	fname
!    write(44,*)'ZONE I=',201,',F=Point'
!
!	do i=0, nsave
!		x=i*psic/nsave
!		y=didpsi(x)
!		write(44,*) x,y
!	enddo
!	close(44)

!	stop

	if( ((eq3_opt==4).or.(eq3_opt==5)).or.((p_opt==4).or.(p_opt==5)) ) then

		!-----------------------------------------------------
		open(44,file='vol.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=volofpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

		!-----------------------------------------------------
		open(44,file='dvoldpsi.plt',status='unknown',action='write')

		write(44,*)'TITLE="solution of GS equation with flow"'
		write(44,*)'Variables =" psi " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*psic/nsave
			y=dvdpsi(x)
			write(44,*) x/psic,y
		enddo
		close(44)

	endif

	!-----------------------------------------------------

	if(bc_type==5) then

		!-----------------------------------------------------
		open(44,file='gtheta.plt',status='unknown',action='write')

		write(44,*)'TITLE="psi(r)"'
		write(44,*)'Variables =" r " "var"'		!, 	fname
		write(44,*)'ZONE I=',nsave,',F=Point'

		do i=1,nsave
			x=i*0.9d0*rmajor/nsave
			y=gtheta_bound_dump(0.5d0,x) !use a random angle
			write(44,*) x/psic,y
		enddo
		close(44)

		!-----------------------------------------------------
		open(44,file='gtheta2.plt',status='unknown',action='write')

		write(44,*)'TITLE="psi(r)"'
		write(44,*)'Variables =" theta " "var"'		!, 	fname
		write(44,*)'ZONE I=', nsave+1,',F=Point'

		do i=0, nsave
			x=i*2.d0*pi/nsave
			y=gtheta_bound_dump(x,0.16d0) !use a random angle
			write(44,*) x,y
		enddo
		close(44)

	endif


	!-----------------------------------------------------
	open(44,file='s.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=sofpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='dsdpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dsdpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='phi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=phiofpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)


	!-----------------------------------------------------
	open(44,file='dphidpsi.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=dphidpsi(x, zone_temp)
		write(44,*) x/psic,y
	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='F.plt',status='unknown',action='write')

	write(44,*)'TITLE="solution of GS equation with flow"'
    write(44,*)'Variables =" psi " "var"'		!, 	fname
    write(44,*)'ZONE I=', nsave+1,',F=Point'

	do i=0, nsave
		x=psic*(xmin + i*(1.d0-xmin)/nsave)
		if(x<0.d0) then
			zone_temp = -1
		else
			zone_temp = 0
		endif
		y=iofpsi(x, zone_temp)*sqrt(mu_mag)
		write(44,*) x/psic,y
	enddo
	close(44)


	open(44,file='d.dat',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=dofpsi(x)
		write(44,*) x/psic,y
	enddo
	close(44)

	open(44,file='p.dat',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=pofpsi(x)
		write(44,*) x/psic,y
	enddo
	close(44)


	open(44,file='b.dat',status='unknown',action='write')

	do i=0,101
		x=i/101.d0*psic
		y=bzero(x)*rmajor
		write(44,*) x/psic,y
	enddo
	close(44)

	continue

end subroutine psiout

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine FINESSE_input
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	integer i
	real (kind=dkind) x,y
	real(kind=dkind) :: gtheta_bound_dump

return

	!-----------------------------------------------------
	open(44,file='p_of_psi.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = pofpsi(x)
		write(44,*) 1.d0-x/psic, y

	enddo
	close(44)

	!-----------------------------------------------------
	open(44,file='F_of_psi.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = bzero(x)*rmajor
		write(44,*) 1.d0-x/psic, y

	enddo
	close(44)


	return

	!-----------------------------------------------------
	open(44,file='Lambda1_1.txt',status='unknown',action='write')
	open(45,file='Lambda1_0.txt',status='unknown',action='write')

	do i = 0,201

		x = i/201.d0*psic
		y = phiofpsi(x)**2*hofpsi(x)
		write(44,*) x/psic, y

		x = (1.d0 - i/201.d0)*psic
		y = phiofpsi(x)**2*hofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='Lambda2_1.txt',status='unknown',action='write')
	open(45,file='Lambda2_0.txt',status='unknown',action='write')

	do i = 0,201

		x = i/201.d0*psic
		y = gamma/(gamma-1.d0)*phiofpsi(x)**(2.d0*gamma)*sofpsi(x)
		write(44,*) x/psic, y

		x = (1.d0 - i/201.d0)*psic
		y = gamma/(gamma-1.d0)*phiofpsi(x)**(2.d0*gamma)*sofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='Lambda3_1.txt',status='unknown',action='write')
	open(45,file='Lambda3_0.txt',status='unknown',action='write')

	do i = 0,201

		x = i/201.d0*psic
		y = -sqrt(0.50d0) * iofpsi(x)
!		y = 0.50d0 * phiofpsi(x)**2 * omegaofpsi(x)**2
		write(44,*) x/psic, y

		x = (1.d0 - i/201.d0)*psic
		y = -sqrt(0.50d0) * iofpsi(x)
!		y = 0.50d0 * phiofpsi(x)**2 * omegaofpsi(x)**2
		write(45,*) 1.d0-x/psic, y

	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='Lambda4_1.txt',status='unknown',action='write')
	open(45,file='Lambda4_0.txt',status='unknown',action='write')

	do i = 0,201

		x = i/201.d0*psic
		y = sqrt(0.5d0) * phiofpsi(x)*omegaofpsi(x)
!		y = -iofpsi(x) / (phiofpsi(x) * omegaofpsi(x))
		write(44,*) x/psic, y

		x = (1.d0 - i/201.d0)*psic
		y = sqrt(0.5d0) * phiofpsi(x)*omegaofpsi(x)
!		y = -iofpsi(x) / (phiofpsi(x) * omegaofpsi(x))
		write(45,*) 1.d0-x/psic, y

	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(44,file='Lambda5_1.txt',status='unknown',action='write')
	open(45,file='Lambda5_0.txt',status='unknown',action='write')

	do i = 0,201

		x = i/201.d0*psic
		y = phiofpsi(x)**2 !* G_gravity * M_gravity
		write(44,*) x/psic, y

		x = (1.d0 - i/201.d0)*psic
		y = phiofpsi(x)**2 !* G_gravity * M_gravity
		write(45,*) 1.d0-x/psic, y

	enddo
	close(44)
	close(45)

	!-----------------------------------------------------
	open(45,file='new_lambda1_0.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = hofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(45)

	!-----------------------------------------------------
	open(45,file='new_lambda2_0.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = Sofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(45)

	!-----------------------------------------------------
	open(45,file='new_lambda3_0.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = omegaofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(45)

	!-----------------------------------------------------
	open(45,file='new_lambda4_0.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = Iofpsi(x)/phiofpsi(x)
		write(45,*) 1.d0-x/psic, y

	enddo
	close(45)

	!-----------------------------------------------------
	open(45,file='new_lambda5_0.txt',status='unknown',action='write')

	do i = 0,201

		x = (1.d0 - i/201.d0)*psic
		y = phiofpsi(x)**2
		write(45,*) 1.d0-x/psic, y

	enddo
	close(45)


		end subroutine FINESSE_input

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psi_derivative(i,j,nx,nz,psi,dpsidx,dpsidz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx_a, dz_a

	implicit none

	integer, intent(in) :: i,j,nx,nz
	real(kind=dkind), intent(inout) :: dpsidx,dpsidz
	real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi


	! d Psi / d x
    if(i == 1) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif left"
       dpsidx = (psi(i+1,j) - psi(i,j))/dx_a(1)
    else if(i == nx) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif right"
       dpsidx = (psi(i,j) - psi(i-1,j))/dx_a(nx)
    else
       ! use a centered difference
!       dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
	   dpsidx = ( dx_a(i-1)**2*psi(i+1,j) +  &
					(dx_a(i)**2-dx_a(i-1)**2)*psi(i,j) -  &
					dx_a(i)**2*psi(i-1,j) ) /  &
					( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
    end if
    ! d Psi / d z
    if(j == 1) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif left"
       dpsidz = (psi(i,j+1) - psi(i,j))/dz_a(j)
    else if(j == nz) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif right"
       dpsidz = (psi(i,j) - psi(i,j-1))/dz_a(j)
    else
       ! use a centered difference
!       dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
		dpsidz = ( dz_a(j-1)**2*psi(i,j+1) +  &
						(dz_a(j)**2-dz_a(j-1)**2)*psi(i,j) -  &
						dz_a(j)**2*psi(i,j-1) ) /  &
						( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )


    end if

	continue

	end subroutine psi_derivative

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine psi_derivative_new(i,j,nx,nz,psi,dpsidx,dpsidz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dx_a, dz_a

	implicit none

	integer, intent(in) :: i,j,nx,nz
	real(kind=dkind), intent(inout) :: dpsidx,dpsidz
	real (kind=dkind), dimension(-1:1,-1:1), intent(in) :: psi


	! d Psi / d x
    if(i == 1) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif left"
       dpsidx = (psi(1,0) - psi(0,0))/dx_a(1)
    else if(i == nx) then
       ! use a 1-sided difference
       ! print *, "One sided x-dif right"
       dpsidx = (psi(0,0) - psi(-1,0))/dx_a(nx)
    else
       ! use a centered difference
!       dpsidx = 0.5d0*(psi(i+1,j) - psi(i-1,j))/dx
	   dpsidx = ( dx_a(i-1)**2*psi(1,0) +  &
					(dx_a(i)**2-dx_a(i-1)**2)*psi(0,0) -  &
					dx_a(i)**2*psi(-1,0) ) /  &
					( dx_a(i)*dx_a(i-1)*(dx_a(i)+dx_a(i-1)) )
    end if
    ! d Psi / d z
    if(j == 1) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif left"
       dpsidz = (psi(0,1) - psi(0,0))/dz_a(j)
    else if(j == nz) then
       ! use a 1-sided difference
       ! print *, "One sided z-dif right"
       dpsidz = (psi(0,0) - psi(0,1))/dz_a(j)
    else
       ! use a centered difference
!       dpsidz = 0.5d0*(psi(i,j+1) - psi(i,j-1))/dz
		dpsidz = ( dz_a(j-1)**2*psi(0,1) +  &
						(dz_a(j)**2-dz_a(j-1)**2)*psi(0,0) -  &
						dz_a(j)**2*psi(0,-1) ) /  &
						( dz_a(j)*dz_a(j-1)*(dz_a(j)+dz_a(j-1)) )


    end if

	continue

end subroutine psi_derivative_new

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine getqstar(psi,jphi,p,bpol,bphi,beta,rho,vr,vphi,vz,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, dx, dz, inv_aspect_ratio
	use triangularity, only : volume, area, R_edge_min, R_edge_max,  &
										z_edge_min, z_edge_max, shape_ellipticity
	use interpolating_functions, only : interp_setup
	use pseudo_IMSL, only : dbsval_safe

	implicit none

	integer, intent(in) :: n
	real (kind=dkind), intent(in), dimension(1:n,1:n) :: jphi, p, psi, bpol,bphi, beta, rho,  &
																					   vr,vphi,vz
	real (kind=dkind) :: b_phi_zero
	real (kind=dkind) :: surf,qstar,dA,ex,ez	!,curr
	real (kind=dkind) :: nu,vol,bigR,epsilon,rminor,  &
						 ptot, Itot, beta_LDX, beta_int, bphiav,  &
						 I_Sauter, I_NCLASS, I_Sauter_phi, I_NCLASS_phi,  &
						B_den, mass, kin_en,  &
						B_p2, big_L_i, small_l_i
	real(kind=dkind) :: bs_loc
	integer :: i,j
	real(kind=dkind) :: R_l, R_r, angle
	real(kind=dkind) :: J_boot_interp(8,0:enq+1+ord_loc)
	real(kind=dkind), dimension(1:2) :: length, Bpint
	real(kind=dkind) :: psi_integral(2)

	! 7/7/2008 total BS current added:
	! set up J_BS interpolation first

	!2/4/2009 added diagnostic for plasma mass

	! Sauter calculation
	if (write_all.and.((bootstrap_option==0).or.(bootstrap_option==2))) then

		J_boot_interp(1,0) = 0.d0
		J_boot_interp(2,0) = J_boot(2,enq) - ( J_boot(2,enq-1)-J_boot(2,enq) ) /  &
							( psival(enq-1)-psival(enq) ) * psival(enq)

		J_boot_interp(1,enq+1) = psic
		J_boot_interp(2,enq+1) = 0.d0

		do i = 1,enq

			J_boot_interp(1,enq+1-i) = psival(i)
			J_boot_interp(2,enq+1-i) = J_boot(2,i)

		enddo

		call interp_setup(enq+2,ord_loc,J_boot_interp(1,0:enq+1),J_boot_interp(2,0:enq+1),  &
								J_boot_interp(3,0:enq+1+ord_loc),J_boot_interp(4,0:enq+1))

	endif

	! NCLASS calculation
	if (write_all.and.((bootstrap_option==1).or.(bootstrap_option==2))) then

		J_boot_interp(5,0) = 0.d0
		J_boot_interp(6,0) = J_boot(6,enq) - ( J_boot(6,enq-1)-J_boot(6,enq) ) /  &
							( psival(enq-1)-psival(enq) ) * psival(enq)

		J_boot_interp(5,enq+1) = psic
		J_boot_interp(6,enq+1) = 0.d0

		do i = 1,enq

			J_boot_interp(5,enq+1-i) = psival(i)
			J_boot_interp(6,enq+1-i) = J_boot(6,i)

		enddo

		call interp_setup(enq+2,ord_loc,J_boot_interp(5,0:enq+1),J_boot_interp(6,0:enq+1),  &
								J_boot_interp(7,0:enq+1+ord_loc),J_boot_interp(8,0:enq+1))

	endif

	! first, get the current and the cross section area


	dA = dx*dz

	curr = 0.d0
	qstar = 0.d0
	surf = 0.d0
	betator = 0.d0
	nu = 0.d0
	vol = 0.d0
	Itot = 0.d0
	beta_LDX = 0.d0
	beta_int = 0.d0
	bphiav = 0.d0
	bpolav = 0.d0
	I_Sauter = 0.d0
	I_NCLASS = 0.d0
	mass = 0.d0
	kin_en = 0.d0
	B_p2 = 0.d0

	if((Broot==4).or.(Broot==5)) then
		B_den = bzero(psic)
	else
		B_den = bzero(0.d0)
	endif

	do j=1,n
		do i=1,n

		    if(sort_grid(i,j,0)>=1) then

				if((psi(i,j)/psic)>=fraction) then	!inner region
!				if(dabs(psi(i,j)/psic)>=fraction) then	!inner region

!				bigR = x_coord(i)

!				bigR = 0.5d0*(x_coord(i)+x_coord(i+1))

!				if(grid_type/=0) dA = dx_a(i)*dz_a(j)

				bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)

				if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
										(dz_a(j) + dz_a(j-1))


!				if(bigR>rmajor+a_elps/3.d0) cycle	!trucco

				surf = surf+dA
				vol = vol + bigR*dA
				curr = curr +jphi(i,j)*dA
				betator = betator+p(i,j)*bigR*dA
				betapol = betapol + bpol(i,j)*bigR*dA
				ptot = ptot + p(i,j)*bigR*dA
				beta_LDX = beta_LDX + bpol(i,j)**2*bigR*dA
				beta_int = beta_int + beta(i,j)**2*bigR*dA
				bphiav = bphiav + bphi(i,j)*dA
				bpolav = bpolav + bpol(i,j)*dA
				mass = mass  + rho(i,j)*bigR*dA
				kin_en = kin_en  + rho(i,j)*(vr(i,j)**2+vphi(i,j)**2+vz(i,j)**2)*bigR*dA
				B_p2 = B_p2 + bpol(i,j)**2*bigR*dA

				if (write_all.and.((bootstrap_option==0).or.(bootstrap_option==2)))  then

					bs_loc = dA * dbsval_safe(psi(i,j),ord_loc,J_boot_interp(3,0:enq+1+ord_loc),enq+2,  &
								J_boot_interp(4,0:enq+1),J_boot_interp(2,0:enq+1))  &
								/ sqrt(bphi(i,j)**2+bpol(i,j)**2)

					I_Sauter = I_Sauter +  bs_loc

					I_Sauter_phi = I_Sauter_phi +  bs_loc * bphi(i,j)/sqrt(bphi(i,j)**2+bpol(i,j)**2)

				endif

				if (write_all.and.((bootstrap_option==1).or.(bootstrap_option==2)))  then

					bs_loc = dA * dbsval_safe(psi(i,j),ord_loc,J_boot_interp(7,0:enq+1+ord_loc),enq+2,  &
								J_boot_interp(8,0:enq+1),J_boot_interp(6,0:enq+1))  &
								/ sqrt(bphi(i,j)**2+bpol(i,j)**2)

					I_NCLASS = I_NCLASS +  bs_loc

					I_NCLASS_phi = I_NCLASS_phi + bs_loc * bphi(i,j)/sqrt(bphi(i,j)**2+bpol(i,j)**2)

				endif

				if (p(i,j)>pedge) then

!					Itot = Itot + dpdpsi(psi(i,j)) * bigR*dA
!					Itot = Itot + dpdpsi(psi(i,j))*bigR**2 * bigR*dA

				endif

				endif

			endif

		enddo
	enddo

	qstar = 2.d0*surf*bzero(0.d0)/(mu_mag*rmajor*curr)
	betastar = betator/vol * 2.d0*mu_mag/bzero(0.d0)**2 *(betator/vol)
	betator = betator/vol * 2.d0*mu_mag/B_den**2
	betapol = betapol/vol
	betapol = betator*B_den**2/betapol**2
	epsilon = (dsqrt(surf/pi)/rmajor)
	nu = betator*qstar**2/ epsilon
	beta_LDX = 2.d0*ptot*mu_mag/beta_LDX
	ptot = ptot*surf/vol
!	Itot = Itot / mu_mag
	bphiav = bphiav/surf
	bpolav = bpolav/surf
	big_L_i = 2.d0*pi*B_p2/(curr**2*mu_mag)
	small_l_i = 2.d0*big_L_i/(mu_mag*rmajor)

	call radius_theta(0.d0,rminor,R_r,ez) !ez not used here
	call radius_theta(pi,rminor,R_l,ez)

	inv_aspect_ratio = (R_r-R_l)/rmajor/2.d0

	R_edge_max = -1.d0
	z_edge_max = -1.d0
	R_edge_min = 1.d6*rmajor
	z_edge_min = 1.d6*rmajor

	do i = 1,r_of_theta_points

		angle = 2.d0*pi*i/r_of_theta_points
		call radius_theta(angle, rminor, ex, ez)

		if(ex>R_edge_max) R_edge_max = ex
		if(ex<R_edge_min) R_edge_min = ex
		if(ez>z_edge_max) z_edge_max = ez
		if(ez<z_edge_min) z_edge_min = ez

	enddo

	shape_ellipticity = (z_edge_max-z_edge_min)/(R_edge_max-R_edge_min)

	if((Broot==0).and.(jump_option<=-5)) call J_jump(psi,bpol,n,length,Bpint,psi_integral)

	print*, 'qstar = ',qstar
!	print*, 'nu = ', nu
	print*, 'plasma current = ',curr
	print*, 'betator = ',betator
	print*, 'betapol = ',betapol
	print*, 'epsilon = ',epsilon

	open(99,file='qstar.out',status='unknown',action='write')

	write(99,*) 'qstar = ',qstar
	write(99,*) 'surf = ',surf
	write(99,*) 'volume = ',vol
	write(99,*) 'plasma current = ',curr
	write(99,*) 'nu = ', nu
	write(99,*) 'betator = ',betator
	write(99,*) 'betapol = ',betapol
	write(99,*) 'average epsilon = ',epsilon
	write(99,*) 'inverse aspect ratio = ',inv_aspect_ratio
	write(99,*) 'total pressure * area = ', ptot
!	write(99,*) 'total current = ', Itot
	write(99,*) 'beta_LDX = ', beta_LDX
	write(99,*) 'B_phi_average = ', bphiav
	write(99,*) 'B_pol_average = ', bpolav
	write(99,*) 'mass = ', mass
	write(99,*) 'kinetic energy = ', kin_en
	write(99,*) 'L_i = ', big_L_i
	write(99,*) 'small_l_i = ', small_l_i
	write(99,*) '               '
	write(99,*) 'I_BS_Sauter = ', I_Sauter
	write(99,*) 'I_BS_NCLASS = ', I_NCLASS
	write(99,*) 'I_BS_Sauter_phi = ', I_Sauter_phi
	write(99,*) 'I_BS_NCLASS_phi = ', I_NCLASS_phi
	write(99,*) 'f_BS_Sauter = ', I_Sauter_phi/curr
	write(99,*) 'f_BS_NCLASS = ', I_NCLASS_phi/curr

	if((Broot==0).and.(jump_option<=-5)) then

		write(99,*) 'Bpol line integral (in, out):', Bpint(1), Bpint(2)
		write(99,*) 'length line integral (in, out):', length(1), length(2)
		write(99,*) 'Bpol average line integral (in, out):', Bpint(1)/length(1), Bpint(2)/length(2)
		write(99,*) 'integral locations:', psi_integral(1), psi_integral(2)

	endif


	close(99)

	if((broot==4).or.(broot==5)) then

		open(99,file='B_av.out',status='unknown',action='write')
		write(99,*) bphiav
		write(99,*) bpolav
		close(99)

	endif

	volume = vol
	area = surf

	return

end subroutine getqstar


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine J_jump(psi,bpol,n,length,Bpint,psival)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! calculates the integral of Bp dl and dl on the two sides of the discontinuity
! (for jump andcurrent spike calculation)

	use constant, only : dkind, dx, dz
	use interpolating_functions, only : lin_interp_2D

	implicit none

	integer, intent(in) :: n
	real (kind=dkind), intent(in), dimension(1:n,1:n) :: psi, bpol
	real(kind=dkind), dimension(1:2) :: length, Bpint
	real(kind=dkind) :: psival(2)
	real(kind=dkind) :: xx(2), zz(2)
	real(kind=dkind) :: dl
	real(kind=dkind) :: xQ, zQ, Bploc
	real(kind=dkind), dimension(1:4) :: psi4, bp4
	real(kind=dkind) :: psi_cmax, psi_cmin
	integer :: i, j, k, h

	length = 0.d0
	Bpint = 0.d0

	psival(1) = psi_degen + 0.01 !(inner)
	psival(2) = psi_degen - 0.01 !(outer)

	psival(1) = psi_degen + sqrt(1.d0/delta_Bern)*psic*2.d0 !(inner)
	psival(2) = psi_degen - sqrt(1.d0/delta_Bern)*psic*2.d0 !(outer)

	do j=1,n
	do i=1,n

		if(sort_grid(i,j,0)<=0) cycle

		psi4(1) = psi(i,j)
		psi4(2) = psi(i+1,j)
		psi4(3) = psi(i+1,j+1)
		psi4(4) = psi(i,j+1)

		psi_cmax = maxval(psi4)
		psi_cmin = minval(psi4)

		do k = 1, 2

			h = 1

			if((psi_cmax-psival(k))*(psi_cmin-psival(k))<0.d0) then
			! useful cell

				! check sides for intersectino points

				if((psi4(1)-psival(k))*(psi4(2)-psival(k))<0.d0) then
					xx(h) = x_coord(i) + (psival(k)-psi4(1))/(psi4(2)-psi4(1))*dx
					zz(h) = z_coord(j)
					h = h+1
				endif
				if((psi4(2)-psival(k))*(psi4(3)-psival(k))<0.d0) then
					xx(h) = x_coord(i+1)
					zz(h) = z_coord(j) + (psival(k)-psi4(2))/(psi4(3)-psi4(2))*dz
					h = h+1
				endif
				if((psi4(3)-psival(k))*(psi4(4)-psival(k))<0.d0) then
					xx(h) = x_coord(i+1) + (psival(k)-psi4(3))/(psi4(4)-psi4(3))*(-dx)
					zz(h) = z_coord(j+1)
					h = h+1
				endif
				if((psi4(4)-psival(k))*(psi4(1)-psival(k))<0.d0) then
					xx(h) = x_coord(i)
					zz(h) = z_coord(j+1) + (psival(k)-psi4(4))/(psi4(1)-psi4(4))*(-dz)
					h = h+1
				endif

				if(h<3) then
				! something dumb happened
					continue
					print*, 'error in J_jump', i, j
				endif

				! midpoint of segment
				xQ = sum(xx)/2.d0
				zQ = sum(zz)/2.d0

				dl = sqrt((xx(1)-xx(2))**2+(zz(1)-zz(2))**2)

				bp4(1) = bpol(i,j)
				bp4(2) = bpol(i+1,j)
				bp4(3) = bpol(i+1,j+1)
				bp4(4) = bpol(i,j+1)

				call lin_interp_2D(bp4,xQ,zQ,Bploc)

				length(k) = length(k) + dl
				Bpint(k) = Bpint(k) + Bploc*dl

			endif

		enddo

	enddo
	enddo

end subroutine J_jump

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine write_MARS_input(psi_all,rho,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! writes the input file for MARS runs (see routine in MARS_reader.for)
! remember, enq<->npsi, nsurf<->ntheta
! MARS does not require the 2pi point, so only write from 0 to 2pi-d_theta

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER

	real(kind=dkind), dimension(1:n,1:n) :: psi_all, psi, rho
	real(kind=dkind), dimension(:,:), allocatable :: rho1_out, rho2_out
	real(kind=dkind), dimension(:), allocatable :: xknot_rho, zknot_rho

	real(kind=dkind) :: psimax, delta, rguessl, rguessr, rguessl0, rguessr0
	real(kind=dkind) :: psiloc, thetaloc, rloc, bigrloc, zloc
	integer :: i, j
	integer :: nsurfdim
	real(kind=dkind) :: xmax, zmax
	integer :: ord_rho = 2	! use linear interpolation to avoid problems at the edge
	real(kind=dkind) :: rho_axis
	integer :: i_degen

	integer :: itemax = 200
	real(kind=dkind) :: xtoll = 1.d-14
	real(kind=dkind) :: ftoll = 1.d-14
	integer :: error

	! first, allocate the arrays

	allocate(rho1_out(1:enq,1:nsurf))
    allocate(bscoef_rho(1:n,1:n))

	if (broot==0) then
		allocate(rho2_out(1:enq,1:nsurf))
	    allocate(bscoef_rho2(1:n,1:n))
	endif

	allocate(bscoef_psi(1:n,1:n))

	allocate(xknot_mag(1:n+ord_loc))
	allocate(zknot_mag(1:n+ord_loc))

	allocate(xknot_rho(1:n+ord_rho))
	allocate(zknot_rho(1:n+ord_rho))

	xmax_mag = xmax
	zmax_mag = zmax

	!before starting, clean up psi to avoid problems with the interpolation

	do j=1,n
	do i=1,n

		if(sort_grid(i,j,0) < 0) then

			psi(i,j) = -psic/10.d0

		else

			psi(i,j) = psi_all(i,j)

		endif

	enddo
	enddo

	!first set up the interpolations

	call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
	call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

	call DBSNAK(n,x_coord(1:n),ord_rho,xknot_rho)
	call DBSNAK(n,z_coord(1:n),ord_rho,zknot_rho)

	call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
					ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)

	if(broot==0) then

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),rho1,n,  &
			 ord_rho,ord_rho,xknot_rho,zknot_rho,bscoef_rho)

	    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),rho2,n,  &
		     ord_rho,ord_rho,xknot_rho,zknot_rho,bscoef_rho2)

	else

		call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),rho,n,  &
			 ord_rho,ord_rho,xknot_rho,zknot_rho,bscoef_rho)

	endif

	! then set up the psi values

	allocate(psival(1:enq))

	psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_psi)

	psival(1) = psimax

	if(broot==0) then

		! psi_degen needs to be repeated, space psi accordingly

	else

		delta = psimax/(enq-1.d0)

		do i = 2,enq-1
			psival(i) = psimax - delta*(i-1.d0)
		enddo

	endif

	psival(enq) = 0.d0

	! then get the magnetic surfaces and q values (surface by surface)

	nsurfdim = nsurf+1+ord_surf

	allocate(rtab(enq,2,nsurfdim))
	allocate(bigrtab(enq,2,nsurfdim))
	allocate(ztab(enq,2,nsurfdim))
	allocate(thetatab(2,nsurfdim))

	! this sets up the 1D grid

	do j = 1, nsurf+1
		thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
	enddo

	! first fill obvious points

	! center
	i = 1

	do j = 1, nsurf+1

		rtab(i,1,j) = 0.d0
		bigrtab(i,1,j) = xmax
		ztab(i,1,j) = zmax

	enddo

	! edge
	i = enq

	do j = 1, nsurf+1

		call radius_theta(thetatab(1,j), rloc, bigrloc, zloc)
		rtab(i,1,j) = rloc
		bigrtab(i,1,j) = bigrloc
		ztab(i,1,j) = zloc

	enddo

	! then fill the rest

	rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
	rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

	do i = 2, enq-1

		surf_index = i

		rguessl = rguessl0
		rguessr = rguessr0

		psiloc = psival(i)

		! this cycle determines the (R,Z) coordinates of each magnetic surface
		do j = 1, nsurf+1

			thetaloc = thetatab(1,j)
			cloc = cos(thetaloc)
			sloc = sin(thetaloc)
			psiloc_mag = psiloc

			call secroot(findsurf,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

			if(error>0) then
				continue
			endif

			rtab(i,1,j) = rloc
			bigrtab(i,1,j) = xmax + rloc*cloc
			ztab(i,1,j) = zmax + rloc*sloc

			rguessr = rloc*1.1
			rguessl = rloc*0.9

		enddo

	enddo

	! then fill rho

	! put a constant value in the center

	rho_axis = DBS2VL(xmax,zmax,ord_rho,ord_rho,xknot_rho, &
		zknot_rho,n,n,bscoef_rho)

	i = 1

	do j = 1, nsurf

		rho1_out(i,j) = rho_axis

	enddo

	! then fill the rest with interpolated values

	do j = 1, nsurf
	do i = 2, enq

		bigrloc = bigrtab(i,1,j)
		zloc = ztab(i,1,j)

		rho1_out(i,j) = DBS2VL(bigrloc,zloc,ord_rho,ord_rho,xknot_rho, &
		zknot_rho,n,n,bscoef_rho)

	enddo
	enddo

	! --------------------------- now write the input file ---------------------------

	open(69, file='FLOW_MARS_input.txt')

	! first the indexes

	if(broot==0) then
		write(69,*) enq, nsurf, i_degen
	else
		write(69,*) enq, nsurf, 0
	endif

	! then the free functions

	do i = 1, enq

		psiloc = psival(i)
		write(69,*) psiloc, Iofpsi(psiloc)*sqrt(mu_mag),Sofpsi(psiloc),omegaofpsi(psiloc),phiofpsi(psiloc)

	enddo

	! last, write R, Z and rho
	! differentiate for transonic case later

	do i = 1, enq
	do j = 1, nsurf

		write(69,*) bigrtab(i,1,j), ztab(i,1,j), rho1_out(i,j)

	enddo
	enddo

	close(69)

	! deallocate stuff for backcompatibility

	deallocate(bscoef_psi)
	deallocate(xknot_mag)
	deallocate(zknot_mag)
	deallocate(xknot_rho)
	deallocate(zknot_rho)
	deallocate(psival)
	deallocate(rtab)
	deallocate(bigrtab)
	deallocate(ztab)
	deallocate(thetatab)
	deallocate(rho1_out)
    deallocate(bscoef_rho)

	if (broot==0) then
		deallocate(rho2_out)
	    deallocate(bscoef_rho2)
	endif

	continue

end subroutine write_MARS_input

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine magnetic_inversion(psi_all,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER

	real(kind=dkind), dimension(1:n,1:n) :: psi_all, psi

	real(kind=dkind) :: psimax, delta, rguessl, rguessr, rguessl0, rguessr0
	real(kind=dkind) :: psiloc, thetaloc, rloc, bigrloc, zloc
	integer :: i, j
	integer :: nsurfdim
	real(kind=dkind) :: xmax, zmax

	integer :: itemax = 200
	real(kind=dkind) :: xtoll = 1.d-14
	real(kind=dkind) :: ftoll = 1.d-14
	integer :: error

	! first, allocate the arrays

	allocate(bscoef_psi(1:n,1:n))

	allocate(xknot_mag(1:n+ord_loc))
	allocate(zknot_mag(1:n+ord_loc))

	xmax_mag = xmax
	zmax_mag = zmax

	!before starting, clean up psi to avoid problems with the interpolation

	do j=1,n
	do i=1,n

		if(sort_grid(i,j,0) < 0) then

			psi(i,j) = -psic/10.d0

		else

			psi(i,j) = psi_all(i,j)

		endif

	enddo
	enddo

	if(MARS_output) then
		continue
	else
		! we might as well leave the number of points assigned by input
!		enq = (n-1)/4
!		enq = (n-1)*3
!		enq = enq/8
	endif

!	enq = n/6

	!first set up the interpolations

	call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
	call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

	call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
					ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)

	! then set up the psi values

	allocate(psival(0:enq+1))

	psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_psi)

	delta = psimax/enq

	psival(0) = psimax

	psival(1) = psimax * ( 0.99d0 + 2.5d-3*log((n-1.d0)/128.d0)/log(2.d0) )

	if(psival(1)>=psimax) psival(1) = psimax*(1.d0-1.d-4)

	do i = 2,enq
		psival(i) = psimax - delta*(i-1.d0)
	enddo

	psival(enq+1) = 0.d0

	! then get the magnetic surfaces and q values (surface by surface)

	nsurfdim = nsurf+1+ord_surf

	allocate(rtab(enq,2,nsurfdim))
	allocate(bigrtab(enq,2,nsurfdim))
	allocate(ztab(enq,2,nsurfdim))
	allocate(thetatab(2,nsurfdim))

	! this sets up the 1D grid

	do j = 1, nsurf+1
		thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
	enddo

	rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
	rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

	do i = 1, enq

		surf_index = i

		rguessl = rguessl0
		rguessr = rguessr0

		psiloc = psival(i)

		! this cycle determines the (R,Z) coordinates of each magnetic surface
		do j = 1, nsurf+1

			thetaloc = thetatab(1,j)
			cloc = cos(thetaloc)
			sloc = sin(thetaloc)
			psiloc_mag = psiloc

			call secroot(findsurf,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

			if(error>0) then
				continue
			endif

			rtab(i,1,j) = rloc
			bigrtab(i,1,j) = xmax + rloc*cloc
			ztab(i,1,j) = zmax + rloc*sloc

			rguessr = rloc*1.1
			rguessl = rloc*0.9

		enddo

	enddo

	! now write the grid

	open(67, file='small_r_of_psi_theta.txt')
	open(68, file='R_of_psi_theta.txt')
	open(69, file='Z_of_psi_theta.txt')

	! also write the boundary (i = enq+1)
	do i = 1, enq+1
	do j = 1, nsurf+1

		if(i==enq+1) then

			call radius_theta(thetatab(1,j), rloc, bigrloc, zloc)
			write(67,*) psival(i), thetatab(1,j), rloc
			write(68,*) psival(i), thetatab(1,j), bigrloc
			write(69,*) psival(i), thetatab(1,j), zloc

		else

			write(67,*) psival(i), thetatab(1,j), rtab(i,1,j)
			write(68,*) psival(i), thetatab(1,j), bigrtab(i,1,j)
			write(69,*) psival(i), thetatab(1,j), ztab(i,1,j)

		endif

	enddo
	enddo

	close(67)
	close(68)
	close(69)

	! deallocate stuff for backcompatibility

	deallocate(bscoef_psi)
	deallocate(xknot_mag)
	deallocate(zknot_mag)
	deallocate(psival)
	deallocate(rtab)
	deallocate(bigrtab)
	deallocate(ztab)
	deallocate(thetatab)

	continue

end subroutine magnetic_inversion

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
subroutine magnetic_output(psi_all,bpol,bphi,rho,csp,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

    use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER

    real(kind=dkind), dimension(1:n,1:n) :: psi_all, bpol, bphi, psi, csp, rho

    integer :: enqbis
    real(kind=dkind), dimension(:), allocatable :: Rknot
    real(kind=dkind) :: psimax, delta, rguessl, rguessr, rguessl0, rguessr0
    real(kind=dkind) :: psiloc, thetaloc, rloc,  rlocprim, bigrloc, zloc,  &
         bpolloc, bphiloc, B_max, B_min, rholoc
    integer :: i, j,k
    integer :: nsurfdim
    integer :: n_int
    real(kind=dkind), dimension(:), allocatable :: w_int, t_int, integrand
    real(kind=dkind), dimension(:), allocatable :: w_int_boot, t_int_boot
    real(kind=dkind) :: a_q, b_q, c_q
    real(kind=dkind) :: xmax, zmax
    real(kind=dkind), dimension(:), allocatable :: L_ratio
    real(kind=dkind), dimension(:,:), allocatable :: Csp_tab, bscoef_csp
    real(kind=dkind) :: L_pol, L_tor

    integer :: itemax = 200
    real(kind=dkind) :: xtoll = 1.d-6
    real(kind=dkind) :: ftoll = 1.d-6
    real(kind=dkind) :: tol ! for max |B|
    integer :: error

    real(kind=dkind), dimension(:,:), allocatable :: viscosity_logs
    real(kind=dkind), dimension(:,:), allocatable :: poloidal_viscosity
    real(kind=dkind) :: rho_0, B_0	! for normalization

    ! first, allocate the arrays that would have been automatic,
    ! if not for the complaints of the LINUX compiler

    allocate(bscoef_psi(1:n,1:n))
    allocate(bscoef_bpol(1:n,1:n))
    allocate(bscoef_bphi(1:n,1:n))
    allocate(bscoef_csp(1:n,1:n))
    allocate(bscoef_rho(1:n,1:n))

    allocate(xknot_mag(1:n+ord_loc))
    allocate(zknot_mag(1:n+ord_loc))

    xmax_mag = xmax
    zmax_mag = zmax

    rho_0 = dofpsi(0.d0)
    B_0 = abs(bzero(0.d0))

    !before starting, clean up psi to avoid problems with the interpolation

    do j=1,n
       do i=1,n

          if(sort_grid(i,j,0) < 0) then

             psi(i,j) = -psic/10.d0

          else

             psi(i,j) = psi_all(i,j)

          endif

       enddo
    enddo

	if(MARS_output) then
		continue
	else
		! we might as well leave the number of points assigned by input
!		enq = (n-1)/4
!		enq = (n-1)*3
!		enq = enq/8
	endif


!	enq = n/6

    !first set up the interpolations

    call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
    call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bpol,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bpol)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),bphi,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_bphi)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),Csp,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_csp)

    call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),rho,n,  &
         ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_rho)

    ! then set up the psi values

    allocate(psival(0:enq+1))
    allocate(qval(enq))
    allocate(Rleft(enq))
    allocate(Rright(enq))
    allocate(L_ratio(enq))

    psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
         zknot_mag,n,n,bscoef_psi)

    delta = psimax/enq

    psival(0) = psimax

    psival(1) = psimax * ( 0.99d0 + 2.5d-3*log((n-1.d0)/128.d0)/log(2.d0) )

    if(psival(1)>=psimax) psival(1) = psimax*(1.d0-1.d-4)

    do i = 2,enq
       psival(i) = psimax - delta*(i-1.d0)
    enddo

    psival(enq+1) = 0.d0

    ! then get the magnetic surfaces and q values (surface by surface)

    nsurfdim = nsurf+1+ord_surf

    n_int = nsurf/2

    allocate(rtab(enq,2,nsurfdim))
    allocate(bigrtab(enq,2,nsurfdim))
    allocate(ztab(enq,2,nsurfdim))
    allocate(thetatab(2,nsurfdim))
    allocate(cross_section(enq,nsurfdim))
    allocate(Csp_tab(enq,nsurfdim))

    allocate(viscosity_logs(2,nsurfdim))
    allocate(poloidal_viscosity(enq,n_int))

    ! if the fraction of trapped particles is also required, we need to start the setup from here

    if(write_all) then

       ibreak_modB = enq + modB_ord
       allocate(modB_coef(4,ibreak_modB))

       tol = sqrt(epsilon(1.d0))

    endif

    ! if the complete output is desired,
    ! here we also set up the bootstrap current calculation,
    ! starting from the effective trapped particle fraction

    if(write_all) then

       if(bootstrap_option==0) then
          allocate(J_boot(5,1:enq))
          !		else	if(bootstrap_option==1) then
       elseif(bootstrap_option==2) then
          allocate(J_boot(7,1:enq))
       endif

       allocate(el_resistivity(1:3,1:enq))
       el_resistivity= 0.d0

       allocate(eff_trap(0:enq+1))
       allocate(surf_length(enq))
       allocate(B2_hat_ave(enq))
       allocate(boot_ne(enq))
       allocate(boot_pe(enq))
       allocate(boot_Te(enq))
       allocate(boot_ni(enq))
       allocate(boot_pi(enq))
       allocate(boot_Ti(enq))
       allocate(boot_neprim(enq))
       allocate(boot_peprim(enq))
       allocate(boot_Teprim(enq))
       allocate(boot_niprim(enq))
       allocate(boot_piprim(enq))
       allocate(boot_Tiprim(enq))
       allocate(inv_asp_ratio(enq))
       allocate(boot_Bmin(enq))
       allocate(boot_Bp_at_Bmin(enq))
       allocate(boot_R_of_Bmin(enq))
       allocate(boot_Rcenter(enq))
       allocate(JparB_ave(1:enq))

       if(bootstrap_option>=1) then

          allocate(boot_tor_flux(0:enq+1))
          allocate(boot_tor_rho(0:enq+1))
          allocate(boot_fhat(0:enq+1))
          allocate(boot_grho_ov_B2(0:enq+1))
          allocate(B2_ave(0:enq+1))
          allocate(Bm2_ave(0:enq+1))

       endif

    endif

    ! this sets up the 1D grid and the integration routine

    do j = 1, nsurf+1
       thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
    enddo

    call DBSNAK(nsurf+1, thetatab(1,1:nsurf+1),  &
         ord_surf,thetatab(2,1:nsurfdim))

    allocate(w_int(1:n_int))
    allocate(t_int(1:n_int))

    call set_weights(n_int,0.d0,2.d0*pi,w_int,t_int)

    if(write_all) then

       allocate(w_int_boot(1:n_int))
       allocate(t_int_boot(1:n_int))

       call set_weights(n_int,0.d0,1.d0,w_int_boot,t_int_boot)

    endif

    rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
    rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

    open(69,file='poloidal_viscosity.plt',action='write')

    write(69,*)'TITLE="poloidal viscosity for solution of GS equation with flow"'
    write(69,*)'Variables ="theta", "r", "R [m] ","z [m]", "viscosity"'
    write(69,*)'ZONE I=',n_int,',J=',enq,',F=Point'

    do i = 1, enq

       surf_index = i

       rguessl = rguessl0
       rguessr = rguessr0
!!$		rguessr = (x_coord(n)-xmax)*0.9d0

       ! check the right limit
!!$		do
!!$			psiloc = DBS2VL(xmax+rguessr,zmax,ord_loc,ord_loc,xknot, &
!!$											zknot,n,n,bscoef_psi)
!!$
!!$			if(psiloc<psival(i)) exit
!!$
!!$			rguessr = rguessr * 1.01d0
!!$
!!$		enddo

       psiloc = psival(i)

       ! this cycle determines the (R,Z) coordinates of each magnetic surface
       do j = 1, nsurf+1

          thetaloc = thetatab(1,j)
          cloc = cos(thetaloc)
          sloc = sin(thetaloc)
          psiloc_mag = psiloc

          call secroot(findsurf,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

          if(error>0) then
             continue
          endif

          rtab(i,1,j) = rloc
          bigrtab(i,1,j) = xmax + rloc*cloc
          ztab(i,1,j) = zmax + rloc*sloc

          rguessl = rloc*0.5d0
          rguessr = rguessr0
!!$			rguessr = rloc*1.05d0

          ! 7/6/2010 poloidal viscosity terms

          bigrloc = bigrtab(i,1,j)
          zloc = ztab(i,1,j)

          rholoc = DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_rho)

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bphi)

          !			print*, rholoc, rho_0
          !			print*, bpolloc, bphiloc, B_0
          viscosity_logs(1,j) = log(rholoc/rho_0) - 1.5d0*log(sqrt(bpolloc**2+bphiloc**2)/B_0)

       enddo

       ! now interpolate the coordinates and minor radius

!!$		call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
!!$			 rtab(2,1:nsurf+1),ord_surf,  &
!!$			 thetatab(2,1:nsurfdim),  &
!!$			 rtab(2,1:nsurf+1))

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            rtab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            rtab(i,2,1:nsurf+1))

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            bigrtab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            bigrtab(i,2,1:nsurf+1))

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            ztab(i,1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            ztab(i,2,1:nsurf+1))

       ! 7/6/2010 poloidal viscosity terms

       call DBSINT(nsurf+1, thetatab(1,1:nsurf+1),  &
            viscosity_logs(1,1:nsurf+1),ord_surf,  &
            thetatab(2,1:nsurfdim),  &
            viscosity_logs(2,1:nsurf+1))

       ! end of poloidal viscosity terms

       ! next, set up the integrand (Jacobian included!)

       ! also include the poloidal viscosity calculation in the same loop
       ! it is actually simpler to immediately save the poloidal viscosity to a tecplot file
       ! (even though it is not elegant at all...)

       allocate (integrand(1:n_int))

       do k = 1, n_int

          thetaloc = t_int(k)

          rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bphi)

          integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
               rloc * bphiloc / (bigrloc * bpolloc)  &
               / (2.d0*pi)

          poloidal_viscosity(i,k) = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, viscosity_logs(2,1:nsurf+1) )

          poloidal_viscosity(i,k) = poloidal_viscosity(i,k) / sqrt(1.d0+(rlocprim/rloc)**2)

          write(69,1234) thetaloc, rloc, bigrloc, zloc, poloidal_viscosity(i,k)

       enddo

       ! the calculation of q is complete, now let's save some numbers

       call integrate(n_int,integrand,w_int,qval(i))

       !----------------added 10/12/2009: length ratio calculation------------

       ! first: poloidal length
       do k = 1, n_int

          thetaloc = t_int(k)

          rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2)

       enddo

       call integrate(n_int,integrand,w_int,L_pol)

       ! second: toroidal length
       do k = 1, n_int

          thetaloc = t_int(k)

          rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
               nsurf+1, rtab(i,2,1:nsurf+1) )

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bphi)

          integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
               rloc * bphiloc / bpolloc

       enddo

       call integrate(n_int,integrand,w_int,L_tor)

       L_ratio(i) = L_tor/L_pol

       !----------------end of length ratio calculation------------

       deallocate (integrand)

       Rleft(i) = dbsval(pi,ord_surf, thetatab(2,:),  &
            nsurf+1, bigrtab(i,2,1:nsurf+1) )

       Rright(i) = dbsval(0.d0,ord_surf, thetatab(2,:),  &
            nsurf+1, bigrtab(i,2,1:nsurf+1) )

       ! if the trapped particle fraction is required, calculate |B|_MAX

       if(write_all) then

          modB_coef(1,enq+1-i) = psiloc
          B_max =  -brent(1.d-6,pi,2.d0*pi-1.d-6,modBval,tol,thetaloc)*abs(bzero(0.d0))
          ! for convenience, the field is normalized to B_phi_0 in modBval
          modB_coef(2,enq+1-i) = B_max
          continue

          ! also get B_min and details for the bootstrap current calculation

          B_min =  brent(1.d-6,pi,2.d0*pi-1.d-6,modBval_min,tol,thetaloc)*abs(b_phi_zero)

          bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, bigrtab(i,2,1:nsurf+1) )

          zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
               nsurf+1, ztab(i,2,1:nsurf+1) )

          bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
               zknot_mag,n,n,bscoef_bpol)

          boot_Bmin(i) = B_min
          boot_Bp_at_Bmin(i) = bpolloc
          boot_R_of_Bmin(i) = bigrloc

       endif

       if(write_all) then

          call bootstrap_setup

       endif

    enddo

    close(69)

1234 format(5(e16.9, 2x))

    !--------------------11/13/2009: poloidal cross section calculation--------------------
    ! note the approximation in the major radius!

    do j = 1, nsurf+1

       thetaloc = thetatab(1,j)

       do i = 1, enq

          if(i==1) then

             cross_section(i,j) = 2.d0*pi * rtab(i,1,j) *  &
                  0.5d0*(bigrtab(i,1,j)-xmax_mag)

          elseif(i>1) then

             cross_section(i,j) = 2.d0*pi * (rtab(i,1,j)-rtab(i-1,1,j)) *  &
                  0.5d0*(bigrtab(i,1,j)+bigrtab(i-1,1,j))

          endif

       enddo

    enddo

    !--------------------end of poloidal cross section calculation--------------------

    !--------------------03/26/2010: poloidal sound speed calculation--------------------

    do j = 1, nsurf+1

       thetaloc = thetatab(1,j)

       do i = 1, enq

          if(i==1) then

             Csp_tab(i,j) = 0.d0

          elseif(i>1) then

             Csp_tab(i,j) = DBS2VL(bigrtab(i,1,j),ztab(i,1,j),ord_loc,ord_loc,xknot_mag, &
                  zknot_mag,n,n,bscoef_csp)

          endif

       enddo

    enddo

    !--------------------end of poloidal cross section calculation--------------------

    if(write_all) then

       ! sets up the average for bootstrap fraction calculation
       call surf_ave_setup(bpol, n, -1.d0)

    endif

    open(4, file='qofR.dat', action='write', status='unknown')

    do i = enq,1,-1
       write(4,*) Rleft(i), qval(i)
    enddo

    do i = 1,enq
       write(4,*) Rright(i), qval(i)
    enddo

    close(4)

    open(4, file='qofpsi.dat', action='write', status='unknown')

    do i = 1,enq
       write(4,*) psival(i), qval(i)
    enddo

    close(4)

    ! finally, let's calculate the magnetic shear r q'/q

    allocate(Rknot(1:enq+ord_surf))
    allocate(shear(2,1:enq))

    call DBSNAK(enq,Rright,ord_surf,Rknot)

    call DBSINT(enq, Rright, qval, ord_surf,  &
         Rknot, shear(2,:))

    do k = 1,enq

       shear(1,k) = dbsder(1,Rright(k), ord_surf, Rknot,  &
            enq, shear(2,:) )

       rloc = Rright(k) - xmax

       shear(1,k) = shear(1,k) * rloc / qval(k)

    enddo

    open(4, file='shearofR.dat', action='write', status='unknown')

    do i = enq,1,-1
       write(4,*) Rleft(i), shear(1,i)
    enddo

    do i = 1,enq
       write(4,*) Rright(i), shear(1,i)
    enddo

    close(4)

    open(4, file='shearofpsi.dat', action='write', status='unknown')

    do i = 1,enq
       write(4,*) psival(i), shear(1,i)
    enddo

    close(4)

    ! last, save some tecplot files

    open(4, file='magnetic_R.plt', action='write', status='unknown')

    write(4,*)'TITLE="solution of GS equation with flow, magnetic postprocessing"'
    write(4,*)'Variables =" R [m] ","q", "shear", "L_ratio"'		!, 	fname
    write(4,*)'ZONE I=',2*enq,',F=Point'

    do i = enq,1,-1
       write(4,88) Rleft(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    do i = 1,enq
       write(4,88) Rright(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    close(4)

    open(4, file='magnetic_psi.plt', action='write', status='unknown')

    write(4,*)'TITLE="solution of GS equation with flow, magnetic postprocessing"'
    write(4,*)'Variables =" psi ","q", "shear", "L_ratio"'		!, 	fname
    write(4,*)'ZONE I=',enq,',F=Point'

    do i = 1,enq
       write(4,88) psival(i), qval(i), shear(1,i), L_ratio(i)
    enddo

    close(4)

    ! write the cross section output

    open(69, file='poloidal_cross_section.plt', action='write', status='unknown')

    write(69,*)'TITLE="solution of GS equation with flow, poloidal cross section"'
    write(69,*)'Variables =" psi ","theta", "area", "AxV"'		!, 	fname
    write(69,*)'ZONE I=',enq,',J=',nsurf+1,',F=Point'

    do j = 1, nsurf+1
       do i = 1,enq
          write(69,66) (psival(i)+psival(i-1))/2.d0, thetatab(1,j), cross_section(i,j), cross_section(i,j)*csp_tab(i,j)
       enddo
    enddo

    close(69)

    ! insert here the trapped particle calculation setup

    if(write_all) then

       call DBSNAK(enq, modB_coef(1,1:enq), modB_ord, modB_coef(3,:))

       call DBSINT(enq, modB_coef(1,1:enq), modB_coef(2,1:enq), modB_ord,  &
            modB_coef(3,:), modB_coef(4,1:enq))

    endif

    ! insert here the q profile interpolation

    ! first, extrapolate q0 with a parabola

    a_q = (qval(1)-qval(2))/(2.d0*(psival(1)-psival(2))*(psival(1)-psimax))
    b_q = -2.d0*a_q*psimax
    c_q = qval(1) - a_q*psival(1)**2 - b_q*psival(1)

    q_c = a_q*psimax**2 + b_q*psimax + c_q

    ! then, extrapolate q_edge with a linear extrapolation

    qe = qval(enq) + (qval(enq)-qval(enq-1))/(psival(enq)-psival(enq-1))*(-psival(enq)) !since psi_edge=0

    ! now set up the interpolation

    enqbis = enq+2 ! added first and last point, axis and edge

    ibreak_q = enq+2+q_ord	!enqbis + q_ord
    allocate(q_coef(4,ibreak_q))

    q_coef(1,1) = 0.d0
    q_coef(2,1) = qe

    do i=1,enq

       q_coef(1,enqbis-i) = psival(i)
       q_coef(2,enqbis-i) = qval(i)

    enddo

    q_coef(1,enqbis) = psimax
    q_coef(2,enqbis) = q_c

    call DBSNAK(enqbis, q_coef(1,1:enqbis), q_ord, q_coef(3,:))

    call DBSINT(enqbis, q_coef(1,1:enqbis), q_coef(2,1:enqbis), q_ord,  &
         q_coef(3,:), q_coef(4,1:enqbis))

    if(write_all) then

       call DBSNAK(enq, modB_coef(1,1:enq), modB_ord, modB_coef(3,:))

       call DBSINT(enq, modB_coef(1,1:enq), modB_coef(2,1:enq), modB_ord,  &
            modB_coef(3,:), modB_coef(4,1:enq))

    endif

    ! now let's clean up all the data

    deallocate(Rknot)
    deallocate(shear)

    deallocate(w_int)
    deallocate(t_int)
    deallocate(L_ratio)

    if(write_all) then

       continue

    else

       deallocate(psival)
       deallocate(qval)
       deallocate(w_int_boot)
       deallocate(t_int_boot)

    endif

    !	deallocate(bscoef_psi)
    !	deallocate(bscoef_bpol)
    !	deallocate(bscoef_bphi)

    continue

66  format(3(e26.17,3x))
88  format(4(e26.17,3x))
    !________________________

  contains

    subroutine bootstrap_setup

      implicit none

      real(kind=dkind), dimension(5,n_int) :: integrand_boot
      real(kind=dkind) :: lambda, btotloc, int_step
      real(kind=dkind) :: Rmin_loc, Rmax_loc

      integer :: it, jt

      ! step 0, compute the length of each magnetic surface (one at a time!)

      Rmin_loc = rmajor
      Rmax_loc = 0.d0

      do it = 1, n_int

         thetaloc = t_int(it)

         rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
              nsurf+1, rtab(i,2,1:nsurf+1) )

         rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
              nsurf+1, rtab(i,2,1:nsurf+1) )

         integrand_boot(1,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
              rloc / (2.d0*pi)

      enddo

      call integrate(n_int,integrand_boot(1,:),w_int,surf_length(i))

      ! first, compute the average in the denominator and the effective trapped particle fraction integral;
      ! then do directly the d_lambda integral

      do jt = 1, n_int

         lambda = t_int_boot(jt)

         do it = 1, n_int

            thetaloc = t_int(it)

            rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, rtab(i,2,1:nsurf+1) )

            rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
                 nsurf+1, rtab(i,2,1:nsurf+1) )

            bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, bigrtab(i,2,1:nsurf+1) )

            zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
                 nsurf+1, ztab(i,2,1:nsurf+1) )

            bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
                 zknot_mag,n,n,bscoef_bpol)

            bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
                 zknot_mag,n,n,bscoef_bphi)

            btotloc = min(sqrt(bpolloc**2+bphiloc**2),B_max)

            ! this integrates the denominator
            integrand_boot(1,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                 sqrt(1.d0-lambda*btotloc/B_max) *  &
                 rloc / (2.d0*pi)

            ! this integrates (B/B_max)**2
            integrand_boot(2,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                 (btotloc/B_max)**2 *  &
                 rloc / (2.d0*pi)

            if(bootstrap_option>=1) then
               ! also calculate other averages

               ! this integrates B**2
               integrand_boot(4,it) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
                    btotloc**2 *  &
                    rloc / (2.d0*pi)

               ! this integrates 1/B**2
               integrand_boot(5,it) = sqrt(1.d0+(rlocprim/rloc)**2) /  &
                    btotloc**2 *  &
                    rloc / (2.d0*pi)

            endif

            ! also check for the aspect ratio at this stage

            if(jt==1) then
               ! just do it once!

               Rmax_loc = max(Rmax_loc,bigrloc)

               Rmin_loc = min(bigrloc,Rmin_loc)

            endif

         enddo

         call integrate(n_int,integrand_boot(1,:),w_int,int_step)
         ! denominator

         call integrate(n_int,integrand_boot(2,:),w_int,B2_hat_ave(i))
         ! B**2

         integrand_boot(3,jt) = lambda / int_step

         if(bootstrap_option>=1) then

            call integrate(n_int,integrand_boot(4,:),w_int,B2_ave(i))
            ! B**2

            call integrate(n_int,integrand_boot(5,:),w_int,Bm2_ave(i))
            ! B**2

         endif

      enddo

      call integrate(n_int,integrand_boot(3,:),w_int_boot,int_step)

      eff_trap(i) = 1.d0 - 0.75d0 * B2_hat_ave(i)* int_step

      if(eff_trap(i)<0.d0) eff_trap(i) = 0.d0
      if(eff_trap(i)>1.d0) eff_trap(i) = 1.d0

      ! set the "aspect ratio"
      inv_asp_ratio(i) = (Rmax_loc-Rmin_loc)/(Rmax_loc+Rmin_loc)

      ! set the "center"
      boot_rcenter(i) = (Rmax_loc+Rmin_loc)/2.d0

      ! the minimum field, the poloidal field at the same location and
      ! also the location are saved in the magnetic output section

      ! then proceed with the various physical quantities (density, temperature, etc.)
      ! the case with more complicated dependencies will be differentiated later

      if(static_equi) then
         ! all depends on psi only

         boot_ne(i) = dofpsi(psival(i))/mass ! in m^-3
         ! ion mass to turn the plasma mass in a number density
         boot_pe(i) = pofpsi(psival(i))*pe_ov_p
         ! this is p_e = pe_ov_p * (total p)

         boot_Te(i) = boot_pe(i)/boot_ne(i)/eV ! in eV
         ! this is T_e = pe_ov_p * (total T)

         boot_ni(i) = dofpsi(psival(i))/mass ! in m^-3
         ! ion mass to turn the plasma mass in a number density
         boot_pi(i) = pofpsi(psival(i))*(1.d0-pe_ov_p)
         ! this is p_i = (total p) - p_e

         boot_Ti(i) = boot_pi(i)/boot_ni(i)/eV ! in eV
         ! this is p_i = (total p) - p_e

         boot_neprim(i) = dddpsi(psival(i))/mass	! * psic?
         boot_peprim(i) = dpdpsi(psival(i)) * pe_ov_p	! * psic?

         boot_Teprim(i) = ( boot_peprim(i)/boot_ne(i) -  &
              boot_neprim(i)*boot_pe(i)/boot_ne(i)**2 )/eV

         !		boot_neprim(i) = boot_neprim(i)/mass

         boot_niprim(i) = dddpsi(psival(i))/mass	! * psic?
         boot_piprim(i) = dpdpsi(psival(i)) * (1.d0-pe_ov_p)	! * psic?

         boot_Tiprim(i) = ( boot_piprim(i)/boot_ni(i) -  &
              boot_niprim(i)*boot_pi(i)/boot_ni(i)**2 )/eV

         !		boot_niprim(i) = boot_niprim(i)/mass

      endif

      continue

      return

    end subroutine bootstrap_setup




  end subroutine magnetic_output


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine plasma_boundary(psi,xmax,zmax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! for free-boundary equilibria, determine the plasma edge (psi=0)

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBSNAK, DBSINT, DBSDER

	real(kind=dkind), dimension(1:n,1:n) :: psi

	real(kind=dkind) :: psimax, rguessl, rguessr, rguessl0, rguessr0
	real(kind=dkind) :: psiloc, thetaloc, rloc
	integer :: i, j,k
	integer :: nsurfdim
	real(kind=dkind) :: xmax, zmax

	integer :: itemax = 200
	real(kind=dkind) :: xtoll = 1.d-6
	real(kind=dkind) :: ftoll = 1.d-6
	integer :: error

	! first, allocate the arrays that would have been automatic,
	! if not for the complaints of the LINUX compiler

	allocate(bscoef_psi(1:n,1:n))

	allocate(xknot_mag(1:n+ord_loc))
	allocate(zknot_mag(1:n+ord_loc))

	xmax_mag = xmax
	zmax_mag = zmax

	!first set up the interpolations

	call DBSNAK(n,x_coord(1:n),ord_loc,xknot_mag)
	call DBSNAK(n,z_coord(1:n),ord_loc,zknot_mag)

	call DBS2IN(n,x_coord(1:n),n,z_coord(1:n),psi,n,  &
					ord_loc,ord_loc,xknot_mag,zknot_mag,bscoef_psi)


	psimax = DBS2VL(xmax,zmax,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_psi)

	! get the magnetic surface corresponding to the plasma edge

	nsurfdim = nsurf+1+ord_surf

	allocate(thetatab(2,nsurfdim))

	! this sets up the 1D grid and the integration routine

	do j = 1, nsurf+1
		thetatab(1,j) = (j-1.d0)*pi*2.d0/nsurf
	enddo

	call DBSNAK(nsurf+1, thetatab(1,1:nsurf+1),  &
		ord_surf,thetatab(2,1:nsurfdim))


	rguessr0 = max(x_coord(n)-xmax, xmax-x_coord(1), z_coord(n)-zmax, zmax-z_coord(1))
	rguessl0 = (x_coord(n)-x_coord(1)+z_coord(n)-z_coord(1))/5.d2

	open(33,file='r_edge_of_theta.dat')

	rguessl = rguessl0
	rguessr = rguessr0

	psiloc = 0.d0
	psiloc_mag = psiloc
	rloc=0.1d0

	! this cycle determines the (R,Z) coordinates of the magnetic surface
	! note that the radius is calculated from the geometric axis (differently from the magnetic_output routine)
	do j = 1, nsurf+1

		thetaloc = thetatab(1,j)
		cloc = cos(thetaloc)
		sloc = sin(thetaloc)

		error = 2

		do while(error>0)

			call secroot(findsurf2,rguessr,rguessl,rloc,itemax,xtoll,ftoll,error)

			if(error>0) then
				continue
				rguessr = rguessr + rloc*1.d-2
			endif

			if(rguessr>2.d0*rloc) then
				print*, 'problem in r_edge at theta =', thetaloc
				exit
			endif

		enddo

		write(33,*) thetaloc, rloc

		rguessl = rloc*0.5d0
		rguessr = rloc*0.9d0

	enddo

	close(33)

	deallocate(bscoef_psi)
	deallocate(xknot_mag)
	deallocate(zknot_mag)
	deallocate(thetatab)



end subroutine plasma_boundary




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function findsurf(r) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2VL

	real(kind=dkind) :: r
	real(kind=dkind) :: answer
	real(kind=dkind) :: bigRtemp, bigZtemp

	bigRtemp = xmax_mag + r*cloc
	bigZtemp = zmax_mag + r*sloc

	if(bigRtemp>x_coord(n)) bigRtemp = x_coord(n)
	if(bigRtemp<x_coord(1)) bigRtemp = x_coord(1)

	if(bigZtemp>z_coord(n)) bigZtemp = z_coord(n)
	if(bigZtemp<z_coord(1)) bigZtemp = z_coord(1)

	answer = DBS2VL(bigRtemp,bigZtemp, &
								ord_loc,ord_loc,xknot_mag,zknot_mag,n,n,bscoef_psi) &
										- psiloc_mag

	return

end function findsurf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function findsurf2(r) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2VL

	real(kind=dkind) :: r
	real(kind=dkind) :: answer
	real(kind=dkind) :: bigRtemp, bigZtemp

	bigRtemp = rmajor + r*cloc
	bigZtemp = r*sloc

	if(bigRtemp>x_coord(n)) bigRtemp = x_coord(n)
	if(bigRtemp<x_coord(1)) bigRtemp = x_coord(1)

	if(bigZtemp>z_coord(n)) bigZtemp = z_coord(n)
	if(bigZtemp<z_coord(1)) bigZtemp = z_coord(1)

	answer = DBS2VL(bigRtemp,bigZtemp, &
								ord_loc,ord_loc,xknot_mag,zknot_mag,n,n,bscoef_psi) &
										- psiloc_mag

	return

end function findsurf2



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		function modBval(th) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			use pseudo_IMSL, only : DBS2VL

			real(kind=dkind) :: th
			real(kind=dkind) :: answer, Bp
			real(kind=dkind) :: bigrloc, zloc, rlocprim, bpolloc, bphiloc


			bigrloc = dbsval(th,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(surf_index,2,1:nsurf+1) )

			zloc = dbsval(th,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(surf_index,2,1:nsurf+1) )

			rlocprim = dbsder(1,th, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(surf_index,2,1:nsurf+1) )

			bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bpol)

			bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bphi)

			answer = -sqrt((bpolloc**2+bphiloc**2)/bzero(0.d0)**2)
			! define this way to "always" have an answer of order 1
			! and negative to exploit minimum search

			continue

			return

		end function modBval

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		function modBval_min(th) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			real(kind=dkind) :: th
			real(kind=dkind) :: answer

			answer = -modBval(th)

			continue

			return

		end function modBval_min


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine bootstrap_setup_ave(psi,rho, p, Br, Bz, T, b_phi, J_par,  &
					base, nn, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: psi,rho, p, Br, Bz, T, b_phi,  &
							J_par, base

	real(kind=dkind), dimension(1:nn,1:nn) :: thing
	real(kind=dkind), dimension(1:enq) :: boot_ptot, boot_ptotprim

	real(kind=dkind) :: base_option

	integer :: i,j



	! new average routine
	! this routine computes the averages of density (same for ions and electrons)
	! and of total pressure

	call get_new_avs(boot_ptot, boot_ne, base, p(1:nn,1:nn), rho(1:nn,1:nn), nn, base_option,  &
										boot_ptotprim, boot_neprim)

	if(bootstrap_option>=1) then

		call get_toroidal_flux(psi(1:nn,1:nn),b_phi(1:nn,1:nn))
		call averages_NCLASS(nn,base, base_option)

	endif

	! arrange the results for ions and elctrons
	! start from the density

	do i = 1,enq

		boot_ni(i) = boot_ne(i)
		boot_niprim(i) = boot_neprim(i)

	enddo

	! then take care of the pressures

	do i = 1,enq

		boot_pe(i) = boot_ptot(i)*pe_ov_p
		boot_pi(i) = boot_ptot(i)-boot_pe(i)

		boot_peprim(i) = boot_ptotprim(i)*pe_ov_p
		boot_piprim(i) = boot_ptotprim(i)-boot_peprim(i)

	enddo

	! finally, take care of the temperatures (no further interpolations needed)

	do i = 1,enq

		boot_Te(i) = boot_pe(i)/boot_ne(i)*mass/eV ! in eV
		boot_Ti(i) = boot_pi(i)/boot_ni(i)*mass/eV ! in eV

		boot_Teprim(i) = ( boot_peprim(i)/boot_ne(i) -  &
								boot_neprim(i)*boot_pe(i)/boot_ne(i)**2 )*mass/eV

		boot_Tiprim(i) = ( boot_piprim(i)/boot_ni(i) -  &
								boot_niprim(i)*boot_pi(i)/boot_ni(i)**2 )*mass/eV

	enddo

	! last thing, change the units of the densities

	do i = 1,enq

		boot_ne(i) = boot_ne(i)/mass ! now in m^-3
		boot_ni(i) = boot_ni(i)/mass ! now in m^-3

		boot_neprim(i) = boot_neprim(i)/mass ! now in m^-3
		boot_niprim(i) = boot_niprim(i)/mass ! now in m^-3

	enddo

	! now all the quantities for the bootstrap current calculation are available
	! still compute <J_par B>, to have the bootstrap fraction

	do j = 1,nn
	do i = 1,nn

		thing(i,j) = J_par(i,j) * sqrt( Br(i,j)**2 + Bz(i,j)**2 + b_phi(i,j)**2 )

	enddo
	enddo

	call get_thing_ave(JparB_ave, base, thing, nn, base_option)

	continue

	return

	end subroutine bootstrap_setup_ave

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_toroidal_flux(psi,b_phi)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind), dimension(1:n,1:n), intent(in) :: psi, b_phi
	real(kind=dkind) :: dA, bigR
	real(kind=dkind) :: bphiloc
	integer :: i,j,k
	real(kind=dkind), dimension(2,0:enq+1) :: temp_array !just to avoid problems with interpolation


	dA = dx*dz

	boot_tor_flux = 0.d0

	do j=1,n
		do i=1,n

			if(sort_grid(i,j,0)>=1) then

				do k=1,enq+1

					if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
											(dz_a(j) + dz_a(j-1))

					bphiloc = b_phi(i,j)

					if(psi(i,j)>psival(k)) then	!inner region for kth surface

						bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)
						boot_tor_flux(k) = boot_tor_flux(k) + bphiloc * dA

					endif

				enddo
	
			endif

		enddo
	enddo

	do k = 0,enq+1

		boot_tor_rho(k) = sqrt(boot_tor_flux(k)/boot_tor_flux(enq+1)) * a_elps
		!WARNING: a_elps = (R_max-R_min)/2 is defined elsewhere (readinput)

	enddo

	do k = 0,enq+1

		temp_array(1,enq+1-k) = psival(k)
		temp_array(2,enq+1-k) = boot_tor_rho(k)

	enddo


	allocate(bscoef_tor_rho(2,enq+2+q_ord))

	! interpolation setup has been moved to a separate function
	call interp_setup(enq+2, q_ord, &
		temp_array(1,0:enq+1),temp_array(2,0:enq+1), &
		bscoef_tor_rho(1,1:enq+2+q_ord),bscoef_tor_rho(2,1:enq+2))

	continue

	return

end subroutine get_toroidal_flux

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine averages_NCLASS(nn,base, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base
	real(kind=dkind) :: base_option
	real(kind=dkind) :: cutoff = 1.d-12

	real(kind=dkind), dimension(:,:), allocatable :: integrand
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, bigrlocprim, zlocprim,  &
								baseloc, bphiloc, der_val, bpolloc

	integer :: i, j, k


	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif


	! then take care of the derivative part
	! the integration procedure is kept here to avoid additional array copying

	allocate (integrand(1:2,1:n_ave_int))
	! "1" is Fhat, "2" is grad rho/B^2

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_boot_base)

			bphiloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bphi)

			bpolloc =  DBS2VL(bigrloc,zloc,ord_loc,ord_loc,xknot_mag, &
											zknot_mag,n,n,bscoef_bpol)

			!------------------normal and derivative------------------

			der_val = dbsder(1,psival(i), q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
					enq+2, bscoef_tor_rho(2,1:enq+2) )
			! note: this is (d rho / d psi), which is equal to 1/(d psi / d rho)

!			if(abs(der_val)<cutoff) der_val = sign(cutoff, der_val)

			integrand(1,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_val * baseloc**base_option  &
									*bphiloc / (2.d0*pi) !2pi? no  mu_mag!

			integrand(2,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * (der_val*bpolloc*bigrloc)**2/  &
									(bphiloc**2+bpolloc**2) *  &
									baseloc**base_option  &
									/ (2.d0*pi)


		enddo

		! then calculate the integrals

		call integrate(n_ave_int, integrand(1,:), w_ave_int, boot_fhat(i))

		boot_fhat(i) = boot_fhat(i) / surf_ave_base(i) * rmajor ! * pi !WARNING!!!

		call integrate(n_ave_int, integrand(2,:), w_ave_int, boot_grho_ov_B2(i))

		boot_grho_ov_B2(i) = boot_grho_ov_B2(i) / surf_ave_base(i)


		continue

	enddo



	deallocate(integrand)

	continue


	return


end subroutine averages_NCLASS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_bootstrap_NCLASS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: ierror
	real(kind=dkind) :: grad_T(2), Tloc(2), densloc(2,1), pprimloc(2,1)
	real(kind=dkind) :: den_cutoff
	real(kind=dkind), dimension(2) :: amu !atomic mass units
	real(kind=dkind) :: der_cutoff = 1.d-8
	real(kind=dkind) :: psiloc, der_val, der_sec_val, eps_loc
	real(kind=dkind) :: gradphi, gradphi_sq, EdotB, ngrth
	real(kind=dkind) :: PS_mom(3)
	integer :: i, k

	real(kind=dkind) :: external_force_fex_iz(3,mx_mi,mx_mz) = 0.d0
	integer, dimension(1:2) :: isotope_status = 0
	integer, dimension(1:2) :: charge_status
	integer :: k_potato = 0
	integer :: m_s = 2

	! see NCLASS documentation for the meaning of the following
	real(kind=dkind) :: p_etap, p_exjb
	real(kind=dkind) :: calm_i(3,3,mx_mi)
	real(kind=dkind) :: caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),  &
					capn_ii(3,3,mx_mi,mx_mi)

	real(kind=dkind) :: boot_p_coeff(mx_ms), boot_T_coeff(mx_ms)

	real(kind=dkind) :: dn_s(mx_ms), gfl_s(5,mx_ms), qfl_s(5,mx_ms), sqz_s(mx_ms),  &
								upar_s(3,3,mx_ms), utheta_s(3,3,mx_ms), vn_s(mx_ms), veb_s(mx_ms), &
								qeb_s(mx_ms), xi_s(mx_ms), ymu_s(3,3,mx_ms)
	real(kind=dkind) :: chip_ss(mx_ms,mx_ms), chit_ss(mx_ms,mx_ms),  &
								dp_ss(mx_ms,mx_ms), dt_ss(mx_ms,mx_ms)



	do i = 1,enq

		psiloc = psival(i)
		eps_loc = inv_asp_ratio(i)

		!------------ set geometry moments for PS current--------------
		! I have no clue of what this is (yet)

		do k=1,3
			PS_mom(k)=0.d0
		enddo

		if(p_eps.gt.0.0) then

			do k=1,3
				PS_mom(k)=k*((1.d0-sqrt(1.d0-eps_loc**2))/eps_loc)**(2.d0*k)  &
					*(1.d0 + k * sqrt(1.d0-eps_loc**2))/((1.d0-eps_loc**2)**1.5d0  &
					*(qval(i)*rmajor)**2)
			enddo

		!---------------end PS current stuff----------------

		endif

		amu(1) = 2.d0
		amu(2) = amu(1)*me_ov_mi !?

		charge_status(1) = 1
		charge_status(2) = -1

		den_cutoff = dofpsi(0.d0) / mass

		EdotB = 0.d0

		der_val = dbsder(1,psiloc, q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
						enq+2, bscoef_tor_rho(2,1:enq+2) )
				! note: this is (d rho / d psi), which is equal to 1/(d psi / d rho)

		if(abs(der_val)<der_cutoff) der_val = sign(der_cutoff, der_val)

		der_sec_val = dbsder(2,psiloc, q_ord, bscoef_tor_rho(1,1:enq+2+q_ord),  &
						enq+2, bscoef_tor_rho(2,1:enq+2) )
				! note: this is (d2 rho / d psi2), which is NOT equal to 1/(d2 psi / d rho2)

		if(abs(der_sec_val)<der_cutoff) der_sec_val = sign(der_cutoff, der_sec_val)


		gradphi = omegaofpsi(psiloc) / der_val

		gradphi_sq = ( domegadpsi(psiloc) - omegaofpsi(psiloc) * der_sec_val/der_val)  &
							/ der_val**2


		Tloc(1) = boot_Ti(i) / 1.d3 ! this is in keV
		Tloc(2) = boot_Te(i) / 1.d3 ! this is in keV

		grad_T(1) = boot_Tiprim(i) / der_val / 1.d3 ! d T / d rho
		grad_T(2) = boot_Teprim(i) / der_val / 1.d3 ! d T / d rho

		densloc(1,1) = boot_ni(i)
		densloc(2,1) = boot_ne(i)

		pprimloc(1,1) = boot_piprim(i) * mass/eV/1.d3 / der_val
		pprimloc(2,1) = boot_peprim(i) * mass/eV/1.d3 / der_val
		pprimloc(1,1) = boot_piprim(i) /eV/1.d3 / der_val	!?
		pprimloc(2,1) = boot_peprim(i) /eV/1.d3 / der_val	!?
		! in whatever ridiculous units these are supposed to be

		ngrth = 1.d0/qval(i)/rmajor	!what the heck is this?

		call NCLASS(2, k_potato, 2, 1, den_cutoff, 0.d0, 0.d0,  &
							B2_ave(i), Bm2_ave(i), EdotB, boot_fhat(i), PS_mom, eff_trap(i), boot_grho_ov_B2(i),  &
							gradphi, gradphi_sq, ngrth, amu(1:2), grad_T, Tloc,  &
							densloc, external_force_fex_iz, pprimloc, m_s, isotope_status, charge_status, J_boot(6,i),  &
							p_etap,p_exjb,calm_i,caln_ii,capm_ii,capn_ii,  &
							boot_p_coeff, boot_T_coeff ,dn_s,gfl_s,qfl_s,sqz_s,upar_s,  &
							utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,  &
							chit_ss,dp_ss,dt_ss, ierror)

		J_boot(6,i) = -J_boot(6,i) ! derivatives have the wrong sign!
		J_boot(7,i) = J_boot(6,i) / JparB_ave(i)

		el_resistivity(3,i) = p_etap ! p_etap is the resistivity

	enddo

! m_i = 2: ions and electrons?
! start without potatoes
! put neutral and charged D, 0 density for the neutrals
! <E dot B> = 0 (ideal MHD)

	continue

end subroutine get_bootstrap_NCLASS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_bootstrap2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: i
	real(kind=dkind) :: a2_el, a2_ion, a1
	real(kind=dkind) :: Coul_log_e, Coul_log_i ! Coulomb logarithm log(Lambda)
	real(kind=dkind) :: taue, vloc
	real(kind=dkind) :: p_nondim, pprim_nd, e_mass
	real(kind=dkind) :: charge_cgs = 4.8032d-10 ! to be fixed
	real(kind=dkind) :: nustar, nustarboh, nuee
	real(kind=dkind) :: ft, ft2, ft3
	real(kind=dkind) :: Z_eff = 1.d0 ! this would be Z_effective
	real(kind=dkind), dimension(:,:), allocatable :: nustar_array
	real(kind=dkind) :: sigma_neo, sigma_spitz
	real(kind=dkind) :: Z_i, n_i
	real(kind=dkind) :: el31, el32, el31_0, el32_0, boot_alpha, boot_alpha_0, el34


	e_mass = mass * me_ov_mi


	allocate(nustar_array(2,enq))

	do i = 1, enq

		Z_i = Z_eff !is this so?

!		Z_i = 1.d0+1.5d0*(abs(psival(i)/psic))**0.25d0

		Z_i = Z_eff_0 + delta_Z_eff*(abs(psival(i)/psic))**Z_eff_exp

		ft = eff_trap(i)
!		ft = eff_trap(i)/(1.d0+sqrt(nustar*Z_eff) + 0.25d0*nustar/Z_eff)

         call sigmaneo(sigma_neo,sigma_spitz,nustar,ft,boot_ne(i),boot_Te(i),Z_i, &
               abs(qval(i)),boot_Rcenter(i),inv_asp_ratio(i))

		el_resistivity(1,i) = 1.d0/sigma_spitz
		el_resistivity(2,i) = 1.d0/sigma_neo

		nustar_array(1,i) = nustar

        nustar = nustar/Z_i

		nustar_array(2,i) = nustar

! WARNING! NOTABLE ASSUMPTION!!
! assumes z_imp=6
         n_i=(6.d0-Z_i)/(6.d0-1.d0) * boot_ne(i)
         call bootstrap_coeffs(ft,abs(qval(i)),boot_rcenter(i),inv_asp_ratio(i),boot_Te(i),boot_ne(i),boot_Ti(i), &
               n_i,Z_i,1.d0,el31_0,el32_0,boot_alpha_0,el31,el32,el34,boot_alpha)


!     For RJBSOS(K,4), use the standard definition. However one mixes p' from chease inputs and ne', Te', Ti' from
!     experimental data. If ne*Te+ni*Ti exper. is not near p_chease, it can give strange profiles, as some cancellation do not appear.
!     Thus develop p'=ne' Te + ne Te' + ni' Ti + ni Ti', and assume ni'/ni = ne'/ne and take ni*Ti=p-pe:
!     => for (K,1-3) use developed formula as well as for nue*=0 case
!!$         RJBSOS(K,1) = - TMF(K)* &
!!$     &     (ZL31_0*(ZA1+ZALFA_0*(1.d0/ZRPEOP-1.d0)*ZA2I) + ZL32_0 * ZA2E)

		! zero nustar case, in physical units
         J_boot(1,i) = rmajor*bzero(psival(i)) * (boot_pi(i)+boot_pe(i)) * &
           (el31_0*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31_0+el32_0)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(1.d0+boot_alpha)*boot_Tiprim(i)/boot_Ti(i))

		! finite nustar case, in physical units
         J_boot(2,i) = rmajor*bzero(psival(i)) * (boot_pi(i)+boot_pe(i)) * &
           (el31*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31+el32)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(el31+el34*boot_alpha)*boot_Tiprim(i)/boot_Ti(i))

		! bootstrap fraction = bootstrap / <J_par B>
		J_boot(3,i) = J_boot(2,i) / JparB_ave(i)

		! zero nustar case CHEASE
		J_boot(4,i) = rmajor*bzero(psival(i)) * boot_ne(i)*boot_Te(i)/pe_ov_p* eV*mu_mag/abs(b_phi_zero)**3 * & !1.602d-19*4.d-07*CPI/B0EXP**2 * &
           (el31_0*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31_0+el32_0)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(1.d0+boot_alpha_0)*el31_0*boot_Tiprim(i)/boot_Ti(i))

		! finite nustar case CHEASE
		J_boot(5,i) = rmajor*bzero(psival(i)) * boot_ne(i)*boot_Te(i)/pe_ov_p* eV*mu_mag/abs(b_phi_zero)**3 *  & !1.602d-19*4.d-07*CPI/B0EXP**2 * &
           (el31*boot_neprim(i)/boot_ne(i) + pe_ov_p*(el31+el32)*boot_Teprim(i)/boot_Te(i) + &
           (1.d0-pe_ov_p)*(el31+el34*boot_alpha)*boot_Tiprim(i)/boot_Ti(i))	!ratio 1.1901e-005?


		continue

	enddo


	deallocate(nustar_array)

	return

	end subroutine get_bootstrap2

!-----------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bootstrap_coeffs(pft,pq,pR,peps,pte,pne,pti,pni,pzeff,pzion, &
      pl31_0,pl32_0,palfa_0,pl31,pl32,pl34,palfa)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

!
!     WARNING: in MKSA
!
!     Compute Bootstrap coefficients using formulas from O. Sauter et al, Phys. Plasmas 7 (1999) 2834.
!
!     Assumes to compute on a single flux surface with:
! Inputs:
!     pft   : trapped fraction
!     pq    : safety factor
!     pR    : Geometrical center of given flux surface in [m]
!     peps  : Inverse aspect ratio of given flux surface
!     pte   : Electron temperature [eV]
!     pne   : Electron density [1/m**3]
!     pti   : Ion temperature [eV]
!     pni   : Main ion density [1/m**3]
!     pzeff : Effective charge (used for Z in electronic terms)
!     pzion : Main ion charge
! Outputs:
!     pl31_0  : L31 coefficient assuming nuestar=0
!     pl32_0  : L32 coefficient assuming nuestar=0
!     palfa_0 : Alfa coefficient assuming nuestar=0
!     pl31    : L31 coefficient
!     pl32    : L32 coefficient
!     pl34    : L34 coefficient (L34 for nuestar=0 is identical to L31_0)
!     palfa   : Alfa coefficient
!

  implicit none

  real(kind=dkind), intent(in) :: pft,pq,pR,peps,pte,pne,pti,pni,pzeff,pzion
  real(kind=dkind), intent(out) :: pl31_0,pl32_0,palfa_0,pl31,pl32,pl34,palfa

  real(kind=dkind)  :: znuestar, znuistar, zlnlam_e, zlnlam_i, zdummy
!-----------------------------------------------------------------------
   
!     basic parameters

  zlnlam_e = 17.d0
  IF ((pne.gt.0.d0).and.(pte.gt.0.d0)) then
     zlnlam_e = 31.3d0 - log(sqrt(pne)/pte)
     zlnlam_e = max(10.d0,zlnlam_e)
  ENDIF
  zlnlam_i = 17.d0
  IF ((pni.gt.0.d0).and.(pti.gt.0.d0)) then
     zlnlam_i = 30.d0 - log(pzion**3.d0*sqrt(pni)/abs(pti)**1.5d0)
     zlnlam_i = max(10.d0,zlnlam_i)
  ENDIF
  znuestar = 6.921d-18 * pq*pR*pne*pzeff*zlnlam_e / (pte*pte*peps**1.5d0)
  znuistar = 4.900d-18 * pq*pR*pni*pzion**4*zlnlam_i / (pti*pti*peps**1.5d0)

!     Compute coefficients for nustar=0

  call final_bootstrap_coeffs(pl31_0,pl32_0,zdummy,palfa_0,pft,pzeff)

!     finite nustar
  call final_bootstrap_coeffs(pl31,pl32,pl34,palfa,pft,pzeff,znuestar,znuistar)

	continue

  return

  end subroutine bootstrap_coeffs

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine final_bootstrap_coeffs(L31, L32, L34, ALFA, ft, Zeff, nuestar, nuistar)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

    real(kind=dkind), intent(in)  :: ft
    real(kind=dkind), OPTIONAL, intent(in)  :: Zeff, nuestar, nuistar
    real(kind=dkind), intent(out) :: L31, L32, L34, ALFA

    real(kind=dkind)  :: ZZ, znuestar, znuistar, zsqnuest, zeffp1, zsqnui, znui2ft6
    real(kind=dkind)  :: zft31eff, zft32ee_eff, zft32ei_eff, zft34eff, zalfa0
!-----------------------------------------------------------------------
!
    ZZ = 2.d0

    if ( PRESENT(zeff) ) ZZ = zeff

    znuestar = 0.d0

    if ( PRESENT(nuestar) ) znuestar = nuestar

    znuistar = 0.d0

    if ( PRESENT(nuistar) ) znuistar = nuistar

    zsqnuest = sqrt(znuestar)

!  effective trapped fractions

    zsqnuest = sqrt(znuestar)

    zft31eff = ft / (1d0+(1.d0-0.1d0*ft)*zsqnuest &
           + 0.5d0*(1.d0-ft)*znuestar/ZZ)

    zft32ee_eff = ft / (1.d0 + 0.26d0*(1.d0-ft)*zsqnuest &
           + 0.18d0*(1.d0-0.37d0*ft)*znuestar/sqrt(ZZ))

    zft32ei_eff = ft / (1.d0 + (1.d0+0.6d0*ft)*zsqnuest &
           + 0.85d0*(1.d0-0.37d0*ft)*znuestar*(1.d0+ZZ))

    zft34eff = ft / (1.d0+(1.d0-0.1d0*ft)*zsqnuest &
           + 0.5d0*(1.d0-0.5d0*ft)*znuestar/ZZ)

    zalfa0 = - 1.17d0*(1.d0-ft) / (1.d0-0.22d0*ft-0.19d0*ft**2)

!coefficients

    zeffp1 = ZZ+1.d0

    L31 = zft31eff * ( (1.d0+1.4d0/zeffp1) &
          - zft31eff* (1.9d0/zeffp1 - zft31eff * (0.3d0/zeffp1 + 0.2d0/zeffp1 * zft31eff)))

    L32 = (0.05d0+0.62d0*ZZ)/ZZ/(1.d0+0.44d0*ZZ)*(zft32ee_eff-zft32ee_eff**4) &
          +  zft32ee_eff**2*(1.d0-1.2d0*zft32ee_eff+0.2d0*zft32ee_eff**2) &
                           /(1.d0+0.22d0*ZZ) &
          - (0.56d0+1.93d0*ZZ)/ZZ/(1.d0+0.44d0*ZZ)*(zft32ei_eff-zft32ei_eff**4) &
          +  zft32ei_eff**2*(1.d0-0.55d0*zft32ei_eff-0.45d0*zft32ei_eff**2) &
                           * 4.95d0/(1.d0+2.48d0*ZZ) &
          + 1.2d0 / (1.d0+0.5d0*ZZ) * (zft32ee_eff**4-zft32ei_eff**4)

    L34 = zft34eff * ( (1.d0+1.4d0/zeffp1) &
          - zft34eff* (1.9d0/zeffp1 - zft34eff * (0.3d0/zeffp1 + 0.2d0/zeffp1 * zft34eff)))

    zsqnui = sqrt(znuistar)
    znui2ft6 = znuistar**2 * ft**6

    ALFA = ((zalfa0 + 0.25d0*(1.d0-ft**2)*zsqnui) &
                             / (1.d0+0.5d0*zsqnui) + 0.315d0*znui2ft6) &
          / (1.d0 + 0.15d0*znui2ft6)

	continue

    return

  end subroutine final_bootstrap_coeffs


  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine sigmaneo(signeo,sigsptz,nuestar,ft,ne,te,zeff,q,R,eps)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! IMPORTANT: ADAPTED FROM THE CHEASE CODE (THANKS TO O. SAUTER)

    ! all are scalar input variables
    ! te in [eV], ne in [m**-3], R in [m]

	implicit none

    real(kind=dkind), intent(out) :: signeo, sigsptz, nuestar
    real(kind=dkind), intent(in)  :: ft, ne, te
    real(kind=dkind), OPTIONAL, intent(in) :: zeff
    real(kind=dkind), OPTIONAL, intent(in) :: q
    real(kind=dkind), OPTIONAL, intent(in) :: R
    real(kind=dkind), OPTIONAL, intent(in) :: eps

    real(kind=dkind)  :: z_zeff, zNZ, zlnL, zft33eff

    z_zeff = 2.d0
    IF ( PRESENT(zeff) ) z_zeff = zeff

    zNZ = 0.58d0 + 0.74d0 / (0.76d0 + z_zeff)
    zlnL = 17.d0

    IF (ne.gt.0.d0 .and. te.gt.0.d0) THEN
       zlnL = 31.3d0 - log(sqrt(ne)/te)
    ENDIF

    sigsptz = 1.9012d4 * abs(te)**1.5d0 / (z_zeff * zNZ * zlnL)

    nuestar = 0.01d0

    IF (PRESENT(q) .AND. PRESENT(R) .AND. PRESENT(eps)) &
          nuestar = 6.921d-18 * q * R * ne * z_zeff * zlnL / (te*te * eps**1.5d0)

    zft33eff = ft / &
          (1.d0+(0.55d0-0.1d0*ft)*sqrt(nuestar) &
           + 0.45d0*(1.d0-ft)*nuestar/z_zeff**1.5d0)

    signeo = sigsptz * (1.d0 - zft33eff*(1.d0+0.36d0/z_zeff &
                                 - zft33eff*(0.59d0/z_zeff - 0.23d0/z_zeff*zft33eff)))

  end subroutine sigmaneo


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_JparB(JparB_ave, bpol, bphi, J_par, nn)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: bpol, bphi, J_par
	real(kind=dkind), dimension(1:enq) :: JparB_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing

	integer :: i,j

	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(bpol, nn, -1.d0)

	endif

	! then create the integrand thing, and give it to the averaging routine

	do j = 1,nn
	do i = 1,nn

		thing(i,j) = J_par(i,j) * sqrt( bpol(i,j)**2 + bphi(i,j)**2 )

	enddo
	enddo

	call surf_ave(bpol, thing, nn, -1.d0, JparB_ave)

	continue

	return

end subroutine get_JparB

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_thing_ave(thing_ave, base, thing, nn, base_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBS2DR

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base, thing
	real(kind=dkind), dimension(1:enq) :: thing_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing_to_pass
	real(kind=dkind) :: base_option

	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif

	! then call the averaging routine

	call surf_ave(base, thing, nn, base_option, thing_ave)

	return

end subroutine get_thing_ave

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine surf_ave_setup(ave_base, nn, ave_option)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: ave_base
	real(kind=dkind), dimension(:), allocatable :: integrand
	real(kind=dkind) :: ave_option !integral of ave_base^ave_option

	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, baseloc

	integer :: i, k

	! first set up the various arrays (note that the weighting stuff will stay)

	allocate(w_ave_int(1:n_ave_int))
	allocate(t_ave_int(1:n_ave_int))

	call set_weights(n_ave_int,0.d0,2.d0*pi,w_ave_int,t_ave_int)

	allocate (integrand(1:n_ave_int))

	allocate(surf_ave_base(1:enq))

	allocate(bscoef_boot_base(1:nn,1:nn))

	! then set up the interpolation of the incoming quantity

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_base,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_boot_base)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_boot_base)

			integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * baseloc**ave_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integral

		call integrate(n_ave_int, integrand, w_ave_int, surf_ave_base(i))

		continue

	enddo



	deallocate(integrand)

	continue

	return

end subroutine surf_ave_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine surf_ave(ave_base, ave_thing, nn, ave_option, ave_result)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use pseudo_IMSL, only : DBS2IN, DBS2VL

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: ave_base, ave_thing
	real(kind=dkind), dimension(:), allocatable :: integrand
	real(kind=dkind), dimension(1:enq) :: ave_result
	real(kind=dkind) :: ave_option !integral of ave_thing * ave_base^ave_option

	real(kind=dkind), dimension(1:nn,1:nn) :: bscoef_thing
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, baseloc, thingloc

	integer :: i, k

	! first set up the various arrays (note that the weighting stuff is given)

	allocate (integrand(1:n_ave_int))

	! then set up the interpolation of the incoming quantities

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_base,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_boot_base)

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),ave_thing,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_thing)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_boot_base)

			thingloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,n,n,bscoef_thing)

			integrand(k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * thingloc * baseloc**ave_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integral

		call integrate(n_ave_int, integrand, w_ave_int, ave_result(i))

		ave_result(i) = ave_result(i) / surf_ave_base(i)

		continue

	enddo



	deallocate(integrand)

	continue

	return

end subroutine surf_ave




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_new_avs(p_ave, dens_ave, base, pres, dens, nn, base_option,  &
										p_prim_ave, dens_prim_ave)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! note: only electron pressure and density get here; ion's are proportional to them

	use pseudo_IMSL, only : DBS2IN, DBS2VL, DBS2DR

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn), intent(in) :: base, pres, dens
	real(kind=dkind), dimension(1:enq) :: p_ave, dens_ave
	real(kind=dkind), dimension(1:nn,1:nn) :: thing_to_pass
	real(kind=dkind), dimension(1:enq) :: p_prim_ave, dens_prim_ave
	real(kind=dkind) :: base_option

	real(kind=dkind), dimension(:,:), allocatable :: integrand
	real(kind=dkind), dimension(1:nn,1:nn) :: bscoef_pres, bscoef_dens
	real(kind=dkind) :: thetaloc, rloc, bigrloc, zloc, rlocprim, bigrlocprim, zlocprim,  &
								baseloc, thingloc, pres_U, pres_L, dens_U, dens_L,  &
								psi_U, psi_L, cos_n, sin_n, der_dens, der_pres

	real(kind=dkind), dimension(1:2) :: U_point, L_point
	! points on the normal to the magnetic surface

	real(kind=dkind) :: dist

	integer :: i, j, k


	! first calculate the denominator in the average, if needed

	if (allocated(surf_ave_base)) then

		continue

	else

		call surf_ave_setup(base, nn, base_option)

	endif

	! then call the averaging routine for pressure and density

	call surf_ave(base, pres, nn, base_option, p_ave)
	call surf_ave(base, dens, nn, base_option, dens_ave)

	! then take care of the derivative part
	! the integration procedure is kept here to avoid additional array copying

	dist = sqrt(dx**2+dz**2)/10.d0

	allocate (integrand(1:2,1:n_ave_int))
	! "1" is pressure, "2" is density

	! then set up the interpolation of the incoming quantities

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),pres,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_pres)

	call DBS2IN(nn,x_coord(1:nn),nn,z_coord(1:nn),dens,nn,  &
					ord_ave,ord_ave,xknot_mag,zknot_mag,bscoef_dens)

	! cycle over the magnetic surfaces

	do i=1,enq

		! first set up the integrand (Jacobian included)

		do k = 1, n_ave_int

			thetaloc = t_ave_int(k)

			rloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zloc = dbsval(thetaloc,ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			rlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, rtab(i,2,1:nsurf+1) )

			bigrlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, bigrtab(i,2,1:nsurf+1) )

			zlocprim = dbsder(1,thetaloc, ord_surf, thetatab(2,:),  &
					nsurf+1, ztab(i,2,1:nsurf+1) )

			baseloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_boot_base)

			!------------------normal and derivative------------------
			cos_n = zlocprim/sqrt(zlocprim**2+bigrlocprim**2)
			sin_n = -bigrlocprim/sqrt(zlocprim**2+bigrlocprim**2)

			U_point(1) = bigrloc + dist*cos_n
			U_point(2) = zloc + dist*sin_n

			L_point(1) = bigrloc - dist*cos_n
			L_point(2) = zloc - dist*sin_n

			psi_L = DBS2VL(L_point(1),L_point(2),ord_loc,ord_loc,xknot_mag, &
											zknot_mag,nn,nn,bscoef_psi)

			psi_U = DBS2VL(U_point(1),U_point(2),ord_loc,ord_loc,xknot_mag, &
											zknot_mag,nn,nn,bscoef_psi)

			pres_L = DBS2VL(L_point(1),L_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_pres)

			pres_U = DBS2VL(U_point(1),U_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_pres)

			dens_L = DBS2VL(L_point(1),L_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_dens)

			dens_U = DBS2VL(U_point(1),U_point(2),ord_ave,ord_ave,xknot_mag, &
											zknot_mag,nn,nn,bscoef_dens)

			der_pres = (pres_U-pres_L)/(psi_U-psi_L)

			der_dens = (dens_U-dens_L)/(psi_U-psi_L)


!			thingloc =  DBS2VL(bigrloc,zloc,ord_ave,ord_ave,xknot_mag, &
!											zknot_mag,n,n,bscoef_thing)

			integrand(1,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_pres * baseloc**base_option  &
									/ (2.d0*pi)

			integrand(2,k) = sqrt(1.d0+(rlocprim/rloc)**2) *  &
									rloc * der_dens * baseloc**base_option  &
									/ (2.d0*pi)

		enddo

		! then calculate the integrals

		call integrate(n_ave_int, integrand(1,:), w_ave_int, p_prim_ave(i))

		p_prim_ave(i) = p_prim_ave(i) / surf_ave_base(i)

		call integrate(n_ave_int, integrand(2,:), w_ave_int, dens_prim_ave(i))

		dens_prim_ave(i) = dens_prim_ave(i) / surf_ave_base(i)

		continue

	enddo



	deallocate(integrand)

	continue


	return

end subroutine get_new_avs


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine bootstrap_cleanup
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


	deallocate(psival)
	deallocate(qval)

	deallocate(J_boot)
	deallocate(eff_trap)
	deallocate(surf_length)
	if(allocated(B2_ave)) deallocate(B2_ave)
	if(allocated(Bm2_ave)) deallocate(Bm2_ave)
	deallocate(B2_hat_ave)
	deallocate(boot_ne)
	deallocate(boot_pe)
	deallocate(boot_Te)
	deallocate(boot_neprim)
	deallocate(boot_peprim)
	deallocate(boot_Teprim)
	deallocate(boot_ni)
	deallocate(boot_pi)
	deallocate(boot_Ti)
	deallocate(boot_niprim)
	deallocate(boot_piprim)
	deallocate(boot_Tiprim)
	if(allocated(boot_tor_flux)) deallocate(boot_tor_flux)
	if(allocated(boot_tor_rho)) deallocate(boot_tor_rho)
	if(allocated(boot_fhat)) deallocate(boot_fhat)
	if(allocated(boot_grho_ov_B2)) deallocate(boot_grho_ov_B2)
	deallocate(inv_asp_ratio)
	deallocate(boot_Bmin)
	deallocate(boot_Bp_at_Bmin)
	deallocate(boot_R_of_Bmin)
	deallocate(boot_Rcenter)
	deallocate(el_resistivity)

	deallocate(Rleft)
	deallocate(Rright)

	deallocate(rtab)
	deallocate(bigrtab)
	deallocate(ztab)
	deallocate(thetatab)

	deallocate(xknot_mag)
	deallocate(zknot_mag)

	deallocate(surf_ave_base)
	deallocate(JparB_ave)

	deallocate(bscoef_boot_base)

end subroutine bootstrap_cleanup


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine pvgamma(nn,psi)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	integer :: nn
	real(kind=dkind),dimension(1:nn,1:nn) :: psi
	integer, parameter :: nsurf=50
	real(kind=dkind), dimension(1:4,1:nsurf) :: pvg
	integer :: i,j,k
	real(kind=dkind) :: dA, bigR

	pvg = 0.d0

	do i=1,nsurf
		pvg(1,i) = psi_out + (psi_in-psi_out)/(nsurf*1.d0)*(i-1.d0)
	enddo

	do i=1,nsurf
		pvg(2,i) = pofpsi(pvg(1,i))
	enddo

	do j=1,nn
		do i=1,nn

		    if(sort_grid(i,j,0)==1) then

				if((psi(i,j)/psic)>=fraction) then	!inner region
!				if(dabs(psi(i,j)/psic)>=fraction) then	!inner region

				bigR = x_coord(i) + dx_a(i) -  dx_a(i-1)

				if(grid_type/=0) dA = 0.25d0*(dx_a(i) + dx_a(i-1)) *  &
										(dz_a(j) + dz_a(j-1))

				do k=1,nsurf

					if(psi(i,j)>=pvg(1,k)) pvg(3,k) = pvg(3,k) + bigR*dA

				enddo

				endif

			endif

		enddo
	enddo

	do i=1,nsurf
		pvg(4,i) = pvg(2,i)*pvg(3,i)**gamma
	enddo

	open(33,file='p_and_V.dat',action='write',status='unknown')

	do i = 1,nsurf

		write(33,88) pvg(1,i),pvg(2,i),pvg(3,i),pvg(4,i)

	enddo

	close(33)

88	format(4(e23.14, 3x))

	continue

	return

	end subroutine pvgamma

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine grid_delta(aa,f,df)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	real(kind=dkind)  aa, f, df

	if(grid_type==1) then

		f = g_ratio * (1.d0-aa**((n_temp-1)/2)) - 1.d0 + aa**(n_temp-1)

		df = (n_temp-1.d0) * aa**(n_temp-2) - g_ratio*((n_temp-1.d0)/2.d0) * aa**((n_temp-3)/2)

	elseif(grid_type==2) then

		f = n_temp*g_ratio/(1.d0-g_ratio) * (1.d0-aa)*aa**(n_temp-1) - (1.d0-aa**n_temp)

!		df = n_temp*aa**(n_temp-1) * (1.d0 + g_ratio/(1.d0-g_ratio) * (&
!				(n_temp-3.d0)/aa**3.d0 - (n_temp-2.d0)/aa**2 ) )

		df = n_temp*aa**(n_temp-1) * (1.d0 + g_ratio/(1.d0-g_ratio) * (&
				(n_temp-1.d0)/g_ratio - n_temp ) )

	endif

	continue

	return

	end subroutine grid_delta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine grid_delta_2(aa,f,df)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	implicit none

	real(kind=dkind)  aa, f, df

	f = aa**n_temp - 1.d0 - g_ratio * (aa-1.d0)

	df = n_temp * aa**(n_temp-1)-g_ratio

	continue

	return

	end subroutine grid_delta_2



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_position(i,j,bound,truebound,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		use constant, only : dx, dz

		integer :: i,j,n
		logical :: bound, truebound
		real(kind=dkind) :: ex,ez


		truebound = .false.

! i,j
		if(sort_grid(i,j,0)==0) then
		! the point is actually on the boundary

			truebound = .true.
			return

		endif

		if((sort_grid(i,j,0)==1).or.(sort_grid(i,j,0)==2)) then

			bound = .false.
			return

		endif

! i,j+1
		if(sort_grid(i,j+1,0)==1) then

			bound = .true.
			return

		endif

! i,j-1
		if(sort_grid(i,j-1,0)==1) then

			bound = .true.
			return

		endif

! i+1,j
		if(sort_grid(i+1,j,0)==1) then

			bound = .true.
			return

		endif

! i-1,j
		if(sort_grid(i-1,j,0)==1) then

			bound = .true.
			return

		endif

		bound = .false.

		continue

		return

	end subroutine check_position

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine point1(i,j,xQ,zQ,Qind)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! finds point1 (lower-left corner of element enclosing xQ,zQ) and puts its i,j coordinates in Qind

		use constant, only : x_coord, z_coord, dx_a, dz_a

		integer :: i,j
		real(kind=dkind) :: xQ,zQ
		integer, dimension(1:2) :: Qind
		integer :: ii,jj

		ii = 0
		jj = 0

		do

!			if( (abs((rmajor - 0.5d0*x_size + (ii-1)*dx)-xQ)<=dx*(1+1.d-9)).and.  &
!				(abs((rmajor - 0.5d0*x_size + ii*dx)-xQ)<=dx*(1+1.d-9))		)then
			if( (abs(x_coord(ii)-xQ)<=dx_a(ii)*(1.d0+1.d-9)).and.  &
				(abs(x_coord(ii+1)-xQ)<=dx_a(ii)*(1.d0+1.d-9))		)then

				Qind(1) = ii
				exit

			endif

			ii = ii + 1

			if((ii>i+2).and.(grid_type==0)) then
				print*, 'error in point1, i', i, j, ii
!				pause 'error in point1, i'
				continue
			endif

		enddo

		do

!			if( (abs((-0.5d0*z_size + (jj-1)*dz)-zQ)<=dz*(1+1.d-9)).and.  &
!				(abs((-0.5d0*z_size + (jj-1)*dz)-zQ)<=dz*(1+1.d-9))	) then
			if( (abs(z_coord(jj)-zQ)<=dz_a(jj)*(1.d0+1.d-9)).and.  &
				(abs(z_coord(jj+1)-zQ)<=dz_a(jj)*(1.d0+1.d-9))	) then

				Qind(2) = jj
				exit

			endif

			jj = jj + 1

			if((jj>j+2).and.(grid_type==0)) then
!				pause 'error in point1, j'
				print*, 'error in point1, j', i, j, jj
				continue
			endif

		enddo

		continue

		return

	end subroutine point1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_point1(i,j,inner)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		integer :: i,j
		logical :: inner

		inner = .true.

! i,j
		if(sort_grid(i,j,0)<=0) then

!			print*, 'error, point 1 is external!'

			inner = .false.
			return

		endif

! i,j+1
		if(sort_grid(i,j+1,0)<=0) then

			inner = .false.
			return

		endif

! i+1,j+1
		if(sort_grid(i+1,j+1,0)<=0) then

			inner = .false.
			return

		endif

! i+1,j
		if(sort_grid(i+1,j,0)<=0) then

			inner = .false.
			return

		endif

		continue

		return

	end subroutine check_point1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_points(i,j,n,dx,dz,inner)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		integer :: i,j,n
		logical :: inner
		real(kind=dkind) :: ex,ez,dx,dz
		integer :: how_many

		how_many = 0

! i,j
		if(sort_grid(i,j,0)<=0) then

			how_many = how_many + 1

		endif

! i,j+1
		if(sort_grid(i,j+1,0)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

! i+1,j+1
		if(sort_grid(i+1,j+1,0)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

! i+1,j
		if(sort_grid(i+1,j,0)<=0) then

			how_many = how_many + 1

			if(how_many>1) then

				inner = .false.
				return

			endif

		endif

		inner = .true.

		continue

		return

	end subroutine check_points

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function tri_area(p1,p2,p3) result(area)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	real(kind=dkind), dimension(1:2) :: p1, p2, p3
	real(kind=dkind) :: area

	area = abs( p1(1)*p2(2) + p2(1)*p3(2) + p3(1)*p1(2)  &
					- p2(1)*p1(2)- p3(1)*p2(2)- p1(1)*p3(2) )/2.d0

	continue

	return

end function tri_area

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine find_pmax(x,f,df)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		implicit none

		real(kind=dkind) :: x,f,df

		f = x**(gammahut-1.d0) *  &
				(1.d0-x)**(aux_fact-1.d0)*  &
				(gammahut-(gammahut+aux_fact)*x)

		df = (1.d0-x)**(aux_fact-2.d0) *  &
				x**(gammahut-2.d0) *  &
				( (gammahut-1.d0)*gammahut -  &
				 2.d0*gammahut*(aux_fact+gammahut-1.d0)*x +  &
				 (aux_fact+gammahut)*(aux_fact+gammahut-1.d0)*x**2 )

		continue

		return

	end subroutine find_pmax

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_interface(psi,n,inorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! This routine find the curve psi=0 inside the computational domain for free-boundary equilibrium.
! The code will get here only for tri_type==13 (obsolete) or tri_type==-2

	use exp_data, only : r_data_psi, xknot_psi, zknot_psi, psi_bscoef
	use triangularity, only : theta_temp, r_ord, r_data
	use pseudo_IMSL, only : DBSNAK, DBSINT

	implicit none

	integer :: i, k, nth, alloc_stat, error
	real(kind=dkind) :: rg1, rg2, smallr_in, smallr_out, r_orp
	integer :: n
	real(kind=dkind), dimension(1:n,1:n) :: psi
	real(kind=dkind) :: inorm
	real(kind=dkind) :: fnew
	real(kind=dkind) :: psi_sol
	logical :: look_for_r

	external psi_sol

	r_orp = 0.05d0

	if(tri_type==13) then

		! first update fraction
		if (n<=inter_switch) then

			fnew = psi_e/psic
			fraction = r_orp*fnew + (1.d0-r_orp)*fraction

		else

			fnew = fraction * psic
			psi_e = r_orp*fnew + (1.d0-r_orp)*psi_e

		endif

	endif

	if (n<inter_switch) return

	inorm = 0.d0

	call psi_interp_setup(psi,n)

	nth = theta_points2

	allocate (r_data_psi(nth+r_ord,3),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in r_data_psi"
		 pause
		 stop
	endif

	do i=1,nth

		theta_temp = (i-1.d0)/(nth-1.d0)*2.d0*pi

		k = 0

		smallr_in = dbsval(theta_temp, r_ord, r_data(1:theta_points2+r_ord,6),  &
			theta_points2, r_cscoef(2,1:theta_points2) )

		smallr_out = dbsval(theta_temp, r_ord, r_data(1:theta_points1+r_ord,3),  &
			theta_points1, r_cscoef(1,1:theta_points1) )

!!$		if(i==1) then
!!$			rg1 = 0.d0
!!$			rg2 = smallr1 * 1.05d0
!!$		else
!!$			rg1 = 0.d0
!!$			rg2 = r_data_psi(i-1,2)*1.09d0
!!$			if(rg2>=smallr2) rg2 = 0.99d0*smallr2
!!$		endif

		rg1 = 0.d0
		rg2 = smallr_in * .5d0

		r_data_psi(i,1) = theta_temp

		look_for_r = .true.

		do while(look_for_r)

			call secroot(psi_sol,rg1,rg2,r_data_psi(i,2),  &
						100,1.d-9,1.d-9, error)

			if(error==0) then

				look_for_r = .false.

			elseif(error==1) then
			! secroot failed because the root was not bracketed, try again
	!!$			rg2 = smallr2 * 0.99**k
	!!$			k = k+1
				rg2 = rg2 * 1.01d0

	!			if((rg2<0.5d0*smallr_in).or.(rg2>0.995d0*smallr_out)) then
				if(rg2>0.995d0*smallr_out) then
	!				r_data_psi(i,2) = smallr1
					r_data_psi(i,2) = smallr_out ! Changed to outer radius on February 3 2022
					print*, 'problem in update_interface, theta = ', theta_temp
					! no need to say this should NOT happen
					look_for_r = .false.
				endif

			endif

		enddo

71		continue

	enddo

	! relaxation process

	do i = 1,nth

		smallr_in = dbsval(r_data_psi(i,1), r_ord, r_data(1:theta_points2+r_ord,6),  &
			theta_points2, r_cscoef(2,1:theta_points2) )

		inorm = max(inorm,dabs( (r_orp*r_data_psi(i,2) + (1.d0-r_orp)*smallr_in)/smallr_in )-1.d0 )

		r_data_psi(i,2) = r_orp*r_data_psi(i,2) + (1.d0-r_orp)*smallr_in

	enddo

	! end of relaxation process

	do i=1,theta_points2 + r_ord
	! just in case

		r_data(i,4) = 0.d0
		r_data(i,5) = 0.d0
		r_data(i,6) = 0.d0

	enddo

	r_cscoef(2,:) = 0.d0

	do i=1,theta_points2

		r_data(i,4) = r_data_psi(i,1)
		r_data(i,5) = r_data_psi(i,2)

	enddo

	call DBSNAK(theta_points2, r_data(1:theta_points2,4),  &
					r_ord,r_data(1:theta_points2+r_ord,6))

	call DBSINT(theta_points2, r_data(1:theta_points2,4),  &
		 r_data(1:theta_points2,5), r_ord,  &
		 r_data(1:theta_points2+r_ord,6),  &
		 r_cscoef(2,1:theta_points2))

	deallocate(r_data_psi)
	deallocate(psi_bscoef)
	deallocate(xknot_psi)
	deallocate(zknot_psi)

	print*, 'interface error:', inorm
	print*,'      '

	if(tri_type==13) then
		! update the grid
		call set_sort_grid(n,n)
	endif
	! Note that this won't be necessary for tri_type=-2, because update_sort_grid is called after update_interface

	continue

	return

end subroutine update_interface


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine initialize_r_free_boundary(psi,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! This routine find the curve psi=0 inside the computational domain for free-boundary equilibrium.
! The location of the curve is initialized here. If there is no closed psi=0 surface, the routine will return
! the plasma boundary as the correspoinding point on the psi=0 curve. Note that contrary to the original
! routine we will always assume that psi=0 is the plasma edge.
! The external radius is also defined here. The way this is done is taken from initialize_bc_equations in FLOW2.
! The code will get here only for tri_type==-2.

	use exp_data, only : r_data_psi, xknot_psi, zknot_psi, psi_bscoef
	use triangularity, only : theta_temp, r_ord, r_data
	use pseudo_IMSL, only : DBSNAK, DBSINT

	implicit none

	integer :: i, k, nth, alloc_stat, error
	real(kind=dkind) :: rg1, rg2, smallr1, smallr2
	integer :: n
	real(kind=dkind), dimension(1:n,1:n) :: psi
	real(kind=dkind) :: fnew
	real(kind=dkind) :: psi_sol

	external psi_sol

	! First initialize the outer boundary

!		if(bc_type==17) then
	! This needs to happen here to avoid cross references

		! big_Psi = psi on the boundary + toroidal flow term assigned through omega_0 and d(Psi)

		! We need to determine the limiting angles corresponding to the corners of the domain first.

!!$		theta_one = atan(b_elps/a_elps)
!!$		theta_two = pi - theta_one
!!$		theta_three = pi + theta_one
!!$		theta_four = Pi + theta_two



	call psi_interp_setup(psi,n)

	nth = theta_points2

	allocate (r_data_psi(nth+r_ord,3),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in r_data_psi"
		 pause
		 stop
	endif

	do i=1,nth

		theta_temp = (i-1.d0)/(nth-1.d0)*2.d0*pi

		k = 0

		smallr1 = x_size/10.d0 ! start with a small value

		smallr2 = dbsval(theta_temp, r_ord, r_data(1:theta_points1+r_ord,3),  &
			theta_points1, r_cscoef(1,1:theta_points1) )

		rg1 = 0.d0
		rg2 = smallr1 * .5d0

		r_data_psi(i,1) = theta_temp

		call secroot(psi_sol,rg1,rg2,r_data_psi(i,2),  &
					100,1.d-9,1.d-9, error)

		do while(error==1)
		! secroot failed because the root was not bracketed, try again
!!$			rg2 = smallr2 * 0.99**k
!!$			k = k+1
			rg2 = rg2 * 1.01d0

			if(rg2>0.995d0*smallr2) then
				r_data_psi(i,2) = smallr2
				print*, 'problem in initialize_r_free_boundary, theta = ', theta_temp
				! no need to say this should NOT happen
				error = 0
			else
				call secroot(psi_sol,rg1,rg2,r_data_psi(i,2),  &
							100,1.d-9,1.d-9, error)
			endif

		enddo

		continue

	enddo

	do i=1,theta_points2 + r_ord
	! just in case

		r_data(i,4) = 0.d0
		r_data(i,5) = 0.d0
		r_data(i,6) = 0.d0

	enddo

	r_cscoef(2,:) = 0.d0

	do i=1,theta_points2

		r_data(i,4) = r_data_psi(i,1)
		r_data(i,5) = r_data_psi(i,2)

	enddo

	call DBSNAK(theta_points2, r_data(1:theta_points2,4),  &
					r_ord,r_data(1:theta_points2+r_ord,6))

	call DBSINT(theta_points2, r_data(1:theta_points2,4),  &
		 r_data(1:theta_points2,5), r_ord,  &
		 r_data(1:theta_points2+r_ord,6),  &
		 r_cscoef(2,1:theta_points2))

	deallocate(r_data_psi)
	deallocate(psi_bscoef)
	deallocate(xknot_psi)
	deallocate(zknot_psi)

	continue

	return

end subroutine initialize_r_free_boundary

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psi_interp_setup(psi,n)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : x_coord, z_coord
	use exp_data, only :  s_ord, nx_FLOW, nz_FLOW,  &
							xknot_psi, zknot_psi, psi_bscoef

	use pseudo_IMSL, only : DBS2IN, DBS2GD, DBSNAK

	implicit none

	integer :: alloc_stat
	integer :: i, j
	integer :: n
	real(kind=dkind), dimension(1:n,1:n) :: psi

	nx_FLOW = n
	nz_FLOW = n

!-----------------------------------------
! interpolation set up

	allocate (psi_bscoef(1:nx_FLOW,1:nz_FLOW),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in psi_bscoef"
		 pause
		 stop
	endif

	allocate (xknot_psi(1:nx_FLOW+s_ord),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in xknot_psi"
		 pause
		 stop
	endif

	allocate (zknot_psi(1:nz_FLOW+s_ord),stat = alloc_stat)

	if(alloc_stat > 0) then
		 print *, "Allocation Error in zknot_psi"
		 pause
		 stop
	endif

	call DBSNAK(nx_FLOW,x_coord,s_ord,xknot_psi)
	call DBSNAK(nz_FLOW,z_coord,s_ord,zknot_psi)

	! (this 2 define the nodes)

! end of interpolation set up
!-----------------------------------------
!-----------------------------------------
! set psi

	call DBS2IN(nx_FLOW,x_coord,nz_FLOW,z_coord,psi,nx_FLOW,  &
				s_ord,s_ord,xknot_psi,zknot_psi,psi_bscoef(:,:) )

!-----------------------------------------

	continue

	return

end subroutine psi_interp_setup

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function dist2(x1,z1,x2,z2) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	real(kind=dkind) ::	x1,z1,x2,z2
	real(kind=dkind) :: answer

	answer = (x1-x2)**2 + (z1-z2)**2

	continue

end function dist2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_sort_grid(psi,nx,nz,inorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The code will get here to update which points are in the main plasma and which ones are in the open field line region.
! For now, only the bc_type==7 option is considered ( March 12 2021).

	implicit none

!	real(kind=dkind), parameter :: small = 1.d-2

	integer :: nx,nz
	real(kind=dkind), intent(in) :: psi(1:nx,1:nz)
	real(kind=dkind) :: inorm
	integer :: i,j,p
	real(kind=dkind) :: ex, ez, rminor, dx, dz, dummy
	real(kind=dkind) :: small = 1.d-2

!	sort_grid = -1
!	January 21 2022: commented the previous line.
!	STILL NEED TO CHECK THAT THIS DOES NOT BREAK TRI_TYPE=-1, -2

! February 4 2022: logical variables should not be needed anymore.
!!$	if((tri_type==-2).and.(initialize_r_for_tri_type_m2)) then
!!$		call initialize_r_free_boundary(psi(:,:),nx)
!!$		initialize_r_for_tri_type_m2 = .false.
!!$		tri_type_m2_ready = .true.
!!$	endif

!!$	if((tri_type==-2).and.(tri_type_m2_ready)) then
	if(tri_type==-2) then
		! Proceed like in old tri_type==13, interpolate psi to find the boundary
		call update_interface(psi(:,:),nx,inorm)
	endif

	do j=1,nz
	do i=1,nx

		if(tri_type==13) then
			cycle
		endif

		if(tri_type==-1) then
			call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
		elseif(tri_type==-2) then
			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
		endif

		if((ex*ex + ez*ez) > rminor**2) then

			sort_grid(i,j,0) = -1

		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then

			sort_grid(i,j,0) = 0

		else

			if(tri_type==-1) then

				 if((((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
						.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then

					sort_grid(i,j,0) = 1 ! External zone

				else

					sort_grid(i,j,0) = 2 ! Internal zone

				endif

			elseif(tri_type==-2) then

				if(p==2) then
					sort_grid(i,j,0) = 2
					! this sets inner zone (plasma) to "2"
				else
					sort_grid(i,j,0) = 1
				endif

			endif

		endif

	enddo
	enddo

	if(tri_type==11) then

		do j=1,nz
		do i=1,nx

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
			sort_grid(i,j,1) = p

		enddo
		enddo

	elseif(tri_type==13) then

		do j=1,nz
		do i=1,nx

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
			if((p==2).and.(sort_grid(i,j,0)==1)) sort_grid(i,j,0) = 2
			! this sets inner zone (plasma) to "2"

		enddo
		enddo

	endif

	open(33,file='grid.plt')

	write(33,*)'TITLE="grid"'
	write(33,*)'Variables =" R [m] ","z [m]", "boh"'
	write(33,*)'ZONE I=',nx,',J=',nz,',F=Point'

	do j = 1,nz
	do i = 1,nx

		write(33,*) x_coord(i), z_coord(j),  sort_grid(i,j,0)

	enddo
	enddo

	close(33)

	continue

end subroutine update_sort_grid



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine update_sort_grid_old(psi,nx,nz,inorm)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! The code will get here to update which points are in the main plasma and which ones are in the open field line region.
! For now, only the bc_type==7 option is considered ( March 12 2021).
! Copy for records: routine heavily edited on February 4 2022.

	implicit none

!	real(kind=dkind), parameter :: small = 1.d-2

	integer :: nx,nz
	real(kind=dkind), intent(in) :: psi(1:nx,1:nz)
	real(kind=dkind) :: inorm
	integer :: i,j,p
	real(kind=dkind) :: ex, ez, rminor, dx, dz, dummy
	real(kind=dkind) :: small = 1.d-2

!	sort_grid = -1
!	January 21 2022: commented the previous line.
!	STILL NEED TO CHECK THAT THIS DOES NOT BREAK TRI_TYPE=-1, -2

	if((tri_type==-2).and.(initialize_r_for_tri_type_m2)) then
		call initialize_r_free_boundary(psi(:,:),nx)
		initialize_r_for_tri_type_m2 = .false.
		tri_type_m2_ready = .true.
	endif

	if((tri_type==-2).and.(tri_type_m2_ready)) then
		! Proceed like in old tri_type==13, interpolate psi to find the boundary
		call update_interface(psi(:,:),nx,inorm)
	endif

	do j=1,nz
	do i=1,nx

		if(tri_type==13) then
			cycle
		endif

		if((tri_type==-1).or.((tri_type==-2).and.(.not.(tri_type_m2_ready)))) then
			call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)
		elseif((tri_type==-2).and.(tri_type_m2_ready)) then
			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
		endif

		if((ex*ex + ez*ez) > rminor**2) then

			sort_grid(i,j,0) = -1

		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then

			sort_grid(i,j,0) = 0

		else

			if((tri_type==-1).or.((tri_type==-2).and.(.not.(tri_type_m2_ready)))) then

				 if((((bc_type==7).or.(bc_type==8)).and.(psi(i,j)<0.d0))  &
						.or. ((bc_type==7).and.(sqrt((x_coord(i)-rmajor)**2+z_coord(j)**2)>0.4d0*z_size))) then

					sort_grid(i,j,0) = 1 ! External zone

				else

					sort_grid(i,j,0) = 2 ! Internal zone

				endif

			elseif((tri_type==-2).and.(tri_type_m2_ready)) then

				if(p==2) then
					sort_grid(i,j,0) = 2
					! this sets inner zone (plasma) to "2"
				else
					sort_grid(i,j,0) = 1
				endif

			endif

		endif

!!$		call radius(i,j,nx,nz,ex,ez,rminor,dx,dz)

!!$		if((ex*ex + ez*ez) > rminor**2) then

!!$			sort_grid(i,j,0) = -1

!!$		elseif( (ex*ex + ez*ez) >= (rminor - small*sqrt(dx_a(i)**2+dz_a(j)**2) )**2 ) then
!		elseif((ex*ex + ez*ez) == rminor**2) then

!!$			sort_grid(i,j,0) = 0

!!$		else

!!$			sort_grid(i,j,0) = 1

!!$		end if

	enddo
	enddo

	if(tri_type==11) then

		do j=1,nz
		do i=1,nx

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
			sort_grid(i,j,1) = p

		enddo
		enddo

	elseif(tri_type==13) then

		do j=1,nz
		do i=1,nx

			call radius_1_3(x_coord(i),z_coord(j),ex,ez,dummy,rminor,p,dummy,dummy,dummy)
			if((p==2).and.(sort_grid(i,j,0)==1)) sort_grid(i,j,0) = 2
			! this sets inner zone (plasma) to "2"

		enddo
		enddo

	endif

	open(33,file='grid.plt')

	write(33,*)'TITLE="grid"'
	write(33,*)'Variables =" R [m] ","z [m]", "boh"'
	write(33,*)'ZONE I=',nx,',J=',nz,',F=Point'

	do j = 1,nz
	do i = 1,nx

		write(33,*) x_coord(i), z_coord(j),  sort_grid(i,j,0)

	enddo
	enddo

	close(33)

	continue

end subroutine update_sort_grid_old



!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine get_minor_radius
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this calculates the minor radius for bootstrap calculation with NCLASS.
! the routine assumes |delta|<1; if needed, that will be patched later

	real(kind=dkind) :: rloc, xloc, zloc, thetaloc
	real(kind=dkind) :: Rmin, Rmax
	integer :: i
	integer :: n_angles = 200

	Rmin = 2.d0*rcenter
	Rmax = rcenter - 2.d2*x_size

	do i = 1,n_angles

		thetaloc = i*2.d0*pi/n_angles

		call radius_theta(thetaloc,rloc,xloc,zloc)

		Rmin = min(Rmin, xloc)
		Rmax = max(Rmax, xloc)

	enddo

	a_elps = (Rmax-Rmin)/2.d0

	continue

	end subroutine get_minor_radius

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius(i,j,nx,nz,ex,ez,r,dxx,dzz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the minor radius r of the plasma
	! as a function of the angle AND of the triangularity option

		use triangularity
		use constant

		integer, intent(in) :: i,j,nx,nz
		real (kind=dkind) :: dxx,dzz,bigR,Z, diff
		real (kind=dkind), intent(out) :: ex,ez,r
		integer k
		real (kind=dkind) :: r_in,r_out,angle

!		ex = ((i-1)*dx - 0.5d0*x_size)
!		ez = ((j-1)*dz - 0.5d0*z_size)

		if( (i==0).or.(j==0).or.(i>nx).or.(j>nz) ) then
		! somehow the code got here with an index exceeding the grid dimension
		! assign this to be an external point

			ex = 1.d0
			ez = 1.d0
			r = 1.d-2
			return

		endif

		ex = x_coord(i)-rmajor
		ez = z_coord(j)

		if (ex==0.d0) then
			angle = pi/2.d0 * dsign(1.d0,ez)
		else
			angle = datan2(ez,ex)
		endif

		r = 0.d0
		bigR = 0.d0
		Z = 0.d0

		! Edited tri_type=-1 and -2 on February 4 2022
		if(tri_type==-2) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle+2.d0*pi

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

			r = r_in


		elseif(tri_type==-1)then
			ex = dabs(((i-1.0d0)/(nx - 1.0d0) - 0.5d0)*x_size/a_elps)
			ez = dabs(((j-1.0d0)/(nz - 1.0d0) - 0.5d0)*z_size/b_elps)
			ex=dmax1(ex,ez)
			ez=0.d0
			r=1.d0

		elseif(tri_type==0) then

!			ex = ((i-1.0d0)/(nx - 1.0d0) - 0.5d0)*x_size/a_elps
!			ez = ((j-1.0d0)/(nz - 1.0d0) - 0.5d0)*z_size/b_elps
			ex = (x_coord(i)-rmajor)/a_elps
			ez = z_coord(j)/b_elps

			r = 1.d0

		elseif(tri_type==1) then

			if((dsin(angle)>=0.d0)) then
				do k=0,n_tri
					r = r + rcoeff_u(k)*dcos(k*angle)
				enddo
			else
          				do k=0,n_tri
					r = r + rcoeff_d(k)*dcos(k*angle)
				enddo
			endif

		elseif(tri_type==2) then

			if((dsin(angle)>=0.d0)) then
				r = smallr0*dsqrt( ellipt**2*(dsin(angle))**2  &
						+(dcos(angle+asin_d_up*dsin(angle)))**2 )
			else
          		r = smallr0*dsqrt( ellipt**2*(dsin(angle))**2  &
						+(dcos(angle+asin_d_down*dsin(angle)))**2 )
			endif

		elseif(tri_type==3) then

			if((dsin(angle)>=0.d0)) then
				r = smallr0/dsqrt( dsin(angle)**2/ellipt**2  &
						+dcos(angle)**2-delta_u_o4*(dcos(3*angle)-dcos(angle)) )
			else
          		r = smallr0/dsqrt( dsin(angle)**2/ellipt**2  &
						+dcos(angle)**2-delta_d_o4*(dcos(3*angle)-dcos(angle)) )
			endif

		elseif(tri_type==4) then

!			if (angle>pi/2) angle=angle-2.d0*pi
			if (angle>1.7451968605014) angle=angle-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,angle,r)

        elseif(tri_type==5) then                       

       		if((dsin(angle)>=0.d0)) then
				do k=0,n_tri
					r = r + rcos_u(k)*dcos(k*angle)+rsin_u(k)*dsin(k*angle)
				enddo
			else
                do k=0,n_tri
					r = r + rcos_d(k)*dcos(k*angle)+rsin_d(k)*dsin(k*angle)
				enddo
			endif

	    elseif(tri_type==6) then                       

			if (angle>pi/2.d0) angle=angle-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,angle,r)

		elseif(tri_type==8) then

			if(angle>=r_data(theta_points,1)) then

				angle = angle - 2.d0*pi

			elseif(angle<r_data(1,1)) then

				angle = angle + 2.d0*pi

			endif

!			if(angle==r_data(theta_points,1)) angle = angle-2.d0*pi

!			if(abs(abs(angle)-pi)<1.d-1) angle=abs(angle)-1.d-1

!!$?			if( (i<nx/2).and.(j==(nz+1)/2) )  angle = pi

!!$?			if( (i<nx/2).and.(j==(nz+1)/2) ) then
!!$?				continue
!!$?			endif

!			r = dcsval(angle, angle_points, r_data(:,3),  &
!					r_cscoef)
			r = dbsval(angle, r_ord, r_data(:,3),  &
						theta_points, r_cscoef(1,1:theta_points) )

!			if( (abs(ez)<2.1d0*dz).and.(ex<-0.1523) ) r=1.d-4	! shape_bookmark	! -.42d0

        elseif(tri_type==9) then

!			bigR = rmajor*dsqrt(1.d0+2.d0*a_elps/rmajor*dcos(angle))
!			Z = rmajor*a_elps*dsin(angle)/bigR
!			r = dsqrt((bigR-rmajor)**2+Z**2)

			bigR = rmajor - 0.5d0*x_size + (i-1)*dx
			Z = -0.5d0*z_size + (j-1)*dz

			diff = bigR**2 * Z**2 + 0.25d0 * ( bigR**2 - rmajor**2 )**2  &
				- a_elps**2 * rmajor**2	

			if(diff<0.d0) then

				ex = .5d0
				ez = .5d0
				r = 1.d0

			elseif(diff==0.d0) then

				ex = 1.d0
				ez = 0.d0
				r = 1.d0

			else

				ex = .5d0
				ez = .5d0
				r = .5d0

			endif

		elseif(tri_type==11) then

			! LDX-like geometry

			if(angle<0.d0) angle = angle+2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
							theta_points1, r_cscoef(1,1:theta_points1) )

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
							theta_points2, r_cscoef(2,1:theta_points2) )

			if( (ex**2+ez**2)<r_in**2 ) then

				ex = 1.d0
				ez = 1.d0
				r = 0.5d0

			else

				r = r_out

			endif

		elseif(tri_type==13) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle+2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
							theta_points1, r_cscoef(1,1:theta_points1) )

			r = r_out

		endif

		continue

	end subroutine radius

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_1_3(x,z,ex,ez,thet,r,zone,alph,rprim,rsec)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! the code used to get here only for tri_type=11
	! but now this is used to get rprim for bcs for tri_type=8, too (LG 4/14/2006)
	! now, here "4" is the zone inside the inner boundary ("hole")
	! and "1" the zone outside the outer boundary ("external vacuum")

		use triangularity
		use constant

		real (kind=dkind), intent(in) :: x,z
		real (kind=dkind), intent(out) :: r
		real (kind=dkind) :: r_in,r_out,ex,ez, thet, angle
		integer :: zone
		real(kind=dkind) :: alph, rprim, rsec

		ex = x - rmajor
		ez = z

		if (ex==0.d0) then
			thet = pi/2.d0 * dsign(1.d0,ez)
			angle = thet
		else
			thet = datan2(ez,ex)
			angle = thet
		endif

		r = 0.d0
		zone = 1
		alph = 1.d0

		if(tri_type==8) then

			! r of theta

			if(angle>=r_data(theta_points,1)) then

				angle = angle - 2.d0*pi

			elseif(angle<r_data(1,1)) then

				angle = angle + 2.d0*pi

			endif

			r = dbsval(angle, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )

			rprim = dbsder(1,angle, r_ord, r_data(:,3),  &
					theta_points, r_cscoef(1,1:theta_points) )

			rsec = dbsder(2,angle, r_ord, r_data(:,3),  &
					theta_points, r_cscoef(1,1:theta_points) )

		elseif(tri_type==11) then

			! LDX-like geometry

			if(thet<0.d0) thet = thet + 2.d0*pi
			if(angle<0.d0) angle = angle + 2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

			if( (ex**2+ez**2)<r_in**2 ) then

				! "move" the point to the outside

				if( (ex**2+ez**2)==0.d0) then

					ex = r_in
					ez = r_in
					! just to have the point "outside"

				else

					alph = sqrt(2*r_in**2/(ex**2+ez**2) - 1.d0)

				endif

!				ex = alph * ex
!				ez = alph * ez
				zone = 4
				r = r_in
				rprim = dbsder(1,angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
						theta_points2, r_cscoef(2,1:theta_points2) )
				rsec = dbsder(2,angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
						theta_points2, r_cscoef(2,1:theta_points2) )

			elseif( (ex**2+ez**2)>r_out**2 ) then

				zone = 1
				r = r_out
				rprim = dbsder(1,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
						theta_points1, r_cscoef(1,1:theta_points1) )
				rsec = dbsder(2,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
						theta_points1, r_cscoef(1,1:theta_points1) )

			else

				if( abs( (ex**2+ez**2)-r_in**2 ) < abs( (ex**2+ez**2)-r_out**2 ) ) then

					! the inner point is closer to the inner than the outer boundary
					zone = 3
					r = r_in
					rprim = dbsder(1,angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
							theta_points2, r_cscoef(2,1:theta_points2) )
					rsec = dbsder(2,angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
							theta_points2, r_cscoef(2,1:theta_points2) )

				else

					zone = 2
					r = r_out
					rprim = dbsder(1,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
							theta_points1, r_cscoef(1,1:theta_points1) )
					rsec = dbsder(2,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
							theta_points1, r_cscoef(1,1:theta_points1) )

				endif

			endif

		elseif((tri_type==13).or.(tri_type==-2)) then

			! plasma - vacuum case

			if(angle<0.d0) angle = angle + 2.d0*pi

			r_out = dbsval(angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
				theta_points1, r_cscoef(1,1:theta_points1) )

			r_in = dbsval(angle, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

			if( (ex**2+ez**2)>r_in**2 ) then

				zone = 1

			else

				zone = 2

			endif

			! the actual radius is always the external radius,
			! the internal one is a separation!
			r = r_out
			rprim = dbsder(1,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
					theta_points1, r_cscoef(1,1:theta_points1) )
			rsec = dbsder(2,angle, r_ord, r_data(1:theta_points1+r_ord,3),  &
					theta_points1, r_cscoef(1,1:theta_points1) )

		else

			print*, 'error in radius_1_3: tri_type=', tri_type
			pause
			stop

		endif

		if(alph/=1.d0) then
			continue
		endif

		thet = angle

		continue

	end subroutine radius_1_3

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_theta(angle,r,x,z)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the minor radius r of the plasma
	! as a function of the angle AND of the triangularity option
	! 1/4/2007: added x and z for eqdsk output

		use triangularity
		use constant

		real (kind=dkind), intent(out) :: r,x,z
		integer k
		real(kind=dkind) :: angle

		r = 0.d0
		theta = angle

		if(tri_type==0) then

			r = a_elps*b_elps/sqrt( (a_elps*sin(angle))**2 +  &
						(b_elps*cos(angle))**2 )

		elseif(tri_type==1) then

			if((dsin(theta)>=0.d0)) then
				do k=0,n_tri
					r = r + rcoeff_u(k)*dcos(k*theta)
				enddo
			else
          		do k=0,n_tri
					r = r + rcoeff_d(k)*dcos(k*theta)
				enddo
			endif

		elseif(tri_type==2) then

			if((dsin(theta)>=0.d0)) then
				r = smallr0*dsqrt( ellipt**2*(dsin(theta))**2  &
						+(dcos(theta+asin_d_up*dsin(theta)))**2 )
			else
          		r = smallr0*dsqrt( ellipt**2*(dsin(theta))**2  &
						+(dcos(theta+asin_d_down*dsin(theta)))**2 )
			endif

		elseif(tri_type==3) then

			if((dsin(theta)>=0.d0)) then
				r = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
						+dcos(theta)**2-delta_u_o4*(dcos(3*theta)-dcos(theta)) )
			else
          		r = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
						+dcos(theta)**2-delta_d_o4*(dcos(3*theta)-dcos(theta)) )
			endif

		elseif(tri_type==4) then

!			if (theta>pi/2) theta=theta-2.d0*pi
			if (theta>1.7451968605014) theta=theta-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,theta,r)

		elseif(tri_type==-1) then
			! square boundary, this will not be needed

        elseif(tri_type==5) then                       

       		if((dsin(theta)>=0.d0)) then
				do k=0,n_tri
					r = r + rcos_u(k)*dcos(k*theta)+rsin_u(k)*dsin(k*theta)
				enddo
			else
                do k=0,n_tri
					r = r + rcos_d(k)*dcos(k*theta)+rsin_d(k)*dsin(k*theta)
				enddo
			endif

	    elseif(tri_type==6) then                       

			if (theta>pi/2.d0) theta=theta-2.d0*pi

			call splint(theta_dat,rminor_dat,d2rminor_dat,theta_points,theta,r)

		elseif(tri_type==8) then

			if(theta>=r_data(theta_points,1)) then

				theta = theta - 2.d0*pi

			elseif(theta<r_data(1,1)) then

				theta = theta + 2.d0*pi

			endif

!			if(abs(abs(theta)-pi)<1.d-1) theta=abs(theta)-1.d-1

!			if(theta<0.d0) theta = theta + 2.d0*pi

!			r = dcsval(theta, theta_points, r_data(:,3),  &
!					r_cscoef)
			r = dbsval(theta, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )

!			if( (abs(ez)<2.1d0*dz).and.(ex<-0.1523) ) r=1.d-4	! shape_bookmark	! -.42d0

        elseif(tri_type==9) then

			! Solovev solution, this should not be needed

		elseif(tri_type==11) then

			! LDX-like geometry
			! this will not be needed

		elseif(tri_type==13) then

			! plasma - vacuum case
			! this time we want the plasma radius!

			r = dbsval(theta, r_ord, r_data(1:theta_points2+r_ord,6),  &
				theta_points2, r_cscoef(2,1:theta_points2) )

		endif

		x = rmajor + r*cos(angle)
		z = r*sin(angle)

		continue

	end subroutine radius_theta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine radius_prim_theta(angle,rp)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine computes the first derivative 
	! of the minor radius r of the plasma as a function of the angle
	! only for numerical input

		use triangularity
		use constant

		real (kind=dkind), intent(out) :: rp
		real(kind=dkind) :: angle

		rp = 0.d0
		theta = angle

		if(tri_type==8) then
			continue
		else
			return
		endif

		if(theta>=r_data(theta_points,1)) then

			theta = theta - 2.d0*pi

		elseif(theta<r_data(1,1)) then

			theta = theta + 2.d0*pi

		endif

			rp = dbsder(1,theta, r_ord, r_data(:,3),  &
				theta_points, r_cscoef(1,1:theta_points) )

		continue

	end subroutine radius_prim_theta



!-----------------------------------------stuff that used to be in root-----------------------------------------



! Given a function fx defined on the interval from x1 -> x2, subdivide
! the interval into n equally spaced segments, and search for zero
! crossings of the function.  nb is input as the maximum number of
! roots sought, and is reset to the number of bracketing pairs
! xb1(1:nb), xb2(1:nb) that are found.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine zbrak(fx,x1,x2,m,xb1,xb2,nb)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dkind
  implicit none

  integer, intent(in) :: m
  integer, intent(inout) :: nb
  real (kind=dkind) :: x1,x2,fx
  external fx
  real (kind=dkind), dimension(1:nb), intent(out) :: xb1,xb2
  integer :: i,nbb
  real (kind=dkind) :: dx,fc,fp,xp,xc

  if(x2 <= x1) then
     print *, "[zbrak]: x1 = ",x1," and x2 = ",x2
     nb = 0
     stop
  end if

  nbb=0
  xp=x1
  dx=(x2-x1)/m ! Determine the spacing appropriate to the mesh
  fp=fx(xp)

  lp: do i=1,m ! Loop over all intervals
     xc=xp+dx
     fc=fx(xc)
     ! if a sign change occurs then record values for the bounds
     if(fc*fp.lt.0.) then 
        nbb=nbb+1
        xb1(nbb)=xp
        xb2(nbb)=xc
        if(nbb.eq.nb) exit lp
     endif
     ! Shift the x position and function value
     xp=xc
     fp=fc
  end do lp

  nb=nbb
  return
end subroutine zbrak

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function rtbis(func,x1,x2,xacc) result (soln)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dkind
  implicit none

  integer, parameter :: jmax = 400
  real (kind=dkind) :: x1,x2,xacc,func
  real (kind=dkind) :: soln
  external func
  integer :: j
  real (kind=dkind) :: dx,f,fmid,xmid

  fmid = func(x2)
  f = func(x1)
  if(f*fmid .ge. 0.0d0) then
     print *, '[rtbis]: root must be bracketed'
     soln = x1
!	 call b_plot(x1,x2,1000)
     return
  end if

  if(f .lt. 0.0d0)then
     soln = x1
     dx = x2 - x1
  else
     soln = x2
     dx = x1 - x2
  endif

  do j = 1, jmax
     dx = 0.5d0*dx
     xmid = soln + dx
     fmid = func(xmid)
     if(fmid .le. 0.0d0) soln = xmid
     if(dabs(dx) .lt. xacc .or. fmid .eq. 0.0d0) return
  end do

  print *, '[rtbis]: too many bisections'

end function rtbis

!------------------------------------------------------------------

! Using a combination of Newton-Raphson and Bisection, find the root
! of a function between x1 and x2.  The root, returned as the value 
! soln, will be refined until its accuracy is known within +/-xacc.
! "funcd" is a user supplied subroutine which returns both the 
! value of the function and its first derivative.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function rtsafe(funcd,x1,x2,xacc,maxit) result (soln)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant, only : dkind
  implicit none

  integer :: maxit
  real (kind=dkind) :: x1,x2,xacc
  real (kind=dkind) :: soln
  external funcd
  integer :: j
  real (kind=dkind) :: df,dxx,dxold,f,fh,fl,temp,xh,xl

  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  if((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) then
     print *, '[rtsafe]: root is not bracketed'
!	 call b_plot(x1,x2,1000)
     soln = 0.0 ! set it to something
     ! In fact since I will be looking for a density which 
     ! I know to be positive this will work fine for an error indicator
	 return
  end if

  if(fl.eq.0.)then
     soln=x1
     return
  else if(fh.eq.0.)then
     soln=x2
     return
  else if(fl.lt.0.)then ! Orient the search so that f(xl) < 0
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  endif

  soln=.5*(x1+x2)  ! Initialize the guess for the root, 
  dxold=abs(x2-x1) ! the "stepsize before last", 
  dxx=dxold         ! and the last step.

  call funcd(soln,f,df) ! Initialize f and df for first guess

  lp: do j=1,maxit ! Loop over the allowed iterations
     ! If the Newton step is out of range or not decreasing 
     ! fast enough, use bisection 
     if(((soln-xh)*df-f)*((soln-xl)*df-f) .ge. 0.d0 &
          .or. abs(2.d0*f) .gt. abs(dxold*df)) then
        dxold=dxx
        dxx=0.5d0*(xh-xl)
        soln=xl+dxx
        if(xl .eq. soln) return
     else
        dxold=dxx
        dxx=f/df
        temp=soln
        soln=soln-dxx
        if(temp .eq. soln) return
     endif
     ! Return if the convergence criterion is met
     if(abs(dxx) .lt. xacc) return

     ! Otherwise, evaluate f and df at our new guess
     call funcd(soln,f,df)

     ! Shrink the bracket enclosing the root
     if(f.lt.0.) then
        xl=soln
     else
        xh=soln
     endif
  end do lp

  print *, '[rtsafe]: maximum number of iterations exceeded'
!  call b_plot(x1,x2,1000)

  return

end function rtsafe

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine secroot(fx,x1,x2,x,itemax,xtoll,ftoll,error)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	implicit none

	real(kind=dkind) :: x1,x2,x
	integer :: i
	integer :: itemax
	real(kind=dkind) :: xtoll, ftoll
	real(kind=dkind) :: f1,f2,fx
	real(kind=dkind) :: m, x3, f3
	integer :: error

!	external fx

	error = 0

	f1 = fx(x1)
	f2 = fx(x2)

	if(f1*f2>0.d0) then

!		pause 'error in secroot, root is not bracketed'
		error = 1
		return

	endif

	if(f2<f1) then
	! (f2>0 by hypothesis)

		x3 = x2
		x2 = x1
		x1 = x3

		f3 = f2
		f2 = f1
		f1 = f3

	endif




	sloop : do i=1,itemax

		m = (f2-f1)/(x2-x1)
		x3 = x1 - f1/m
		f3 = fx(x3)

		if(f3>0.d0) then

			x2 = x3
			f2 = f3

		else

			x1 = x3
			f1 = f3

		endif

		if( (abs(f2-f1)<ftoll).or.(abs(x2-x1)<xtoll).or.  &
			 (abs(f2)<ftoll).or.(abs(f1)<ftoll) ) then
		! converged!

			exit sloop

		endif

	enddo sloop

	if(abs(f2)<abs(f1)) then

		x = x2

	else

		x = x1

	endif

end subroutine secroot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      function derfc (x) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!!$***BEGIN PROLOGUE  DERFC
!!$***PURPOSE  Compute the complementary error function.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C8A, L5A1E
!!$***TYPE      DOUBLE PRECISION (ERFC-S, DERFC-D)
!!$***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
!!$             SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$ DERFC(X) calculates the double precision complementary error function
!!$ for double precision argument X.
!!$
!!$ Series for ERF        on the interval  0.          to  1.00000E+00
!!$                                        with weighted Error   1.28E-32
!!$                                         log weighted Error  31.89
!!$                               significant figures required  31.05
!!$                                    decimal places required  32.55
!!$
!!$ Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00
!!$                                        with weighted Error   2.67E-32
!!$                                         log weighted Error  31.57
!!$                               significant figures required  30.31
!!$                                    decimal places required  32.42
!!$
!!$ Series for ERFC       on the interval  0.          to  2.50000E-01
!!$                                        with weighted error   1.53E-31
!!$                                         log weighted error  30.82
!!$                               significant figures required  29.47
!!$                                    decimal places required  31.70
!!$
!!$***REFERENCES  (NONE)
!!$***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770701  DATE WRITTEN
!!$   890531  Changed all specific intrinsics to generic.  (WRB)
!!$   890531  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$   920618  Removed space from variable names.  (RWC, WRB)
!!$***END PROLOGUE  DERFC

	use constant, only : dkind

	implicit none

      real(kind=dkind) :: X, ERFCS(21), ERFCCS(59), ERC2CS(49), SQEPS,  &
       SQRTPI, XMAX, TXMAX, XSML, Y, answer	! D1MACH, DCSEVL,
	real (kind=dkind) :: eta
	integer ::  nterf, nterfc, nterc2
      LOGICAL FIRST
      SAVE ERFCS, ERC2CS, ERFCCS, SQRTPI, NTERF,  &
      NTERFC, NTERC2, XSML, XMAX, SQEPS, FIRST

!	  real(kind=dkind) :: initds

      DATA ERFCS(  1) / -.49046121234691808039984544033376D-1     /
      DATA ERFCS(  2) / -.14226120510371364237824741899631D+0     /
      DATA ERFCS(  3) / +.10035582187599795575754676712933D-1     /
      DATA ERFCS(  4) / -.5768764699767484765082702550917D-3     /
      DATA ERFCS(  5) / +.27419931252196061034422160791471D-4     /
      DATA ERFCS(  6) / -.11043175507344507604135381295905D-5     /
      DATA ERFCS(  7) / +.38488755420345036949961311498174D-7     /
      DATA ERFCS(  8) / -.11808582533875466969631751801581D-8     /
      DATA ERFCS(  9) / +.32334215826050909646402930953354D-10    /
      DATA ERFCS( 10) / -.79910159470045487581607374708595D-12    /
      DATA ERFCS( 11) / +.17990725113961455611967245486634D-13    /
      DATA ERFCS( 12) / -.37186354878186926382316828209493D-15    /
      DATA ERFCS( 13) / +.71035990037142529711689908394666D-17    /
      DATA ERFCS( 14) / -.12612455119155225832495424853333D-18    /
      DATA ERFCS( 15) / +.20916406941769294369170500266666D-20    /
      DATA ERFCS( 16) / -.32539731029314072982364160000000D-22    /
      DATA ERFCS( 17) / +.47668672097976748332373333333333D-24    /
      DATA ERFCS( 18) / -.65980120782851343155199999999999D-26    /
      DATA ERFCS( 19) / +.86550114699637626197333333333333D-28    /
      DATA ERFCS( 20) / -.10788925177498064213333333333333D-29    /
      DATA ERFCS( 21) / +.12811883993017002666666666666666D-31    /
      DATA ERC2CS(  1) / -.6960134660230950112739150826197D-1      /
      DATA ERC2CS(  2) / -.4110133936262089348982212084666D-1      /
      DATA ERC2CS(  3) / +.3914495866689626881561143705244D-2      /
      DATA ERC2CS(  4) / -.4906395650548979161280935450774D-3      /
      DATA ERC2CS(  5) / +.7157479001377036380760894141825D-4      /
      DATA ERC2CS(  6) / -.1153071634131232833808232847912D-4      /
      DATA ERC2CS(  7) / +.1994670590201997635052314867709D-5      /
      DATA ERC2CS(  8) / -.3642666471599222873936118430711D-6      /
      DATA ERC2CS(  9) / +.6944372610005012589931277214633D-7      /
      DATA ERC2CS( 10) / -.137122090210436601953460514121D-7      /
      DATA ERC2CS( 11) / +.2788389661007137131963860348087D-8      /
      DATA ERC2CS( 12) / -.5814164724331161551864791050316D-9      /
      DATA ERC2CS( 13) / +.1238920491752753181180168817950D-9      /
      DATA ERC2CS( 14) / -.2690639145306743432390424937889D-10     /
      DATA ERC2CS( 15) / +.5942614350847910982444709683840D-11     /
      DATA ERC2CS( 16) / -.1332386735758119579287754420570D-11     /
      DATA ERC2CS( 17) / +.3028046806177132017173697243304D-12     /
      DATA ERC2CS( 18) / -.6966648814941032588795867588954D-13     /
      DATA ERC2CS( 19) / +.1620854541053922969812893227628D-13     /
      DATA ERC2CS( 20) / -.3809934465250491999876913057729D-14     /
      DATA ERC2CS( 21) / +.9040487815978831149368971012975D-15     /
      DATA ERC2CS( 22) / -.2164006195089607347809812047003D-15     /
      DATA ERC2CS( 23) / +.5222102233995854984607980244172D-16     /
      DATA ERC2CS( 24) / -.1269729602364555336372415527780D-16     /
      DATA ERC2CS( 25) / +.3109145504276197583836227412951D-17     /
      DATA ERC2CS( 26) / -.7663762920320385524009566714811D-18     /
      DATA ERC2CS( 27) / +.1900819251362745202536929733290D-18     /
      DATA ERC2CS( 28) / -.4742207279069039545225655999965D-19     /
      DATA ERC2CS( 29) / +.1189649200076528382880683078451D-19     /
      DATA ERC2CS( 30) / -.3000035590325780256845271313066D-20     /
      DATA ERC2CS( 31) / +.7602993453043246173019385277098D-21     /
      DATA ERC2CS( 32) / -.1935909447606872881569811049130D-21     /
      DATA ERC2CS( 33) / +.4951399124773337881000042386773D-22     /
      DATA ERC2CS( 34) / -.1271807481336371879608621989888D-22     /
      DATA ERC2CS( 35) / +.3280049600469513043315841652053D-23     /
      DATA ERC2CS( 36) / -.8492320176822896568924792422399D-24     /
      DATA ERC2CS( 37) / +.2206917892807560223519879987199D-24     /
      DATA ERC2CS( 38) / -.5755617245696528498312819507199D-25     /
      DATA ERC2CS( 39) / +.1506191533639234250354144051199D-25     /
      DATA ERC2CS( 40) / -.395450295901879695310428569599D-26     /
      DATA ERC2CS( 41) / +.1041529704151500979984645051733D-26     /
      DATA ERC2CS( 42) / -.275148779527876507945017890133D-27     /
      DATA ERC2CS( 43) / +.7290058205497557408997703680000D-28     /
      DATA ERC2CS( 44) / -.193693964591594780407750109866D-28     /
      DATA ERC2CS( 45) / +.5160357112051487298370054826666D-29     /
      DATA ERC2CS( 46) / -.137841932219309409938964480000D-29     /
      DATA ERC2CS( 47) / +.3691326793107069042251093333333D-30     /
      DATA ERC2CS( 48) / -.990938959062436542065322666666D-31     /
      DATA ERC2CS( 49) / +.2666491705195388413323946666666D-31     /
      DATA ERFCCS(  1) / +.715179310202924774503697709496D-1        /
      DATA ERFCCS(  2) / -.265324343376067157558893386681D-1        /
      DATA ERFCCS(  3) / +.171115397792085588332699194606D-2        /
      DATA ERFCCS(  4) / -.163751663458517884163746404749D-3        /
      DATA ERFCCS(  5) / +.198712935005520364995974806758D-4        /
      DATA ERFCCS(  6) / -.284371241276655508750175183152D-5        /
      DATA ERFCCS(  7) / +.460616130896313036969379968464D-6        /
      DATA ERFCCS(  8) / -.822775302587920842057766536366D-7        /
      DATA ERFCCS(  9) / +.159214187277090112989358340826D-7        /
      DATA ERFCCS( 10) / -.329507136225284321486631665072D-8        /
      DATA ERFCCS( 11) / +.722343976040055546581261153890D-9        /
      DATA ERFCCS( 12) / -.166485581339872959344695966886D-9        /
      DATA ERFCCS( 13) / +.401039258823766482077671768814D-10       /
      DATA ERFCCS( 14) / -.100481621442573113272170176283D-10       /
      DATA ERFCCS( 15) / +.260827591330033380859341009439D-11       /
      DATA ERFCCS( 16) / -.699111056040402486557697812476D-12       /
      DATA ERFCCS( 17) / +.192949233326170708624205749803D-12       /
      DATA ERFCCS( 18) / -.547013118875433106490125085271D-13       /
      DATA ERFCCS( 19) / +.158966330976269744839084032762D-13       /
      DATA ERFCCS( 20) / -.472689398019755483920369584290D-14       /
      DATA ERFCCS( 21) / +.143587337678498478672873997840D-14       /
      DATA ERFCCS( 22) / -.444951056181735839417250062829D-15       /
      DATA ERFCCS( 23) / +.140481088476823343737305537466D-15       /
      DATA ERFCCS( 24) / -.451381838776421089625963281623D-16       /
      DATA ERFCCS( 25) / +.147452154104513307787018713262D-16       /
      DATA ERFCCS( 26) / -.489262140694577615436841552532D-17       /
      DATA ERFCCS( 27) / +.164761214141064673895301522827D-17       /
      DATA ERFCCS( 28) / -.562681717632940809299928521323D-18       /
      DATA ERFCCS( 29) / +.194744338223207851429197867821D-18       /
      DATA ERFCCS( 30) / -.682630564294842072956664144723D-19       /
      DATA ERFCCS( 31) / +.242198888729864924018301125438D-19       /
      DATA ERFCCS( 32) / -.869341413350307042563800861857D-20       /
      DATA ERFCCS( 33) / +.315518034622808557122363401262D-20       /
      DATA ERFCCS( 34) / -.115737232404960874261239486742D-20       /
      DATA ERFCCS( 35) / +.428894716160565394623737097442D-21       /
      DATA ERFCCS( 36) / -.160503074205761685005737770964D-21       /
      DATA ERFCCS( 37) / +.606329875745380264495069923027D-22       /
      DATA ERFCCS( 38) / -.231140425169795849098840801367D-22       /
      DATA ERFCCS( 39) / +.888877854066188552554702955697D-23       /
      DATA ERFCCS( 40) / -.344726057665137652230718495566D-23       /
      DATA ERFCCS( 41) / +.134786546020696506827582774181D-23       /
      DATA ERFCCS( 42) / -.531179407112502173645873201807D-24       /
      DATA ERFCCS( 43) / +.210934105861978316828954734537D-24       /
      DATA ERFCCS( 44) / -.843836558792378911598133256738D-25       /
      DATA ERFCCS( 45) / +.339998252494520890627359576337D-25       /
      DATA ERFCCS( 46) / -.137945238807324209002238377110D-25       /
      DATA ERFCCS( 47) / +.563449031183325261513392634811D-26       /
      DATA ERFCCS( 48) / -.231649043447706544823427752700D-26       /
      DATA ERFCCS( 49) / +.958446284460181015263158381226D-27       /
      DATA ERFCCS( 50) / -.399072288033010972624224850193D-27       /
      DATA ERFCCS( 51) / +.167212922594447736017228709669D-27       /
      DATA ERFCCS( 52) / -.704599152276601385638803782587D-28       /
      DATA ERFCCS( 53) / +.297976840286420635412357989444D-28       /
      DATA ERFCCS( 54) / -.126252246646061929722422632994D-28       /
      DATA ERFCCS( 55) / +.539543870454248793985299653154D-29       /
      DATA ERFCCS( 56) / -.238099288253145918675346190062D-29       /
      DATA ERFCCS( 57) / +.109905283010276157359726683750D-29       /
      DATA ERFCCS( 58) / -.486771374164496572732518677435D-30       /
      DATA ERFCCS( 59) / +.152587726411035756763200828211D-30       /
      DATA SQRTPI / 1.77245385090551602729816748334115D0 /
      DATA FIRST /.TRUE./
!!$***FIRST EXECUTABLE STATEMENT  DERFC
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NTERF = INITDS (ERFCS, 21, ETA)
         NTERFC = INITDS (ERFCCS, 59, ETA)
         NTERC2 = INITDS (ERC2CS, 49, ETA)
!!$
         XSML = -SQRT(-LOG(SQRTPI*D1MACH(3)))
         TXMAX = SQRT(-LOG(SQRTPI*D1MACH(1)))
         XMAX = TXMAX - 0.5D0*LOG(TXMAX)/TXMAX - 0.01D0
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
!!$
      IF (X.GT.XSML) GO TO 20
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
!!$
      answer = 2.0D0
      RETURN
!!$
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0D0) GO TO 30
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
!!$
      IF (Y.LT.SQEPS) answer = 1.0D0 - 2.0D0*X/SQRTPI
      IF (Y.GE.SQEPS) answer = 1.0D0 - X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,  &
       ERFCS, NTERF))
      RETURN
!!$
!!$ ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
!!$
 30   Y = Y*Y
      IF (Y.LE.4.D0) answer = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (  &
       (8.D0/Y-5.D0)/3.D0, ERC2CS, NTERC2) )
      IF (Y.GT.4.D0) answer = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (  &
       8.D0/Y-1.D0, ERFCCS, NTERFC) )
      IF (X.LT.0.D0) answer = 2.0D0 - answer
      RETURN
!!$
 40  print*, 'SLATEC', 'DERFC', 'X SO BIG ERFC UNDERFLOWS', 1, 1
      answer = 0.D0
      return
!!$
      end function derfc



 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     function initds (OS, NOS, ETA) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$***BEGIN PROLOGUE  INITDS
!!$***PURPOSE  Determine the number of terms needed in an orthogonal
!!$            polynomial series so that it meets a specified accuracy.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C3A2
!!$***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
!!$***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!!$             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$  Initialize the orthogonal series, represented by the array OS, so
!!$  that INITDS is the number of terms needed to insure the error is no
!!$  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!!$  machine precision.
!!$
!!$             Input Arguments --
!!$   OS     double precision array of NOS coefficients in an orthogonal
!!$          series.
!!$   NOS    number of coefficients in OS.
!!$   ETA    single precision scalar containing requested accuracy of
!!$          series.
!!$
!!$***REFERENCES  (NONE)
!!$***ROUTINES CALLED  XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770601  DATE WRITTEN
!!$   890531  Changed all specific intrinsics to generic.  (WRB)
!!$   890831  Modified array declarations.  (WRB)
!!$   891115  Modified error message.  (WRB)
!!$   891115  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$***END PROLOGUE  INITDS

	use constant, only : dkind

	implicit none

    real(kind=dkind) OS(*), eta, err
	real(kind=dkind) :: answer
	integer :: nos
	integer :: i, ii

!!$***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1)  print*, 'SLATEC', 'INITDS',  &
        'Number of coefficients is less than 1', 2, 1
!!$
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
!!$
   20 IF (I .EQ. NOS) print*, 'SLATEC', 'INITDS',  &
        'Chebyshev series too short for specified accuracy', 1, 1
      answer = I
!!$
      RETURN
      END function initds




 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      function DCSEVL (X, CS, N) result(answer)
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$***BEGIN PROLOGUE  DCSEVL
!!$***PURPOSE  Evaluate a Chebyshev series.
!!$***LIBRARY   SLATEC (FNLIB)
!!$***CATEGORY  C3A2
!!$***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
!!$***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!!$***AUTHOR  Fullerton, W., (LANL)
!!$***DESCRIPTION
!!$
!!$  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!!$  a method presented in the paper by Broucke referenced below.
!!$
!!$       Input Arguments --
!!$  X    value at which the series is to be evaluated.
!!$  CS   array of N terms of a Chebyshev series.  In evaluating
!!$       CS, only half the first coefficient is summed.
!!$  N    number of terms in array CS.
!!$
!!$***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!!$                 Chebyshev series, Algorithm 446, Communications of
!!$                 the A.C.M. 16, (1973) pp. 254-256.
!!$               L. Fox and I. B. Parker, Chebyshev Polynomials in
!!$                 Numerical Analysis, Oxford University Press, 1968,
!!$                 page 56.
!!$***ROUTINES CALLED  D1MACH, XERMSG
!!$***REVISION HISTORY  (YYMMDD)
!!$   770401  DATE WRITTEN
!!$   890831  Modified array declarations.  (WRB)
!!$   890831  REVISION DATE from Version 3.2
!!$   891214  Prologue converted to Version 4.0 format.  (BAB)
!!$   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!!$   900329  Prologued revised extensively and code rewritten to allow
!!$           X to be slightly outside interval (-1,+1).  (WRB)
!!$   920501  Reformatted the REFERENCES section.  (WRB)
!!$***END PROLOGUE  DCSEVL

	use constant, only : dkind

	implicit none

    real(kind=dkind) :: B0, B1, B2, CS(*), ONEPL, TWOX, X	!, D1MACH
	real(kind=dkind) :: answer
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./

	  integer :: n, i, ni

!!$***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) print*, 'SLATEC', 'DCSEVL',  'NUMBER OF TERMS .LE. 0', 2, 2
      IF (N .GT. 1000) print*, 'SLATEC', 'DCSEVL', 'NUMBER OF TERMS .GT. 1000', 3, 2
      IF (ABS(X) .GT. ONEPL) print*, 'SLATEC', 'DCSEVL', 'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1
!!$
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!!$
      answer = 0.5D0*(B0-B2)
!!$
      RETURN

      end function dcsevl



  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     function D1MACH(I)
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

      INTEGER I
	  real(kind=dkind) :: D1MACH, DMACH(5)

!!$
!!$  DOUBLE-PRECISION MACHINE CONSTANTS
!!$  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!!$  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!!$  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!!$  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!!$  D1MACH( 5) = LOG10(B)
!!$
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!!$  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
!!$  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
!!$  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
!!$  MANY MACHINES YET.
!!$  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!!$  ON THE NEXT LINE
      DATA SC/0/
!!$  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!!$  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!!$          mail netlib@research.bell-labs.com
!!$          send old1mach from blas
!!$  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!!$
!!$     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!!$      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!!$      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!!$      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!!$      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!!$      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
!!$
!!$     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!!$     32-BIT INTEGERS.
!!$      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!!$      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!!$      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!!$      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!!$      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
!!$
!!$     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!!$      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!!$      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!!$      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!!$      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!!$      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
!!$
!!$     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532.AND. SMALL(2) .EQ. -448790528) THEN
!!$           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532 .AND. SMALL(1) .EQ. -448790528) THEN
!!$           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935.AND. SMALL(2) .EQ. 10752) THEN
!!$               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943.AND. SMALL(2) .EQ. 704643072) THEN
!!$               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684.AND. SMALL(2) .EQ. -448790528) THEN
!!$           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074.AND. SMALL(2) .EQ. 58688) THEN
!!$           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
!!$                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
!!$    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/' appropriate for your machine.')

!!$* /* Standard C source for D1MACH -- remove the * in column 1 */
!!$*#include <stdio.h>
!!$*#include <float.h>
!!$*#include <math.h>
!!$*double d1mach_(long *i)
!!$*{
!!$*	switch(*i){
!!$*	  case 1: return DBL_MIN;
!!$*	  case 2: return DBL_MAX;
!!$*	  case 3: return DBL_EPSILON/FLT_RADIX;
!!$*	  case 4: return DBL_EPSILON;
!!$*	  case 5: return log10((double)FLT_RADIX);
!!$*	  }
!!$*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
!!$*	exit(1); return 0; /* some compilers demand return values */
!!$*}
      end function D1MACH

      SUBROUTINE I1MCRY(A, A1, B, C, D)
!!$**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END SUBROUTINE I1MCRY



! Test out the solver for the bernoulli function by finding num 
! different solutions
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine b_plot(min,max,num)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	real (kind=dkind), intent(in) :: min,max
	integer, intent(in) :: num
	real (kind=dkind) :: rho,drho,b,db
	integer :: i

	drho = (max - min)/(num - 1)

	open (unit=13,file='b_plot.dat',status='unknown',action='write')

	open (unit=11,file='db_plot.dat',status='unknown',action='write')

	rho = min
	do i = 1,num
	call newt_rap(rho,b,db)

	write (13,*) rho,b
	write (11,*) rho,db

	rho = rho + drho
	end do

	close(13)
	close(11)

end subroutine b_plot



! Test out the solver for the bernoulli function by finding num 
! different solutions
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine b_plot_gauss(min,max,num)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind

	real (kind=dkind), intent(in) :: min,max
	integer, intent(in) :: num
	real (kind=dkind) :: rho,drho,b,db,b_g,db_g
	integer :: i

	drho = (max - min)/(num - 1)

	open (unit=13,file='b_plot_gauss.dat',status='unknown',action='write')

	open (unit=11,file='db_plot_gauss.dat',status='unknown',action='write')

	rho = min
	do i = 1,num
	call newt_rap(rho,b,db)
	call newt_rap_gauss(rho,b_g,db_g)

	write (13,*) rho,b,b_g
	write (11,*) rho,db,db_g

	rho = rho + drho
	end do

	close(13)
	close(11)

end subroutine b_plot_gauss


!-----------------------------------------end of stuff that used to be in root-----------------------------------------




end module solver

