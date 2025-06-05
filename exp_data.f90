module exp_data

use constant

integer, parameter :: data_dim_lim = 2000
integer :: ibreak_p, ibreak_fi, ibreak_th, ibreak_b0, ibreak_v, ibreak_pp, ibreak_d,  &
			ibreak_ppar, ibreak_pperp, ibreak_pb

real(kind=dkind) p_iso_data(data_dim_lim,3)	!numerical data for p(psi)
real(kind=dkind) p_par_data(data_dim_lim,3)	!numerical data for p_par(psi)
real(kind=dkind) p_perp_data(data_dim_lim,3)	!numerical data for p_perp(psi)
real(kind=dkind) omega_data(data_dim_lim,3)	!numerical data for omega(psi)
real(kind=dkind) F_data(data_dim_lim,3)	!numerical data for F(psi)=R0*B0(psi)
real(kind=dkind) d_data(data_dim_lim,3)	!numerical data for D(psi)
real(kind=dkind) vol_data(data_dim_lim,3)	!numerical data for p(psi) (LDX)
real(kind=dkind) psiprim_data(data_dim_lim,3)	!numerical data for bc_type=5 (LDX)
real(kind=dkind) psib_data(data_dim_lim,3)	!numerical data for bc_type=7, 8 (LDX, free-boundary)

logical ::	numerical_n , numerical_p_iso , numerical_p_par ,  &
			numerical_p_perp , numerical_F , numerical_omega	!whether to use numerical data

logical :: numerical_psiprim

logical :: numerical_mphi, numerical_mtheta

real(kind=dkind) mphi_data(data_dim_lim,3)	!numerical data for Mphi(psi)
real(kind=dkind) mtheta_data(data_dim_lim,3)	!numerical data for Mtheta(psi)

real(kind=dkind) X0_data(data_dim_lim,3)	!numerical data for X0(psi)
real(kind=dkind) X1_data(data_dim_lim,3)	!numerical data for X1(psi)

real(kind=dkind), dimension(:,:), allocatable :: psi_bscoef

real(kind=dkind), dimension(:,:), allocatable :: p_iso_cscoef, F_cscoef,   &
												 mphi_cscoef, mtheta_cscoef,  &
												 X0_cscoef, X1_cscoef, vol_cscoef,  &
												 psiprim_cscoef, d_cscoef,  &
												 ppar_cscoef, pperp_cscoef,  &
												 psib_cscoef
! for IMSL spline routine

integer :: p_iso_ord = 3	! interpolation order
integer :: F_ord = 3		! interpolation order
integer :: mphi_ord = 3		! interpolation order
integer :: mtheta_ord = 3	! interpolation order
integer :: vol_ord = 3		! interpolation order
integer :: psiprim_ord = 3	! interpolation order
integer :: d_ord = 3		! interpolation order
integer :: ppar_ord = 3	! interpolation order
integer :: pperp_ord = 3	! interpolation order
integer :: psib_ord = 2		! interpolation order
integer :: s_ord = 3	! interpolation order for psi
!integer :: p_iso_ord = 3	! interpolation order


integer :: omega_option		!1=mphi, 2=omega

real(kind=dkind), dimension(:,:), allocatable :: r_data_psi
! to update the interface for the separatrix case (tri_type==13)

real(kind=dkind), dimension(:), allocatable :: xknot_psi, zknot_psi
integer :: nx_FLOW, nz_FLOW

end module exp_data