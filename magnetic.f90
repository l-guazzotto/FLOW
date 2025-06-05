module magnetic

  use constant
  implicit none

  ! A scale for Bphi, fvaccuum and pcenter
  real (kind=dkind) :: b_phi_zero

  integer :: F_opt  ! option for b_zero(psi)

  ! Using the usual notation, $ F = B_\phi * R $
  real (kind=dkind) :: fvaccuum

  ! fcenter is F = B_phi*R at the Magnetic Axis
  real (kind=dkind) :: fcenter

  ! If F_opt = 1
  ! B_phi * R = F = fvaccuum + (fcenter - fvaccuum) * (psi/psic)^kappa
  real (kind=dkind) :: kappa

  ! If F_opt = 2
  ! B_phi * R = F = ( fvaccuum**2-eta_P*2.d0*mu_mag*rmajor**2*pofpsi(psi) )**0.5
  real(kind=dkind) :: eta_P

  ! stuff for super-Alfvénic equilibria begins here...
  real (kind=dkind) :: b_vert
  real (kind=dkind) :: b_edge

  real (kind=dkind) :: ffrac		!for bzero(psi)
  real (kind=dkind) :: amult		!for bzero(psi)
  real (kind=dkind) :: pmult		!for bzero(psi)	
  real (kind=dkind) :: pfact		!for bzero(psi)
  real (kind=dkind) :: pexp			!for bzero(psi)
  real (kind=dkind) :: p_eps		!for bzero(psi)
  ! ... and ends here

  ! mu_mag is $\mu$ for the magnetic field.
  ! In mks this is $\mu_0 = 4.0d-7*pi$.
  ! In cgs this is $4\pi$
  ! Riccardo prefers to set this to 1
  real (kind=dkind) :: mu_mag		!

  real (kind=dkind) :: psi_degen = 0.d0
  ! 12/20/2008 psi_degen moved here for consistency

  ! Psi is defined to equal 0 at the boundary 
  real (kind=dkind) :: psic = .1d0 !0.062d0 ! $\psi$ at the magnetic axis
  real (kind=dkind) :: psic_13 = .81d0 !0.062d0 ! $\psi$ at the magnetic axis
  real (kind=dkind) :: delta_psi = 1.d0/30.d0 ! for numerical input

  real (kind=dkind) :: fraction 
  ! this is to set inner and outer regions:
  ! psi/psic >  frac -> inner region,
  ! psi/psic <  frac -> outer region

  real(kind=dkind) :: psi_ext
  ! for boundary condition on external boundary

  real(kind=dkind) :: psi_in, psi_out
  ! for boundary conditions for LDX

  real(kind=dkind) :: psi_pmax
  ! psi of maximum p in LDX profiles (eq3_opt==9)

  real(kind=dkind) :: B_pmax
  ! B at maximum p in LDX profiles (eq3_opt==9)

  real(kind=dkind) :: I_coil = -1.5d6
  !current in LDX coil

  real(kind=dkind) :: delta_coil = -5.d-2
  !coil shift with respect to the geometric center

  real(kind=dkind) :: bpol0_fix = 1.d-1
  real(kind=dkind) :: bpol0_temp = 1.d-1
  real(kind=dkind) :: psi_B_min = 0.3038d0 !(psi/psic)
  real(kind=dkind) :: ana_fac0 = 1.d-1
  real(kind=dkind) :: ana_fac = 1.d-1
  !for fiddling with anisotropy in the LDX case

  real(kind=dkind) :: psi_e
  ! position of the boundary in the separatrix case (tri_type==13)

  integer :: inter_switch

  real(kind=dkind), dimension(:,:), allocatable :: modB_coef
  ! added on 11/27/2007 for trapped particle calculation
  ! IMPORTANT NOTE. This is used differently from other interpolations:
  ! all data is in the same array (therefore defined with size 4,:)

  integer :: modB_ord = 2 ! interpolation order for |B|
  integer :: ibreak_modB

  real(kind=dkind) :: curr
  ! toroidal current

  real(kind=dkind) :: betator, betapol
  ! toroidal and poloidal beta

  real(kind=dkind) :: betastar
  ! beta star for eqdsk output (never heard of this thing before)

  real(kind=dkind), dimension(:,:), allocatable :: q_coef
  ! added on 12/29/2007 for q profile interpolation for eqdsk
  ! IMPORTANT NOTE. This is used differently from other interpolations:
  ! all data is in the same array (therefore defined with size (4,:))

  real(kind=dkind), dimension (:), allocatable :: qval
  ! moved here to be used in the bootstrap current calculation

  integer :: q_ord = 3 ! interpolation order for |B|
  integer :: ibreak_q

  real(kind=dkind) :: q_c, qe
  ! axis and edge safety factor

  real(kind=dkind) :: mu_RFP
  ! mu0 for RFP B_0 profile

  real(kind=dkind), dimension(:), allocatable :: psival

!____________________________________
! what follows is stuff for the magnetic output,
! moved here on January 10 2008 to humor the LINUX compiler

	integer, parameter :: ord_surf = 3
	integer, parameter :: ord_loc = 2
	integer :: nsurf = 0

	integer :: enq = 0

	real(kind=dkind), dimension(:,:,:), allocatable :: rtab, bigrtab, ztab
	real(kind=dkind), allocatable, dimension(:,:) ::thetatab, shear,  &
			bscoef_psi, bscoef_bpol,bscoef_bphi, bscoef_rho, bscoef_rho2
	real(kind=dkind), allocatable, dimension(:) :: xknot_mag, zknot_mag

	real(kind=dkind), dimension(:,:), allocatable :: bscoef_boot_base

	real(kind=dkind) :: xmax_mag, zmax_mag ! position of the magnetic axis
	real(kind=dkind) :: cloc, sloc, psiloc_mag

	real(kind=dkind), dimension (:), allocatable :: Rleft,Rright
	! position of the various magnetic surfaces along the midplane

	integer :: surf_index
	! to pass an index to a one-argument funtion

	real(kind=dkind), dimension(:,:), allocatable :: cross_section
	! cross section for poloidal flow

! end of magnetic output stuff
!_____________________________________

! stuff for bootstrap current calculation
! started on January 21 2008

  integer :: bootstrap_option = 2
  ! 0->Sauter, 1->NCLASS, 2->both

  real(kind=dkind), dimension(:), allocatable :: eff_trap
  ! effective trapped particle fraction (as a function of psi)

  real(kind=dkind), dimension(:), allocatable :: surf_length
  ! length of magnetic surface, for averaging

  real(kind=dkind), dimension(:), allocatable :: B2_ave
  ! average of B**2
  real(kind=dkind), dimension(:), allocatable :: Bm2_ave
  ! average of 1/B**2
  real(kind=dkind), dimension(:), allocatable :: B2_hat_ave
  ! average of (B/Bmax)**2

  real(kind=dkind), dimension(:), allocatable :: boot_ne, boot_pe, boot_Te, boot_ni, boot_pi, boot_Ti
  real(kind=dkind), dimension(:), allocatable :: boot_neprim, boot_peprim, boot_Teprim,  &
																	  boot_niprim, boot_piprim, boot_Tiprim
  real(kind=dkind), dimension(:), allocatable :: boot_tor_flux, boot_tor_rho, boot_fhat, boot_grho_ov_B2
  ! various physical quantities (names are intuitive)
  ! boot_tor_rho = normalized sqrt(toroidal flux) [m] (rho from now on)
  ! boot_fhat = something with F and (d rho / d psi)
  ! boot_ghro_ov_B2 = <(grad rho/B)^2>

  real(kind=dkind), allocatable, dimension(:,:) :: bscoef_tor_rho

  real(kind=dkind), dimension(:), allocatable :: inv_asp_ratio
  ! "aspect ratio" of each magnetic surface (no elongation!)

  real(kind=dkind), dimension(:), allocatable :: boot_Rcenter
  ! geometrical center of each magnetic surface (no elongation!)

  real(kind=dkind), dimension(:), allocatable :: boot_Bmin
  ! minimum value of |B| on each magnetic surface

  real(kind=dkind), dimension(:), allocatable :: boot_Bp_at_Bmin
  ! value of B_poloidal at the location of the minimum
  ! value of |B| on each magnetic surface

  real(kind=dkind), dimension(:), allocatable :: boot_R_of_Bmin
  ! major radius at the location of the minimum
  ! value of |B| on each magnetic surface

  real(kind=dkind), dimension(:,:), allocatable :: J_boot
  ! bootstrap current stuff

  real(kind=dkind), dimension(:,:), allocatable :: el_resistivity

  real(kind=dkind), dimension(:), allocatable :: surf_ave_base
  ! average on a magnetic surface of the "base"/"weight" quantity used for averaging

  real(kind=dkind), dimension(:), allocatable :: w_ave_int, t_ave_int
  ! weights and locations for integrations to compute surface averages

  integer :: n_ave_int = 100

  integer :: ord_ave = 2
  ! order of interpolation for averaging over a magnetic surface

  real(kind=dkind), dimension(:), allocatable :: JparB_ave
  ! average over a magnetic surface of J_par B

  real(kind=dkind) :: bpolav
  ! average poloidal field (moved here on February 6 2008 for bootstrap calculation)

  real(kind=dkind) :: Z_eff_0, delta_Z_eff, Z_eff_exp
  ! for Z_effective in bootstrap calculation: Z_eff = Z_eff_0 + delta_Z_eff * (psi/psic)**Z_eff_exp

end module magnetic
