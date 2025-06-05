module constant

  implicit none

  ! Compiler specific parameters
!  integer, parameter :: skind = 4 ! 4 byte real.
!  integer, parameter :: dkind = 8 ! 8 byte real.

  integer, parameter :: skind = kind(1.0e0)
  integer, parameter :: dkind = kind(1.0d0)

  ! Essentially Global parameters for the torus
  
  ! RMAJOR is the major radius of the Tokamak
  real (kind=dkind) :: rmajor

  ! rminor is the minor radius of the plasma
  real (kind=dkind) :: rminor

  ! x_size is the size of the computational domain in the x-direction
  real (kind=dkind) :: x_size

  ! z_size is the size of the computational domain in the z-direction
  real (kind=dkind) :: z_size

  integer, dimension(:,:,:), allocatable :: sort_grid
!  real(kind=dkind), dimension(:,:,:), allocatable :: sort_grid
  !distinguishes between internal and external grid points;
  !for LDX, second index marks the zone

  ! Pi
  real (kind=dkind), parameter :: pi = &
       3.14159265358979323846d0

  ! The problem to be solved:
  ! 1=MHD, 2=CGL, 3=KINETIC WITH TOROIDAL FLOW ONLY,
  ! 4=LDX WITH TOROIDAL FLOW ONLY
  integer :: eq_type
  integer :: eq3_opt
  ! if eq_type==3 this gives:
  ! eq3_opt=1	-> T_par(psi) is assigned
  ! eq3_opt=2	-> P_par(psi) = D(psi)*T_par(psi) is assigned
  ! eq3_opt=3	-> P_par(psi), theta(psi) are assigned (for LDX)
  ! eq3_opt=4	-> P_par(psi), theta(psi) are assigned from V(psi) (for LDX)
  ! eq3_opt=5	-> as 4,  but P is forced to 0 in a range close to psimax

  ! Are we using numerical data?
  logical :: numerical
  integer :: numerical_opt	! if so, how?
  ! if((numerical) .AND.(eq_type==3)) -> eq3_opt=2
  integer :: data_dim, vol_dim, psiprim_dim, psib_dim
  ! number of points in numerical input

  ! Broot is a control parameter for the Bernoulli solver
  ! Broot = 0 -> A maximum mach theta max is searched for which gives
  !              a subsonic solution interior of a supersonic soln.
  ! Broot = 1 -> Choose the subsonic branch (The heavy solution)
  ! Broot = 2 -> Choose the supersonic branch (The light solution)
  ! Broot = 3 -> Super-Alfvénic equilibrium
  ! Broot = 4 -> RFP equilibrium with subsonic solution
  ! Broot = 5 -> RFP equilibrium with transonic solution
  integer :: Broot

  ! The ion mass (this better be set equal to 1 if mu_mag is 1)
  real (kind=dkind) :: mass
  ! The elementary charge
  real (kind=dkind) :: eV

  real(kind=dkind) :: me_ov_mi
  ! ratio between electron and ion mass

  real(kind=dkind) :: pe_ov_p = 0.5d0
  ! ratio between electron and total pressure

  integer :: n	   ! number of points of the finest grid
  integer :: n_min ! number of points of the coarsest grid
  integer :: min_it ! minimum number of iterations for each grid
  integer :: max_it ! maximum number of iterations for each grid
  logical :: accelerate ! whether to use Chebyshev acceleration
  real (kind=dkind) :: fix_orp
  ! if no acceleration is used, the over relaxation parameter is equal to fix_orp
  real (kind=dkind) :: max_orp ! maximum over relaxation parameter if bc_type=2
  integer :: bc_type	! type of BC's:
  ! 1 -> psi=0 for all outer points
  ! 2 -> psi=0 on the geometric boundary (1D linear interpolation used) OBSOLETE!!
  ! 3 -> psi=0 on the geometric boundary (2D linear interpolation used)
  ! 4 -> 1 for n<bc_switch, 3 for n>=bc_switch
  ! 7 -> for free-boundary (psi=0 on plasma boundary, psi assigned on fixed boundary)
  ! 8 -> like 7, but no need for a magnetic axis in the plasma

  integer :: bc_setup_option
  ! >0 for old method (solution of system)
  ! <0 for new method (minimization of distance)

  integer :: bc_switch
  ! when to switch bc for bc_type=4, 7, 8

  real (kind=dkind) :: fix_orp1 = 1.1d-5
  real (kind=dkind) :: fix_orp0

  real(kind=dkind), dimension(:), allocatable :: x_coord, z_coord
  ! these arrays are set once for each grid

  real(kind=dkind) :: rcenter
  ! geometrical center of the plasma

  real(kind=dkind) :: zcenter
  ! geometrical center of the plasma

  real(kind=dkind) :: dx, dz

  real(kind=dkind), dimension(:), allocatable :: dx_a, dz_a
  ! these arrays are set once for each grid

  real(kind=dkind), dimension(:), allocatable :: dx_ext, dz_ext
  ! arbitrary grids are assigned in these arrays from external input

  integer :: grid_type
  ! for non-uniform grid distribution:
  ! 0=uniform grid
  ! 1=linear increment grid in R

  logical :: write_all
  ! whether to write all output

  logical :: write_all_bin
  !whether to save the output in binary form

!  logical :: FLOS_input
  ! whether to write the input for FLOS
  ! January 9 2008: the FLOS project has been abandoned, all the FLOS stuff is hence removed

  logical :: restart
  ! whether to proceed a previous run or start from scratch
  ! NOTE: for the time being, the grid assignement is left in inputfile.dat,
  ! and the input is not repeated: make sure that inputfile.dat is consistent!

  integer :: r_of_theta_points = 200
  ! how many points to use when saving the plasma shape

  integer :: mir, miz
  ! the index position of the magnetic axis

  real(kind=dkind) :: inv_aspect_ratio
  ! the inverse aspect ratio, only used for output purposes

  logical :: input_EQDSK
  ! added on January 07 2008 for EQDSK input

  character (len=64)  :: EQDSK_file
  ! the name of the eqdsk file to be read as input


!_________________________NCLASS stuff_________________________

  integer, parameter :: mx_ms=40
  integer, parameter :: mx_mi=9
  integer, parameter :: mx_mz=18

!_________________________end of NCLASS stuff_________________________

!_________________________gravity stuff___________________________

  integer :: gravity_type
  ! 0 -> no gravity; 1 -> point mass at origin; 2-> constant gravity in the R direction

  real(kind=dkind) :: G_gravity
  ! gravitational constant (6.67300d-11 in SI units)

  real(kind=dkind) :: M_gravity
  ! mass originating the gravity field

   logical :: Kepler_Omega
   ! whether to use (constant at R0) Keplerian frequency


!_________________________end of gravity stuff___________________________

 integer :: jump_option = 1 ! for transonic equilibria:
! -3 -> like -2, with "consistent" calculation of |grad psi| in half grid points
! -2 -> use new formula with grad Ptot in the RHS
! -1 -> use no van Leer limiter in ngs_solve
! 0 -> use old ngs_solve
! 1-> impose jump condition
! 2-> as 1, with lambda0 fixed in time
! 3 -> as 2, but lambda0 also constant in space
! 4 ->as in 3, but changes sign on the two sides of the transonic surface
 real(kind=dkind) :: lambda0_fix
 integer :: dir_switch	! for jump_option>=-2
 real(kind=dkind) :: delta_Bern_fact = 1.d0

	logical :: MARS_output
	! whether to save file for MARS

end module constant
