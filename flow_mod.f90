module FLOW_mod
  use constant
  use solver
  use exp_data
  use triangularity
  use magnetic
  use p_d_profile

  implicit none

  ! Interface to allow some input parameters to be passed in functionally rather than via file IO
  ! As well as allowing FLOW output file data to be passed directly to this module and subsequently to Python code
  ! (Interface is necessary for optional & keyword arguments in FORTRAN)
  interface 
          subroutine readinput(m,input_n, input_p_iso, input_p_par, input_p_perp,&
                          input_F, input_omega, input_mtheta, input_psiprim,&
                          input_a_elps, input_b_elps, input_k_ellipt, input_delta_up, input_delta_down)
                  integer :: m
                  real, dimension(:,:), intent(in), optional :: input_n, input_p_iso, input_p_par, input_p_perp,&
                          input_F, input_omega, input_mtheta, input_psiprim
                  real, dimension(:), intent(in), optional :: input_a_elps, input_b_elps, input_k_ellipt,&
                          input_delta_up, input_delta_down
          end subroutine
        
          subroutine radial_plot(psi, truepsi, nx, nz, fname, iz, fname_output_data)
                  import :: skind
                  import :: dkind
                  integer, intent(in) :: nx,nz,iz
                  real (kind=skind), dimension(1:nx,1:nz), intent(in) :: psi
                  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: truepsi
                  character (len=*), intent(in) ::fname
                  real (kind=dkind), dimension(1:nx*nz,1:3), intent(out), optional :: fname_output_data
          end subroutine

          !subroutine psiout(phi_output_data)
          !        !integer, parameter :: nsave = 201
          !        real (kind=dkind), dimension(1:nsave+1,1:2), intent(out) :: phi_output_data
          !end subroutine
  end interface

  integer :: m	!= 1 + n/2
  real (kind=dkind), dimension(:,:), allocatable :: psi, rho, residual,b_phi
  real (kind=dkind), dimension(:,:), allocatable :: vr,vphi,vz,br,bz,p,ppar,pperp
  real (kind=dkind), dimension(:,:), allocatable :: temp,tpar,tperp,mtor	!bphi,
  real (kind=dkind), dimension(:,:), allocatable :: malfven2 ! Alfvem mode Mach number squared
  real (kind=dkind), dimension(:,:), allocatable :: mslow ! slow mode Mach number
  real (kind=dkind), dimension(:,:), allocatable :: mpol ! poloidal sonic Mach number
  real (kind=dkind), dimension(:,:), allocatable :: mcusp ! cusp Mach number
  real (kind=dkind), dimension(:,:), allocatable :: beta,betapar,betaperp ! plasma beta
  real (kind=skind), dimension(:,:), allocatable :: out ! For outputing data
  real (kind=dkind), dimension(:,:), allocatable :: outd ! For outputing data
  real (kind=dkind), dimension(:,:), allocatable :: j_phi ! toroidal current
  real (kind=dkind), dimension(:,:), allocatable :: j_phi_new ! toroidal current (from derivatives)
  real (kind=dkind), dimension(:,:), allocatable :: bpol ! poloidal component of the field
  real (kind=dkind), dimension(:,:), allocatable :: rbtheta ! for cylindrical checks
  real (kind=dkind), dimension(:,:), allocatable :: trapped ! fraction of trapped particles
  real (kind=dkind), dimension(:,:), allocatable :: j_par ! parallel current
  real (kind=dkind), dimension(:,:), allocatable :: j_x ! R current
  real (kind=dkind), dimension(:,:), allocatable :: j_z ! Z current
  real (kind=dkind), dimension(:,:), allocatable :: cs ! sound speed
  real (kind=dkind), dimension(:,:), allocatable :: csp ! poloidal sound speed
  real (kind=dkind), dimension(:,:), allocatable :: hyper ! discriminant for hyperbolic region
  real (kind=dkind), dimension(:,:), allocatable :: br_2,bz_2,br_3,bz_3, br_gg, bz_gg
  real (kind=dkind), dimension(:,:,:), allocatable :: gg

  integer :: alloc_stat
  integer :: ir, iz
  real (kind = dkind) :: mx
  integer i
  real (kind=dkind) :: x,y
  real (kind=dkind) :: vb,vs
  
  ! Variables to go into readinput() if assigned values through the Python script
  real, dimension(:,:), allocatable :: n_density
  real, dimension(:,:), allocatable :: p_iso
  real, dimension(:,:), allocatable :: p_par
  real, dimension(:,:), allocatable :: p_perp
  real, dimension(:,:), allocatable :: F
  real, dimension(:,:), allocatable :: omega
  real, dimension(:,:), allocatable :: mtheta
  ! Deprecated: real, dimension(:,:), allocatable :: psiprim
  
  ! More variables to go into readinput() like above, but these ones are scalars
  ! Note: yes, they're not really scalars, but f2py only knows FORTRAN 90, and allocatable scalars weren't a thing yet.
  ! Hence this dumb thing were these inputs are actually rank-1 arrays, but I'll only be using the first value
  real, dimension(:), allocatable :: input_a_elpss, input_b_elpss, input_k_ellipts, input_delta_ups, input_delta_downs
  
  ! Just holds the viscosity data for output
  real (kind=dkind), dimension(:,:), allocatable :: pol_viscosity_data
  ! As well as other output data for SHS and Braginskii viscosity calculation
  real (kind=dkind), dimension(:,:), allocatable :: psi_data, temp_data, rho_data, mpol_data
  real (kind=dkind), dimension(:,:), allocatable :: vr_data, vphi_data, vz_data, vp_data
  real (kind=dkind), dimension(:,:), allocatable :: br_data, bphi_data, bz_data, bpol_data, bfield_data
  ! data above from radial_plot (data_dump.f90), phi_data from psiout (trans_solve.f90)
  real (kind=dkind), dimension(:,:), allocatable :: phi_data

  contains
    subroutine run_FLOW()

  ! Could be useful for debug purposes
  !if (allocated(mtheta)) then
  !       print *, "mtheta present: "
  !       print *, mtheta
  !else
  !       print *, "mtheta not supplied functionally, therefore using datafiles instead"
  !endif
  
  
  call readinput(m,input_n=n_density,input_p_iso=p_iso,input_p_par=p_par,input_p_perp=p_perp,&
        input_F=F,input_omega=omega,input_mtheta=mtheta,&
        input_a_elps=input_a_elpss, input_b_elps=input_b_elpss, input_k_ellipt=input_k_ellipts,&
        input_delta_up=input_delta_ups, input_delta_down=input_delta_downs)

! allocation section
! this section is necessary because n is an input value and not a parameter

! Allocation section for Python output data holders
  allocate(psi_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(temp_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(rho_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(mpol_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(vr_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(vphi_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(vz_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(vp_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif
  
  allocate(br_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(bphi_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(bz_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  allocate(bpol_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif
  
  allocate(bfield_data(1:n*n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

  ! WARNING: hardcoding the array size for phi(psi) here,
  ! because nsave is just a hard coded value in psiout subroutine anyway
  ! This is still terrible coding practice and should be fixed in the future
  allocate(phi_data(1:201+1,1:2),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi_data"
  endif

! Initializations for Python output data holders
  psi_data(1:n*n,1:3) = 0.d0
  temp_data(1:n*n,1:3) = 0.d0
  rho_data(1:n*n,1:3) = 0.d0
  mpol_data(1:n*n,1:3) = 0.d0
  vr_data(1:n*n,1:3) = 0.d0
  vphi_data(1:n*n,1:3) = 0.d0
  vz_data(1:n*n,1:3) = 0.d0
  vp_data(1:n*n,1:3) = 0.d0
  br_data(1:n*n,1:3) = 0.d0
  bphi_data(1:n*n,1:3) = 0.d0
  bz_data(1:n*n,1:3) = 0.d0
  bpol_data(1:n*n,1:3) = 0.d0
  bfield_data(1:n*n,1:3) = 0.d0
  phi_data(1:201+1,1:2) = 0.d0

!  allocate(psi(0:n+1,0:n+1),stat = alloc_stat)
  allocate(psi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in psi"
  end if
  psi(1:n,1:n) = 0.d0

  allocate(rho(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in rho"
  end if
  rho(1:n,1:n) = 0.d0

  allocate(residual(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in residual"
  end if
  residual(1:n,1:n) = 0.0

  allocate(b_phi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in b_phi"
  end if
  b_phi(1:n,1:n) = 0.0

! end of allocation section

  if((tri_type==11).and.((eq3_opt==7).or.(eq3_opt==8).or.(eq3_opt==9)).or.  &
			((p_opt==7).or.(p_opt==8).or.(p_opt==9))) then
	psi_pres = psi_in
	psic = psi_pres
  endif

!  if(Broot==3) call psiout
   call psiout(phi_data)
   call FINESSE_input


!----------------------------------------------------------------------------

  call set_triangularity(asin_d_up, asin_d_down,   &
						theta_dat,rminor_dat,0.d0,0.d0,d2rminor_dat,  &
						rcos_u,rcos_d,rsin_u,rsin_d,  &
						bigR_cos,bigR_sin,Z_cos,Z_sin,'R_Z.dat'  )

	if((write_all).and.(tri_type/=-1).and.(tri_type/=0)  &
		.and.(tri_type/=9).and.(tri_type/=11)) call get_minor_radius

!----------------------------------------------------------------------------

	open(111,file='residual.dat',status='unknown',action='write')

  call mgrid(psi,rho,residual,b_phi)

  close(111)

  ! find the max of psi and rho
  ! -----------------------------------------------------
  mx = 0.0d0
  do iz = 1, n
     do ir = 1, n
		if(sort_grid(ir,iz,0)<0) cycle
        if(((bc_type==7).and.(psi(ir,iz)) > mx).or.((bc_type/=7).and.(dabs(psi(ir,iz)) > mx))) then
           mx = psi(ir,iz)
           mir = ir
           miz = iz
        end if
     end do
  end do
  print *, "Psi_max = ",mx," at (",mir,",",miz,")"," (",  &
       x_coord(mir),z_coord(miz),")"
  ! -----------------------------------------------------

! output allocation section

  allocate(vr(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in vr"
  end if
  vr(1:n,1:n) = 0.0

  allocate(vphi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in vphi"
  end if
  vphi(1:n,1:n) = 0.0

  allocate(vz(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in vz"
  end if
  vz(1:n,1:n) = 0.0

  allocate(br(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in br"
  end if
  br(1:n,1:n) = 0.0

  allocate(bz(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in bz"
  end if
  bz(1:n,1:n) = 0.0

  allocate(br_2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in br_2"
  end if
  br_2(1:n,1:n) = 0.0

  allocate(bz_2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in bz_2"
  end if
  bz_2(1:n,1:n) = 0.0

  allocate(br_3(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in  br_3"
  end if
  br_3(1:n,1:n) = 0.0

  allocate(bz_3(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in bz_3"
  end if
  bz_3(1:n,1:n) = 0.0

  allocate(p(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in p"
  end if
  p(1:n,1:n) = 0.0

  allocate(ppar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in ppar"
  end if
  ppar(1:n,1:n) = 0.0

  allocate(pperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error pperp"
  end if
  pperp(1:n,1:n) = 0.0

  allocate(temp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in temp"
  end if
  temp(1:n,1:n) = 0.0

  allocate(tpar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in tpar"
  end if
  tpar(1:n,1:n) = 0.0

  allocate(tperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in tperp"
  end if
  tperp(1:n,1:n) = 0.0

  allocate(mtor(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in mtor"
  end if
  mtor(1:n,1:n) = 0.0

  allocate(malfven2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in malfven2"
  end if
  malfven2(1:n,1:n) = 0.0

  allocate(mslow(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in mslow"
  end if
  mslow(1:n,1:n) = 0.0

  allocate(mcusp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in mcusp"
  end if
  mcusp(1:n,1:n) = 0.0

  allocate(mpol(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in mpol"
  end if
  mpol(1:n,1:n) = 0.0

  allocate(beta(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in beta"
  end if
  beta(1:n,1:n) = 0.0

  allocate(betapar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in betapar"
  end if
  betapar(1:n,1:n) = 0.0

  allocate(betaperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in betaperp"
  end if
  betaperp(1:n,1:n) = 0.0

  allocate(j_phi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in j_phi"
  end if
  j_phi(1:n,1:n) = 0.0

  allocate(j_phi_new(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in j_phi_new"
  end if
  j_phi_new(1:n,1:n) = 0.0

  allocate(j_par(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in j_par"
  end if
  j_par(1:n,1:n) = 0.0

  allocate(j_x(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in j_x"
  end if
  j_x(1:n,1:n) = 0.0

  allocate(j_z(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in j_z"
  end if
  j_z(1:n,1:n) = 0.0

  allocate(rbtheta(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in rbtheta"
  end if
  rbtheta(1:n,1:n) = 0.0

  allocate(trapped(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in trapped"
  end if
  trapped(1:n,1:n) = 0.0

  allocate(cs(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in cs"
  end if
  cs(1:n,1:n) = 0.d0

  allocate(csp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in csp"
  end if
  csp(1:n,1:n) = 0.d0

  allocate(hyper(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in hyper"
  end if
  hyper(1:n,1:n) = 0.d0

  allocate(br_gg(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in br_gg"
  end if
  br_gg(1:n,1:n) = 0.0

  allocate(bz_gg(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in bz_gg"
  end if
  bz_gg(1:n,1:n) = 0.0

  allocate(gg(1:n,1:n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in gg"
  end if
  gg(1:n,1:n,1:3) = 0.0


  allocate(out(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in out"
  end if
  out(1:n,1:n) = 0.0

  allocate(outd(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error in outd"
  end if
  outd(1:n,1:n) = 0.d0

! end of allocation section

  out(1:n,1:n) = psi(1:n,1:n)

  if(tri_type==9) then

	  call radial_plot_2(psi,n,n,"psi")

  else

	  call radial_plot(out,psi,n,n,"psi",m,psi_data)

  endif

  call radial_plot(out,psi,n,n,"psi_clean",m)

!!$  call q_data(out,n,n,"psi")

  out(1:n,1:n) = rho(1:n,1:n)/mass
  call radial_plot(out,psi,n,n,"rho",m,rho_data)

  out(1:n,1:n) = residual(1:n,1:n)
  call radial_plot(out,psi,n,n,"residual",m)

	call r_of_theta

	call write_restart_data(psi,n,n,"psi")
	call write_restart_data(rho,n,n,"rho")

	call psiout(phi_data)

  call physical_var(psi,rho,n,n,vr,vphi,vz,br,b_phi,bz,p,ppar,pperp,  &
					malfven2,mslow,mcusp,mpol,beta,betapar,betaperp,j_phi,temp,  &
					tpar,tperp,mtor,j_par,j_x,j_z,cs,csp,hyper,br_2,bz_2,br_3,bz_3,gg,  &
					br_gg, bz_gg, j_phi_new)

	call write_restart_data(b_phi,n,n,"b_phi")

	if(write_all_bin) then

		call write_restart_data(vr,n,n,"vr")
		call write_restart_data(vphi,n,n,"vphi")
		call write_restart_data(vz,n,n,"vz")

		call write_restart_data(br,n,n,"br")
		call write_restart_data(bz,n,n,"bz")

		call write_restart_data(p,n,n,"pres")

	endif

  if(tri_type==11) call get_rbtheta(rbtheta,br,bz,n)

  !-------------------------------------------------------------------

  out(1:n,1:n) = vphi(1:n,1:n)
  call radial_plot(out,psi,n,n,"vphi",m,vphi_data)



  out(1:n,1:n) = br(1:n,1:n)
  call radial_plot(out,psi,n,n,"br",m,br_data)

!!$  out(1:n,1:n) = br_2(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"br_2",m)

!!$  out(1:n,1:n) = br_3(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"br_3",m)

  out(1:n,1:n) = b_phi(1:n,1:n)
  call radial_plot(out,psi,n,n,"bphi",m,bphi_data)

!!$  call q_data(out,n,n,"bphi")

  out(1:n,1:n) = bz(1:n,1:n)
  call radial_plot(out,psi,n,n,"bz",m,bz_data)

!!$  out(1:n,1:n) = bz_2(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"bz_2",m)

!!$  out(1:n,1:n) = bz_3(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"bz_3",m)

  ! Calculate bpoloidal
  allocate(bpol(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  bpol(1:n,1:n) = 0.d0

  do iz = 1, n
     do ir = 1, n
		bpol(ir,iz) = dsqrt(br(ir,iz)**2 + bz(ir,iz)**2)
        out(ir,iz) = bpol(ir,iz)
     end do
  end do
  call radial_plot(out,psi,n,n,"bpol",m,bpol_data)

!!$  call q_data(out,n,n,"bpol")

  ! Calculate b total
  do iz = 1, n
     do ir = 1, n
        out(ir,iz) = dsqrt(br(ir,iz)**2 + bz(ir,iz)**2+b_phi(ir,iz)**2)
     end do
  end do
  call radial_plot(out,psi,n,n,"b_field",m,bfield_data)

  if (eq_type == 1) then

	out(1:n,1:n) = p(1:n,1:n)
	call radial_plot(out,psi,n,n,"pres",m)
!	if(tri_type==11) call q_data(out,n,n,"pres")

  elseif ((eq_type==2).OR.(eq_type==3)) then

    out(1:n,1:n) = ppar(1:n,1:n)
    call radial_plot(out,psi,n,n,"p_par",m)

	out(1:n,1:n) = pperp(1:n,1:n)
    call radial_plot(out,psi,n,n,"p_perp",m)
  endif

  out(1:n,1:n) = temp(1:n,1:n)
  call radial_plot(out,psi,n,n,"temp",m,temp_data)

  if ((eq_type==2).OR.(eq_type==3)) then

    out(1:n,1:n) = tpar(1:n,1:n)
    call radial_plot(out,psi,n,n,"T_par",m)

	out(1:n,1:n) = tperp(1:n,1:n)
    call radial_plot(out,psi,n,n,"T_perp",m)

  endif

  if(eq_type==3) then

	continue

  else

	  vb=0.d0
	  vs=0.d0

	  do iz=1,n
		do ir=1,n

			vb=max(vr(ir,iz),vb)
			vs=min(vr(ir,iz),vs)

		enddo
	  enddo


	  if( (vb+abs(vs))>0.d0) then

		  out(1:n,1:n) = vz(1:n,1:n)
		  call radial_plot(out,psi,n,n,"vz",m,vz_data)

		  out(1:n,1:n) = vr(1:n,1:n)
		  call radial_plot(out,psi,n,n,"vr",m,vr_data)

		  ! Calculate vpoloidal
		  do iz = 1, n
			 do ir = 1, n
				out(ir,iz) = dsqrt(vr(ir,iz)**2 + vz(ir,iz)**2)
			 end do
		  end do
		  call radial_plot(out,psi,n,n,"vp",m,vp_data)

		  out(1:n,1:n) = malfven2(1:n,1:n)
		  call radial_plot(out,psi,n,n,"malfven2",m)

	  out(1:n,1:n) = mpol(1:n,1:n)
	  call radial_plot(out,psi,n,n,"mpol",m,mpol_data)

!	  out(1:n,1:n) = mcusp(1:n,1:n)
!	  call radial_plot(out,psi,n,n,"mcusp",m)

	  endif

  endif

  out(1:n,1:n) = mtor(1:n,1:n)
  call radial_plot(out,psi,n,n,"mtor",m)

  out(1:n,1:n) = beta(1:n,1:n)
  call radial_plot(out,psi,n,n,"beta",m)

  if ((eq_type==2).OR.(eq_type==3)) then

    out(1:n,1:n) = betapar(1:n,1:n)
    call radial_plot(out,psi,n,n,"betapar",m)

    out(1:n,1:n) = betaperp(1:n,1:n)
    call radial_plot(out,psi,n,n,"betaperp",m)

  endif

  out(1:n,1:n) = j_phi(1:n,1:n)
  call radial_plot(out,psi,n,n,"j_phi",m)

  out(1:n,1:n) = j_phi_new(1:n,1:n)
  call radial_plot(out,psi,n,n,"j_phi_new",m)

  out(1:n,1:n) = j_par(1:n,1:n)
  call radial_plot(out,psi,n,n,"j_par",m)

  out(1:n,1:n) = j_x(1:n,1:n)
  call radial_plot(out,psi,n,n,"j_x",m)

  out(1:n,1:n) = j_z(1:n,1:n)
  call radial_plot(out,psi,n,n,"j_z",m)

  if(tri_type==11) then

	out(1:n,1:n) = rbtheta(1:n,1:n)
	call radial_plot(out,psi,n,n,"rbtheta",m)

    call get_r3btheta(rbtheta,br,bz,n)

	out(1:n,1:n) = rbtheta(1:n,1:n)
	call radial_plot(out,psi,n,n,"r3btheta",m)

  endif

  out(1:n,1:n) = cs(1:n,1:n)
  call radial_plot(out,psi,n,n,"cs",m)

  out(1:n,1:n) = csp(1:n,1:n)
  call radial_plot(out,psi,n,n,"csp",m)

!  out(1:n,1:n) = hyper(1:n,1:n)
!  call radial_plot(out,psi,n,n,"hyper",m)

	do ir = 1,n
	do iz = 1, n
		out(ir,iz) = p(ir,iz) + (br(ir,iz)**2 + bz(ir,iz)**2+b_phi(ir,iz)**2)/(2.d0*mu_mag)
	enddo
	enddo
  call radial_plot(out,psi,n,n,"ptot",m)

!	call radial_plot_gg(gg,br_gg,bz_gg,n,n,m)

  call geom(n,n)

!  if(Broot==3) then

!	  call bandpsi(psi,br,b_phi,bz,n,n)


	if(tri_type==11) call pvgamma(n,psi)


  call grid_output(n,n)

  open(11,file='psimax.dat')
  write(11,*) x_coord(mir), z_coord(miz)
  close(11)

	if((bc_type==7).or.(bc_type==8)) call plasma_boundary(psi,x_coord(mir), z_coord(miz))

	if(MARS_output) then
	
!		call write_MARS_input(psi,rho,x_coord(mir), z_coord(miz))

		call write_full_grid(p(1:n,1:n),psi,n,n,"pres",pedge)
		call write_full_grid(rho(1:n,1:n),psi,n,n,"rho",dedge)

	endif

!	call magnetic_inversion(psi,x_coord(mir), z_coord(miz))

  if(write_all) then
  ! this is differentiated because get_trapped will not run if write_all is false

	! WARNING: SOMETHING WENT WRONG SOMEWHERE IN MAGNETIC OUTPUT
	! GETTING ERRORS WITH FILEDS = 0
	! AND DBSINT RECEIVNG ALLEGEDLY NON-INCREASING INPUT
	! SKIP TRAPPED PARTICLES ETC.

	  call magnetic_output(psi,bpol,b_phi,rho,csp,x_coord(mir), z_coord(miz),pol_viscosity_data)

	  call get_trapped(n, n, psi, bpol, b_phi, trapped)
	  out(1:n,1:n) = trapped(1:n,1:n)
	  call radial_plot(out,psi,n,n,"trapped",m)


!!$	  if((static_equi).and.(bootstrap_option==0)) then
!!$
!!$		  call get_JparB(JparB_ave, bpol, b_phi, J_par, n)
!!$
!!$	  else
!!$
!!$		call bootstrap_setup_ave(psi,rho, p, Br, Bz, Temp, b_phi, J_par, bpol, n, -1.d0)
!!$
!!$	  endif



!!$	  call surf_ave(bpol, p(1:n,1:n), n, -1.d0, boot_pi)
!!$
!!$		outd(1:n,1:n) = 1.d0
!!$
!!$  	  call surf_ave(bpol, outd(1:n,1:n), n, -1.d0, boot_ni)
!!$  	  call surf_ave(bpol, temp(1:n,1:n), n, -1.d0, boot_Ti)

!!$	  call get_bootstrap2
!!$	  call get_bootstrap_NCLASS
!!$	  call bootstrap_output

  endif


	if(eq_type==1) then

	  call getqstar(psi,j_phi,p,bpol,b_phi,beta,rho,vr,vphi,vz,n)

	else

		do iz=1,n
		do ir=1,n
			outd(ir,iz) = (ppar(ir,iz) + 2.d0*pperp(ir,iz))/3.d0
		enddo
		enddo

		call getqstar(psi,j_phi,outd,bpol,b_phi,beta,rho,vr,vphi,vz,n)

	endif

	if(write_all) call bootstrap_cleanup
	! 7/7/2008 this moved here  to allow for I_BS calculation in getqstar


	call r_of_theta
        
!	call eqdskwrite3(n,n,psi,p)
!	call eqdskwrite4(n,n,psi,p)



!------------------------------------------------------------------------------
! deallocation section

!111 deallocate(psi) [commented out because I don't think this was intentional]
  deallocate(psi)
  deallocate(rho)
  deallocate(residual)
  deallocate(b_phi)
  deallocate(vr)
  deallocate(vphi)
  deallocate(vz)
  deallocate(br)
  deallocate(bz)
  deallocate(bpol)
  deallocate(p)
  deallocate(ppar)
  deallocate(pperp)
  deallocate(temp)
  deallocate(tpar)
  deallocate(tperp)
  deallocate(mtor)
  deallocate(malfven2)
  deallocate(mslow)
  deallocate(mcusp)
  deallocate(beta)
  deallocate(betapar)
  deallocate(betaperp)
  deallocate(j_phi)
  deallocate(j_par)
  deallocate(j_x)
  deallocate(j_z)
  if(allocated(rbtheta)) deallocate(rbtheta)
  deallocate(trapped)
  deallocate(out)
  deallocate(outd)

  !call allocation_clean_up

! end of deallocation section




  !Uncomment this if you want to pause at the end of FLOW's execution
  !pause

  continue

end subroutine run_FLOW


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine allocation_clean_up
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use solver
	use triangularity
	use magnetic

	implicit none

	if(allocated(ind_bound)) deallocate(ind_bound)
	if(allocated(coord_bound)) deallocate(coord_bound)

	if(allocated(x_coord)) deallocate(x_coord)
	if(allocated(z_coord)) deallocate(z_coord)

	if(allocated(dx_a)) deallocate(dx_a)
	if(allocated(dz_a)) deallocate(dz_a)

	if(allocated(sort_grid)) deallocate(sort_grid)

	if(allocated(d_cscoef)) deallocate(d_cscoef)
	if(allocated(p_iso_cscoef)) deallocate(p_iso_cscoef)
	if(allocated(ppar_cscoef)) deallocate(ppar_cscoef)
	if(allocated(pperp_cscoef)) deallocate(pperp_cscoef)
	if(allocated(F_cscoef)) deallocate(F_cscoef)
	if(allocated(mphi_cscoef)) deallocate(mphi_cscoef)
	if(allocated(mtheta_cscoef)) deallocate(mtheta_cscoef)

	if(allocated(dx_ext)) deallocate(dx_ext)
	if(allocated(dz_ext)) deallocate(dz_ext)

	if(allocated(q_coef)) deallocate(q_coef)

	if(allocated(r_data)) deallocate(r_data)
	if(allocated(r_cscoef)) deallocate(r_cscoef)

        if(allocated(bscoef_rho)) deallocate(bscoef_rho)
	if(allocated(bscoef_psi)) deallocate(bscoef_psi)
	if(allocated(bscoef_bpol)) deallocate(bscoef_bpol)
	if(allocated(bscoef_bphi)) deallocate(bscoef_bphi)
	if(allocated(cross_section)) deallocate(cross_section)

        ! Added by Ian
        if(allocated(w_ave_int)) deallocate(w_ave_int)
	if(allocated(t_ave_int)) deallocate(t_ave_int)
        if(allocated(br_2)) deallocate(br_2)
        if(allocated(bz_2)) deallocate(bz_2)
        if(allocated(br_3)) deallocate(br_3)
        if(allocated(bz_3)) deallocate(bz_3)
        if(allocated(mpol)) deallocate(mpol)
        if(allocated(j_phi_new)) deallocate(j_phi_new)
        if(allocated(cs)) deallocate(cs)
        if(allocated(csp)) deallocate(csp)
        if(allocated(hyper)) deallocate(hyper)
        if(allocated(br_gg)) deallocate(br_gg)
        if(allocated(bz_gg)) deallocate(bz_gg)
        if(allocated(gg)) deallocate(gg)

        if(allocated(psi_data)) deallocate(psi_data)
        if(allocated(temp_data)) deallocate(temp_data)
        if(allocated(rho_data)) deallocate(rho_data)
        if(allocated(mpol_data)) deallocate(mpol_data)
        if(allocated(vr_data)) deallocate(vr_data)
        if(allocated(vphi_data)) deallocate(vphi_data)
        if(allocated(vz_data)) deallocate(vz_data)
        if(allocated(vp_data)) deallocate(vp_data)
        if(allocated(br_data)) deallocate(br_data)
        if(allocated(bphi_data)) deallocate(bphi_data)
        if(allocated(bz_data)) deallocate(bz_data)
        if(allocated(bpol_data)) deallocate(bpol_data)
        if(allocated(bfield_data)) deallocate(bfield_data)
        if(allocated(phi_data)) deallocate(phi_data)

!	if(allocated()) deallocate()

end subroutine allocation_clean_up
end module FLOW_mod
