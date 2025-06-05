
program FLOW_mac

  use constant
  use solver
  use exp_data
  use triangularity
  use magnetic

  implicit none

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


  call readinput(m)

! allocation section
! this section is necessary because n is an input value and not a parameter

!  allocate(psi(0:n+1,0:n+1),stat = alloc_stat)
  allocate(psi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  psi(1:n,1:n) = 0.d0

  allocate(rho(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  rho(1:n,1:n) = 0.d0

  allocate(residual(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  residual(1:n,1:n) = 0.0

  allocate(b_phi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  b_phi(1:n,1:n) = 0.0

! end of allocation section

  if((tri_type==11).and.((eq3_opt==7).or.(eq3_opt==8).or.(eq3_opt==9)).or.  &
			((p_opt==7).or.(p_opt==8).or.(p_opt==9))) then
	psi_pres = psi_in
	psic = psi_pres
  endif

!  if(Broot==3) call psiout
   call psiout
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
     print *, "Allocation Error"
  end if
  vr(1:n,1:n) = 0.0

  allocate(vphi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  vphi(1:n,1:n) = 0.0

  allocate(vz(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  vz(1:n,1:n) = 0.0

  allocate(br(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  br(1:n,1:n) = 0.0

  allocate(bz(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  bz(1:n,1:n) = 0.0

  allocate(br_2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  br_2(1:n,1:n) = 0.0

  allocate(bz_2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  bz_2(1:n,1:n) = 0.0

  allocate(br_3(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  br_3(1:n,1:n) = 0.0

  allocate(bz_3(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  bz_3(1:n,1:n) = 0.0

  allocate(p(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  p(1:n,1:n) = 0.0

  allocate(ppar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  ppar(1:n,1:n) = 0.0

  allocate(pperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  pperp(1:n,1:n) = 0.0

  allocate(temp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  temp(1:n,1:n) = 0.0

  allocate(tpar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  tpar(1:n,1:n) = 0.0

  allocate(tperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  tperp(1:n,1:n) = 0.0

  allocate(mtor(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  mtor(1:n,1:n) = 0.0

  allocate(malfven2(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  malfven2(1:n,1:n) = 0.0

  allocate(mslow(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  mslow(1:n,1:n) = 0.0

  allocate(mcusp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  mcusp(1:n,1:n) = 0.0

  allocate(mpol(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  mpol(1:n,1:n) = 0.0

  allocate(beta(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  beta(1:n,1:n) = 0.0

  allocate(betapar(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  betapar(1:n,1:n) = 0.0

  allocate(betaperp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  betaperp(1:n,1:n) = 0.0

  allocate(j_phi(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  j_phi(1:n,1:n) = 0.0

  allocate(j_phi_new(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  j_phi_new(1:n,1:n) = 0.0

  allocate(j_par(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  j_par(1:n,1:n) = 0.0

  allocate(j_x(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  j_x(1:n,1:n) = 0.0

  allocate(j_z(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  j_z(1:n,1:n) = 0.0

  allocate(rbtheta(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  rbtheta(1:n,1:n) = 0.0

  allocate(trapped(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  trapped(1:n,1:n) = 0.0

  allocate(cs(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  cs(1:n,1:n) = 0.d0

  allocate(csp(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  csp(1:n,1:n) = 0.d0

  allocate(hyper(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  hyper(1:n,1:n) = 0.d0

  allocate(br_gg(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  br_gg(1:n,1:n) = 0.0

  allocate(bz_gg(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  bz_gg(1:n,1:n) = 0.0

  allocate(gg(1:n,1:n,1:3),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  gg(1:n,1:n,1:3) = 0.0


  allocate(out(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  out(1:n,1:n) = 0.0

  allocate(outd(1:n,1:n),stat = alloc_stat)
  if(alloc_stat > 0) then
     print *, "Allocation Error"
  end if
  outd(1:n,1:n) = 0.d0

! end of allocation section

  out(1:n,1:n) = psi(1:n,1:n)

  if(tri_type==9) then

	  call radial_plot_2(psi,n,n,"psi")

  else

	  call radial_plot(out,psi,n,n,"psi",m)

  endif

  call radial_plot(out,psi,n,n,"psi_clean",m)

!!$  call q_data(out,n,n,"psi")

  out(1:n,1:n) = rho(1:n,1:n)/mass
  call radial_plot(out,psi,n,n,"rho",m)

  out(1:n,1:n) = residual(1:n,1:n)
  call radial_plot(out,psi,n,n,"residual",m)

	call r_of_theta

	call write_restart_data(psi,n,n,"psi")
	call write_restart_data(rho,n,n,"rho")

	call psiout

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
  call radial_plot(out,psi,n,n,"vphi",m)



  out(1:n,1:n) = br(1:n,1:n)
  call radial_plot(out,psi,n,n,"br",m)

!!$  out(1:n,1:n) = br_2(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"br_2",m)

!!$  out(1:n,1:n) = br_3(1:n,1:n)
!!$  call radial_plot(out,psi,n,n,"br_3",m)

  out(1:n,1:n) = b_phi(1:n,1:n)
  call radial_plot(out,psi,n,n,"bphi",m)

!!$  call q_data(out,n,n,"bphi")

  out(1:n,1:n) = bz(1:n,1:n)
  call radial_plot(out,psi,n,n,"bz",m)

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
  call radial_plot(out,psi,n,n,"bpol",m)

!!$  call q_data(out,n,n,"bpol")

  ! Calculate b total
  do iz = 1, n
     do ir = 1, n
        out(ir,iz) = dsqrt(br(ir,iz)**2 + bz(ir,iz)**2+b_phi(ir,iz)**2)
     end do
  end do
  call radial_plot(out,psi,n,n,"b_field",m)

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
  call radial_plot(out,psi,n,n,"temp",m)

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
		  call radial_plot(out,psi,n,n,"vz",m)

		  out(1:n,1:n) = vr(1:n,1:n)
		  call radial_plot(out,psi,n,n,"vr",m)

		  ! Calculate vpoloidal
		  do iz = 1, n
			 do ir = 1, n
				out(ir,iz) = dsqrt(vr(ir,iz)**2 + vz(ir,iz)**2)
			 end do
		  end do
		  call radial_plot(out,psi,n,n,"vp",m)

		  out(1:n,1:n) = malfven2(1:n,1:n)
		  call radial_plot(out,psi,n,n,"malfven2",m)

	  out(1:n,1:n) = mpol(1:n,1:n)
	  call radial_plot(out,psi,n,n,"mpol",m)

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
		call write_MARS_input(psi,rho,x_coord(mir), z_coord(miz))
	endif

!	call magnetic_inversion(psi,x_coord(mir), z_coord(miz))

  if(write_all) then
  ! this is differentiated because get_trapped will not run if write_all is false

	! WARNING: SOMETHING WENT WRONG SOMEWHERE IN MAGNETIC OUTPUT
	! GETTING ERRORS WITH FILEDS = 0
	! AND DBSINT RECEIVNG ALLEGEDLY NON-INCREASING INPUT
	! SKIP TRAPPED PARTICLES ETC.

	  call magnetic_output(psi,bpol,b_phi,rho,csp,x_coord(mir), z_coord(miz))

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

111  deallocate(psi)
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

  call allocation_clean_up

! end of deallocation section




  pause

  continue

end program FLOW_mac

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

	if(allocated(bscoef_psi)) deallocate(bscoef_psi)
	if(allocated(bscoef_bpol)) deallocate(bscoef_bpol)
	if(allocated(bscoef_bphi)) deallocate(bscoef_bphi)
	if(allocated(cross_section)) deallocate(cross_section)

!	if(allocated()) deallocate()

end subroutine allocation_clean_up





