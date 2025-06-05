!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine radial_plot(psi,truepsi,nx,nz,fname,iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  use magnetic
  implicit none

  integer, intent(in) :: nx,nz,iz
  real (kind=skind), dimension(1:nx,1:nz), intent(in) :: psi
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: truepsi
  character (len=*), intent(in) :: fname
  integer :: ix
  real (kind=dkind) :: x,z,ex,ez
  integer i,j

	if((write_all).or.(fname=='psi_step').or.(fname=='psi')   &
			.or.(fname=='residual_step')) then

		continue

	else

		return

	endif



    open (unit=17,file=fname//'.plt')

    write(17,*)'TITLE="solution of GS equation with flow"'
    write(17,*)'Variables =" x ","y", "var"'		!, 	fname
!    write(17,*)'Variables =" R [m] ","z [m]", "var"'		!, 	fname
    write(17,*)'ZONE I=',nx,',J=',nz,',F=Point'


!  open (unit=44,file=fname//'.out',status='unknown', &
 !      form='formatted',action='write')


	do j=1,nz

		z = z_coord(j)

		do i=1,nx

			if(rmajor>1.d8) then
			! likely to be an astrophysical case, normalize to rmajor

				x = x_coord(i)/rmajor

			elseif(2.d0*a_elps/(rmajor+rcenter)<1.d-3) then
			! almost cylindrical case, reset x coordinate

				x = x_coord(i)-rcenter

			else
			!standard toroidal case

				x = x_coord(i)

			endif


!		    if ((sort_grid(i,j,0)==1).or.(fname=='psi').or.(fname=='psi_step')  &
		    if ( ((sort_grid(i,j,0)==2) .and.(tri_type==13)).or.  &
					((sort_grid(i,j,0)>=1) .and.(tri_type/=13)).or.(fname=='psi_step')  &
					.or.(fname=='residual_step').or.(fname=='rho_step').or.(fname=='psi').or.(fname=='psi_eqdsk1')  &
					.or.(fname=='psi_eqdsk2').or.(fname=='psi_eqdsk3')  &
					.or.(((bc_type==7).or.(bc_type==8)).and.(fname/='psi_clean'))) then

				if((truepsi(i,j)>=fraction*psic).or.(fname=='psi').or.(fname=='psi_step')  &
					.or.(fname=='residual_step').or.(fname=='rho_step').or.(fname=='psi_chease' ).or.(fname=='psi_eqdsk1')  &
					.or.(fname=='psi_eqdsk2').or.(fname=='psi_eqdsk3').or.(bc_type==7).or.(bc_type==8)) then
					write(17,88) x,z,psi(i,j)
	!				write(44,88) x,z,psi(i,j)
				else
					write(17,88) x,z,0.d0
				endif

			else
				write(17,88) x,z,0.d0
!				write(44,88) x,z,0.d0
			endif

		end do
	end do

	close(17)
 !   close(44)
88   format(3(e12.6,3x))



  open (unit=33,file=fname//'.txt',status='unknown', &
       form='formatted',action='write')


!  write (33,'(a)') "# radius                 value "
  do ix = 1, nx
	 x = x_coord(ix)
     write (33,'(e20.14,''     '',e20.14)') x, psi(ix,iz)
  end do

  close(33)


end subroutine radial_plot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine radial_plot_2(psi,nx,nz,fname)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi
  character (len=*), intent(in) :: fname
  integer :: ix
  real (kind=dkind) :: x,z,ex,ez
  integer i,j



    open (unit=17,file=fname//'.plt')

    write(17,*)'TITLE="solution of GS equation with flow"'
    write(17,*)'Variables =" R [m] ","z [m]", "var"'		!, 	fname
    write(17,*)'ZONE I=',nx,',J=',nz,',F=Point'


!  open (unit=44,file=fname//'.out',status='unknown', &
!       form='formatted',action='write')

	do j=1,nz

		z = z_coord(j)

		do i=1,nx

			x = x_coord(i)


		    if(sort_grid(i,j,0)==1) then

				write(17,88) x,z,psi(i,j)
!				write(44,88) x,z,psi(i,j)
			else
				write(17,88) x,z,0.d0
!				write(44,88) x,z,0.d0
			endif

		end do
	end do

	close(17)
!    close(44)

	continue

88   format(3(e26.17,3x))



end subroutine radial_plot_2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bootstrap_output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use magnetic

	implicit none

	integer :: i

	! -------------------output vs. psi-------------------

	open(69,file='bootstrap.plt',action='write')

	if(bootstrap_option==0) then

		write(69,*)'TITLE="Bootstrap current for solution of GS equation with flow"'
		write(69,123)'Variables ="psi ","J_BS_0", "J_BS", "f_BS", "<J_par B>", "J_BS_0_CHEASE", "J_BS_CHEASE"'
		write(69,*)'ZONE I=',enq,',F=Point'

!	elseif(bootstrap_option==1) then

	elseif(bootstrap_option==2) then

		write(69,*)'TITLE="Bootstrap current for solution of GS equation with flow"'
		write(69,123)'Variables ="psi ","J_BS_0", "J_BS", "f_BS", "<J_par B>", "J_BS_0_CHEASE",  "J_BS_CHEASE", "J_BS_NCLASS", "f_BS_NCLASS"'
		write(69,*)'ZONE I=',enq,',F=Point'

	endif

	if(bootstrap_option==0) then

		do i =1,enq

			write(69,1234) psival(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i), J_boot(4,i), J_boot(5,i)

		enddo

	! elseif(bootstrap_option==1) then
	elseif(bootstrap_option==2) then

		do i =1,enq

			write(69,2345) psival(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i),  &
								J_boot(4,i), J_boot(5,i), J_boot(6,i), J_boot(7,i)

		enddo

	endif

	close(69)

	! -------------------end of output vs. psi-------------------

	! -------------------output vs. R-------------------

	open(4, file='bootstrap_R.plt', action='write', status='unknown')

	if(bootstrap_option==0) then

		write(4,*)'TITLE="Bootstrap current for solution of GS equation with flow"'
		write(4,123)'Variables ="R ","J_BS_0", "J_BS", "f_BS", "<J_par B>", "J_BS_0_CHEASE",  "J_BS_CHEASE"'
		write(4,*)'ZONE I=',2*enq,',F=Point'

!	elseif(bootstrap_option==1) then

	elseif(bootstrap_option==2) then

		write(4,*)'TITLE="Bootstrap current for solution of GS equation with flow"'
		write(4,123)'Variables ="R ","J_BS_0", "J_BS", "f_BS", "<J_par B>", "J_BS_0_CHEASE",  "J_BS_CHEASE", "J_BS_NCLASS", "f_BS_NCLASS"'
		write(4,*)'ZONE I=',2*enq,',F=Point'

	endif


	if(bootstrap_option==0) then

		do i = enq,1,-1
			write(4,1234) Rleft(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i), J_boot(4,i), J_boot(5,i)
		enddo

		do i = 1,enq
			write(4,1234) Rright(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i), J_boot(4,i), J_boot(5,i)
		enddo

!	elseif(bootstrap_option==1) then

	elseif(bootstrap_option==2) then

		do i = enq,1,-1
			write(4,2345) Rleft(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i),  &
								J_boot(4,i), J_boot(5,i), J_boot(6,i), J_boot(7,i)
		enddo

		do i = 1,enq
			write(4,2345) Rright(i), J_boot(1,i), J_boot(2,i), J_boot(3,i), JparB_ave(i),  &
								J_boot(4,i), J_boot(5,i), J_boot(6,i), J_boot(7,i)
		enddo

	endif

	close(4)

	!----------------------- conductivity output -----------------------------
	! NOTE: the same output is written regardless of bootstrap_option:
	! depending on bootstrap_option, some entries may be 0

	! -------------------output vs. psi-------------------

	open(4, file='neoclass_res.plt', action='write', status='unknown')


	write(4,*)'TITLE="Neoclassical resistivity for solution of GS equation with flow"'
	write(4,123)'Variables ="psi ","sigma Spitzer", "sigma_Sauter", "sigma_NCLASS"'
	write(4,*)'ZONE I=',enq,',F=Point'


	do i = 1,enq
		write(4,1234) psival(i), el_resistivity(1,i), el_resistivity(2,i), el_resistivity(3,i)
	enddo

	close(4)

	! -------------------end of output vs. psi-------------------

	! -------------------output vs. R-------------------

	open(4, file='neoclass_res_R.plt', action='write', status='unknown')


	write(4,*)'TITLE="Neoclassical resistivity for solution of GS equation with flow"'
	write(4,123)'Variables ="R ","sigma Spitzer", "sigma_Sauter", "sigma_NCLASS"'
	write(4,*)'ZONE I=',2*enq,',F=Point'


	do i = enq,1,-1
		write(4,1234) Rleft(i), el_resistivity(1,i), el_resistivity(2,i), el_resistivity(3,i)
	enddo

	do i = 1,enq
		write(4,1234) Rright(i), el_resistivity(1,i), el_resistivity(2,i), el_resistivity(3,i)
	enddo

	close(4)


	continue

123 format(A130)
1234 format(7(e16.10, 3x))
2345 format(9(e16.10, 3x))

end subroutine bootstrap_output

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine short_output(j_phi,p,ppar,pperp,b_phi,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
!  use constant, only : eq_type	!,a_elps,b_elps
  use magnetic, only : mu_mag,b_phi_zero
  use triangularity

  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: j_phi,p,ppar,pperp,b_phi
  real (kind=dkind) :: x,ex,ez,curr,b0,p_av,betator, ptot
  integer i,j,k,ncenter
  real (kind=dkind) :: counter,vol
  logical :: b_in = .false.


  area = dx*dz
  curr = 0.d0
  p_av = 0.d0
  ptot = 0.d0
  b0 = 0.d0
  counter = 0.
  vol = 0.

  ncenter = (nx-1)/2+1
  j = 0

!  do while (NOT(b_in))

	j = j+1

!	if ((b_phi(ncenter,j)>1.d-10)) then
!		b0 = b_phi(ncenter,j)
!		b_in=.true.
!	endif

	if ((b_phi(ncenter,j)>1.d-10)) then
		b_in=.true.
	endif
	if(b_in) b0 = b_phi(ncenter,j+1)

!  enddo

  if (dabs((b0-b_phi_zero)/(b0+b_phi_zero))>1.d-2) b0 = b_phi_zero


	do i=1,nx

!	if ((b_phi(i,(nz-1)/2+1)>1.d-10).AND.NOT(b_in)) then
!		b0 = b_phi(i,(nz-1)/2+1)
!		b_in=.true.
!	endif
!	if ((b_phi(i,(nz-1)/2+1)<1.d-10).AND.(b_in)) then
!		b0 = b0 + b_phi(i-1,(nz-1)/2+1)
!		b_in=.false.
!	endif

		x = 0.5d0*(x_coord(i)+x_coord(i+1))

		if(grid_type/=0) area = dx_a(i)*dz_a(j)

		do j=1,nz


!!$          ex = ((i-1)*dx - 0.5d0*x_size)
!!$          ez = ((j-1)*dz - 0.5d0*z_size)
!!$		  theta = datan2(ez,ex)
!!$		  rminor = 0.d0
!!$
!!$		  if(tri_type==1) then
!!$
!!$			if((dsin(theta)>=0.d0)) then
!!$				do k=0,n_tri
!!$					rminor = rminor + rcoeff_u(k)*dcos(k*theta)
!!$				enddo
!!$			else
!!$          		do k=0,n_tri
!!$					rminor = rminor + rcoeff_d(k)*dcos(k*theta)
!!$				enddo	
!!$			endif
!!$
!!$		  elseif(tri_type==2) then
!!$
!!$			if((dsin(theta)>=0.d0)) then
!!$				rminor = smallr0*dsqrt( ellipt**2*dsin(theta)**2  &
!!$						+dcos(theta+asin_d_up*dsin(theta))**2 )
!!$			else
!!$          		rminor = smallr0*dsqrt( ellipt**2*dsin(theta)**2  &
!!$						+dcos(theta+asin_d_down*dsin(theta))**2 )
!!$			endif
!!$
!!$		  elseif(tri_type==3) then
!!$
!!$			if((dsin(theta)>=0.d0)) then
!!$				rminor = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
!!$						+dcos(theta)**2-delta_u_o4*(dcos(3*theta)-dcos(theta)) )
!!$			else
!!$          		rminor = smallr0/dsqrt( dsin(theta)**2/ellipt**2  &
!!$						+dcos(theta)**2-delta_d_o4*(dcos(3*theta)-dcos(theta)) )
!!$			endif
!!$
!!$		  endif
			
!		  if((ex*ex + ez*ez) <= 1.0d0) then
!		  if((ex*ex + ez*ez) <= rminor**2) then


			if(sort_grid(i,j,0)<=0) then
!			if(sort_grid(i,j,0)==-1) then
				 cycle
			end if

			counter = counter+1.d0
			vol = vol + area*x

			curr = curr + j_phi(i,j)*area	!dabs(j_phi(i,j))*area

			if (eq_type == 1) then
				p_av = p_av + p(i,j)*area*x
			else
				p_av = p_av + (ppar(i,j)+2.d0*pperp(i,j))*area*x/3.d0
			endif

			continue

!		  endif

	enddo
	enddo

	ptot = p_av

!	p_av = p_av/counter
	p_av = p_av/vol

!	if(b_in) 	b0 = b0 + b_phi(nx,(nz-1)/2+1)
		
!	b0 = b0/2.d0
	betator = 2.d0*mu_mag*p_av/b0**2

    open (unit=17,file='output.dat')
	write(17,88) curr
	write(17,88) betator
	write(17,*) 'total pressure', ptot
	close(17)

	continue

88 format (e20.14)

end subroutine short_output

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine bandpsi(psi,br,b_phi,bz,nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: psi
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: br,b_phi,bz
  integer :: ix
  real (kind=dkind) :: x,z,ex,ez,b
  integer i,j



    open (unit=17,file='bofpsi.dat')

	do i=1,nx

		x = x_coord(i)

		do j=1,nz

!			z = -0.5d0*z_size + (j-1)*dz
			z = z_coord(j)

		    if(sort_grid(i,j,0)==1) then

				b=dsqrt(br(i,j)**2+b_phi(i,j)**2+bz(i,j)**2)

				write(17,88) x,z,psi(i,j),b

			endif

		end do
	end do

	close(17)


88   format(4(e20.12,3x))


end subroutine bandpsi

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine q_data(var,nx,nz,fname)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=skind), dimension(1:nx,1:nz), intent(in) :: var
  real (kind=skind), dimension(1:nx,1:nz) :: out
  character (len=*), intent(in) :: fname
  real (kind=dkind) :: x,z,ex,ez
  integer i,j


    

	do i=1,nx
    
		x = x_coord(i)

		do j=1,nz

!			z = -0.5d0*z_size + (j-1)*dz
			z = z_coord(j)

		    if(sort_grid(i,j,0)==1) then
				out(i,j) = var(i,j)
			else
				out(i,j) = 0.d0
			endif

		end do
	end do


    open (unit=17,file=fname//'only.dat')

	do j=1,nz

		if(nx==33) then
			write(17,33) out(j,1:nx)
		elseif(nx==65) then
			write(17,65) out(j,1:nx)
		elseif(nx==129) then
			write(17,129) out(j,1:nx)
		elseif(nx==257) then
			write(17,257) out(j,1:nx)
		elseif(nx==513) then
			write(17,513) out(j,1:nx)
		elseif(nx==1025) then
			write(17,1025) out(j,1:nx)
		elseif(nx==2051) then
			write(17,2051) out(j,1:nx)
		endif

	enddo

	close(17)

	continue

33   format(33(f14.8, 3x))
65   format(65(f14.8, 3x))
129   format(129(f14.8, 3x))
257   format(257(f14.8, 3x))
513   format(513(f14.8, 3x))
1025   format(1025(f14.8, 3x))
2051   format(2051(f14.8, 3x))


	
	return
 
end subroutine q_data

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine geom(nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind) :: x,z


 

  open(22,file='geom.dat', status='unknown', action='write')

!  x = rmajor - 0.5d0*x_size 
	x = x_coord(1)

  write(22,*) x

  x = x_coord(nx)

  write(22,*) x

!  z = -0.5d0*z_size
  z = z_coord(1)

  write(22,*) z

!  z = -0.5d0*z_size + (nz-1)*dz
  z = z_coord(nz)

  write(22,*) z

  close(22)

  continue

  return

end subroutine geom

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine grid_output(nx,nz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind) :: x,z,ex,ez
  integer i,j


    open (unit=69,file='FLOW_n.dat')

	write(69,*) nx
	write(69,*) nz

	close(69)


    open (unit=619,file='FLOW_xgrid.dat')

	do i=1,nx

		write(619,88) x_coord(i)

	end do

	close(619)


    open (unit=619,file='FLOW_zgrid.dat')

	do j=1,nz

!		z = -0.5d0*z_size + (j-1)*dz
		z = z_coord(j)
		write(619,88) z

	end do

	close(619)


88   format(e26.17,3x)



end subroutine grid_output


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  function gtheta_bound_dump(th,smallr) result(answer)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! this function computes grad psi on the inner boundary
  ! for the LDX case

	use constant
	use magnetic
	use exp_data
	use pseudo_IMSL

	real(kind=dkind), intent(in) :: smallr, th
	real(kind=dkind) :: answer
	real(kind=dkind) :: bigR, z, K_ell, E_ell, k2, k, R02
	real(kind=dkind) :: elle, ellf

	R02 = rmajor**2
	bigR = rmajor + smallr*cos(th)
	z = smallr*sin(th)

	k2 = 4.d0*rmajor*bigR/((rmajor+bigR)**2+z**2)
	k = 2.d0*sqrt(rmajor*bigR/((rmajor+bigR)**2+z**2))

	K_ell = ellf(0.5d0*pi,k)
	E_ell = elle(0.5d0*pi,k)

!	answer = mu_mag/(pi*k)*I_coil*sqrt(rmajor*bigR)*  &
!			((1.d0-0.5*k2)*K_ell-E_ell)

	answer = mu_mag*I_coil/(2.d0*pi)/sqrt((bigR+rmajor)**2+z**2) *  &
			( (bigR*cos(th)+z*sin(th))*K_ell + (bigR*cos(th)-z*sin(th))*  &
				(R02+bigR**2+z**2)/((bigR-rmajor)**2+z**2)*E_ell )


	R0 = rmajor + delta_coil
	R02 = R0**2

	bigR = rmajor + smallr*cos(th)
	z = smallr*sin(th)

	k2 = 4.d0*R0*bigR/((R0+bigR)**2+z**2)
	k = 2.d0*sqrt(R0*bigR/((R0+bigR)**2+z**2))

	K_ell = ellf(0.5d0*pi,k)
	E_ell = elle(0.5d0*pi,k)

!	answer = mu_mag/(pi*k)*I_coil*sqrt(rmajor*bigR)*  &
!			((1.d0-0.5*k2)*K_ell-E_ell)

	answer = mu_mag*I_coil/(2.d0*pi)/sqrt((bigR+R0)**2+z**2) *  &
			( (bigR*cos(th)+z*sin(th))*K_ell + (bigR*cos(th)-z*sin(th))*  &
				(R02+bigR**2+z**2)/((bigR-R0)**2+z**2)*E_ell )

	answer = dbsval(th, psiprim_ord, psiprim_data(1:ibreak_pp,3),  &
	psiprim_dim, psiprim_cscoef(1,1:psiprim_dim) )


	continue

  end function gtheta_bound_dump

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_rbtheta(rbtheta,br,bz,nn)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use solver

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn) :: rbtheta, br, bz
	integer :: i, j, zone
	real(kind=dkind) :: x, z, r, th
	real(kind=dkind) :: alph, rprim, rsec, ex, ez

	do i=1,nn

		x = x_coord(i)

		do j=1,nn

			z = z_coord(j)

		    if (sort_grid(i,j,0)==1) then

				call radius_1_3(x,z,ex,ez,th,r,zone,alph,rprim,rsec)

				r = sqrt(ex**2+ez**2)

				rbtheta(i,j) = r*abs( bz(i,j)*cos(th) - br(i,j)*sin(th) )

			else

				rbtheta(i,j) = 0.d0

			endif

		enddo

	enddo

	continue

	return


end subroutine get_rbtheta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine get_r3btheta(rbtheta,br,bz,nn)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use solver

	implicit none

	integer :: nn
	real(kind=dkind), dimension(1:nn,1:nn) :: rbtheta, br, bz
	integer :: i, j, zone
	real(kind=dkind) :: x, z, r, th
	real(kind=dkind) :: alph, rprim, rsec, ex, ez

	do i=1,nn

		x = x_coord(i)

		do j=1,nn

			z = z_coord(j)

		    if (sort_grid(i,j,0)==1) then

				call radius_1_3(x,z,ex,ez,th,r,zone,alph,rprim,rsec)

				r = (sqrt(ex**2+ez**2))**3

				rbtheta(i,j) = r*abs( bz(i,j)*cos(th) - br(i,j)*sin(th) )

			else

				rbtheta(i,j) = 0.d0

			endif

		enddo

	enddo

	continue

	return


end subroutine get_r3btheta


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine r_of_theta
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  implicit none

  real (kind=dkind) :: angle, smallr, xdum, zdum
  integer :: i
  integer :: imax

	imax = r_of_theta_points

	open(33, file='r_of_theta.dat', status='unknown', action = 'write')

	do i=0,imax

		angle = 2.d0*pi*i/imax

		call radius_theta(angle,smallr,xdum,zdum)

!		write(33,288) angle, smallr
		write(33,289) angle, smallr

	enddo

	close(33)

	if(tri_type==8) then

		open(33, file='rprim_of_theta.dat', status='unknown', action = 'write')

		do i=0,imax

			angle = 2.d0*pi*i/imax

			call radius_prim_theta(angle,smallr)

			write(33,289) angle, smallr

		enddo

		close(33)

	endif

288   format(2(e26.17,3x))
289	  format(2(f21.18, 3x))

	continue

	return

end subroutine r_of_theta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine write_restart_data(thing,nx,nz,fname)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  use magnetic
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: thing
  character (len=*), intent(in) :: fname
  integer i,j

    open (unit=17, file=fname//'bin.out', form='unformatted', status='unknown', action='write')

	do 33 j=1,nz
	do 33 i=1,nx

	33 write(17) thing(i,j)

	close(17)

	continue

	return

end subroutine write_restart_data


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_restart_data(thing,nx,nz,fname)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  use magnetic
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz) :: thing
  character (len=*), intent(in) :: fname
  integer i,j

    open (unit=17, file=fname//'bin.out', form='unformatted', status='old', action='read')

	do 33 j=1,nz
	do 33 i=1,nx

	33 read(17) thing(i,j)

	close(17)

	continue

	return

end subroutine read_restart_data

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine eqdskwrite3(nx,nz,psi,p)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine gives a first bare-bones version
	! for an interface between FLOW and other codes (e.g. NOVA, PEST, MARS, CHEASE...)

		use constant
		use pseudo_IMSL, only : dbsval, dbsval_safe
		use magnetic
		use solver, only : bzero, pofpsi, dpdpsi, dbzerodpsi, radius_theta

		implicit none

        integer nx,nz,i,j,iomap,iodat,imax,jmax,jmid
!        real(kind=dkind), dimension(1:nz,1:nx):: psi_eq
        real(kind=dkind), dimension(1:nx,1:nz), intent(in):: psi,p
        real(kind=dkind) :: psimax,psi_max2,r0_eq,x0_eq,xma,zma
        real(kind=dkind), dimension(1:nx):: p_eq,g_eq
        real(kind=dkind) :: dummy
        real(kind=dkind) :: ex,ez,r,rloc2,scal,psiloc, angle
		character  cdate*8
		integer :: idummy
		real(kind=dkind) :: RC0P = 0.D0 ! no clue of what this is
		integer, parameter :: ilimiter = 5
		real(kind=dkind), dimension(ilimiter) :: zr, zz
		real(kind=dkind), dimension(r_of_theta_points) :: redge, zedge
		real(kind=dkind), dimension(nz) :: dummy_array

		! check b_phi_zero

		psimax = psi(mir,miz)


        iomap = 22
        iodat = 44

!!$        open(iomap,file='eqdsk',form='unformatted')
        open(iomap,file='eqdsk_FLOW.out',form='formatted')


!        FIRST LINE: 48 CHAR, DUMMY INTEGER, NRBOX, NZBOX

         cdate='20050802'
!   CALL DATE(ZDATE)
         idummy = 3
!CJEM         write(iomap,9380) 'FROM CHEASE, ALL IN MKSA UNITS     ',
         write(iomap,9380) ' FLOW    00/00/00      #000000,00000ms', &
			       cdate,idummy,nx,nz

!        2ND LINE: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZBOXMID

         write(iomap,9381) x_size,z_size,rmajor,x_coord(1) ,z_coord(miz)!*rmajor
		 ! in FLOW the grid is in physical units



!C        3RD LINE: RMAG, ZMAG, PSIMAG, PSIEDGE, B0
!C        ZMAG HAS BEEN SHifTED TO ZERO IN PSIBOX
!CYQL     KEEP ZMAG UNSHifTED

!         ZMAG0 = RZMAG
 
         write(iomap,9381) x_coord(mir),z_coord(miz),  psimax,RC0P,bzero(0.d0)

!C        4TH LINE: PLASMA CURRENT, PSIAX1, PSIAX2, RAXIS1, RAXIS2

!         if (R0EXP.EQ.1 .AND. B0EXP.EQ.1) ZMU0 = 1.0

         write(iomap,9381) curr,  &	!? *rmajor*b_phi_zero/mu_mag  what kind of units is the current in? is it total or toroidal?
            psimax,RC0P,x_coord(mir),RC0P

!C        5TH LINE: ZAXIS1, ZAXIS2, PSI_SEP, R_XPOINT, Z_XPOINT

         write(iomap,9381) z_coord(miz),RC0P,RC0P,RC0P,RC0P

!!$!C     EQUISTANT PSI-MESH FOR PROFILES IN S=SQRT(1.-PSI/PSIMIN)
!!$
!!$         ZDPSI = ABS(SPSIM) / FLOAT(NRBOX-1)
!!$         DO I=1,NRBOX
!!$           ZSTEMP(I) = SPSIM + FLOAT(I-1)*ZDPSI
!!$         endDO
!!$
!!$!C        6TH ENTRY: T(PSI) (OR G)
!!$
!!$         CALL SPLINE(CSM,TMF,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,TMF,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = T0
!!$         ZCOF = R0EXP*B0EXP
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = rmajor*bzero(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!!$!C        7TH ENTRY: PRESSURE
!!$
!!$         CALL SPLINE(CSM,CPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPR,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = CP0
!!$         ZCOF = B0EXP*B0EXP / ZMU0
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nz-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = pofpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!!$!C        8TH ENTRY: TT' (OR GG')
!!$
!!$         CALL SPLINE(CSM,TTP,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,TTP,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = DTTP0
!!$         ZCOF = B0EXP
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = rmajor**2*bzero(psiloc)*dbzerodpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!C        9TH ENTRY: P'

!!$         CALL SPLINE(CSM,CPPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPPR,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = DPDP0
!!$         ZCOF = B0EXP / ZMU0 / R0EXP / R0EXP
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = dpdpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!C        10TH ENTRY: PSI(I,J)

!!$         ZCOF = B0EXP * R0EXP**2
!!$         write(iomap,9381) ((EQDSPSI(I,J)*ZCOF,I=1,NRBOX),  &
!!$            J=1,NZBOX)

		write(iomap,9381) ((psi(i,j), i=1,nx), j=1,nz)	! check: psi!?

!C        11TH ENTRY: Q PROFILE

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = dbsval_safe(psiloc,q_ord, q_coef(3,:),  &
					ibreak_q-q_ord, q_coef(4,1:ibreak_q-q_ord), q_coef(2,1:ibreak_q-q_ord) )	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!		Bmax = dbsval(psiloc,modB_ord, modB_coef(3,:),  &
!					enq, modB_coef(4,1:enq) )
!!$         CALL SPLINE(CSM,QPSI,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,QPSI,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = Q0
!!$         write(iomap,9381) (ZTEMP(I),I=1,NRBOX)

!C        12TH ENTRY: (R,Z) OF PLASMA BOUNDARY AND DUMMY LIMITER POSITION

!C     LIMITER ARRAY
!CJEM         ILIMITER = 4
         
         zr(1) = x_coord(1)
         zz(1) = z_coord(1)
         zr(2) = x_coord(nx/2)
         zz(2) = z_coord(1)
         zr(3) = x_coord(nx/2)
         zz(3) =  z_coord(nz)
         zr(4) = x_coord(nx)
         zz(4) =  z_coord(nz)
         zr(5) = x_coord(nx)
         zz(5) = z_coord(1)

         write(iomap,1204) r_of_theta_points, ilimiter
 !!$        ZCOF = R0EXP

		do i = 1,r_of_theta_points

			angle = 2.d0*pi*i/r_of_theta_points
			call radius_theta(angle,dummy,redge(i),zedge(i))

		enddo

         write(iomap,9381) (redge(i),zedge(i), i=1,r_of_theta_points)
         write(iomap,9381) (zr(i),zz(i),i=1,ilimiter)

!C        LAST LINES: SOME EQUILIBRIUM VALUES
!C        (NO HEADER)

         CALL OUTMKSA(iomap,1,psimax)

         close(unit=iomap,status='keep')



1204    format(5I5)
9380    format(A40,A8,3I4)
9381    format(1P,5E16.9)

    end subroutine eqdskwrite3


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    subroutine eqdskwrite4(nx,nz,psi,p)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	! this subroutine gives a first bare-bones version
	! for an interface between FLOW and other codes (e.g. NOVA, PEST, MARS, CHEASE...)

		use constant
		use pseudo_IMSL, only : dbsval, dbsval_safe
		use magnetic
		use solver, only : bzero, pofpsi, dpdpsi, dbzerodpsi, radius_theta

		implicit none

        integer nx,nz,i,j,iomap,iodat,imax,jmax,jmid
!        real(kind=dkind), dimension(1:nz,1:nx):: psi_eq
        real(kind=dkind), dimension(1:nx,1:nz), intent(in):: psi,p
        real(kind=dkind) :: psimax,psi_max2,r0_eq,x0_eq,xma,zma
        real(kind=dkind), dimension(1:nx):: p_eq,g_eq
        real(kind=dkind) :: dummy
        real(kind=dkind) :: ex,ez,r,rloc2,scal,psiloc, angle
		character  cdate*8
		integer :: idummy
		real(kind=dkind) :: RC0P = 0.D0 ! no clue of what this is
		integer, parameter :: ilimiter = 5
		real(kind=dkind), dimension(ilimiter) :: zr, zz
		real(kind=dkind), dimension(r_of_theta_points) :: redge, zedge
		real(kind=dkind), dimension(nz) :: dummy_array

		! check b_phi_zero

		psimax = psi(mir,miz)

        iomap = 22
        iodat = 44

!!$        open(iomap,file='eqdsk',form='unformatted')
        open(iomap,file='eqdsk_FLOW_for_input.out',form='formatted')


!        FIRST LINE: 48 CHAR, DUMMY INTEGER, NRBOX, NZBOX

         cdate='20050802'
!   CALL DATE(ZDATE)
         idummy = 3
!CJEM         write(iomap,9380) 'FROM CHEASE, ALL IN MKSA UNITS     ',
         write(iomap,9380) ' FLOW    00/00/00      #000000,00000ms', &
			       cdate,idummy,nx,nz

!        2ND LINE: RBOXLEN, ZBOXLEN, R0, RBOXLFT, ZBOXMID

         write(iomap,9381) x_size,z_size,rmajor,x_coord(1)
		 ! in FLOW the grid is in physical units



!C        3RD LINE: RMAG, ZMAG, PSIMAG, PSIEDGE, B0
!C        ZMAG HAS BEEN SHifTED TO ZERO IN PSIBOX
!CYQL     KEEP ZMAG UNSHifTED

!         ZMAG0 = RZMAG
 
         write(iomap,9381) x_coord(mir),z_coord(miz),  psimax,RC0P,bzero(0.d0)

!C        4TH LINE: PLASMA CURRENT, PSIAX1, PSIAX2, RAXIS1, RAXIS2

!         if (R0EXP.EQ.1 .AND. B0EXP.EQ.1) ZMU0 = 1.0

         write(iomap,9381) -curr

!C        5TH LINE: ZAXIS1, ZAXIS2, PSI_SEP, R_XPOINT, Z_XPOINT

         write(iomap,9381) z_coord(miz)


		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = rmajor*bzero(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!!$!C        7TH ENTRY: PRESSURE
!!$
!!$         CALL SPLINE(CSM,CPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPR,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = CP0
!!$         ZCOF = B0EXP*B0EXP / ZMU0
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nz-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = pofpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!!$!C        8TH ENTRY: TT' (OR GG')
!!$
!!$         CALL SPLINE(CSM,TTP,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,TTP,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = DTTP0
!!$         ZCOF = B0EXP
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = rmajor**2*bzero(psiloc)*dbzerodpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!C        9TH ENTRY: P'

!!$         CALL SPLINE(CSM,CPPR,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,CPPR,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = DPDP0
!!$         ZCOF = B0EXP / ZMU0 / R0EXP / R0EXP
!!$         write(iomap,9381) (ZTEMP(I)*ZCOF,I=1,NRBOX)

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = dpdpsi(psiloc)	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!C        10TH ENTRY: PSI(I,J)

!!$         ZCOF = B0EXP * R0EXP**2
!!$         write(iomap,9381) ((EQDSPSI(I,J)*ZCOF,I=1,NRBOX),  &
!!$            J=1,NZBOX)

		write(iomap,9381) ((psi(i,j), i=1,nx), j=1,nz)	! check: psi!?

!C        11TH ENTRY: Q PROFILE

		do i=1,nx
			psiloc = psimax/(nx-1.d0)*(i-1.d0)
			dummy_array(nx+1-i) = dbsval_safe(psiloc,q_ord, q_coef(3,:),  &
					ibreak_q-q_ord, q_coef(4,1:ibreak_q-q_ord), q_coef(2,1:ibreak_q-q_ord) )	! check: psiloc!?
		enddo

		write(iomap,9381) (dummy_array(i), i=1,nx)

!		Bmax = dbsval(psiloc,modB_ord, modB_coef(3,:),  &
!					enq, modB_coef(4,1:enq) )
!!$         CALL SPLINE(CSM,QPSI,NPSI+1,ZD2TMP,ZWORK,ZWORK1)
!!$         CALL PPSPLN(NRBOX,ZSTEMP,NPSI,CSM,QPSI,ZD2TMP,ZTEMP)
!!$         ZTEMP(1) = Q0
!!$         write(iomap,9381) (ZTEMP(I),I=1,NRBOX)

!C        12TH ENTRY: (R,Z) OF PLASMA BOUNDARY AND DUMMY LIMITER POSITION

!C     LIMITER ARRAY
!CJEM         ILIMITER = 4

         
         zr(1) = x_coord(1)
         zz(1) = z_coord(1)
         zr(2) = x_coord(nx/2)
         zz(2) = z_coord(1)
         zr(3) = x_coord(nx/2)
         zz(3) =  z_coord(nz)
         zr(4) = x_coord(nx)
         zz(4) =  z_coord(nz)
         zr(5) = x_coord(nx)
         zz(5) = z_coord(1)

         write(iomap,1204) r_of_theta_points !, ilimiter
 !!$        ZCOF = R0EXP

		do i = 1,r_of_theta_points

			angle = 2.d0*pi*i/r_of_theta_points
			call radius_theta(angle,dummy,redge(i),zedge(i))

		enddo

         write(iomap,9381) (redge(i),zedge(i), i=1,r_of_theta_points)
!         write(iomap,9381) (zr(i),zz(i),i=1,ilimiter)

!C        LAST LINES: SOME EQUILIBRIUM VALUES
!C        (NO HEADER)

         CALL OUTMKSA(iomap,1,psimax)

         close(unit=iomap,status='keep')



1204    format(5I5)
9380    format(A40,A8,3I4)
9381    format(1P,5E16.9)

    end subroutine eqdskwrite4

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     SUBROUTINE OUTMKSA(iunit,kopt,psimax)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	 ! adapted from CHEASE
!                                        AUTHORS:
!                                        O.SAUTER,  CRPP-EPFL

!        KUNIT: DISK UNIT ON WHICH TO write
!        KOPT = 1: NO HEADER
!        KOPT = 2: WITH HEADER

	use constant
	use magnetic
	use triangularity
                      
	implicit none

	integer :: iunit, kopt
	real(kind=dkind) :: psimax




	 write(iunit,*)

	 if (kopt .eq. 2) then
	   write(iunit,9100)
	 endif

	write(iunit,9101) rmajor,' R0 [M]'
	write(iunit,9101) b_phi_zero,' B0 [T]'

	 write(iunit,9102) x_coord(mir)/rmajor,' R of magnetic axis --> [M]   ', x_coord(mir)
	 write(iunit,9102) z_coord(miz)/rmajor,' Z of magnetic axis --> [M]   ', z_coord(miz)

	 write(iunit,9102) psimax/rmajor**2/b_phi_zero, ' psi-axis --> [T M**2] ', psimax

	 write(iunit,9102) 2.d0*pi*psimax/rmajor**2/b_phi_zero,' 2*pi*psi-axis -->     ',  &
					2.d0*pi*psimax
	 write(iunit,9102) curr*mu_mag/rmajor/b_phi_zero, ' toroidal current --> [A] ', curr

	 write(iunit,9101) betapol ,' poloidal beta'
	 write(iunit,9101) betator,' toroidal beta' !not sure of what this is in CHEASE
	 write(iunit,9101) betastar,' beta* (sqrt(<P**2>))'
	 write(iunit,9101) betator,' beta_exp=<P>*2*mu0/B0**2'
	 write(iunit,9101) 0.d0,' li (not computed in FLOW)'
	 write(iunit,9101) q_c,' q_0'
	 write(iunit,9101) qe,' q_edge'

	 write(iunit,9102) 1.d0/inv_aspect_ratio,' aspect ratio ; a/R= ', inv_aspect_ratio
	 write(iunit,9101) shape_ellipticity,' b/a '
	 write(iunit,9102) volume/rmajor**3, ' volume -> ',volume
	 write(iunit,9102) area/rmajor**2,   ' area   -> ', area
	 write(iunit,9102) 1.d0, ' LENGTH -> not implemented, yet',1.d0
!	 write(iunit,9102)RLENG(NPSI1), ' LENGTH -> ',RLENG(NPSI1)*R0EXP

	 write(iunit,9102) R_edge_min/rmajor,' R_min -> R_min [m] ', R_edge_min
	 write(iunit,9102) R_edge_max/rmajor,' R_max -> R_max [m] ', R_edge_max
	 write(iunit,9102) z_edge_min/rmajor,' z_min -> z_min [m] ', z_edge_min
	 write(iunit,9102) z_edge_max/rmajor,' ZMAX -> ZMAX [m] ', z_edge_max
	 write(iunit,9102) 0.5d0*(R_edge_min+R_edge_max)/rmajor,' R_geom -> R_geom [m] ',  &
								0.5d0*(R_edge_min+R_edge_max)
	 write(iunit,9102) 0.5d0*(R_edge_max-R_edge_min)/rmajor,' minor radius -> a [m] ',  &
								0.5d0*(R_edge_max-R_edge_min)/rmajor

         return
 9100    format(/,1X,'*************************************',  &
              //,1X,'SOME QUANTITIES AND THEIR MKSA VALUES',  &
              //,1X,'*************************************',/)
 9101    format(1PE18.8,A)
 9102    format(1PE18.8,A,E18.8)
 9103    format(1PE18.8,2A)

end subroutine outmksa


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine psi_boundary_plot(nx,nz,psi)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant, only : dkind, x_coord, z_coord

	implicit none

	integer :: nx, nz
	real(kind=dkind) :: psi(1:nx,1:nz)
	integer :: imin, imax, jmin, jmax
	integer :: i, j

	imin = 3
	imax = nx-2
	jmin = 5
	jmax = nz-7

	open(69, file = 'psi_boundary.txt', action = 'write')

	! left side
	i = imax
	do j = jmin, jmax

		write(69,*) x_coord(i), z_coord(j), psi(i,j)

	enddo

	! top
	j = jmax
	do i = imax-1, imin, -1

		write(69,*) x_coord(i), z_coord(j), psi(i,j)

	enddo

	! right side
	i = imin
	do j = jmax-1, jmin, -1

		write(69,*) x_coord(i), z_coord(j), psi(i,j)

	enddo

	! bottom
	j = jmin
	do i = imin+1, imax-1

		write(69,*) x_coord(i), z_coord(j), psi(i,j)

	enddo

	close(69)

	continue

	return

end subroutine psi_boundary_plot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine radial_plot_gg(gg,br_gg,bz_gg,nx,nz,iz)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  use magnetic
  implicit none

  integer, intent(in) :: nx,nz,iz
  real (kind=dkind), dimension(1:nx,1:nz,1:3), intent(in) :: gg
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: br_gg, bz_gg
  integer :: ix
  real (kind=dkind) :: x,z,ex,ez
  integer i,j


    open (unit=61,file='gg_r.plt')
    open (unit=62,file='gg_z.plt')
    open (unit=63,file='div_gg.plt')
    open (unit=64,file='br_gg.plt')
    open (unit=65,file='bz_gg.plt')

  open (unit=31,file='gg_r.txt',status='unknown', &
       form='formatted',action='write')
  open (unit=32,file='gg_z.txt',status='unknown', &
       form='formatted',action='write')
  open (unit=33,file='div_gg.txt',status='unknown', &
       form='formatted',action='write')
  open (unit=34,file='br_gg.txt',status='unknown', &
       form='formatted',action='write')
  open (unit=35,file='bz_gg.txt',status='unknown', &
       form='formatted',action='write')


	do i = 61,65

		write(i,*)'TITLE="solution of GS equation with flow"'
		write(i,*)'Variables =" R [m] ","z [m]", "var"'
		write(i,*)'ZONE I=',nx,',J=',nz,',F=Point'

	enddo


	do j=1,nz

		z = z_coord(j)

		do i=1,nx

			x = x_coord(i)

			if(sort_grid(i,j,0)>0) then

				write(61,88) x-dx/2, z, gg(i,j,1)
				write(62,88) x, z-dz/2, gg(i,j,2)
				write(63,88) x, z, gg(i,j,3)
				write(64,88) x, z-dz/2, br_gg(i,j)
				write(65,88) x-dx/2, z, bz_gg(i,j)

				if(j==iz) then

					write(31,'(e20.14,''     '',e20.14)') x-dx/2, gg(i,j,1)
					write(32,'(e20.14,''     '',e20.14)') x, gg(i,j,2)
					write(33,'(e20.14,''     '',e20.14)') x, gg(i,j,3)
					write(34,'(e20.14,''     '',e20.14)') x, br_gg(i,j)
					write(35,'(e20.14,''     '',e20.14)') x-dx/2, bz_gg(i,j)

				endif

			else

				write(61,88) x-dx/2, z, 0.d0
				write(62,88) x, z-dz/2, 0.d0
				write(63,88) x, z, 0.d0
				write(64,88) x, z-dz/2, 0.d0
				write(65,88) x-dx/2, z, 0.d0

				if(j==iz) then

					write(31,'(e20.14,''     '',e20.14)') x-dx/2,  0.d0
					write(32,'(e20.14,''     '',e20.14)') x,  0.d0
					write(33,'(e20.14,''     '',e20.14)') x,  0.d0
					write(34,'(e20.14,''     '',e20.14)') x,  0.d0
					write(35,'(e20.14,''     '',e20.14)') x-dx/2, 0.d0

				endif

			endif

		end do
	end do

	do i=31,35
		close(i)
		close(i+30)
	enddo


88   format(3(e12.6,3x))



end subroutine radial_plot_gg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine write_full_grid(thing,truepsi,nx,nz,fname,edgeval)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  use constant
  use solver
  use magnetic
  implicit none

  integer, intent(in) :: nx,nz
  real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: thing, truepsi
  real(kind=dkind) :: edgeval
  character (len=*), intent(in) :: fname
  real (kind=dkind) :: x,z,ex,ez
  integer i,j

    open (unit=17, file=fname//'_full_bin.out', form='unformatted', status='unknown', action='write')



	do j=1,nz

		z = z_coord(j)

		do i=1,nx

			if(rmajor>1.d8) then
			! likely to be an astrophysical case, normalize to rmajor

				x = x_coord(i)/rmajor

			elseif(2.d0*a_elps/(rmajor+rcenter)<1.d-3) then
			! almost cylindrical case, reset x coordinate

				x = x_coord(i)-rcenter

			else
			!standard toroidal case

				x = x_coord(i)

			endif


		    if ( ((sort_grid(i,j,0)==2) .and.(tri_type==13)).or.  &
					((sort_grid(i,j,0)>=1) .and.(tri_type/=13)).or.(fname=='psi_step')  &
					.or.(fname=='residual_step').or.(fname=='rho_step').or.(fname=='psi').or.(fname=='psi_eqdsk1')  &
					.or.(fname=='psi_eqdsk2').or.(fname=='psi_eqdsk3')  &
					.or.(((bc_type==7).or.(bc_type==8)).and.(fname/='psi_clean'))) then

				if((truepsi(i,j)>=fraction*psic).or.(fname=='psi').or.(fname=='psi_step')  &
					.or.(fname=='residual_step').or.(fname=='rho_step').or.(fname=='psi_chease' ).or.(fname=='psi_eqdsk1')  &
					.or.(fname=='psi_eqdsk2').or.(fname=='psi_eqdsk3').or.(bc_type==7).or.(bc_type==8)) then
					write(17) x,z,thing(i,j)
				else
					write(17) x,z,edgeval
				endif

			else
				write(17) x,z,edgeval
			endif

		end do
	end do





	do 33 j=1,nz
	do 33 i=1,nx

	33 write(17) thing(i,j)

	close(17)

	continue

	return

end subroutine write_full_grid

