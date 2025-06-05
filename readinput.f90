!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine readinput(m,input_n, input_p_iso, input_p_par, input_p_perp,&
                input_F, input_omega, input_mtheta, input_psiprim, input_a_elps,&
                input_b_elps, input_k_ellipt, input_delta_up, input_delta_down)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use array_dimensions
	use constant
	use flow
	use magnetic
	use p_d_profile
	use exp_data
	use triangularity
	use solver

	implicit none

        interface
                subroutine read_numerical(input_n, input_p_iso, input_p_par, input_p_perp,&
                                input_F, input_omega, input_mtheta, input_psiprim)
                        real, dimension(:,:), intent(in), optional :: input_n, input_p_iso, input_p_par, input_p_perp,&
                                input_F, input_omega, input_mtheta, input_psiprim 
                end subroutine
        end interface

	integer :: m
        real, dimension(:,:), intent(in), optional :: input_n, input_p_iso, input_p_par, input_p_perp,&
                input_F, input_omega, input_mtheta, input_psiprim 
        real, dimension(:), intent(in), optional :: input_a_elps, input_b_elps, input_k_ellipt,&
                input_delta_up, input_delta_down

	real(kind=dkind) :: Fc_o_Fv		!Fcenter/Fvacuum
	real(kind=dkind) :: de_o_dc		!Dedge/Dcenter
	real(kind=dkind) :: tpae_o_tpac !Tparedge/Tparcenter 
	real(kind=dkind) ::	pe_o_pc		!Pedge/Pcenter

        
!	real(kind=dkind) :: rtsafe
!	external  rtsafe


	include "inp_out.f90"

	open(unit=5,file='inputfile.dat',status='old')

	read(5,input_which)
	read(5,input_triangularity)

        ! If triangularity inputs are provided functionally, then use those to replace the ones read from the inputfile
        ! As noted in flow_mod, these inputs are scalars, but FORTRAN 95 only allows allocatable arrays, so this is the
        ! dumb workaround to make this still work. (if there's a better of doing this, please fix this)
        if(present(input_a_elps)) then
                print*, "a_elps provided functionally; overriding information in inputfile."
                a_elps = input_a_elps(1) !only taking first element in array
        endif
        if(present(input_b_elps)) then
                print*, "b_elps provided functionally; overriding information in inputfile."
                b_elps = input_b_elps(1)
        endif
        if(present(input_k_ellipt)) then
                print*, "k_ellipt provided functionally; overriding information in inputfile."
                k_ellipt = input_k_ellipt(1)
        endif
        if(present(input_delta_up)) then
                print*, "delta_up provided functionally; overriding information in inputfile."
                delta_up = input_delta_up(1)
        endif
        if(present(input_delta_down)) then
                print*, "delta_down provided functionally; overriding information in inputfile."
                delta_down = input_delta_down(1)
        endif

	read(5,input_constants)
	read(5,input_flow)
	read(5,input_magnetic)
	read(5,input_p_d_profile)
	read(5,input_numerical)
	read(5,input_solver)
	read(5,input_gravity)
	read(5,input_MARS)

	if(eq_type==1) eq3_opt=1
	! to avoid conflicts with the p_opt option

	
	if((tri_type==11).and.((eq3_opt==9).or.(p_opt==9))) then
		B_pmax = ana_fac
	endif

	read(5,input_Alf)

	if(Broot==3) then

		read(5,input_Alf)

		M0 = mach_alf_0*dsqrt(2.d0/(beta_center*gamma))

	endif

	if(tri_type==11) then

		read(5,input_LDX)
		b_phi_zero = bpolzero(psic)

	else

		psi_in = 0.d0
		psi_out = 0.d0

	endif

	close(5)


	continue

	m	= 1 + n/2

!	if(tri_type/=0) then
!		a_elps=1.d0
!		b_elps=1.d0
!	endif

	fvaccuum = b_phi_zero*rmajor
	fcenter	= Fc_o_Fv*fvaccuum

	dedge = dcenter * de_o_dc
	pcenter = beta_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag
	pedge = pcenter * pe_o_pc
	qperpcenter =  betaperp_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag 
	qperpedge = qperpcenter * qpee_o_qpec
	qparcenter = betapar_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag
	qparedge = qparcenter * qpae_o_qpac
	tparedge = tparcenter * tpae_o_tpac

	if(input_EQDSK) then

		! all that follows is not needed
		continue

	else

		if(numerical) then

				continue

		else

				numerical_n =		.false.
				numerical_p_iso =	.false.
				numerical_p_par =	.false.
				numerical_p_perp =	.false.
				numerical_F =		.false.
				numerical_omega =	.false.
				numerical_mtheta =	.false.
				numerical_psiprim = .false.

		endif

		if( (eq3_opt==4).or.(eq3_opt==5).or.(p_opt==4).or.(p_opt==5) ) then

			numerical = .true.
			numerical_opt = 1
			gammahut = gamma + alpha

		elseif((eq3_opt==6).or.(p_opt==6)) then

			gammahut = 1.d0 - gamma*alpha

		elseif((eq3_opt==7).or.(eq3_opt==8).or.(eq3_opt==9).or.  &
				(p_opt==7).or.(p_opt==8).or.(p_opt==9)) then

			gammahut = gamma*alpha

		endif


		if((eq3_opt==9).or.(p_opt==9)) then

			pcenter = pedge/((psi_out/psi_in)**gammahut *  &
				(1.d0-psi_out/psi_in)**aux_fact)

			if(eq_type==3) then

				qparedge = pedge
				qparcenter = pcenter

				psi_pmax = rtsafe(find_pmax,psi_out*1.1d0,psi_in*0.9d0,1.d-12,200)

				continue

			endif

		endif

		if((bc_type==7).or.(bc_type==8)) call read_numerical_bc

		continue

		! some controls on the input
		if ((eq_type==1).or.(eq_type==3)) then
			continue
		else
			print*,'error in eq_type :  ',eq_type
		endif

		if(eq_type==3) then
			if ((broot==0).or.(broot==2).or.(broot==5)) then
				print*, 'warning: wrong option for broot for eq_type=3'
				print*, 'broot changed from',broot,'to 1'
				print*, '  '
				broot = 1
			endif
			if ((mach_theta_max>0.d0).or.(mach_theta_edge>0.d0)) then
				print*, 'warning: wrong value for Mach_theta for eq_type=3'
				print*, 'Mach_theta set to 0'
				print*, '  '
				mach_theta_edge = 0.d0
				mach_theta_max = 0.d0
			endif

		endif

		if( (bc_type==7).and.(tri_type==13)) then
			psic_13 = psic - psi_e
		endif

		if( (bc_type==7).and.((tri_type==13).or.(tri_type==-1))) then	! WARNING!!!!!
			continue
		elseif( (bc_type>=3).and.(tri_type/=0).and.(tri_type/=8).and.(tri_type/=9)  &
					.and.(tri_type/=18).and.(tri_type/=88) ) then
			print*, 'warning, option not yet implemented for'
			print*, 'bc_type =', bc_type, 'and tri_type =', tri_type
			print*, 'bc_type changed to 1'
			bc_type = 1
		endif

		if(grid_type==3) call read_grid

	endif

	bc_setup_option = sign(1,bc_type)
	bc_type = abs(bc_type)

	
	  if (numerical) then

			if(numerical_opt==1) then

				call read_numerical(input_n=input_n, input_p_iso=input_p_iso, &
                                input_p_par=input_p_par, input_p_perp=input_p_perp, input_F=input_F, &
                                input_omega=input_omega, input_mtheta=input_mtheta, input_psiprim=input_psiprim)

			elseif(numerical_opt>=2) then

				call read_numerical_super

			else

				pause 'error in numerical_opt'
				stop

			endif

			continue

	  endif


	! the eqdsk file contains several data that are duplicated above, and which will be changed
	! the previous part of the reading process has been preserved to leave control over
	! the details of the FLOW input not contained in eqdsk (TO BE EDITED LATER)
	if(input_EQDSK) then
		call readeqdsk(trim(EQDSK_file))
	endif

	! last check, verifies whether there is any flow in the equilibrium

	if(enforce_static==0) then

		call check_flow

	elseif(enforce_static==1) then

		static_equi = .true.

	elseif(enforce_static==2) then

		static_equi = .false.

	endif



end subroutine readinput

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical(input_n, input_p_iso, input_p_par, input_p_perp,&
                input_F, input_omega, input_mtheta, input_psiprim)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
	use flow
	use pseudo_IMSL, only : DBSNAK, DBSINT, dbsval
	implicit none
	integer i
	integer :: d_dim
        real, dimension(:,:), optional, intent(in) :: input_n, input_p_iso, input_p_par, input_p_perp,&
                input_F, input_omega, input_mtheta, input_psiprim

	if(numerical_n) then

                !if data is provided in a functional way:
                if(present(input_n)) then
                        do i = 1, size(input_n,1)
                                d_data(i,1) = input_n(i,1)
                                d_data(i,2) = input_n(i,2)
                        enddo
                        !print *, d_data
                        d_dim = size(input_n,1)

                !otherwise, open the datafiles instead
                else

		        i = 0

		        open(33,file='n.dat', status='old', action='read')

		        do

			        i = i+1
			        read(33,*,end=19) d_data(i,1), d_data(i,2)

		        enddo

	                19	close(33)
        
		        d_dim = i-1
                endif
		ibreak_d = d_dim + d_ord
                
                if (.not.allocated(d_cscoef)) then
		allocate(d_cscoef(1,d_dim))
                endif


		!this sets up the knots
		call DBSNAK(d_dim, d_data(1:d_dim,1),  &
						d_ord,d_data(1:ibreak_d,3))

		call DBSINT(d_dim, d_data(1:d_dim,1),  &
			d_data(1:d_dim,2), d_ord,  &
			d_data(1:ibreak_d,3),  &
			d_cscoef(1,1:d_dim))
	
	endif

!!$	if(input_EQDSK) return
!!$ March 6 2008: this return has been removed to allow for other numerical inputs LG

!----------------------------------------------------------------------------

	if(numerical_p_iso) then

                !if data is provided in a functional way:
                if(present(input_p_iso)) then
                        do i = 1, size(input_p_iso,1)
                                p_iso_data(i,1) = input_p_iso(i,1)
                                p_iso_data(i,2) = input_p_iso(i,2)
                        enddo
                        !print *, p_iso_data
                        data_dim = size(input_p_iso,1)

                !otherwise, open the datafiles instead
                else      

	                i = 0

	                open(33,file='p_iso.dat', status='old', action='read')

	                do 

		                i = i+1
		                read(33,*,end=66) p_iso_data(i,1), p_iso_data(i,2)

	                enddo

                        66	close(33)

	                data_dim = i-1
                endif
	        ibreak_p = data_dim + p_iso_ord

	        allocate(p_iso_cscoef(1,data_dim))

	        !this sets up the knots
	        call DBSNAK(data_dim, p_iso_data(1:data_dim,1),  &
					        p_iso_ord,p_iso_data(1:ibreak_p,3))

	        call DBSINT(data_dim, p_iso_data(1:data_dim,1),  &
		         p_iso_data(1:data_dim,2), p_iso_ord,  &
		         p_iso_data(1:ibreak_p,3),  &
		         p_iso_cscoef(1,1:data_dim))



!!$	open(33,file='p_iso.dat',  &
!!$			status='old',action='read')
!!$
!!$	do i=1,data_dim
!!$		read(33,*) p_iso_data(i,1),p_iso_data(i,2)
!!$		p_iso_data(i,2) = p_iso_data(i,2) !*2.d0
!!$	enddo
!!$
!!$	close(33)
!!$	call spline(p_iso_data(:,1),p_iso_data(:,2),data_dim,1.d40,1.d40,p_iso_data(:,3))

	endif

!----------------------------------------------------------------------------

	if(numerical_p_par) then

                !if data is provided in a functional way:
                if(present(input_p_par)) then
                        do i = 1, size(input_p_par,1)
                                p_par_data(i,1) = input_p_par(i,1)
                                p_par_data(i,2) = input_p_par(i,2)
                        enddo
                        !print *, p_par_data
                        data_dim = size(input_p_par,1)
                
                !otherwise, open the datafiles instead
                else

	                i = 0

	                open(33,file='p_par.dat',  &
				status='old',action='read')

	                do 

		                i = i+1
		                read(33,*,end=67) p_par_data(i,1), p_par_data(i,2)

	                enddo

67	close(33)

	                data_dim = i-1
                endif
	        ibreak_ppar = data_dim + ppar_ord

	        allocate(ppar_cscoef(1,data_dim))

	        !this sets up the knots
	        call DBSNAK(data_dim, p_par_data(1:data_dim,1),  &
					ppar_ord,p_par_data(1:ibreak_ppar,3))

	        call DBSINT(data_dim, p_par_data(1:data_dim,1),  &
		        p_par_data(1:data_dim,2), ppar_ord,  &
		        p_par_data(1:ibreak_ppar,3),  &
		        ppar_cscoef(1,1:data_dim))


        endif



!----------------------------------------------------------------------------

	if(numerical_p_perp) then

                !if data is provided in a functional way:
                if(present(input_p_perp)) then
                        do i = 1, size(input_p_perp,1)
                                p_perp_data(i,1) = input_p_perp(i,1)
                                p_perp_data(i,2) = input_p_perp(i,2)
                        enddo
                        !print *, p_perp_data
                        data_dim = size(input_p_perp,1)
                
                !otherwise, open the datafiles instead
                else

	                i = 0

	                open(33,file='p_perp.dat',  &
				status='old',action='read')

	                do 

		                i = i+1
		                read(33,*,end=68) p_perp_data(i,1), p_perp_data(i,2)

	                enddo

68	close(33)

	                data_dim = i-1
                endif
	        ibreak_pperp = data_dim + pperp_ord

	        allocate(pperp_cscoef(1,data_dim))

	        !this sets up the knots
	        call DBSNAK(data_dim, p_perp_data(1:data_dim,1),  &
					pperp_ord,p_perp_data(1:ibreak_pperp,3))

	        call DBSINT(data_dim, p_perp_data(1:data_dim,1),  &
		        p_perp_data(1:data_dim,2), pperp_ord,  &
		        p_perp_data(1:ibreak_pperp,3),  &
		        pperp_cscoef(1,1:data_dim))

        endif


!----------------------------------------------------------------------------

	if(numerical_F) then

                !if data is provided in a functional way:
                if(present(input_F)) then
                        do i = 1, size(input_F,1)
                                F_data(i,1) = input_F(i,1)
                                F_data(i,2) = input_F(i,2)
                                F_data(i,2) = F_data(i,2)/rmajor
                        enddo
                        !print *, F_data
                        data_dim = size(input_F,1)
                
                !otherwise, open the datafiles instead
                else
                        
                        open(33,file='b0.dat', status='old', action='read')

	                i = 0

	                do
		                i = i+1
		                read(33,*,end=69) F_data(i,1),F_data(i,2)
		                F_data(i,2) = F_data(i,2)/rmajor
	                enddo

                        69	close(33)

	                data_dim = i-1
                endif
	        ibreak_B0 = data_dim + F_ord

	        allocate(F_cscoef(1,data_dim))

	        !this sets up the knots
	        call DBSNAK(data_dim, F_data(1:data_dim,1),  &
					F_ord,F_data(1:ibreak_B0,3))

	        call DBSINT(data_dim, F_data(1:data_dim,1),  &
		        F_data(1:data_dim,2), F_ord,  &
		        F_data(1:ibreak_B0,3),  &
		        F_cscoef(1,1:data_dim))

        endif


!----------------------------------------------------------------------------

	if(numerical_omega) then

!	inquire('omega.dat')

!!$	open(33,file='omega.dat',  &
!!$				status='old',action='read')
!!$
!!$	do i=1,data_dim
!!$		read(33,*) omega_data(i,1),omega_data(i,2)
!!$		omega_data(i,2) = omega_data(i,2)	*1.d0
!!$	enddo
!!$
!!$	close(33)
!!$	call spline(omega_data(:,1),omega_data(:,2),data_dim,1.d40,10.d40,omega_data(:,3))

		if(omega_option==1) then

			numerical_mphi = .true.
			numerical_omega = .false.

		endif

                !NOTE: not accounting for an omega.dat as that seems to be deprecated
                !if data is provided in a functional way:
                if(present(input_omega)) then
                        do i = 1, size(input_omega,1)
                                mphi_data(i,1) = input_omega(i,1)
                                mphi_data(i,2) = input_omega(i,2)
                        enddo
                        !print *, mphi_data
                        data_dim = size(input_omega,1)

                !otherwise, open the datafiles instead
                else 

		        open(33,file='mphi.dat', status='old', action='read')

		        i = 0

		        do

			        i = i+1

			        read(33,*,end=98) mphi_data(i,1),mphi_data(i,2)

		        enddo

		98	close(33)

		        data_dim = i-1
                
                endif

		ibreak_fi = data_dim + mphi_ord

		allocate(mphi_cscoef(1,data_dim))


		!this sets up the knots
		call DBSNAK(data_dim, mphi_data(1:data_dim,1),  &
						mphi_ord,mphi_data(1:ibreak_fi,3))

		call DBSINT(data_dim, mphi_data(1:data_dim,1),  &
			 mphi_data(1:data_dim,2), mphi_ord,  &
			 mphi_data(1:ibreak_fi,3),  &
			 mphi_cscoef(1,1:data_dim))



	endif

!----------------------------------------------------------------------------

	if(numerical_mtheta) then

                !if data is provided in a functional way:
                if(present(input_mtheta)) then
                        do i = 1, size(input_mtheta,1)
                                mtheta_data(i,1) = input_mtheta(i,1)
                                mtheta_data(i,2) = input_mtheta(i,2)
                        enddo
                        !print *, mtheta_data
                        data_dim = size(input_mtheta,1)

                !otherwise, open the datafiles instead
                else 

		        open(33,file='mtheta.dat', status='old', action='read')

		        i = 0

		        do

			        i = i+1

			        read(33,*,end=109) mtheta_data(i,1),mtheta_data(i,2)

		        enddo

		        109	close(33)

		        data_dim = i-1
                endif

!----------------------------------------------------------------------------

		ibreak_th = data_dim + mtheta_ord

		allocate(mtheta_cscoef(1,data_dim))

		! interpolation setup has been moved to a separate function
		call interp_setup(data_dim,mtheta_ord, &
			mtheta_data(1:data_dim,1),mtheta_data(1:data_dim,2), &
			mtheta_data(1:ibreak_th,3),mtheta_cscoef(1,1:data_dim))

		if((Broot==0).or.(Broot==5)) then
		! set maximum value of Mach_theta

			mach_theta_num = 0.d0

			do i=1,data_dim

				if(abs(mtheta_data(i,2))>mach_theta_num) then

				        mach_theta_num = mtheta_data(i,2)

				endif

			enddo

			mach_theta_max = mach_theta_num

		else

			! any two numbers will do, as long as they are equal
			mach_theta_max = 1.d-1
			mach_theta_num = 1.d-1

		endif
                
	endif

!----------------------------------------------------------------------------
!NOTE: This seems to be deprecated now, as I don't know what vol.dat is supposed to be

	if( (eq3_opt==4).or.(eq3_opt==5).or.  &
			(p_opt==4).or.(p_opt==5) ) then

		open(33,file='vol.dat', status='old', action='read')

		i = 0

		do

			i = i+1
			read(33,*,end=189) vol_data(i,1), vol_data(i,2)

		enddo

	189	close(33)

		vol_dim = i-1
		ibreak_v = vol_dim + vol_ord

		allocate(vol_cscoef(1,vol_dim))


		!this sets up the knots
		call DBSNAK(vol_dim, vol_data(1:vol_dim,1),  &
						vol_ord,vol_data(1:ibreak_v,3))

		call DBSINT(vol_dim, vol_data(1:vol_dim,1),  &
			 vol_data(1:vol_dim,2), vol_ord,  &
			 vol_data(1:ibreak_v,3),  &
			 vol_cscoef(1,1:vol_dim))

		if((eq3_opt==4).or.(p_opt==4)) vol1 = dbsval(psifrac, vol_ord, vol_data(1:ibreak_v,3),  &
									vol_dim, vol_cscoef(1,1:vol_dim) )

		if((eq3_opt==5).or.(p_opt==5)) vol1 = dbsval(1.d0, vol_ord, vol_data(1:ibreak_v,3),  &
									vol_dim, vol_cscoef(1,1:vol_dim) )

	endif


!----------------------------------------------------------------------------

	if((numerical_psiprim).and.(bc_type==5)) then

                !if data is provided in a functional way:
                if(present(input_psiprim)) then
                        do i = 1, size(input_psiprim,1)
                                psiprim_data(i,1) = input_psiprim(i,1)
                                psiprim_data(i,2) = input_psiprim(i,2)
                        enddo
                        !print *, psiprim_data
                        psiprim_dim = size(input_psiprim,1)

                !otherwise, open the datafiles instead
                else 
		        i = 0

		        open(33,file='psiprim.dat', status='old', action='read')

		        do

			        i = i+1
			        read(33,*,end=199) psiprim_data(i,1), psiprim_data(i,2)

		        enddo

	199	close(33)

		        psiprim_dim = i-1
                endif
		ibreak_pp = psiprim_dim + psiprim_ord

		allocate(psiprim_cscoef(1,psiprim_dim))


		!this sets up the knots
		call DBSNAK(psiprim_dim, psiprim_data(1:psiprim_dim,1),  &
						psiprim_ord,psiprim_data(1:ibreak_pp,3))

		call DBSINT(psiprim_dim, psiprim_data(1:psiprim_dim,1),  &
			 psiprim_data(1:psiprim_dim,2), psiprim_ord,  &
			 psiprim_data(1:ibreak_pp,3),  &
			 psiprim_cscoef(1,1:psiprim_dim))

	endif

!----------------------------------------------------------------------------


	continue

end subroutine read_numerical

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical_bc
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
	use pseudo_IMSL, only : DBSNAK, DBSINT, dbsval
	implicit none
	integer i


	if((bc_type==7).or.(bc_type==8)) then

		i = 0

		open(33,file='psioftheta.dat', status='old', action='read')

		do

			i = i+1
			read(33,*,end=299) psib_data(i,1), psib_data(i,2)

		enddo

	299	close(33)

		psib_dim = i-1
		ibreak_pb = psib_dim + psib_ord

		allocate(psib_cscoef(1,psib_dim))


		!this sets up the knots
		call DBSNAK(psib_dim, psib_data(1:psib_dim,1),  &
						psib_ord,psib_data(1:ibreak_pb,3))

		call DBSINT(psib_dim, psib_data(1:psib_dim,1),  &
			 psib_data(1:psib_dim,2), psib_ord,  &
			 psib_data(1:ibreak_pb,3),  &
			 psib_cscoef(1,1:psib_dim))

	endif



end subroutine read_numerical_bc


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_numerical_super
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! this one needs to read B0, P, M_theta, M_phi
! data_dim is not assigned, but is the same for all data, 
! so it's derivad from the data

	use constant
	use magnetic
	use exp_data
	use interpolating_functions
	use p_d_profile
!	use IMSL, only : dcsdec, dc2con
	use pseudo_IMSL, only : DBSNAK, DBSINT
	implicit none
	integer i

	integer ::  option = 1 !0=normal, 1=shape preserving

	real(kind=dkind), dimension(:,:), allocatable :: dummy
	integer, dimension(data_dim_lim) :: dumint
	integer :: itmax=2000

	if(input_EQDSK) return

! start from the pressure

	i = 0

	open(33,file='p_iso.dat', status='old', action='read')

	do 

		i = i+1
		read(33,*,end=99) p_iso_data(i,1), p_iso_data(i,2)

	enddo

99	close(33)

	data_dim = i-1

	ibreak_p = data_dim + p_iso_ord
	ibreak_B0 = data_dim + F_ord
	ibreak_fi = data_dim + mphi_ord
	ibreak_th = data_dim + mtheta_ord




	allocate(p_iso_cscoef(1,data_dim))
	allocate(F_cscoef(1,data_dim))
	allocate(mphi_cscoef(1,data_dim))
	allocate(mtheta_cscoef(1,data_dim))


	!this sets up the knots
	call DBSNAK(data_dim, p_iso_data(1:data_dim,1),  &
					p_iso_ord,p_iso_data(1:ibreak_p,3))

	call DBSINT(data_dim, p_iso_data(1:data_dim,1),  &
		 p_iso_data(1:data_dim,2), p_iso_ord,  &
		 p_iso_data(1:ibreak_p,3),  &
		 p_iso_cscoef(1,1:data_dim))



! next B0

	open(33,file='b0.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) F_data(i,1),F_data(i,2)
		F_data(i,2) = F_data(i,2)/rmajor
	enddo

	close(33)

	!this sets up the knots
	call DBSNAK(data_dim, F_data(1:data_dim,1),  &
					F_ord,F_data(1:ibreak_B0,3))

	call DBSINT(data_dim, F_data(1:data_dim,1),  &
		 F_data(1:data_dim,2), F_ord,  &
		 F_data(1:ibreak_B0,3),  &
		 F_cscoef(1,1:data_dim))



! next M_theta

	open(33,file='mtheta.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) mtheta_data(i,1),mtheta_data(i,2)
!		mtheta_data(i,2) = mtheta_data(i,2)*.0d0
	enddo

	close(33)

	!this sets up the knots
	call DBSNAK(data_dim, mtheta_data(1:data_dim,1),  &
					mtheta_ord,mtheta_data(1:ibreak_th,3))

	call DBSINT(data_dim, mtheta_data(1:data_dim,1),  &
		 mtheta_data(1:data_dim,2), mtheta_ord,  &
		 mtheta_data(1:ibreak_th,3),  &
		 mtheta_cscoef(1,1:data_dim))



! next M_phi

	open(33,file='mphi.dat', status='old', action='read')

	do i=1,data_dim
		read(33,*) mphi_data(i,1),mphi_data(i,2)
	enddo

	close(33)


	!this sets up the knots
	call DBSNAK(data_dim, mphi_data(1:data_dim,1),  &
					mphi_ord,mphi_data(1:ibreak_fi,3))

	call DBSINT(data_dim, mphi_data(1:data_dim,1),  &
		 mphi_data(1:data_dim,2), mphi_ord,  &
		 mphi_data(1:ibreak_fi,3),  &
		 mphi_cscoef(1,1:data_dim))


	if((numerical_omega).and.(omega_option==1) ) then

	 	numerical_mphi = .true.
		numerical_omega = .false.

	endif


!	set other flags

!	numerical_p_iso = .true.
!	numerical_F = .true.
!	numerical_mtheta = .true.

 	numerical_n 		=	.false. 			!density
!!$	numerical_p_par 	=	.false. 			!parallel pressure
!!$	numerical_p_perp 	=	.false. 			!perpendicular pressure
!!$	numerical_omega 	=	.false. 			!toroidal rotation

	if(numerical_opt==3) then
	! also read and interpolate X

		print*, 'option has been disabled'
		numerical_opt = 2

!!$		allocate(X0_cscoef(4,data_dim))
!!$		allocate(X1_cscoef(4,data_dim))
!!$
!!$		! first X0
!!$
!!$		open(33,file='X0.dat', status='old', action='read')
!!$
!!$		do i=1,data_dim
!!$			read(33,*) X0_data(i,1),X0_data(i,2)
!!$		enddo
!!$
!!$		close(33)
!!$
!!$		call dcsdec (data_dim, X0_data(1:data_dim,1),  &
!!$					 X0_data(1:data_dim,2), 2, 0.d0, 2,  &
!!$					 0.d0, X0_data(1:data_dim,3), X0_cscoef)
!!$
!!$		! then X1
!!$
!!$		open(33,file='X1.dat', status='old', action='read')
!!$
!!$		do i=1,data_dim
!!$			read(33,*) X1_data(i,1),X1_data(i,2)
!!$		enddo
!!$
!!$		close(33)
!!$
!!$		call dcsdec (data_dim, X1_data(1:data_dim,1),  &
!!$					 X1_data(1:data_dim,2), 2, 0.d0, 2,  &
!!$					 0.d0, X1_data(1:data_dim,3), X1_cscoef)

	endif


	continue

	return

end subroutine read_numerical_super

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine read_grid
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant

	implicit none

	integer :: i, dsize
	real(kind=dkind), dimension(5000) :: xtemp

	i = 0

	open(33,file='input_xgrid.dat', status='old', action='read')

	do

		i = i+1
		read(33,*,end=111) xtemp(i)

	enddo

111	close(33)

	dsize = i-1

	allocate(dx_ext(1:dsize))

	do i=1,dsize

		dx_ext(i) = xtemp(i)

	enddo

	allocate(dz_ext(1:dsize))

	open(33,file='input_zgrid.dat', status='old', action='read')

	do i = 1,dsize

		read(33,*) dz_ext(i)

		if(i>2000) then
			continue
		endif

	enddo

	continue


end subroutine read_grid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine readeqdsk(filename)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	use constant
	use triangularity
	use magnetic
	use exp_data

	implicit none

        interface
                subroutine radial_plot(psi, truepsi, nx, nz, fname, iz, fname_output_data)
                        import :: skind
                        import :: dkind
                        integer, intent(in) :: nx,nz,iz
                        real (kind=skind), dimension(1:nx,1:nz), intent(in) :: psi
                        real (kind=dkind), dimension(1:nx,1:nz), intent(in) :: truepsi
                        character (len=*), intent(in) ::fname
                        real (kind=dkind), dimension(1:nx*nz,1:3), intent(out), optional :: fname_output_data
                end subroutine 
        end interface

	real(kind=dkind), dimension(:,:), allocatable :: psi_eqdsk
	real(kind=skind), dimension(:,:), allocatable :: out
	real(kind=dkind), dimension(:), allocatable :: dummy_arr
	real(kind=dkind), dimension(:), allocatable :: R_input, Z_input
	real(kind=dkind) :: dummy, psimax, rmin, rmax
	real(kind=dkind), dimension(1:20) :: dummy_list
	real(kind=dkind) :: zmin, zmax, rcenter_one,rcenter_two, z_center

	character  cdate*8, dummychar*48
	character(len=*) :: filename

	integer :: i, j
	integer :: iomap
	integer :: idummy, nxloc, nzloc
	integer :: grid_type_save	! this allows saving the incoming data

	grid_type_save = grid_type

	if(grid_type<0) grid_type = -(10 + grid_type)


	iomap = 22

	open(iomap,file=filename,form='formatted',action='read')

	! as far as FLOW is concerned, the first line is garbage

	read(iomap,9380)dummychar,cdate,idummy,nxloc,nzloc

	! the second line contains useful staff: x_size,z_size,rmajor, dummies

	read(iomap,9381) x_size,z_size,rmajor,dummy_list(1), dummy_list(2)
	! in FLOW the grid is in physical units
	rcenter = rmajor ! this would be the default in most codes

	! the third line is essentially garbage, but we can read psimax, just in case

     read(iomap,9381) dummy,dummy, psimax,dummy,b_phi_zero

	! as far as FLOW is concerned, the fourth and fifth lines are garbage

	read(iomap,9381) dummy_list(3),dummy_list(4),dummy_list(5),dummy_list(6),dummy_list(7)	!4th line
	read(iomap,9381) dummy_list(8),dummy_list(9),dummy_list(10),dummy_list(11),dummy_list(12)	!5th line

	! the sixth entry is F(psi), this one is needed by FLOW

	numerical_F = .true.

	ibreak_B0 = nxloc + F_ord

	if(allocated(F_cscoef)) deallocate(F_cscoef)
	! F_cscoef could have been created earlier in the read_numerical part:
	! if so, it needs to be replaced
	allocate(F_cscoef(1,nxloc))
	! F_data is already allocated, that will be changed
	F_data(:,:) = 0.d0

	! need to set up the "x" points for the interpolation

	do i=1,nxloc
		F_data(i,1) = 1.d0/(nxloc-1.d0)*(i-1.d0)
	enddo

	read(iomap,9381) (F_data(nxloc+1-i,2), i=1,nxloc)

	do i=1,nxloc
		F_data(i,2) = F_data(i,2)/rmajor
		! FLOW uses B0, not F
	enddo

	! interpolation setup has been moved to a separate routine
	call interp_setup(nxloc,F_ord, &
		F_data(1:nxloc,1),F_data(1:nxloc,2), &
		F_data(1:nxloc+F_ord,3),F_cscoef(1,1:nxloc))


	! the seventh entry is P(psi), this one is also needed by FLOW

	numerical_p_iso = .true.

	ibreak_p = nxloc + p_iso_ord

	if(allocated(p_iso_cscoef)) deallocate(p_iso_cscoef)
	! p_iso_cscoef could have been created earlier in the read_numerical part:
	! if so, it needs to be replaced
	allocate(p_iso_cscoef(1,nxloc))
	! p_iso_data is already allocated, that will be changed
	p_iso_data(:,:) = 0.d0

	! need to set up the "x" points for the interpolation

	do i=1,nxloc
		p_iso_data(i,1) = 1.d0/(nxloc-1.d0)*(i-1.d0)
	enddo

	read(iomap,9381) (p_iso_data(nxloc+1-i,2), i=1,nxloc)

	! interpolation setup has been moved to a separate routine
	call interp_setup(nxloc,p_iso_ord, &
		p_iso_data(1:nxloc,1),p_iso_data(1:nxloc,2), &
		p_iso_data(1:nxloc+p_iso_ord,3),p_iso_cscoef(1,1:nxloc))

	! entries 8 to 11 are not needed

	allocate(dummy_arr(1:nxloc))
	allocate(psi_eqdsk(1:nxloc,1:nzloc))

	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 8th
	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 9th
	read(iomap,9381) ((psi_eqdsk(i,j), i=1,nxloc), j=1,nzloc)	! 10th
	read(iomap,9381) (dummy_arr(i), i=1,nxloc)	! 11th

	deallocate(dummy_arr)



	! the twelwth entry is the shape of the boundary, this one is needed

	read(iomap,1204) theta_points, idummy

	allocate(R_input(theta_points))
	allocate(Z_input(theta_points))

	read(iomap,9381) (R_input(i),Z_input(i), i=1,theta_points)

	! the next entry is the limiter, not needed by FLOW,
	! and the last lines contain some equilibrium values, also not needed by FLOW

	! read the limiter nevertheless, to check the grid

	allocate(dummy_arr(2*idummy))

	read(iomap,9381) (dummy_arr(i),dummy_arr(i+idummy),i=1,idummy)

	continue

	close(unit=iomap,status='keep')

	rmin = rmajor
	rmax = 0.d0

	do i = 1,theta_points

		if(R_input(i)>rmax) rmax = R_input(i)
		if(R_input(i)<rmin) rmin = R_input(i)

	enddo

	rcenter_one = (rmax+rmin)/2.d0


	rmin = rmajor
	rmax = 0.d0
	zmin = rmajor
	zmax = 0.d0

	do i = 1,idummy

		if(dummy_arr(i)>rmax) rmax = dummy_arr(i)
		if(dummy_arr(i)<rmin) rmin = dummy_arr(i)

		if(dummy_arr(i+idummy)>zmax) zmax = dummy_arr(i+idummy)
		if(dummy_arr(i+idummy)<zmin) zmin = dummy_arr(i+idummy)

	enddo

	rcenter_two = (rmax+rmin)/2.d0
	z_center = (zmax+zmin)/2.d0	! just out of curiosity

	deallocate(dummy_arr)

	! the input eqdsk psi is saved three times, since the definition of the grid is ambiguous

	allocate(out(nxloc,nzloc))

	out(1:nxloc,1:nzloc) = psi_eqdsk(1:nxloc,1:nzloc)

	call set_grid(nxloc,nzloc)

	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk1",nxloc/2)
	call psi_boundary_plot(nxloc,nzloc,psi_eqdsk)


	rcenter = rcenter_two

	call set_grid(nxloc,nzloc)

	out(1:nxloc,1:nzloc) = psi_eqdsk(1:nxloc,1:nzloc)
	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk2",nxloc/2)

	rcenter = rcenter_one	! this is the final rcenter
	call set_grid(nxloc,nzloc)
	call radial_plot(out,psi_eqdsk,nxloc,nzloc,"psi_eqdsk3",nxloc/2)


	deallocate(out)
	deallocate(psi_eqdsk)

	rcenter = rcenter_one	! this is the final rcenter, repeated here for clarity

	call plasma_shape_conversion(R_input,Z_input)


	deallocate(R_input)
	deallocate(Z_input)

	grid_type = grid_type_save

	continue

	return

1204    format(2I5)	!format(5I5)
9380    format(A40,A8,3I4)
9381    format(1P,5E16.9)

end subroutine readeqdsk


!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$subroutine M3D_readeqdsk
!!$!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!$
!!$!!$ ***************************************************************
!!$!!$
!!$!!$  read file in Lao 12-1-95 format for DIIID shot.
!!$!!$  interpolate psi(x,y) to unstructured mesh
!!$!!$  interpolate p(psi) onto unstructured mesh.
!!$!!$   R0 = 1.660 m, a = 0.644 m, Kappa = 1.84, delta = 0.78
!!$!!$   B0 = 2.078 T at R = 1.6955 m, Ip = 1.59 MA
!!$!!$   q(0) = 3.26, q95 = 5.81, qmin = 2.17    , r(qmin)/a = 0.695
!!$!!$   P(0)/ <P> = 5.7, Volume = 22.55 m3, Betan = 2.0
!!$!!$   betat = 2.3 %, betap = 0.97, li = 1.03, betat* = 4.0 %
!!$!!$
!!$!!$  Date: Mon, 4 Dec 1995 16:19:00 -0800 (PST)
!!$!!$  From: LAO@GAV.GAT.COM ( Lang Lao )
!!$!!$  The units of PSIRZ and PRES are (as described in
!!$!!$  the format file) Weber/radian (poloidal flux/two pi)
!!$!!$  and Newton/m2.  FPOL is R BT (poloidal current function
!!$!!$  BT twopi R = mu0 poloidal current).
!!$!!$  J_Toroidal is related to PPRIME and FFPRIM through
!!$!!$  J_T= -(R PPRIME + FFPRIM / R)
!!$!!$
!!$!!$   HS normalization: L = 1m,  (x,y) = L (x~, y~)
!!$!!$   psi~ = psi(Weber) / ( B0(T) L(m)^2 ),  
!!$!!$   p~ = 4 pi mu0 p(Newton/m^2) R0(m) / ( L(m) B0(T)^2 ), 
!!$!!$   mu0 = 1e-7H/m, where .~ is dimensionless
!!$
!!$	use constant, only : dkind
!!$
!!$	implicit none
!!$
!!$#include "param1" 
!!$#include "cli1"
!!$#include "cli2"
!!$#include "grfblk"
!!$#include "eqdsk.h"
!!$
!!$	logical ok
!!$	data B0 / 2.078 /, q0 / 3.26 /, R0 / 1.66 /
!!$	DATA TX/.3955/, TY/.9765/
!!$
!!$	sone = 1.
!!$	pi = 4.*atan(sone)
!!$	R0 = rmajor
!!$
!!$	inquire( file='eqdsk.data', exist=ok )
!!$
!!$	if( .not. ok ) then
!!$
!!$		print*, ' inplasma:d3read: file eqdsk.data not found '
!!$		return
!!$
!!$	else
!!$
!!$		open(9,file='eqdsk.data',status='old')
!!$		read (9,2000) (case(i),i=1,6),idum,nw,nh
!!$		print*,' d3read 1', case(1)
!!$		print*,' d3read cases', (case(i),i=1,6)
!!$		icase = 0
!!$		if( case(1) .eq. '  JSOLVE') icase = 1
!!$		if( case(2) .eq. 'M CHEASE') icase = 5
!!$		if( case(2) .eq. 'M CHEASE') icase = 3
!!$		if( case(1) .eq. '  ITERB ') icase = 2
!!$		if( case(1) .eq. '  EFITD ') icase = 3
!!$
!!$		if( idum .eq. 13 ) icase = 4
!!$
!!$		print*,' d3read icase ',icase
!!$		print *,'d3read: nw, nh ', nw, nh
!!$		read (9,2020) rdim,zdim,rcentr,rleft,zmid
!!$		read (9,2020) rmaxis,zmaxis,simag,sibry,bcentr
!!$		print *,'d3read: rdim, zdim, rleft ', rdim, zdim, rleft
!!$		print *,'d3read: rmaxis, zmaxis, rcentr ',
!!$		&                    rmaxis, zmaxis, rcentr
!!$		read (9,2020) current,simag,xdum1,rmaxis,xdum2
!!$		print *, 'current ', current
!!$		read (9,2020) zmaxis,xdum3,sibry,xdum,xdum
!!$		print *, 'zmaxis ', zmaxis
!!$
!!$		if( icase .eq. 4 ) print *, 'beta,beta0,betaN = ',xdum1,xdum2,xdum3
!!$
!!$		read (9,2020) (fpol(i),i=1,nw)
!!$		print *,'d3read: fpol(); R0 B0' , fpol(1), R0 * B0
!!$		read (9,2020) (pres(i),i=1,nw)
!!$		read (9,2020) (ffprim(i),i=1,nw)
!!$		read (9,2020) (pprime(i),i=1,nw)
!!$
!!$		read (9,2020) ((psirz(i,j),i=1,nw),j=1,nh)
!!$		read (9,2020) (qpsi(i),i=1,nw)
!!$		read (9,2022) nbbbs,limitr
!!$		print *,'d3read: nbbbs,limitr ',nbbbs,limitr
!!$		read (9,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
!!$		print *,'d3read: rbbbs, zbbbs read '
!!$		read (9,2020) (xlim(i),ylim(i),i=1,limitr)
!!$		print *,'d3read: xlim, ylim read '
!!$		nmass = 0
!!$
!!$		close(9)
!!$
!!$	endif
!!$
!!$	write(6,*) 'd3read: nmass=',nmass
!!$	write(0,*) 'd3read: nmass=',nmass
!!$	print*,' d3read a B_T / I =', .5*rdim*bcentr/(current*1.e-6)
!!$
!!$	!!$  shift eqdsk grid vertically to align with vmec grid
!!$
!!$	delz = zmaxis - zmid
!!$
!!$	if( itearing.eq.1 ) delz = 0.
!!$
!!$	zmid = zmid - delz
!!$	zmaxis = zmaxis - delz
!!$
!!$	do l = 1, limitr
!!$		ylim(l) = ylim(l) - delz
!!$	enddo
!!$
!!$	2000 format (6a8,3i4)
!!$	2020 format (5e16.9)
!!$	2022 format (2i5)
!!$
!!$	acen = psirz(nw/2,nh/2)
!!$	print*,' d3read: acen = ',acen
!!$	tc1 = 1.
!!$
!!$	!!$  get psi on separatrix
!!$
!!$	is = (rmaxis - rleft) * (nw - 1.) / rdim + 1
!!$	js = (zmaxis + .5*zdim - zmid) * (nh - 1.) / zdim + 1
!!$	psimaxis = tc1*psirz(is,js)
!!$	psis = 0.0
!!$	rmin = rbbbs(1)
!!$	rmax = rbbbs(1)
!!$
!!$	do i = 1, nbbbs
!!$
!!$		is = (rbbbs(i) - rleft) * (nw - 1.) / rdim + 1
!!$		js = (zbbbs(i) + .5*zdim - zmid) * (nh - 1.) / zdim + 1
!!$		psis = psis + tc1*psirz(is,js)
!!$		if( rmin .gt. rbbbs(i) ) rmin = rbbbs(i)
!!$		if( rmax .lt. rbbbs(i) ) rmax = rbbbs(i)
!!$
!!$	enddo
!!$
!!$	psis = psis / nbbbs
!!$
!!$	do j = 1, nh
!!$	do i = 1, nw
!!$
!!$		psirz(i,j) = psirz(i,j) - psis
!!$
!!$	enddo
!!$	enddo
!!$
!!$	psimaxis = psimaxis - psis
!!$	rminor = .5 * ( rmax - rmin )
!!$	print*,' d3read: psis,psimaxis,rminor = ',psis,psimaxis,rminor
!!$	print*,' d3read: rminor,ratio = ',rminor,rdim/(2.*rminor)
!!$	if( psimaxis .gt. psis ) tc1 = -tc1
!!$
!!$	!!$  get psi on axis and separatrix
!!$
!!$	rminor = .5 * ( rmax - rmin )
!!$	print*,' d3read: psis,psimaxis,rminor = ',psis,psimaxis,rminor
!!$	print*,' d3read: rminor,ratio = ',rminor,rdim/(2.*rminor)
!!$
!!$	!!$  get psi on axis and separatrix
!!$	call inteqdsk( psimaxis, rmaxis, zmaxis )
!!$
!!$	psimin = tc1*psimaxis
!!$	psimax = tc1*psis
!!$
!!$	do j = 1, nh
!!$	do i = 1, nw
!!$		psirz(i,j) = tc1*psirz(i,j)
!!$	enddo
!!$	enddo
!!$
!!$	call intrpEQ( a(1,1), psirz, xlim, ylim, nw, nh )
!!$
!!$	do j = 1, nh
!!$	do i = 1, nw
!!$		psirz(i,j) = rleft + rdim*(i - 1.)/(nw - 1.)
!!$	enddo
!!$	enddo
!!$
!!$
!!$
!!$	do k = 2, ku
!!$	do l = 1, ld
!!$		a(l,k) = a(l,1)
!!$	enddo
!!$	enddo
!!$
!!$	if( icase .eq. 5 ) then
!!$
!!$		do l = 1, nw
!!$			ffprim(l) = -ffprim(l)
!!$			pprime(l) = -pprime(l)
!!$		enddo
!!$
!!$	endif
!!$
!!$	print*,' d3read: pres max/min = ',pmax,pmin
!!$	print*,' d3read: fpol max/min = ',fmax,fmin
!!$	print*,' d3read: ffprim max/min = ',fpmax,fpmin
!!$	print*,' d3read: pprime max/min = ',ppmax,ppmin
!!$	write(6,*)' d3read: pres max/min = ',pmax,pmin
!!$	write(6,*)' d3read: fpol max/min = ',fmax,fmin
!!$	write(6,*)' d3read: ffprim max/min = ',fpmax,fpmin
!!$	write(6,*)' d3read: pprime max/min = ',ppmax,ppmin
!!$	cmu = 4.*pi*1.e-7
!!$	scale = 1.
!!$	!!$     if( itearing .eq. 4 ) scale = 1./eqscale
!!$	!!$     if( icase .eq. 4 ) cmu = 1.
!!$	tc1 = 1.
!!$	if( icase .eq. 5 ) tc1 = -1.
!!$	do k = 1, ku
!!$	do l = 1, ld
!!$	!!$     c(l,k) = -tc1*(cmu*wb(l,k)*(scale*rpls1(l,k))**2 + wa(l,k))*scale
!!$	!!$     c(l,k) = -tc1*(cmu*wb(l,k)*wd(l,1)**2 + wa(l,k))*scale**2
!!$	c(l,k) = -tc1*(cmu*wb(l,1)*wd(l,1)**2 + wa(l,1))*scale
!!$	a(l,k) = a(l,1)/scale*tc1
!!$	!!$     p(l,k) = p(l,k)/scale
!!$	si(l,k) = wc(l,1)
!!$	if( rv(l) .gt. 1. ) then
!!$	p(l,k) = 0.
!!$	c(l,k) = 0.
!!$	endif
!!$	enddo 
!!$	enddo 
!!$	#ifdef PLOT
!!$	call wplot(a,'a0',zero,0)
!!$	call wplot(wa,'ffp',zero,0)
!!$	call wplot(wb,'ppr',zero,0)
!!$	call wplot(si,'si',zero,0)
!!$	call wplot(c,'c',zero,0)
!!$	call wpri (c,'c',zero,0)
!!$	#endif
!!$
!!$	return
!!$
!!$end subroutine M3D_readeqdsk
