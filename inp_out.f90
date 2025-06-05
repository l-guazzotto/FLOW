!	input and open files 

	namelist/input_which/ input_EQDSK, EQDSK_file

	namelist/input_triangularity/ tri_type, n_tri, theta_points_max,  &
								  a_elps, b_elps, k_ellipt, delta_up,  &
								  delta_down

	namelist/input_constants/  rmajor, rcenter, zcenter, x_size , z_size,  &
							grid_type,  &
							eq_type, eq3_opt, data_dim, numerical,   &
							Broot, mass, me_ov_mi, pe_ov_p, eV

	namelist/input_flow/ mach_phi_max, alpha_mphi, mphi_min, &
						 mach_theta_max, mach_theta_edge, t_mth, w_mth, enforce_static

	namelist/input_magnetic/  b_phi_zero, F_opt, Fc_o_Fv, kappa, eta_P, mu_RFP,  &
							  mu_mag, psi_e, inter_switch,  Z_eff_0, delta_Z_eff, Z_eff_exp

	namelist/input_p_d_profile/ gamma, dcenter, de_o_dc, alpha, alpha_rho,  &
								beta_center, betaperp_center, betapar_center,  &
								pe_o_pc, qpee_o_qpec, qpae_o_qpac,  &
								tparcenter, tpae_o_tpac, theteps, p_opt
  
	namelist/input_numerical/  numerical_opt, numerical_n, numerical_p_iso,   &
							   numerical_p_par, numerical_p_perp, numerical_F,  &
							   numerical_omega, omega_option, numerical_mtheta,  &
							   numerical_psiprim

	namelist/input_solver/ n_min, n, min_it, max_it, accelerate, fix_orp,  &
						   bc_type, max_orp, bc_switch, write_all, write_all_bin, restart,  &
						   jump_option, lambda0_fix, dir_switch, delta_Bern_fact, r_orp

	namelist/input_Alf/ mach_alf_0, fraction, ffrac, amult, pmult, pfact, pexp,  &
						p_eps, b_vert, b_edge

	namelist/input_LDX/ psi_in, psi_out, psifrac, alpha_th, aux_fact

	namelist/input_gravity/ gravity_type, G_gravity, M_gravity, Kepler_Omega

	namelist/input_MARS/ MARS_output, enq, nsurf

!------------------------------------------------------------------------

!	open(unit=5,file='inputfile.dat',status='old')
