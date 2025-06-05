module p_d_profile

  use constant, only: dkind
!  use magnetic
  implicit none

  ! The adiabatic exponent, $\gamma$
  real (kind=dkind) :: gamma

  real (kind=dkind) :: dcenter
  real (kind=dkind) :: dedge

  ! The Pressure P = pedge + pcenter * (psi/psic)^alpha
  ! The Density D = dedge + dcenter * (psi/psic)^alpha_rho
  ! pcenter is the Gas Pressure at the Magnetic Axis
  ! dcenter is the Gas Density at the Magnetic Axis
  real (kind=dkind) :: alpha
  real (kind=dkind) :: alpha_rho
  real (kind=dkind) :: alpha_th ! anisotropy for LDX

  real (kind=dkind) :: beta_center
  real (kind=dkind) :: betaperp_center
  real (kind=dkind) :: betapar_center

  real (kind=dkind) :: pcenter
!      pcenter =  beta_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag

  real (kind=dkind) :: pedge

  real (kind=dkind) :: qperpcenter
!	    qperpcenter	=   betaperp_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag

  real (kind=dkind) :: qperpedge

  real(kind=dkind) ::	qpee_o_qpec !qperpedge/qperpcenter
  real(kind=dkind) ::	qpae_o_qpac !qparedge/qparcenter

  real (kind=dkind) :: qparcenter
! qparcenter =  betapar_center*b_phi_zero*b_phi_zero/2.0d0/mu_mag

  real (kind=dkind) :: qparedge

  real (kind=dkind) :: tparcenter	!T_par=P_par*D

  real (kind=dkind) :: tparedge

  real (kind=dkind) :: theteps		!anisotropy

  real(kind=dkind), parameter :: q0 = 2.d0

  real(kind=dkind) :: psifrac
  ! for edge pressure in LDX

  real(kind=dkind) :: K_LDX, vol1
  ! for peak pressure in LDX
  
  real(kind=dkind) :: gammahut
  ! for LDX profile: p*v**gammahut=const.

  real(kind=dkind) :: aux_fact
  ! for pressure slope in LDX

  integer :: p_opt
  ! option for pressure equation


end module p_d_profile
