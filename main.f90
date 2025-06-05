! This program is partially just kept around for testing pursposes as well as backwards compatibility
! Everything else has been moved to flow_mod.f90 (not a great name, I know...)
program FLOW_runner
  use FLOW_mod
  !real (kind=8), dimension(44,100) :: test
  real, dimension(11,2) :: mthetatest
  real (kind=dkind), dimension(:,:), allocatable :: poloidal_viscosity
  !interface
  !  subroutine run_FLOW(mtheta,pol_visc)
  !    real, dimension(11,2), intent(in), optional :: mtheta
  !    real (kind=8), dimension(44,100), intent(out) :: pol_visc
  !  end subroutine
  !end interface
  mtheta = reshape((/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, &
          1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732, 1.732 /), shape(mthetatest))
  !call run_FLOW(mtheta=mthetatest,pol_visc=test)
  call run_FLOW()
  call allocation_clean_up
  !print *, pol_viscosity_data
end program
