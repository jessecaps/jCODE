module monitor_cases

  ! External modules
  use monitor

  implicit none

end module monitor_cases


! ============================================== !
! Setup the monitor for various simulation cases !
! ============================================== !
subroutine monitor_cases_setup(mode, controlIteration)

  ! Internal modules
  use monitor

  ! External modules
  use solver_options

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  if (mode.eq.ADJOINT .or. controlIteration.gt.0) return

  call monitor_shear_setup
  call monitor_rayleigh_taylor_setup
  call monitor_impulse_setup
  call monitor_homogeneous_setup
  call monitor_shock_tube_setup
  call monitor_drag_setup

  return
end subroutine monitor_cases_setup


! =============================================== !
! Monitor a timestep for various simulation cases !
! =============================================== !
subroutine monitor_cases_timestep(controlIteration)

  ! Internal modules
  use monitor

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  if (controlIteration.gt.0) return

  call monitor_shear_timestep
  call monitor_rayleigh_taylor_timestep
  call monitor_impulse_timestep
  call monitor_homogeneous_timestep
  call monitor_shock_tube_timestep
  call monitor_drag_timestep

  return
end subroutine monitor_cases_timestep
