module boundary

  ! External modules
  use precision
  use grid_patch

  implicit none
  
  real(WP) :: defaultInviscidPenaltyAmount, defaultViscousPenaltyAmount

end module boundary


! ========================== !
! Setup the boundary patches !
! ========================== !
subroutine boundary_setup

  ! Internal modules
  use boundary

  ! External modules
  use parser

  implicit none

  ! Read in default SAT parameters
  call parser_read('default inviscid penalty amount', defaultInviscidPenaltyAmount, 1.0_WP)
  call parser_read('default viscous penalty amount', defaultViscousPenaltyAmount, 1.0_WP)

  call sponge_setup
  call farfield_setup
  call outflow_setup
  call inflow_setup
  call slip_setup
  call isothermal_setup
  call adiabatic_setup
  call boundary_strong_setup
  call excitation_setup

  return
end subroutine boundary_setup


! ============================ !
! Cleanup the boundary patches !
! ============================ !
subroutine boundary_cleanup

  ! Internal modules
  use boundary

  implicit none

  call sponge_cleanup
  call farfield_cleanup
  call outflow_cleanup
  call inflow_cleanup
  call slip_cleanup
  call isothermal_cleanup
  call adiabatic_cleanup
  call boundary_strong_cleanup
  call excitation_cleanup

  return
end subroutine boundary_cleanup


! ======================== !
! Add bundary source terms !
! ======================== !
subroutine boundary_sources(mode, source)

  ! Internal modules
  use boundary

  ! External modules
  use solver_options
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Start the boundary timer
  call timing_start('boundary')

  select case (mode)

  case (FORWARD)
     call sponge_forward(source)
     call farfield_forward(source)
     call outflow_forward(source)
     call inflow_forward(source)
     call slip_forward(source)
     call isothermal_forward(source)
     call adiabatic_forward(source)
     call excitation_forward(time, source)

  case (ADJOINT)
     call sponge_adjoint(source)
     call farfield_adjoint(source)
     call outflow_adjoint(source)
     call inflow_adjoint(source)
     call slip_adjoint(source)
     call isothermal_adjoint(source)
     call adiabatic_adjoint(source)

  end select

  ! Stop the boundary timer
  call timing_stop('boundary')

  return
end subroutine boundary_sources


! ================================================ !
! Store the viscous fluxes on the boundary patches !
! ================================================ !
subroutine boundary_store_viscous_fluxes(fluxes)

  ! Internal modules
  use boundary

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns, nDimensions), intent(in) :: fluxes

  call farfield_store_viscous_fluxes(fluxes)
  call slip_store_viscous_fluxes(fluxes)

  return
end subroutine boundary_store_viscous_fluxes
