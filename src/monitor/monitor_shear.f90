module monitor_shear

  ! External modules
  use monitor

  implicit none

  real(WP) :: velocityDifference

end module monitor_shear


! ========================================== !
! Setup the routine to monitor a shear layer !
! ========================================== !
subroutine monitor_shear_setup

  ! Internal modules
  use monitor_shear

  ! External modules
  use parser
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  if (trim(simulationName) .ne. 'shear layer') return

  ! Read in the velocity difference
  call parser_read('shear layer velocity difference', velocityDifference)

  ! Set the monitor names
  call monitor_create('shear_layer', 3)
  call monitor_set_header(1, 'delta_m', 'r')
  call monitor_set_header(2, 'delta_w', 'r')
  call monitor_set_header(3, 'growth_rate', 'r')

  return
end subroutine monitor_shear_setup


! ============================== !
! Compute shear layer statistics !
! ============================== !
subroutine monitor_shear_timestep

  ! Internal modules
  use monitor_shear

  ! External modules
  use solver_options
  use state
  use state_functions

  implicit none

  ! Local variables
  real(WP) :: momentumThickness, vorticityThickness, growthRate

  if (trim(simulationName) .ne. 'shear layer') return

  ! Compute the shear layer statistics
  call compute_momentum_thickness(conservedVariables, velocityDifference, 2, 1,            &
       momentumThickness)
  call compute_vorticity_thickness(conservedVariables(:,1), conservedVariables(:,2),       &
       velocityDifference, 2, vorticityThickness)
  call compute_shear_growth_rate(conservedVariables, velocityDifference, 2, 1, growthRate)

  ! Set the shear layer parameters
  call monitor_select('shear_layer')
  call monitor_set_single_value(1, momentumThickness)
  call monitor_set_single_value(2, vorticityThickness)
  call monitor_set_single_value(3, growthRate)

  return
end subroutine monitor_shear_timestep
