module chemical_controller

  ! External modules
  use controller
  use precision
  use onestep

  implicit none

end module chemical_controller


! =========================== !
! Setup the ignition actuator !
! =========================== !
subroutine chemical_controller_setup

  ! Internal modules
  use chemical_controller

  ! External modules
  use string
  use parser
  use simulation_flags
  use combustion

  implicit none

  if (predictionOnly) return

  ! Ignition actuator requires ignition source to be used
  if (chemistryModel .ne. ONE_STEP)                                                          &
       call die('chemical_controller_setup: chemical actuator requires one-step chemistry!')

  ! Initialize the baseline ignition parameters
  nControlParameters = 1
  allocate(baselineValue(nControlParameters))
  baselineValue = controlParameter

  return
end subroutine chemical_controller_setup


! ============================= !
! Cleanup the ignition actuator !
! ============================= !
subroutine chemical_controller_cleanup

  ! Internal modules
  use chemical_controller

  implicit none

  ! Nothing to do here

  return
end subroutine chemical_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine chemical_controller_update_parameters(actuationAmount)

  ! Internal modules
  use chemical_controller

  ! External modules
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Local variables
  real(WP), parameter :: smallNumber = 1.0E-9_WP

  controlParameter = max( smallNumber, baselineValue(1) +                                    &
       real(gradientDirection, WP) * actuationAmount * controlGradient(1) )

  return
end subroutine chemical_controller_update_parameters


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine chemical_controller_compute_sensitivity

  ! Internal modules
  use chemical_controller

  ! External modules
  use geometry
  use solver_options
  use grid, only : coordinates
  use grid_functions, only : inner_product
  use state, only : adjointVariables

  implicit none

  ! Local variables
  integer :: i
  real(WP), dimension(:,:), allocatable :: F, chemicalSource

  ! Allocate temporary arrays
  allocate(F(nGridPoints, 2))
  allocate(chemicalSource(nGridPoints, nUnknowns))

  ! Compute the source terms
  chemicalSource = 0.0_WP
  call onestep_forward(chemicalSource)

  ! Compute sensitivities of the chemical source terms
  instantaneousSensitivity = 0.0_WP
  do i = nDimensions+2, nUnknowns
     F(:,1) = adjointVariables(:, i)! * controlMollifier(:,1)
     F(:,2) = chemicalSource(:, i) / controlParameter
     instantaneousSensitivity(1) = instantaneousSensitivity(1) + inner_product(F(:,1), F(:,2))
  end do

  ! Clean up
  deallocate(F)
  deallocate(chemicalSource)

  return
end subroutine chemical_controller_compute_sensitivity


! ========================= !
! Set the parameters bounds !
! ========================= !
subroutine chemical_controller_parameter_bound(lowerBound, upperBound, scale)

  ! Internal modules
  use perturbation_controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(out) :: lowerBound, upperBound, scale

  ! Local variables
  real(WP), parameter :: smallNumber = 1.0E-9_WP

  lowerBound = smallNumber
  upperBound = huge(1.0_WP)
  scale = 1.0_WP

  return
end subroutine chemical_controller_parameter_bound
