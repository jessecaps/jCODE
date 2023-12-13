module fuel_controller

  ! External modules
  use controller
  use precision

  implicit none

  integer :: fuelComponent
  real(WP) :: actuatorStart, actuatorDuration

end module fuel_controller


! ======================= !
! Setup the fuel actuator !
! ======================= !
subroutine fuel_controller_setup

  ! Internal modules
  use fuel_controller

  ! External modules
  use string
  use parser
  use solver_options
  use simulation_flags

  implicit none

  if (predictionOnly) return

  ! Make sure species are used
  if (nSpecies .le. 0)                                                                       &
       call die('fuel_controller_setup: fuel actuator requires nSpecies > 0!')

  ! Make sure number of control parameters = 1
  if (nControlParameters .ne. 1)                                                             &
       call die('fuel_controller_setup: number of control parameters must be 1!')

  call parser_read('fuel actuator time duration', actuatorDuration, -1.0_WP)
  if (actuatorDuration .gt. 0.0_WP)                                                          &
       call parser_read('fuel actuator time start', actuatorStart)

  ! Read the fuel component to actuate
  call parser_read('actuator fuel component', fuelComponent)
  if (fuelComponent .lt. 1 .or. fuelComponent .gt. nSpecies)                                 &
       call die('fuel_controller_setup: invalid fuel component!')

  ! Fuel actuator gradient is space-time dependent
  spaceTimeGradient = .true.

  return
end subroutine fuel_controller_setup


! ========================= !
! Cleanup the fuel actuator !
! ========================= !
subroutine fuel_controller_cleanup

  ! Internal modules
  use fuel_controller

  implicit none

  ! Nothing to do

  return
end subroutine fuel_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine fuel_controller_update_parameters(actuationAmount)

  ! Internal modules
  use fuel_controller

  ! External modules
  use math, only : pi
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Local variables
  real(wp) :: timePortion

  if (actuatorDuration .gt. 0.0_WP) then
     timePortion = exp( -0.5_WP * (time - actuatorStart)**2 / actuatorDuration**2 ) /        &
          actuatorDuration / sqrt(2.0_WP * pi)
  else
     timePortion = 1.0_WP
  end if

  iGradientBuffer = iGradientBuffer - 1

  if (iGradientBuffer .eq. size(gradientBuffer, 3)) call load_sensitivity_gradient

  if (controllerPatch%nPatchPoints .gt. 0) then
     controlForcing(:,:) = 0.0_WP
     controlForcing(:,nDimensions+2+fuelComponent) = - timePortion * actuationAmount *       &
          gradientBuffer(:,1,iGradientBuffer)
  end if

  if (iGradientBuffer .eq. 1) iGradientBuffer = size(gradientBuffer, 3) + 1

  return
end subroutine fuel_controller_update_parameters


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine fuel_controller_compute_sensitivity

  ! Internal modules
  use fuel_controller

  ! External modules
  use math, only : pi
  use geometry
  use grid_functions, only : inner_product
  use state, only : adjointVariables
  use time_info

  implicit none

  ! Local variables
  real(WP) :: timePortion
  real(WP), allocatable :: F(:,:)

  iGradientBuffer = iGradientBuffer + 1

  if (actuatorDuration .gt. 0.0_WP) then
     timePortion = exp( -0.5_WP * (time - actuatorStart)**2 /actuatorDuration**2 ) /         &
          actuatorDuration / sqrt(2.0_WP * pi)
  else
     timePortion = 1.0_WP
  end if

  allocate(F(nGridPoints, 1))
  F(:,1) = adjointVariables(:,nDimensions+2+fuelComponent) * controlMollifier(:,1) *         &
       timePortion
  instantaneousSensitivity = inner_product(F, F)
  deallocate(F)

  if (controllerPatch%nPatchPoints .gt. 0) then
     allocate(F(controllerPatch%nPatchPoints, 2))
     call patch_collect(controllerPatch, adjointVariables(:,nDimensions+2+fuelComponent),    &
          F(:,1))
     call patch_collect(controllerPatch, controlMollifier(:,1), F(:,2))
     F(:,2) = F(:,2) * timePortion
     gradientBuffer(:, 1, iGradientBuffer) = product(F, dim = 2)
     deallocate(F)
  end if

  return
end subroutine fuel_controller_compute_sensitivity
