module thermal_controller

  ! External modules
  use controller
  use precision

  implicit none

  real(WP) :: actuatorStart, actuatorDuration

end module thermal_controller


! ========================== !
! Setup the thermal actuator !
! ========================== !
subroutine thermal_controller_setup

  ! Internal modules
  use thermal_controller

  ! External modules
  use string
  use parser
  use simulation_flags

  implicit none

  if (predictionOnly) return

  ! Make sure number of control parameters = 1
  if (nControlParameters .ne. 1)                                                             &
       call die('thermal_controller_setup: number of control parameters must be 1!')

  call parser_read('thermal actuator time duration', actuatorDuration, -1.0_WP)
  if (actuatorDuration .gt. 0.0_WP)                                                          &
       call parser_read('thermal actuator time start', actuatorStart)

  ! Thermal actuator gradient is space-time dependent
  spaceTimeGradient = .true.

  ! Initialize the baseline ignition parameters
  allocate(baselineValue(1))
  baselineValue = 0.0_WP

  return
end subroutine thermal_controller_setup


! ============================ !
! Cleanup the thermal actuator !
! ============================ !
subroutine thermal_controller_cleanup

  ! Internal modules
  use thermal_controller

  implicit none

  ! Nothing to do

  return
end subroutine thermal_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine thermal_controller_update_parameters(actuationAmount)

  ! Internal modules
  use thermal_controller

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
     controlForcing(:,nDimensions+2) = - timePortion * actuationAmount *                     &
          gradientBuffer(:,1,iGradientBuffer)
  end if

  if (iGradientBuffer .eq. 1) iGradientBuffer = size(gradientBuffer, 3) + 1

  return
end subroutine thermal_controller_update_parameters


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine thermal_controller_compute_sensitivity

  ! Internal modules
  use thermal_controller

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
  F(:,1) = adjointVariables(:,nDimensions+2) * controlMollifier(:,1) * timePortion
  instantaneousSensitivity = inner_product(F, F)
  deallocate(F)

  if (controllerPatch%nPatchPoints .gt. 0) then
     allocate(F(controllerPatch%nPatchPoints, 2))
     call patch_collect(controllerPatch, adjointVariables(:,nDimensions+2), F(:,1))
     call patch_collect(controllerPatch, controlMollifier(:,1), F(:,2))
     F(:,2) = F(:,2) * timePortion
     gradientBuffer(:, 1, iGradientBuffer) = product(F, dim = 2)
     deallocate(F)
  end if

  return
end subroutine thermal_controller_compute_sensitivity
