module momentum_controller

  ! External modules
  use controller
  use precision

  implicit none

  integer :: momentumComponent
  real(WP) :: actuatorStart, actuatorDuration

end module momentum_controller


! =========================== !
! Setup the momentum actuator !
! =========================== !
subroutine momentum_controller_setup

  ! Internal modules
  use momentum_controller

  ! External modules
  use string
  use parser
  use simulation_flags
  use geometry

  implicit none

  if (predictionOnly) return

  ! Make sure number of control parameters = 1
  if (nControlParameters .ne. 1)                                                             &
       call die('momentum_controller_setup: number of control parameters must be 1!')

  call parser_read('momentum actuator time duration', actuatorDuration, -1.0_WP)
  if (actuatorDuration .gt. 0.0_WP)                                                          &
       call parser_read('momentum actuator time start', actuatorStart)

  ! Read the momentum component to actuate
  call parser_read('actuator momentum component', momentumComponent)
  if (momentumComponent .lt. 1 .or. momentumComponent .gt. nDimensions)                      &
       call die('momentum_controller_setup: invalid momentum component!')

  ! Momentum actuator gradient is space-time dependent
  spaceTimeGradient = .true.

  return
end subroutine momentum_controller_setup


! ============================= !
! Cleanup the momentum actuator !
! ============================= !
subroutine momentum_controller_cleanup

  ! Internal modules
  use momentum_controller

  implicit none

  ! Nothing to do

  return
end subroutine momentum_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine momentum_controller_update_parameters(actuationAmount)

  ! Internal modules
  use momentum_controller

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
     controlForcing(:,momentumComponent+1) = - timePortion * actuationAmount *               &
          gradientBuffer(:,1,iGradientBuffer)
  end if

  if (iGradientBuffer .eq. 1) iGradientBuffer = size(gradientBuffer, 3) + 1

  return
end subroutine momentum_controller_update_parameters


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine momentum_controller_compute_sensitivity

  ! Internal modules
  use momentum_controller

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
  F(:,1) = adjointVariables(:,momentumComponent + 1) * controlMollifier(:,1) * timePortion
  instantaneousSensitivity = inner_product(F, F)
  deallocate(F)

  if (controllerPatch%nPatchPoints .gt. 0) then
     allocate(F(controllerPatch%nPatchPoints, 2))
     call patch_collect(controllerPatch, adjointVariables(:,momentumComponent + 1), F(:,1))
     call patch_collect(controllerPatch, controlMollifier(:,1), F(:,2))
     F(:,2) = F(:,2) * timePortion
     gradientBuffer(:, 1, iGradientBuffer) = product(F, dim = 2)
     deallocate(F)
  end if

  return
end subroutine momentum_controller_compute_sensitivity
