module ic_controller

  ! External modules
  use math, only : pi
  use controller
  use precision
  use solver_options

  implicit none

end module ic_controller


! ============================= !
! Setup the generic IC actuator !
! ============================= !
subroutine ic_controller_setup

  ! Internal modules
  use ic_controller

  ! External modules
  use string
  use parallel
  use parser
  use equation_of_state
  use random

  implicit none

  if (predictionOnly) return

  ! Make sure number of control parameters = 1
  if (nControlParameters .ne. 1)                                                             &
       call die('ic_controller_setup: number of control parameters must be 1!')

  return
end subroutine ic_controller_setup


! =============================== !
! Cleanup the generic IC actuator !
! =============================== !
subroutine ic_controller_cleanup

  ! Internal modules
  use ic_controller

  implicit none

  ! Nothing to do

  return
end subroutine ic_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine ic_controller_update_parameters(actuationAmount)

  ! Internal modules
  use ic_controller

  ! External modules

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Nothing to do
  
  return
end subroutine ic_controller_update_parameters


! ===================================== !
! Update initial conditions using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine ic_controller_update_ic(actuationAmount)

  ! Internal modules
  use ic_controller

  ! External modules
  use geometry
  use grid
  use state, only : conservedVariables, adjointVariables

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Local variables
  integer  :: i, j, k, l, gridIndex

  if (controllerPatch%nPatchPoints .gt. 0) then
     do l = 1, nUnknowns
        do k = controllerPatch%iStart(3), controllerPatch%iEnd(3)
           do j = controllerPatch%iStart(2), controllerPatch%iEnd(2)
              do i = controllerPatch%iStart(1), controllerPatch%iEnd(1)
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))

                 conservedVariables(gridIndex, l) = conservedVariables(gridIndex, l) +       &
                      real(gradientDirection, WP) * controlMollifier(gridIndex, 1) *         &
                      actuationAmount * adjointVariables(gridIndex, l)
              end do
           end do
        end do
     end do
  end if

  return
end subroutine ic_controller_update_ic


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine ic_controller_compute_sensitivity

  ! Internal modules
  use ic_controller

  implicit none

  instantaneousSensitivity = 0.0_WP

  return
end subroutine ic_controller_compute_sensitivity


! ======================================================== !
! Compute the adjoint sensitivity of the initial condition !
! ======================================================== !
subroutine ic_controller_compute_sensitivity_ic

  ! Internal modules
  use ic_controller

  ! External modules

  use parallel
  use geometry
  use solver_options
  use grid_functions, only : inner_product
  use state
  use time_info

  implicit none

  ! Local variables
  integer  :: i
  real(WP), allocatable :: F(:,:)

  ! Only measure sensitivity of initial condition
  if (timestep .ne. 0) return

  instantaneousSensitivity = 0.0_WP

  allocate(F(nGridPoints, 1))
  do i = 1, nUnknowns
     F(:,1) = adjointVariables(:,i) * controlMollifier(:,1)
     instantaneousSensitivity = instantaneousSensitivity + inner_product(F, F)
  end do
  deallocate(F)

  return
end subroutine ic_controller_compute_sensitivity_ic
