module pressure_gradient

  ! External modules
  use precision

  implicit none

  real(WP) :: pressureGradient(3)
  logical :: usePressureGradient

end module pressure_gradient


! =========================== !
! Setup the pressure gradient !
! =========================== !
subroutine pressure_gradient_setup

  ! Internal modules
  use pressure_gradient

  ! External modules
  use parser

  implicit none

  pressureGradient = 0.0_WP
  call parser_is_defined('pressure gradient', usePressureGradient)
  if (usePressureGradient) call parser_read('pressure gradient', pressureGradient)

  return
end subroutine pressure_gradient_setup


! ============================= !
! Cleanup the pressure gradient !
! ============================= !
subroutine pressure_gradient_cleanup

  ! Internal modules
  use pressure_gradient

  implicit none

  return
end subroutine pressure_gradient_cleanup


! ==================================================== !
! Add the the pressure gradient during the forward run !
! ==================================================== !
subroutine pressure_gradient_forward(source)

  ! Internal modules
  use pressure_gradient

  ! External modules
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  if (.not. usePressureGradient) return

  do i = 1, nGridPoints
     source(i, 2:nDimensions+1) = source(i, 2:nDimensions+1) + pressureGradient(1:nDimensions)
     source(i, nDimensions+2) = source(i, nDimensions+2) +                                   &
          sum(velocity(i,1:nDimensions) * pressureGradient(1:nDimensions))
  end do
 
  return
end subroutine pressure_gradient_forward


! ==================================================== !
! Add the the pressure gradient during the adjoint run !
! ==================================================== !
subroutine pressure_gradient_adjoint(source)

  ! Internal modules
  use pressure_gradient

  ! External modules
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i
  real(WP), allocatable :: localSourceJacobian(:,:), temp(:)

  if (.not. usePressureGradient) return

  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  localSourceJacobian = 0.0_WP 

  do i = 1, nGridPoints
     localSourceJacobian(nDimensions+2,2:nDimensions+1) = pressureGradient(1:nDimensions) *  &
          specificVolume(i,1)
     temp = matmul(transpose(localSourceJacobian), adjointVariables(i,:))
     source(i,:) = source(i,:) - temp
  end do

  deallocate(localSourceJacobian)
  deallocate(temp)
 
  return
end subroutine pressure_gradient_adjoint
