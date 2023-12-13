module drag_force

  ! External modules
  use functional

  implicit none

  real(WP) :: dragDirection(3)

contains

  subroutine verify_drag_force_patch(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die('verify_drag_force_patch: Normal direction is invalid!')

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die('verify_drag_force_patch: &
         &Extends more than 1 grid point along normal direction!')

    return
  end subroutine verify_drag_force_patch

end module drag_force


! ==================================== !
! Setup the drag force cost functional !
! ==================================== !
subroutine drag_force_setup

  ! Internal modules
  use drag_force

  ! External modules
  use parser
  use first_derivative

  implicit none

  ! Verify the patch type
  call verify_drag_force_patch(functionalPatch)

  ! Get the drag direction
  dragDirection = 0.0_WP

  call parser_read('drag direction x', dragDirection(1), 1.0_WP)
  if (nDimensions .ge. 2) call parser_read('drag direction y', dragDirection(2), 0.0_WP)
  if (nDimensions .eq. 3) call parser_read('drag direction z', dragDirection(3), 0.0_WP)

  if (sum(dragDirection ** 2) .le. epsilon(0.0_WP)) call die('drag_force_setup:              &
       &Unable to determine a unit vector for computing drag force!')

  dragDirection = dragDirection / sqrt(sum(dragDirection ** 2))

  return
end subroutine drag_force_setup


! ====================================== !
! Cleanup the drag force cost functional !
! ====================================== !
subroutine drag_force_cleanup

  ! Internal modules
  use drag_force

  implicit none

  ! Nothiing to do

  return
end subroutine drag_force_cleanup


! ====================================== !
! Compute the drag force cost functional !
! ====================================== !
subroutine drag_force_compute

  ! Internal modules
  use drag_force

  ! External modules
  use geometry
  use simulation_flags
  use solver_options
  use grid
  use grid_functions, only : inner_product
  use first_derivative
  use state, only : stressTensor

  implicit none

  ! Local variables
  integer :: i, j
  real(WP) :: normBoundaryFactor
  real(WP), allocatable :: F(:,:)

  i = abs(functionalPatch%normalDirection)
  normBoundaryFactor = 1.0_WP / firstDerivative(i)%normBoundary(1)

  allocate(F(nGridPoints, 1))

  F(:,1) = 0.0_WP
  do j = 1, nDimensions
     if (useViscosity) then
        F(:,1) = F(:,1) + dragDirection(j) *                                                 &
             sum(metrics(:,1+nDimensions*(i-1):nDimensions*i) *                              &
             stressTensor(:,1+nDimensions*(j-1):nDimensions*j), dim = 2)
     end if
  end do
  F(:,1) = normBoundaryFactor * F(:,1)

  instantaneousCostFunctional = inner_product(F(:,1), targetMollifier(:,1))

  deallocate(F)

  return
end subroutine drag_force_compute


! ====================================== !
! Compute the drag force adjoint forcing !
! ====================================== !
subroutine drag_force_adjoint_source(source)

  ! Internal modules
  use drag_force

  ! External modules
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_patch
  use first_derivative
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: direction
  real(WP) :: forcingFactor, normBoundaryFactor
  real(WP), allocatable :: temp1(:,:), temp2(:,:)

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  direction = abs(functionalPatch%normalDirection)
  normBoundaryFactor = 1.0_WP / firstDerivative(direction)%normBoundary(1)

  if (functionalPatch%nPatchPoints .gt. 0) then

     allocate(temp1(nGridPoints, nUnknowns))
     allocate(temp2(nGridPoints, 1))

     ! Hack for TBL:

     temp1 = 0.0_WP

     temp2(:,1) = jacobian(:,1) * metrics(:,1) * dynamicViscosity(:,1)
     call adjoint_first_derivative_project_boundary_and_apply(1, temp2,                      &
          functionalPatch%normalDirection)
     call first_derivative_apply_norm(1, temp2)
     temp1(:,3) = jacobian(:,1) * specificVolume(:,1) * temp2(:,1)

     temp2(:,1) = jacobian(:,1) * metrics(:,5) * dynamicViscosity(:,1)
     call adjoint_first_derivative_project_boundary_and_apply(2, temp2,                      &
          functionalPatch%normalDirection)
     call first_derivative_apply_norm(1, temp2)
     temp1(:,2) = jacobian(:,1) * specificVolume(:,1) * temp2(:,1)

     source = source + forcingFactor * normBoundaryFactor * temp1

     deallocate(temp2)
     deallocate(temp1)

  end if

  return
end subroutine drag_force_adjoint_source
