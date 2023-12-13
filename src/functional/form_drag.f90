module form_drag

  ! External modules
  use functional

  implicit none

  real(WP) :: dragDirection(3), inviscidPenaltyAmount

contains

  subroutine verify_form_drag_patch(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) /= extent((i-1)*2+2)) call die('verify_form_drag_patch: &
         &Extends more than 1 grid point along normal direction!')

    if ((patch%normalDirection .gt. 0 .and. extent((i-1)*2+1) .ne. 1) .or.                   &
         (patch%normalDirection .lt. 0 .and. extent((i-1)*2+2) .ne. globalGridSize(i)))            &
         call die('verify_form_drag_patch: Not aligned with a computational boundary!')

    return
  end subroutine verify_form_drag_patch

end module form_drag


! ======================================== !
! Setup the acoustic noise cost functional !
! ======================================== !
subroutine form_drag_setup

  ! Internal modules
  use form_drag

  ! External modules
  use parser
  use first_derivative

  implicit none

  ! Verify the patch type
  call verify_form_drag_patch(functionalPatch)

  ! Get the drag direction
  dragDirection = 0.0_WP

  call parser_read('drag direction x', dragDirection(1), 1.0_WP)
  if (nDimensions .ge. 2) call parser_read('drag direction y', dragDirection(2), 0.0_WP)
  if (nDimensions .eq. 3) call parser_read('drag direction z', dragDirection(3), 0.0_WP)

  if (sum(dragDirection ** 2) .le. epsilon(0.0_WP)) call die('form_drag_setup: &
       &Unable to determine a unit vector for computing pressure drag!')

  dragDirection = dragDirection / sqrt(sum(dragDirection ** 2))

  ! Inviscid penalty amount
  call parser_read(trim(functionalPatch%name) // 'inviscid_penalty_amount',                  &
       inviscidPenaltyAmount, 2.0_WP)
  inviscidPenaltyAmount = sign(inviscidPenaltyAmount,                                        &
       real(functionalPatch%normalDirection, WP))
  inviscidPenaltyAmount = inviscidPenaltyAmount /                                            &
       firstDerivative(abs(functionalPatch%normalDirection))%normBoundary(1)

  return
end subroutine form_drag_setup


! ========================================== !
! Cleanup the acoustic noise cost functional !
! ========================================== !
subroutine form_drag_cleanup

  ! Internal modules
  use form_drag

  implicit none

  ! Nothiing to do

  return
end subroutine form_drag_cleanup


! ========================================== !
! Compute the acoustic noise cost functional !
! ========================================== !
subroutine form_drag_compute

  ! Internal modules
  use form_drag

  ! External modules
  use geometry
  use solver_options
  use grid
  use grid_functions, only : inner_product
  use first_derivative
  use state, only : pressure

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: normBoundaryFactor
  real(WP), allocatable :: F(:,:)

  i = abs(functionalPatch%normalDirection)
  normBoundaryFactor = 1.0_WP / firstDerivative(i)%normBoundary(1)

  allocate(F(nGridPoints, 2))
  F(:,1) = - (pressure(:,1) - 1.0_WP / ratioOfSpecificHeats)
  F(:,2) = matmul(metrics(:,1+nDimensions*(i-1):nDimensions*i),                              &
       dragDirection(1:nDimensions)) * normBoundaryFactor
  instantaneousCostFunctional = inner_product(F(:,1), F(:,2))
  deallocate(F)

  return
end subroutine form_drag_compute


! ========================================== !
! Compute the acoustic noise adjoint forcing !
! ========================================== !
subroutine form_drag_adjoint_source(source)

  ! Internal modules
  use form_drag

  ! External modules
  use state_jacobian
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use first_derivative
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, direction, gridIndex, patchIndex
  real(WP) :: forcingFactor, normBoundaryFactor, F
  real(WP), allocatable :: localConservedVariables(:), metricsAlongNormalDirection(:),       &
       unitNormal(:), incomingJacobianOfInviscidFlux(:,:)

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  direction = abs(functionalPatch%normalDirection)
  normBoundaryFactor = sign(1.0_WP / firstDerivative(direction)%normBoundary(1),             &
       real(functionalPatch%normalDirection, WP))

  allocate(localConservedVariables(nUnknowns))
  allocate(unitNormal(nDimensions))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           patchIndex = i - functionalPatch%offset(1) +                                      &
                functionalPatch%localSize(1) * (j - 1 - functionalPatch%offset(2) +          &
                functionalPatch%localSize(2) * (k - 1 - functionalPatch%offset(3)))

           localConservedVariables = conservedVariables(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                metrics(gridIndex,1+nDimensions*(direction-1):                               &
                nDimensions*direction)
           unitNormal = metricsAlongNormalDirection /                                        &
                sqrt(sum(metricsAlongNormalDirection ** 2))

           if (useContinuousAdjoint) then

              F = forcingFactor * jacobian(gridIndex, 1) *                                   &
                   dot_product(adjointVariables(gridIndex,2:nDimensions+1) -                 &
                   sign(dragDirection(1:nDimensions),                                        &
                   real(functionalPatch%normalDirection, WP)), unitNormal)

              call compute_incoming_jacobian_of_inviscid_flux(localConservedVariables,       &
                   metricsAlongNormalDirection, - functionalPatch%normalDirection,           &
                   incomingJacobianOfInviscidFlux,                                           &
                   specificVolume = specificVolume(gridIndex, 1),                            &
                   pressure = pressure(gridIndex, 1))

              source(gridIndex,:) = source(gridIndex,:) - inviscidPenaltyAmount * F *        &
                   matmul(transpose(incomingJacobianOfInviscidFlux(2:nDimensions+1,:)),      &
                   unitNormal)

           else

              F = forcingFactor * jacobian(gridIndex, 1) * normBoundaryFactor *              &
                   (ratioOfSpecificHeats - 1.0_WP) *                                         &
                   dot_product(metricsAlongNormalDirection, dragDirection(1:nDimensions))

              source(gridIndex,1) = source(gridIndex,1) +                                    &
                   0.5_WP * sum(velocity(gridIndex,:) ** 2) * F
              source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) -        &
                   velocity(gridIndex,:) * F
              source(gridIndex,nDimensions+2) = source(gridIndex,nDimensions+2) +  F

           end if

        end do
     end do
  end do

  deallocate(localConservedVariables)
  deallocate(unitNormal)
  deallocate(metricsAlongNormalDirection)
  deallocate(incomingJacobianOfInviscidFlux)

  return
end subroutine form_drag_adjoint_source
