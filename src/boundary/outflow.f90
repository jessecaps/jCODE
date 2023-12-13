module outflow

  ! External modules
  use boundary
  use precision

  implicit none

  integer :: nOutflows

  type, private :: t_Outflow
     real(WP) :: inviscidPenaltyAmount
  end type t_Outflow

  type(t_Outflow), allocatable :: outflowData(:)
  type(t_Patch), pointer :: outflowPatch(:)

contains

  subroutine setup_outflow_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use geometry
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Outflow), intent(inout) :: bc

    ! Local variables
    integer :: i, j, k, gridIndex, patchIndex
    real(WP) :: radius, holeRadius, holePosition(3)
    logical :: puncture, inverted
    character(len = str_medium) :: holeShape

    ! Verify the patch type
    call verify_outflow_patch(patch)

    ! Inviscid penalty amount
    call parser_read(trim(patch%name) // ' inviscid penalty amount',                         &
         bc%inviscidPenaltyAmount, defaultInviscidPenaltyAmount)
    bc%inviscidPenaltyAmount = sign(bc%inviscidPenaltyAmount, real(patch%normalDirection, WP))
    if (useUpwinding) then
       bc%inviscidPenaltyAmount = bc%inviscidPenaltyAmount /                                 &
            upwindLeft(abs(patch%normalDirection))%normBoundary(1)
    else
       bc%inviscidPenaltyAmount = bc%inviscidPenaltyAmount /                                 &
            firstDerivative(abs(patch%normalDirection))%normBoundary(1)
    end if

    ! Add a hole to the patch
    call parser_read(trim(patch%name) // ' include hole', puncture, .false.)
    if (puncture) then

       ! Initialize the hole
       allocate(patch%hole(patch%nPatchPoints))
       patch%hole = 0

       ! Hole properties
       call parser_read(trim(patch%name) // ' hole is inverted', inverted, .false.)
       call parser_read(trim(patch%name) // ' hole shape', holeShape)

       select case(trim(holeShape))

       case('circle')

          call parser_read(trim(patch%name) // ' hole radius', holeRadius)
          holePosition = 0.0_WP
          do i = 1, nDimensions
             call parser_read(trim(patch%name) // ' hole position', holePosition(1:nDimensions))
          end do

          do k = patch%iStart(3), patch%iEnd(3)
             do j = patch%iStart(2), patch%iEnd(2)
                do i = patch%iStart(1), patch%iEnd(1)
                   gridIndex = i - gridOffset(1) + localGridSize(1) *                        &
                        (j - 1 - gridOffset(2) + localGridSize(2) *                          &
                        (k - 1 - gridOffset(3)))
                   patchIndex = i - patch%offset(1) + patch%localSize(1) *                   &
                        (j - 1 - patch%offset(2) + patch%localSize(2) *                      &
                        (k - 1 - patch%offset(3)))

                   radius = sqrt(sum((coordinates(gridIndex, 1:nDimensions) -                &
                        holePosition(1:nDimensions)) ** 2))
                   if (.not.inverted .and. radius .le. holeRadius) then
                      patch%hole(patchIndex) = 1
                   else if (inverted .and. radius .gt. holeRadius) then
                      patch%hole(patchIndex) = 1
                   end if

                end do
             end do
          end do

       case default

          call die("setup_outflow_patch: Unknown hole shape on patch " //                    &
               trim(patch%name) // " '" // trim(holeShape) // "' !")

       end select

    end if

    return
  end subroutine setup_outflow_patch


  subroutine cleanup_outflow_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Outflow), intent(inout) :: bc

    return
  end subroutine cleanup_outflow_patch


  subroutine verify_outflow_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    if (patch%patchType .ne. SAT_OUTFLOW)                                                    &
         call die('verify_outflow_patch: Patch type mismatch!')

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_outflow_patch: Normal direction is invalid for '" //               &
         trim(patch%name) // "'!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_outflow_patch: Invalid extent on '" // trim(patch%name) // "'!")
    end do

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die("verify_outflow_patch: '" //     &
         trim(patch%name) //  "' extends more than 1 grid point along normal direction!")

    return
  end subroutine verify_outflow_patch


  subroutine outflow_patch_forward(patch, bc, source)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Outflow), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, incomingDirection, gridIndex, gridIndex2, patchIndex,  &
         ijk(3)
    real(WP) :: localConservedVariables(nUnknowns), localTargetState(nUnknowns),             &
         rho_, velocity_(3), metricsAlongNormalDirection(nDimensions),                       &
         incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns)

    direction = abs(patch%normalDirection)
    incomingDirection = patch%normalDirection
    incomingJacobianOfInviscidFlux = 0.0_WP

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))
             if (allocated(patch%hole)) then
                if (patch%hole(patchIndex) .eq. 1) cycle
             end if

             ! Get the index of the neighboring grid point
             ijk(1) = i; ijk(2) = j; ijk(3) = k
             ijk(direction) = ijk(direction) + sign(1, patch%normalDirection)
             gridIndex2 = ijk(1) - gridOffset(1) + localGridSize(1) *                        &
                  (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                           &
                  (ijk(3) - 1 - gridOffset(3)))

             ! Get the neighboring velocity and clip to ensure outflow
             velocity_(1:nDimensions) = velocity(gridIndex2, 1:nDimensions)
             if (patch%normalDirection .lt. 0) then
                velocity_(direction) = max(velocity_(direction), 0.0_WP)
             else
                velocity_(direction) = min(velocity_(direction), 0.0_WP)
             end if

             ! Set the target state
             rho_ = conservedVariables(gridIndex2, 1)
             if (twoWayCoupling) rho_ = rho_ / volumeFraction(gridIndex2, 1)
             localTargetState(1) = rho_
             localTargetState(2:nDimensions+1) = rho_ * velocity_(1:nDimensions)
             localTargetState(nDimensions+2) = pressure(gridIndex2, 1) /                     &
                  (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * rho_ *                          &
                  sum(velocity_(1:nDimensions)**2)
             do l = 1, nSpecies
                localConservedVariables(nDimensions+2+l) = rho_ *                            &
                     massFraction(gridIndex2, l)
             end do
             
             localConservedVariables = conservedVariables(gridIndex,:)
             if (twoWayCoupling) localConservedVariables = localConservedVariables /         &
                  volumeFraction(gridIndex, 1)
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             call compute_incoming_jacobian_of_inviscid_flux(localTargetState,               &
                  metricsAlongNormalDirection, incomingDirection,                            &
                  incomingJacobianOfInviscidFlux)

             source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *          &
                  jacobian(gridIndex, 1) * matmul(incomingJacobianOfInviscidFlux,            &
                  localConservedVariables - localTargetState)

          end do
       end do
    end do

    return
  end subroutine outflow_patch_forward


  subroutine outflow_patch_adjoint(patch, bc, source)

    ! External modules
    use state_jacobian
    use simulation_flags
    use solver_options
    use geometry
    use grid
    use state
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Outflow), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, incomingDirection, gridIndex, patchIndex, gridIndex2,  &
         ijk(3)
    real(WP) :: rho_, velocity_(3)
    real(WP), allocatable :: localTargetState(:), metricsAlongNormalDirection(:),            &
         metricsAlongDirection2(:), incomingJacobianOfInviscidFlux(:,:),                     &
         localMassFraction(:)

    direction = abs(patch%normalDirection)
    if (useContinuousAdjoint) then
       incomingDirection = -patch%normalDirection
    else
       incomingDirection = +patch%normalDirection
    end if

    allocate(localTargetState(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))
             if (allocated(patch%hole)) then
                if (patch%hole(patchIndex) .eq. 1) cycle
             end if

             ! Get the index of the neighboring grid point
             ijk(1) = i; ijk(2) = j; ijk(3) = k
             ijk(direction) = ijk(direction) + sign(1, patch%normalDirection)
             gridIndex2 = ijk(1) - gridOffset(1) + localGridSize(1) *                        &
                  (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                           &
                  (ijk(3) - 1 - gridOffset(3)))

             ! Get the neighboring velocity and clip to ensure outflow
             velocity_(1:nDimensions) = velocity(gridIndex2, 1:nDimensions)
             if (patch%normalDirection .lt. 0) then
                velocity_(direction) = max(velocity_(direction), 0.0_WP)
             else
                velocity_(direction) = min(velocity_(direction), 0.0_WP)
             end if

             ! Set the target state
             rho_ = conservedVariables(gridIndex2, 1)
             localTargetState(1) = rho_
             localTargetState(2:nDimensions+1) = rho_ * velocity_(1:nDimensions)
             localTargetState(nDimensions+2) = pressure(gridIndex2, 1) /                     &
                  (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * rho_ *                          &
                  sum(velocity_(1:nDimensions)**2)
             do l = 1, nSpecies
                localTargetState(nDimensions+2+l) =                                          &
                     conservedVariables(gridIndex2,nDimensions+2+l)
             end do
             
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             call compute_incoming_jacobian_of_inviscid_flux(localTargetState,               &
                  metricsAlongNormalDirection, incomingDirection,                            &
                  incomingJacobianOfInviscidFlux)

             if (useContinuousAdjoint) then
                source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *       &
                     jacobian(gridIndex, 1) *                                                &
                     matmul(transpose(incomingJacobianOfInviscidFlux),                       &
                     adjointVariables(gridIndex,:))
             else
                source(gridIndex,:) = source(gridIndex,:) + bc%inviscidPenaltyAmount *       &
                     jacobian(gridIndex, 1) *                                                &
                     matmul(transpose(incomingJacobianOfInviscidFlux),                       &
                     adjointVariables(gridIndex,:))
             end if

          end do
       end do
    end do

    deallocate(incomingJacobianOfInviscidFlux)
    deallocate(metricsAlongNormalDirection)
    deallocate(localTargetState)
    if (allocated(metricsAlongDirection2)) deallocate(metricsAlongDirection2)
    if (allocated(localMassFraction)) deallocate(localMassFraction)

    return
  end subroutine outflow_patch_adjoint

end module outflow


! ================================== !
! Setup the outflow boundary patches !
! ================================== !
subroutine outflow_setup

  ! Internal modules
  use outflow

  ! External modules
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use state_functions
  use state

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nOutflows = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_OUTFLOW) then
        nOutflows = nOutflows + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nOutflows .eq. 0) return

  ! Allocate the outflow type
  allocate(outflowData(nOutflows))

  ! Connect the boundary patch
  outflowPatch => patches(j:j+nOutflows-1)

  ! Setup the boundary conditions
  do i = 1, nOutflows
     call setup_outflow_patch(outflowPatch(i), outflowData(i))
  end do

  return
end subroutine outflow_setup


! ==================================== !
! Cleanup the outflow boundary patches !
! ==================================== !
subroutine outflow_cleanup

  ! Internal modules
  use outflow

  ! Local variables
  integer :: i

  if (nOutflows .gt. 0) then
     do i = 1, nOutflows
        call cleanup_outflow_patch(outflowPatch(i), outflowData(i))
     end do
     deallocate(outflowData)
     nullify(outflowPatch)
  end if

  nOutflows = 0

  return
end subroutine outflow_cleanup


! ============================================= !
! Add the outflow source during the forward run !
! ============================================= !
subroutine outflow_forward(source)

  ! Internal modules
  use outflow

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nOutflows
     call outflow_patch_forward(outflowPatch(i), outflowData(i), source)
  end do

  return
end subroutine outflow_forward


! ============================================= !
! Add the outflow source during the adjoint run !
! ============================================= !
subroutine outflow_adjoint(source)

  ! Internal modules
  use outflow

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nOutflows
     call outflow_patch_adjoint(outflowPatch(i), outflowData(i), source)
  end do

  return
end subroutine outflow_adjoint
