module inflow

  ! External modules
  use boundary
  use precision

  implicit none

  integer :: nInflowWalls

  type :: t_InflowWall
     real(WP) :: inviscidPenaltyAmount
     real(WP), allocatable :: targetPenalty(:,:)
  end type t_InflowWall

  type(t_InflowWall), allocatable :: inflowWallData(:)
  type(t_Patch), pointer :: inflowWallPatch(:)

contains

  subroutine setup_inflow_patch(patch, bc)

    ! External modules
    use parser
    use geometry
    use simulation_flags
    use first_derivative
    use equation_of_state
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_InflowWall), intent(inout) :: bc

    ! Local variables
    integer :: i, j, k, l, gridIndex, patchIndex, direction
    real(WP) :: radius, holeRadius, holePosition(3), normalMomentum
    real(WP), allocatable :: localTargetState(:), metricsAlongNormalDirection(:),            &
         targetPressure(:), targetVelocity(:,:), targetMassFraction(:,:)
    logical :: puncture, inverted
    character(len = str_medium) :: holeShape

    ! Verify the patch type
    call verify_inflow_patch(patch)

    ! Get SAT penalty amounts
    call parser_read(trim(patch%name) //  ' inviscid penalty amount',                        &
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
             call parser_read(trim(patch%name) // ' hole position',                          &
                  holePosition(1:nDimensions))
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

          call die("setup_inflow_patch: Unknown hole shape on patch " //               &
               trim(patch%name) // " '" // trim(holeShape) // "' !")

       end select

    end if

    ! Inflow requires a target state
    if (.not. useTargetState) call die('setup_inflow_patch: target state required!')

    ! Store target penalty
    allocate(bc%targetPenalty(patch%nPatchPoints, nUnknowns))
    allocate(localTargetState(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(targetPressure(nGridPoints))
    allocate(targetVelocity(nGridPoints, nDimensions))
    if (nSpecies .gt. 0) then
       allocate(targetMassFraction(nGridPoints, nSpecies))
       call compute_dependent_variables(targetState, velocity = targetVelocity,              &
            pressure = targetPressure, massFraction = targetMassFraction)
    else
       call compute_dependent_variables(targetState, velocity = targetVelocity,        &
            pressure = targetPressure)
    end if
    direction = abs(patch%normalDirection)
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

             localTargetState = targetState(gridIndex,:)
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             normalMomentum = dot_product(localTargetState(2:nDimensions+1),                 &
                  metricsAlongNormalDirection)

             bc%targetPenalty(patchIndex, 1) = normalMomentum
             bc%targetPenalty(patchIndex, 2:nDimensions+1) = normalMomentum *                &
                  targetVelocity(gridIndex,:)
             bc%targetPenalty(patchIndex, nDimensions+2) =                                   &
                  normalMomentum / targetState(gridIndex, 1) *                               &
                  (localTargetState(nDimensions+2) + targetPressure(gridIndex))
             do l = 1, nSpecies 
                bc%targetPenalty(patchIndex, nDimensions+2+l) = normalMomentum *                &
                     targetMassFraction(gridIndex, l)
             end do
          end do
       end do
    end do

    ! Cleanup
    deallocate(targetVelocity, targetPressure, localTargetState, metricsAlongNormalDirection)
    if (nSpecies .gt. 0) then
       deallocate(targetMassFraction)
    end if

    return
  end subroutine setup_inflow_patch


  subroutine cleanup_inflow_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_InflowWall), intent(inout) :: bc

    if (allocated(bc%targetPenalty)) deallocate(bc%targetPenalty)

    return
  end subroutine cleanup_inflow_patch


  subroutine verify_inflow_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_inflow_patch: Normal direction is invalid for '" //          &
         trim(patch%name) // "'!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_inflow_patch: Invalid extent on '" //                     &
            trim(patch%name) // "'!")
    end do

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die("verify_inflow_patch: '" // &
         trim(patch%name) //  "' extends more than 1 grid point along normal direction!")

    return
  end subroutine verify_inflow_patch


  subroutine inflow_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_InflowWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: localConservedVariables(:), metricsAlongNormalDirection(:),     &
         inviscidPenalty(:)
    real(WP) :: normalMomentum

    direction = abs(patch%normalDirection)

    allocate(localConservedVariables(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(inviscidPenalty(nUnknowns))

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

             localConservedVariables = conservedVariables(gridIndex,:)
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             normalMomentum = dot_product(localConservedVariables(2:nDimensions+1),          &
                  metricsAlongNormalDirection)

             inviscidPenalty(1) = normalMomentum
             inviscidPenalty(2:nDimensions+1) = normalMomentum * velocity(gridIndex,:)
             inviscidPenalty(nDimensions+2) =                                                &
                  normalMomentum * specificVolume(gridIndex, 1) *                            &
                  (localConservedVariables(nDimensions+2) + pressure(gridIndex, 1))
             do l = 1, nSpecies 
                inviscidPenalty(nDimensions+2+l) = normalMomentum * massFraction(gridIndex, l)
             end do
             inviscidPenalty = inviscidPenalty - bc%targetPenalty(patchIndex, :)

             source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *          &
                  jacobian(gridIndex, 1) * inviscidPenalty

          end do
       end do
    end do

    deallocate(inviscidPenalty)
    deallocate(metricsAlongNormalDirection)
    deallocate(localConservedVariables)

    return
  end subroutine inflow_patch_forward


  subroutine inflow_patch_adjoint(patch, bc, source)

    ! External modules
    use state_jacobian
    use simulation_flags
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_InflowWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: localConservedVariables(:), metricsAlongNormalDirection(:),     &
         inviscidPenalty(:), deltaPressure(:), deltaInviscidPenalty(:,:)

    direction = abs(patch%normalDirection)

    allocate(localConservedVariables(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(inviscidPenalty(nUnknowns))
    allocate(deltaPressure(nUnknowns))
    allocate(deltaInviscidPenalty(nUnknowns, nUnknowns))

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

             localConservedVariables = conservedVariables(gridIndex,:)
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             deltaPressure(1) = 0.5_WP * sum(velocity(gridIndex,:) ** 2)
             deltaPressure(2:nDimensions+1) = - velocity(gridIndex,:)
             deltaPressure(nDimensions+2) = 1.0_WP
             do l = 1, nSpecies
                deltaPressure(nDimensions+2+l) = 0.0_WP
             end do
             deltaPressure = deltaPressure * (ratioOfSpecificHeats - 1.0_WP)

             call compute_jacobian_of_inviscid_flux(localConservedVariables,                 &
                  metricsAlongNormalDirection, deltaInviscidPenalty,                         &
                  specificVolume = specificVolume(gridIndex, 1),                             &
                  velocity = velocity(gridIndex, :), pressure = pressure(gridIndex, 1),      &
                  massFraction = massFraction(gridIndex, :))

             do l = 1, nDimensions
                deltaInviscidPenalty(l+1,:) = deltaInviscidPenalty(l+1,:) -                  &
                     metricsAlongNormalDirection(l) * deltaPressure
             end do

             if (useContinuousAdjoint) then
                source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *       &
                     jacobian(gridIndex, 1) * matmul(transpose(deltaInviscidPenalty),        &
                     adjointVariables(gridIndex,:))
             else
                source(gridIndex,:) = source(gridIndex,:) + bc%inviscidPenaltyAmount *       &
                     jacobian(gridIndex, 1) * matmul(transpose(deltaInviscidPenalty),        &
                     adjointVariables(gridIndex,:))
             end if

          end do
       end do
    end do

    deallocate(deltaInviscidPenalty)
    deallocate(deltaPressure)
    deallocate(inviscidPenalty)
    deallocate(metricsAlongNormalDirection)
    deallocate(localConservedVariables)

    return
  end subroutine inflow_patch_adjoint

end module inflow


! =================================== !
! Setup the inflow wall patches !
! =================================== !
subroutine inflow_setup

  ! Internal modules
  use inflow

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nInflowWalls = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_INFLOW) then
        nInflowWalls = nInflowWalls + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nInflowWalls .eq. 0) return

  ! Allocate the sponge type
  allocate(inflowWallData(nInflowWalls))

  ! Connect the boundary patch
  inflowWallPatch => patches(j:j+nInflowWalls-1)

  ! Setup the boundary conditions
  do i = 1, nInflowWalls
     call setup_inflow_patch(inflowWallPatch(i), inflowWallData(i))
  end do

  return
end subroutine inflow_setup


! ===================================== !
! Cleanup the inflow wall patches !
! ===================================== !
subroutine inflow_cleanup

  ! Internal modules
  use inflow

  implicit none

  ! Local variables
  integer :: i

  if (nInflowWalls .gt. 0) then
     do i = 1, nInflowWalls
        call cleanup_inflow_patch(inflowWallPatch(i), inflowWallData(i))
     end do
     deallocate(inflowWallData)
     nullify(inflowWallPatch)
  end if

  nInflowWalls = 0

  return
end subroutine inflow_cleanup


! ======================================================= !
! Add the inflow wall source during the forward run !
! ======================================================= !
subroutine inflow_forward(source)

  ! Internal modules
  use inflow

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nInflowWalls
     call inflow_patch_forward(inflowWallPatch(i),                               &
          inflowWallData(i), source)
  end do

  return
end subroutine inflow_forward


! ======================================================= !
! Add the inflow wall source during the adjoint run !
! ======================================================= !
subroutine inflow_adjoint(source)

  ! Internal modules
  use inflow

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nInflowWalls
     call inflow_patch_adjoint(inflowWallPatch(i),                               &
          inflowWallData(i), source)
  end do

  return
end subroutine inflow_adjoint
