module impenetrable

  ! External modules
  use boundary
  use precision

  implicit none

  type :: t_ImpenetrableWall
     real(WP) :: inviscidPenaltyAmount
  end type t_ImpenetrableWall

contains

  subroutine setup_impenetrable_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use geometry
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_ImpenetrableWall), intent(inout) :: bc

    ! Local variables
    integer :: i, j, k, gridIndex, patchIndex
    real(WP) :: radius, holeRadius, holePosition(3)
    logical :: puncture, inverted
    character(len = str_medium) :: holeShape

    ! Verify the patch type
    call verify_impenetrable_patch(patch)

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

          call die("setup_impenetrable_patch: Unknown hole shape on patch " //               &
               trim(patch%name) // " '" // trim(holeShape) // "' !")

       end select

    end if

    return
  end subroutine setup_impenetrable_patch


  subroutine cleanup_impenetrable_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_ImpenetrableWall), intent(inout) :: bc

    ! Nothing to do

    return
  end subroutine cleanup_impenetrable_patch


  subroutine verify_impenetrable_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_impenetrable_patch: Normal direction is invalid for '" //          &
         trim(patch%name) // "'!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_impenetrable_patch: Invalid extent on '" //                     &
            trim(patch%name) // "'!")
    end do

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die("verify_impenetrable_patch: '" // &
         trim(patch%name) //  "' extends more than 1 grid point along normal direction!")

    return
  end subroutine verify_impenetrable_patch


  subroutine impenetrable_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_ImpenetrableWall), intent(inout) :: bc
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

             if (twoWayCoupling) then
                localConservedVariables = conservedVariables(gridIndex,:)                    &
                     / volumeFraction(gridIndex, 1)
             else
                localConservedVariables = conservedVariables(gridIndex,:)
             end if
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
             if (twoWayCoupling) inviscidPenalty = inviscidPenalty *                         &
                  volumeFraction(gridIndex, 1)

             source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *          &
                  jacobian(gridIndex, 1) * inviscidPenalty

          end do
       end do
    end do

    deallocate(inviscidPenalty)
    deallocate(metricsAlongNormalDirection)
    deallocate(localConservedVariables)

    return
  end subroutine impenetrable_patch_forward


  subroutine impenetrable_patch_adjoint(patch, bc, source)

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
    type(t_ImpenetrableWall), intent(inout) :: bc
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
  end subroutine impenetrable_patch_adjoint

end module impenetrable
