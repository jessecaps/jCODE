module farfield

  ! External modules
  use boundary
  use precision

  implicit none

  integer :: nFarFields

  type, private :: t_FarField
     real(WP) :: inviscidPenaltyAmount, viscousPenaltyAmount
     real(WP), allocatable, dimension(:,:,:) :: viscousFluxes, targetViscousFluxes
  end type t_FarField

  type(t_FarField), allocatable :: farFieldData(:)
  type(t_Patch), pointer :: farFieldPatch(:)

contains

  subroutine setup_farfield_patch(patch, bc, targetViscousFluxes)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use geometry
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_FarField), intent(inout) :: bc
    real(WP), intent(in) :: targetViscousFluxes(:,:,:)

    ! Local variables
    integer :: i, j, k, gridIndex, patchIndex
    real(WP) :: radius, holeRadius, holePosition(3)
    logical :: puncture, inverted
    character(len = str_medium) :: holeShape

    ! Verify the patch type
    call verify_farfield_patch(patch)

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

    ! Viscous penalty amount
    if (useViscosity) then
       call parser_read(trim(patch%name) //  ' viscous penalty amount',                      &
            bc%viscousPenaltyAmount, defaultViscousPenaltyAmount)
       bc%viscousPenaltyAmount = sign(bc%viscousPenaltyAmount,                               &
            real(patch%normalDirection, WP))
       bc%viscousPenaltyAmount = bc%viscousPenaltyAmount /                                   &
            firstDerivative(abs(patch%normalDirection))%normBoundary(1)
    else
       bc%viscousPenaltyAmount = 0.0_WP
    end if

    ! Setup memory for storing viscous fluxes on the patch
    if (patch%nPatchPoints .gt. 0 .and. useViscosity) then
       allocate(bc%viscousFluxes(patch%nPatchPoints, nUnknowns, nDimensions))
       allocate(bc%targetViscousFluxes(patch%nPatchPoints, nUnknowns, nDimensions))
    end if

    ! Send the target viscous fluxes to the boundary patch
    if (useViscosity) call patch_collect(patch, targetViscousFluxes, bc%targetViscousFluxes)

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

          call die("setup_farfield_patch: Unknown hole shape on patch " //                   &
               trim(patch%name) // " '" // trim(holeShape) // "' !")

       end select

    end if

    return
  end subroutine setup_farfield_patch


  subroutine cleanup_farfield_patch(patch, bc)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_FarField), intent(inout) :: bc

    if (allocated(bc%viscousFluxes)) deallocate(bc%viscousFluxes)
    if (allocated(bc%targetViscousFluxes)) deallocate(bc%targetViscousFluxes)

    return
  end subroutine cleanup_farfield_patch


  subroutine verify_farfield_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, extent(6)

    if (patch%patchType .ne. SAT_FAR_FIELD)                                                  &
         call die('verify_farfield_patch: Patch type mismatch!')

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .eq. 0)            &
         call die("verify_farfield_patch: Normal direction is invalid for '" //              &
         trim(patch%name) // "'!")

    if (.not. useTargetState)                                                                &
         call die('verify_farfield_patch: No target state available for far-field patch!')

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_farfield_patch: Invalid extent on '" // trim(patch%name) // "'!")
    end do

    i = abs(patch%normalDirection)
    if (extent((i-1)*2+1) .ne. extent((i-1)*2+2)) call die("verify_farfield_patch: '" //     &
         trim(patch%name) //  "' extends more than 1 grid point along normal direction!")

    return
  end subroutine verify_farfield_patch


  subroutine farfield_patch_forward(patch, bc, source)

    ! External modules
    use state_jacobian
    use solver_options
    use geometry
    use grid
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_FarField), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, incomingDirection, gridIndex, patchIndex
    real(WP), allocatable :: localConservedVariables(:), localTargetState(:),                &
         metricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:), vfCorrection

    direction = abs(patch%normalDirection)
    incomingDirection = patch%normalDirection

    allocate(localConservedVariables(nUnknowns))
    allocate(localTargetState(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
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

             localTargetState = targetState(gridIndex,:)
             localConservedVariables = conservedVariables(gridIndex,:)
             if (twoWayCoupling) then
                ! Remove volume fraction from conserved variables
                do l = 1, nUnknowns 
                   localConservedVariables(l) = conservedVariables(gridIndex,l) /            &
                        volumeFraction(gridIndex,1)
                end do
                vfCorrection = volumeFraction(gridIndex,1)
             else
                vfCorrection = 1.0_WP
             end if
             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             call compute_incoming_jacobian_of_inviscid_flux(localTargetState,               &
                  metricsAlongNormalDirection, incomingDirection,                            &
                  incomingJacobianOfInviscidFlux)

             source(gridIndex,:) = source(gridIndex,:) - bc%inviscidPenaltyAmount *          &
                  jacobian(gridIndex, 1) * matmul(incomingJacobianOfInviscidFlux,            &
                  localConservedVariables - localTargetState) * vfCorrection

             if (useViscosity) source(gridIndex,:) = source(gridIndex,:) +                   &
                  bc%viscousPenaltyAmount * jacobian(gridIndex, 1) *                         &
                  vfCorrection * matmul(bc%viscousFluxes(patchIndex,:,:) -                   &
                  bc%targetViscousFluxes(patchIndex,:,:), metricsAlongNormalDirection)

          end do
       end do
    end do

    deallocate(incomingJacobianOfInviscidFlux)
    deallocate(metricsAlongNormalDirection)
    deallocate(localTargetState)

    return
  end subroutine farfield_patch_forward


  subroutine farfield_patch_adjoint(patch, bc, source)

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
    type(t_FarField), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, incomingDirection, gridIndex, patchIndex
    real(WP), allocatable :: localTargetState(:), metricsAlongNormalDirection(:),            &
         metricsAlongDirection2(:), incomingJacobianOfInviscidFlux(:,:),                     &
         localViscousFluxJacobian(:,:), localViscousFluxJacobian2(:,:), localMassFraction(:),&
         localEnthalpyFlux(:), localSpeciesFlux(:,:), temp1(:,:,:), temp2(:,:)

    direction = abs(patch%normalDirection)
    if (useContinuousAdjoint) then
       incomingDirection = -patch%normalDirection
    else
       incomingDirection = +patch%normalDirection
    end if

    allocate(localTargetState(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
    if (useViscosity) then
       allocate(metricsAlongDirection2(nDimensions))
       allocate(localViscousFluxJacobian(nUnknowns, nUnknowns))
       allocate(localViscousFluxJacobian2(nUnknowns-1, nUnknowns-1))
       if (nSpecies .gt. 0) then
          allocate(localMassFraction(nSpecies))
          allocate(localSpeciesFlux(nSpecies,nDimensions))
          if (equationOfState .eq. IDEAL_GAS_MIXTURE) allocate(localEnthalpyFlux(nDimensions))
       end if
       allocate(temp1(nGridPoints, nUnknowns - 1, nDimensions)); temp1 = 0.0_WP
    end if

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

             if (useViscosity) then

                if (nSpecies .gt. 0) then
                   localMassFraction = massFraction(gridIndex,:)
                   localSpeciesFlux = speciesFlux(gridIndex,:,:)
                   if (equationOfState .eq. IDEAL_GAS_MIXTURE)                               &
                        localEnthalpyFlux = enthalpyFlux(gridIndex,:)
                end if
                call compute_first_partial_viscous_jacobian(conservedVariables(gridIndex,:), &
                     metricsAlongNormalDirection, stressTensor(gridIndex,:),                 &
                     heatFlux(gridIndex,:), localEnthalpyFlux, localSpeciesFlux,             &
                     localViscousFluxJacobian, specificVolume(gridIndex,1),                  &
                     velocity(gridIndex,:), temperature(gridIndex,1), localMassFraction)
                source(gridIndex,:) = source(gridIndex,:) - bc%viscousPenaltyAmount *        &
                     jacobian(gridIndex, 1) * matmul(transpose(localViscousFluxJacobian),    &
                     adjointVariables(gridIndex,:))

                if (useContinuousAdjoint) cycle

                do l = 1, nDimensions
                   metricsAlongDirection2 =                                                  &
                        metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)
                   call compute_second_partial_viscous_jacobian(velocity(gridIndex,:),       &
                        dynamicViscosity(gridIndex,1),                                       &
                        secondCoefficientOfViscosity(gridIndex,1),                           &
                        thermalDiffusivity(gridIndex,1), massDiffusivity(gridIndex,:),       &
                        temperature(gridIndex,1), jacobian(gridIndex,1),                     &
                        metricsAlongNormalDirection, metricsAlongDirection2,                 &
                        localViscousFluxJacobian2)
                   temp1(gridIndex,:,l) = - matmul(transpose(localViscousFluxJacobian2),     &
                        adjointVariables(gridIndex,2:nUnknowns))
                end do

             end if
          end do
       end do
    end do

    deallocate(incomingJacobianOfInviscidFlux)
    deallocate(metricsAlongNormalDirection)
    deallocate(localTargetState)
    if (allocated(metricsAlongDirection2)) deallocate(metricsAlongDirection2)
    if (allocated(localViscousFluxJacobian)) deallocate(localViscousFluxJacobian)
    if (allocated(localViscousFluxJacobian2)) deallocate(localViscousFluxJacobian2)
    if (allocated(localMassFraction)) deallocate(localMassFraction)
    if (allocated(localSpeciesFlux)) deallocate(localSpeciesFlux)
    if (allocated(localEnthalpyFlux)) deallocate(localEnthalpyFlux)

    if (useContinuousAdjoint .or. .not. useViscosity) return

    allocate(temp2(nGridPoints, nUnknowns - 1)); temp2 = 0.0_WP

    do i = 1, nDimensions
       call adjoint_first_derivative_apply(i, temp1(:,:,i))
    end do
    temp2 = sum(temp1, dim = 3) !... divergence of the adjoint flux

    deallocate(temp1)  !... no longer needed

    select case (equationOfState)

    case (IDEAL_GAS)

       temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *                 &
            temp2(:,nDimensions+1)
       do i = 1, nDimensions
          temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                    &
               temp2(:,nDimensions+1)
       end do
       do k = 1, nSpecies
          temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k)
       end do

    case (IDEAL_GAS_MIXTURE)

       temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *                 &
            mixtureMolecularWeight(:,1) * temp2(:,nDimensions+1)
       do i = 1, nDimensions
          temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *                    &
               temp2(:,nDimensions+1)
       end do
       do k = 1, nSpecies
          temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k) +        &
               temperature(:,1) * (molecularWeightInverse(nSpecies+1) -                      &
               molecularWeightInverse(k)) * temp2(:,nDimensions+1) /                         &
               ratioOfSpecificHeats
       end do

    end select

    do i = 1, size(temp2, 2)
       temp2(:,i) = bc%viscousPenaltyAmount * jacobian(:,1) * temp2(:,i)
    end do

    source(:,2:nUnknowns) = source(:,2:nUnknowns) + temp2
    source(:,1) = source(:,1) - specificVolume(:,1) * conservedVariables(:,nDimensions+2)    &
         * temp2(:,nDimensions+1) - sum(velocity * temp2(:,1:nDimensions), dim = 2)
    if (nSpecies .gt. 0) source(:,1) = source(:,1) -                                         &
         sum(massFraction * temp2(:,nDimensions+2:nUnknowns-1), dim = 2)
    
    deallocate(temp2)

    return
  end subroutine farfield_patch_adjoint

end module farfield


! ==================================== !
! Setup the far-field boundary patches !
! ==================================== !
subroutine farfield_setup

  ! Internal modules
  use farfield

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
  real(WP), allocatable :: targetViscousFluxes(:,:,:)

  ! Find the number of boundaries of this type
  nFarFields = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_FAR_FIELD) then
        nFarFields = nFarFields + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nFarFields .eq. 0) return

  ! Allocate the far-field type
  allocate(farFieldData(nFarFields))

  ! Connect the boundary patch
  farFieldPatch => patches(j:j+nFarFields-1)

  ! Compute the target viscous fluxes
  if (useViscosity) then
     allocate(targetViscousFluxes(nGridPoints, nUnknowns, nDimensions))
     call update_state(targetState)
     if (equationOfState .eq. IDEAL_GAS) then
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
             targetViscousFluxes, speciesFlux = speciesFlux)
     else
        call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
             targetViscousFluxes, speciesFlux, enthalpyFlux)
     end if
  end if

  ! Setup the boundary conditions
  do i = 1, nFarFields
     call setup_farfield_patch(farFieldPatch(i), farFieldData(i), targetViscousFluxes)
  end do

  ! Cleanup
  if (allocated(targetViscousFluxes)) deallocate(targetViscousFluxes)

  return
end subroutine farfield_setup


! ====================================== !
! Cleanup the far-field boundary patches !
! ====================================== !
subroutine farfield_cleanup

  ! Internal modules
  use farfield

  ! Local variables
  integer :: i

  if (nFarFields .gt. 0) then
     do i = 1, nFarFields
        call cleanup_farfield_patch(farFieldPatch(i), farFieldData(i))
     end do
     deallocate(farFieldData)
     nullify(farFieldPatch)
  end if

  nFarFields = 0

  return
end subroutine farfield_cleanup


! =============================================== !
! Add the far-field source during the forward run !
! =============================================== !
subroutine farfield_forward(source)

  ! Internal modules
  use farfield

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nFarFields
     call farfield_patch_forward(farFieldPatch(i), farFieldData(i), source)
  end do

  return
end subroutine farfield_forward


! =============================================== !
! Add the far-field source during the adjoint run !
! =============================================== !
subroutine farfield_adjoint(source)

  ! Internal modules
  use farfield

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nFarFields
     call farfield_patch_adjoint(farFieldPatch(i), farFieldData(i), source)
  end do

  return
end subroutine farfield_adjoint


! ================================================= !
! Store the viscous fluxes on the far-field patches !
! ================================================= !
subroutine farfield_store_viscous_fluxes(fluxes)

  ! Internal modules
  use farfield

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns, nDimensions), intent(in) :: fluxes

  ! Local variables
  integer :: i

  do i = 1, nFarFields
     if (allocated(farFieldData(i)%viscousFluxes))                                           &
          call patch_collect(farFieldPatch(i), fluxes, farFieldData(i)%viscousFluxes)
  end do

  return
end subroutine farfield_store_viscous_fluxes
