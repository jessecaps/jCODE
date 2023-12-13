module isothermal

  ! External modules
  use boundary
  use impenetrable, only : t_ImpenetrableWall
  use precision

  implicit none

  integer :: nIsothermalWalls

  type, private :: t_IsothermalWall
     type(t_ImpenetrableWall) :: impenetrableWallData
     real(WP) :: viscousPenaltyAmount, neumannCoefficient, targetVelocity(3)
     real(WP), allocatable :: temperature(:), massFraction(:,:), dynamicViscosity(:),        &
          secondCoefficientOfViscosity(:), thermalDiffusivity(:), massDiffusivity(:,:)
     logical :: enforceMassFraction
  end type t_IsothermalWall

  type(t_IsothermalWall), allocatable :: isothermalWallData(:)
  type(t_Patch), pointer :: isothermalWallPatch(:)

contains

  subroutine setup_isothermal_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use first_derivative
    use state
    use equation_of_state
    use impenetrable, only : setup_impenetrable_patch

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_IsothermalWall), intent(inout) :: bc

    ! Local variables
    integer :: i
    real(WP) :: wallTemperature, wallMassFraction
    real(WP), allocatable :: targetTemperature(:), targetMassFraction(:,:)
    character(len = str_medium) :: species
    logical :: useTarget

    ! Verify the patch type
    call verify_isothermal_patch(patch)

    ! Isothermal wall requires impenetrable wall data
    call setup_impenetrable_patch(patch, bc%impenetrableWallData)

    ! Neumann or Derichlet condition for species mass fraction?
    call parser_read(trim(patch%name) // ' enforce wall mass fraction',                      &
         bc%enforceMassFraction, .false.)
    if (nSpecies .eq. 0) bc%enforceMassFraction = .false.

    if (patch%nPatchPoints .gt. 0 .and. useViscosity) then

       allocate(bc%temperature(patch%nPatchPoints))
       allocate(bc%dynamicViscosity(patch%nPatchPoints))
       allocate(bc%secondCoefficientOfViscosity(patch%nPatchPoints))
       allocate(bc%thermalDiffusivity(patch%nPatchPoints))
       if (bc%enforceMassFraction) then
          allocate(bc%massDiffusivity(patch%nPatchPoints, nSpecies))
          allocate(bc%massFraction(patch%nPatchPoints, nSpecies))
       end if

       ! Assign target velocity (default is no-slip)
       call parser_is_defined(trim(patch%name) // ' velocity', useTarget)
       if (useTarget .and. predictionOnly) then
          call parser_read(trim(patch%name) // ' velocity', bc%targetVelocity)
       else
          bc%targetVelocity = 0.0_WP
       end if

       ! Assign the wall temperature & mass fraction
       if (useTargetState) then ! ... get values from target state
          allocate(targetTemperature(nGridPoints))
          if (nSpecies .gt. 0) then
             allocate(targetMassFraction(nGridPoints, nSpecies))
             call compute_dependent_variables(targetState, temperature = targetTemperature,  &
                  massFraction = targetMassFraction)
          else
             call compute_dependent_variables(targetState, temperature = targetTemperature)
          end if
          call patch_collect(patch, targetTemperature, bc%temperature)
          if (bc%enforceMassFraction)                                                        &
               call patch_collect(patch, targetmassFraction, bc%massFraction)
          deallocate(targetTemperature)
          if (allocated(targetMassFraction)) deallocate(targetMassFraction)

       else ! ... read in values from the input file
          call parser_read(trim(patch%name) // ' temperature', wallTemperature,              &
               1.0_WP / (ratioOfSpecificHeats - 1.0_WP))
          bc%temperature(:) = wallTemperature
          if (bc%enforceMassFraction) then
             do i = 1, nSpecies
                write(species, "(A,I1.1)") "mass fraction ", i
                call parser_read(trim(patch%name) // trim(species), wallMassFraction)
                bc%massFraction(:,i) = wallMassFraction
             end do
          end if
       end if

       ! Assign the viscous fluxes
       if (bc%enforceMassFraction) then
          call compute_transport_variables(bc%temperature, bc%dynamicViscosity,              &
               bc%secondCoefficientOfViscosity, bc%thermalDiffusivity, bc%massDiffusivity)
       else
          call compute_transport_variables(bc%temperature, bc%dynamicViscosity,              &
               bc%secondCoefficientOfViscosity, bc%thermalDiffusivity)
       end if

    end if

    ! Viscous penalty amount
    if (useViscosity) then
       call parser_read(trim(patch%name) // ' viscous penalty amount',                       &
            bc%viscousPenaltyAmount, defaultViscousPenaltyAmount)
       bc%viscousPenaltyAmount = bc%viscousPenaltyAmount * 0.25_WP * reynoldsNumberInverse / &
            firstDerivative(abs(patch%normalDirection))%normBoundary(1) *                    &
            max(ratioOfSpecificHeats * prandtlNumberInverse, 5.0_WP / 3.0_WP)

       ! Species penalty amount to enforce Neumann condition on mass fraction
       if (.not. bc%enforceMassFraction) then
          bc%neumannCoefficient = 1.0_WP /                                                   &
               firstDerivative(abs(patch%normalDirection))%normBoundary(1)
          bc%neumannCoefficient = sign(bc%neumannCoefficient,                                &
               real(patch%normalDirection, WP))
       else
          bc%neumannCoefficient = 0.0_WP
       end if
    else
       bc%viscousPenaltyAmount = 0.0_WP
       bc%neumannCoefficient = 0.0_WP
    end if

    return
  end subroutine setup_isothermal_patch


  subroutine cleanup_isothermal_patch(patch, bc)

    ! External modules
    use impenetrable, only : cleanup_impenetrable_patch

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_IsothermalWall), intent(inout) :: bc

    call cleanup_impenetrable_patch(patch, bc%impenetrableWallData)

    if (allocated(bc%temperature)) deallocate(bc%temperature)
    if (allocated(bc%massFraction)) deallocate(bc%massFraction)
    if (allocated(bc%dynamicViscosity)) deallocate(bc%dynamicViscosity)
    if (allocated(bc%secondCoefficientOfViscosity))                                          &
         deallocate(bc%secondCoefficientOfViscosity)
    if (allocated(bc%thermalDiffusivity)) deallocate(bc%thermalDiffusivity)
    if (allocated(bc%massDiffusivity)) deallocate(bc%massDiffusivity)

    return
  end subroutine cleanup_isothermal_patch


  subroutine verify_isothermal_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    if (patch%patchType .ne. SAT_ISOTHERMAL_WALL)                                            &
         call die('verify_isothermal_patch: patch type mismatch for SAT_ISOTHERMAL_WALL')

    return
  end subroutine verify_isothermal_patch


  subroutine isothermal_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use state
    use impenetrable, only : impenetrable_patch_forward

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_IsothermalWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: metricsAlongNormalDirection(:), penaltyAtBoundary(:), temp(:,:)

    ! Enforce the no-penetration condition
    call impenetrable_patch_forward(patch, bc%impenetrableWallData, source)

    if (.not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(penaltyAtBoundary(nUnknowns))
    if (nSpecies .gt. 0 .and. .not. bc%enforceMassFraction)                                  &
         allocate(temp(nUnknowns, nDimensions))

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

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             penaltyAtBoundary(1) = 0.0_WP
             penaltyAtBoundary(2:nDimensions+1) =                                            &
                  conservedVariables(gridIndex,2:nDimensions+1) -                            &
                  conservedVariables(gridIndex, 1) * bc%targetVelocity(1:nDimensions)
             penaltyAtBoundary(nDimensions+2) =                                              &
                  conservedVariables(gridIndex,nDimensions+2) - 0.5_WP *                     &
                  conservedVariables(gridIndex, 1) * sum(bc%targetVelocity(1:nDimensions)**2)
             if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
                penaltyAtBoundary(nDimensions+2) = penaltyAtBoundary(nDimensions+2) -        &
                     conservedVariables(gridIndex,1) * bc%temperature(patchIndex) /          &
                     ratioOfSpecificHeats / mixtureMolecularWeight(gridIndex, 1)
             else
                penaltyAtBoundary(nDimensions+2) = penaltyAtBoundary(nDimensions+2) -        &
                     conservedVariables(gridIndex,1) * bc%temperature(patchIndex) /          &
                     ratioOfSpecificHeats
             end if
             if (bc%enforceMassFraction) then
                do l = 1, nSpecies
                   penaltyAtBoundary(nDimensions+2+l) =                                      &
                        conservedVariables(gridIndex,nDimensions+2+l) -                      &
                        conservedVariables(gridIndex, 1) * bc%massFraction(patchIndex, l)
                end do
             else
                do l = 1, nSpecies
                   penaltyAtBoundary(nDimensions+2+l) = 0.0_WP !... will be applied later
                end do
             end if
             penaltyAtBoundary =  penaltyAtBoundary * jacobian(gridIndex, 1)**2 *            &
                  dot_product(metricsAlongNormalDirection, metricsAlongNormalDirection)

             source(gridIndex,:) = source(gridIndex,:) - bc%viscousPenaltyAmount *           &
                  penaltyAtBoundary

             if (nSpecies .gt. 0 .and. .not. bc%enforceMassFraction) then
                temp = 0.0_WP
                select case (nDimensions)

                case (1)
                   do l = 1, nSpecies
                      temp(l+3,1) = - speciesFlux(gridIndex,l,1)
                   end do

                case (2)
                   do l = 1, nSpecies
                      temp(l+4,1) = - speciesFlux(gridIndex,l,1)
                      temp(l+4,2) = - speciesFlux(gridIndex,l,2)
                   end do

                case (3)
                   do l = 1, nSpecies
                      temp(l+5,1) = - speciesFlux(gridIndex,l,1)
                      temp(l+5,2) = - speciesFlux(gridIndex,l,2)
                      temp(l+5,3) = - speciesFlux(gridIndex,l,3)
                   end do

                end select

                ! Zero flux condition for mass fraction
                source(gridIndex,:) = source(gridIndex,:) + bc%neumannCoefficient *          &
                     jacobian(gridIndex, 1) * matmul(temp, metricsAlongNormalDirection)
             end if

          end do
       end do
    end do

    deallocate(metricsAlongNormalDirection)
    deallocate(penaltyAtBoundary)
    if (nSpecies .gt. 0 .and. .not. bc%enforceMassFraction) deallocate(temp)

    return
  end subroutine isothermal_patch_forward


  subroutine isothermal_patch_adjoint(patch, bc, source)

    ! External modules
    use simulation_flags
    use solver_options
    use geometry
    use grid
    use state
    use impenetrable, only : impenetrable_patch_adjoint
    use first_derivative

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_IsothermalWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: adjointPenalties(:), metricsAlongNormalDirection(:),            &
         metricsAlongDirection2(:), localViscousFluxJacobian(:,:),                           &
         localViscousFluxJacobian2(:,:), localMassFraction(:), localSpeciesFlux(:,:),        &
         temp1(:,:,:), temp2(:,:)

    ! Apply no-penetration conditions
    call impenetrable_patch_adjoint(patch, bc%impenetrableWallData, source)

    if (useContinuousAdjoint .or. .not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(adjointPenalties(nUnknowns))
    allocate(metricsAlongNormalDirection(nDimensions))
    if (nSpecies .gt. 0) then
       allocate(metricsAlongDirection2(nDimensions))
       allocate(localViscousFluxJacobian(nUnknowns, nUnknowns))
       allocate(localViscousFluxJacobian2(nUnknowns-1, nUnknowns-1))
       allocate(localMassFraction(nSpecies))
       allocate(localSpeciesFlux(nSpecies,nDimensions))
    end if
    allocate(temp1(nGridPoints, nUnknowns - 1, nDimensions)); temp1 = 0.0_WP

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

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             adjointPenalties = 0.0_WP
             adjointPenalties(1) = - adjointVariables(gridIndex,nDimensions+2) *             &
                  bc%temperature(patchIndex) / ratioOfSpecificHeats
             if (equationOfState .eq. IDEAL_GAS_MIXTURE) adjointPenalties(1) =               &
                  adjointPenalties(1) * molecularWeightInverse(nSpecies+1)

             if (bc%enforceMassFraction) then
                do l = 1, nSpecies
                   adjointPenalties(1) = adjointPenalties(1) -                               &
                        adjointVariables(gridIndex, nDimensions+2+l) *                       &
                        bc%massFraction(patchIndex, l)
                end do
             else
                !... nothing to do here
             end if

             adjointPenalties(2:nDimensions+2) = adjointVariables(gridIndex,2:nDimensions+2)
             if (bc%enforceMassFraction) then
                do l = 1, nSpecies
                   adjointPenalties(nDimensions+2+l) =                                       &
                        adjointVariables(gridIndex,nDimensions+2+l)
                end do
             else
                !... nothing to do here
             end if

             if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
                do l = 1, nSpecies
                   adjointPenalties(nDimensions+2+l) = adjointPenalties(nDimensions+2+l)     &
                        - adjointVariables(gridIndex,nDimensions+2) *                        &
                        bc%temperature(patchIndex) / ratioOfSpecificHeats *                  &             
                        (molecularWeightInverse(l) - molecularWeightInverse(nSpecies+1))
                end do
             end if

             adjointPenalties = adjointPenalties * jacobian(gridIndex, 1)**2 *               &
                  dot_product(metricsAlongNormalDirection, metricsAlongNormalDirection)

             source(gridIndex,:) = source(gridIndex,:) + bc%viscousPenaltyAmount *           &
                  adjointPenalties(:)

             if (nSpecies .gt. 0 .and. .not. bc%enforceMassFraction) then
                localMassFraction = massFraction(gridIndex,:)
                localSpeciesFlux = speciesFlux(gridIndex,:,:)
             call compute_first_partial_fluxes_jacobian(conservedVariables(gridIndex,:),     &
                  metricsAlongNormalDirection, specificVolume(gridIndex,1),                  &
                  velocity(gridIndex,:), temperature(gridIndex,1), localViscousFluxJacobian, &
                  localMassFraction, localSpeciesFlux, gridIndex)
             source(gridIndex,:) = source(gridIndex,:) - bc%neumannCoefficient *             &
                  jacobian(gridIndex, 1) * matmul(transpose(localViscousFluxJacobian),       &
                  adjointVariables(gridIndex,:))

             do l = 1, nDimensions
                metricsAlongDirection2 = metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)
                call compute_second_partial_fluxes_jacobian(dynamicViscosity(gridIndex,1),   &
                     massDiffusivity(gridIndex,:), jacobian(gridIndex,1),                    &
                     metricsAlongNormalDirection, metricsAlongDirection2,                    &
                     localViscousFluxJacobian2, gridIndex)
                temp1(gridIndex,:,l) = - matmul(transpose(localViscousFluxJacobian2),        &
                     adjointVariables(gridIndex,2:nUnknowns))
             end do
             end if
          end do
       end do
    end do

    deallocate(adjointPenalties)
    deallocate(metricsAlongNormalDirection)

    if (nSpecies .le. 0 .or. bc%enforceMassFraction) return

    deallocate(metricsAlongDirection2)
    deallocate(localViscousFluxJacobian)
    deallocate(localViscousFluxJacobian2)
    deallocate(localMassFraction)
    deallocate(localSpeciesFlux)

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
       temp2(:,i) = bc%neumannCoefficient * jacobian(:,1) * temp2(:,i)
    end do

    source(:,2:nUnknowns) = source(:,2:nUnknowns) + temp2
    source(:,1) = source(:,1) - specificVolume(:,1) * conservedVariables(:,nDimensions+2)    &
         * temp2(:,nDimensions+1) - sum(velocity * temp2(:,1:nDimensions), dim = 2)
    if (nSpecies .gt. 0) source(:,1) = source(:,1) -                                         &
         sum(massFraction * temp2(:,nDimensions+2:nUnknowns-1), dim = 2)

    deallocate(temp2)

    return
  end subroutine isothermal_patch_adjoint


  subroutine compute_first_partial_fluxes_jacobian(conservedVariables, metrics,              &
       specificVolume, velocity, temperature, firstPartialViscousJacobian, massFraction,     &
       speciesFlux, gridIndex)

    ! External modules
    use simulation_flags
    use solver_options

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:), metrics(:), specificVolume, velocity(:),  &
         temperature
    integer, intent(in) :: gridIndex
    real(WP), intent(in), optional :: massFraction(:), speciesFlux(:,:)
    real(WP), intent(out) :: firstPartialViscousJacobian(:,:)

    ! Local variables
    integer :: k
    real(WP) :: phiSquared, contravariantSpeciesFlux(nSpecies),                              &
         deltaViscosity(nDimensions+nSpecies+2)

    ! Zero-out first partial viscous Jacobian
    firstPartialViscousJacobian = 0.0_WP

    select case (nDimensions)

    case (1)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * velocity(1) ** 2
       do k = 1, nSpecies
          contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) !... not normalized
       end do
       deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume /        &
            temperature * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature /      &
            ratioOfSpecificHeats)
       deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /      &
            temperature * velocity(1)
       deltaViscosity(3) = powerLawExponent * ratioOfSpecificHeats * specificVolume /        &
            temperature
       do k = 1, nSpecies
          deltaViscosity(3+k) = 0.0_WP
       end do

       firstPartialViscousJacobian(1,1) = 0.0_WP
       firstPartialViscousJacobian(2,1) = 0.0_WP
       firstPartialViscousJacobian(3,1) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

    case (2)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                               &
            (velocity(1) ** 2 + velocity(2) ** 2)
       do k = 1, nSpecies
          contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                      &
               metrics(2) * speciesFlux(k,2) !... not normalized
       end do
       deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume /        &
            temperature * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature /      &
            ratioOfSpecificHeats)
       deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /      &
            temperature * velocity(1)
       deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /      &
            temperature * velocity(2)
       deltaViscosity(4) = powerLawExponent * ratioOfSpecificHeats * specificVolume /        &
            temperature
       do k = 1, nSpecies
          deltaViscosity(4+k) = 0.0_WP
       end do

       firstPartialViscousJacobian(1,1) = 0.0_WP
       firstPartialViscousJacobian(2,1) = 0.0_WP
       firstPartialViscousJacobian(3,1) = 0.0_WP
       firstPartialViscousJacobian(4,1) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = 0.0_WP
       firstPartialViscousJacobian(4,2) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = 0.0_WP
       firstPartialViscousJacobian(4,3) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,4) = 0.0_WP
       firstPartialViscousJacobian(2,4) = 0.0_WP
       firstPartialViscousJacobian(3,4) = 0.0_WP
       firstPartialViscousJacobian(4,4) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,4) = - deltaViscosity(4) *                         &
               contravariantSpeciesFlux(k)
       end do

    case (3)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                               &
            (velocity(1) ** 2 + velocity(2) ** 2 + velocity(3) ** 2)
       do k = 1, nSpecies
          contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                      &
               metrics(2) * speciesFlux(k,2) +                                               &
               metrics(3) * speciesFlux(k,3) !... not normalized
       end do
       deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume /       &
            temperature * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature /     &
            ratioOfSpecificHeats)
       deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /     &
            temperature * velocity(1)
       deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /     &
            temperature * velocity(2)
       deltaViscosity(4) = - powerLawExponent * ratioOfSpecificHeats * specificVolume /     &
            temperature * velocity(3)
       deltaViscosity(5) = powerLawExponent * ratioOfSpecificHeats * specificVolume /       &
            temperature
       do k = 1, nSpecies
          deltaViscosity(5+k) = 0.0_WP
       end do

       firstPartialViscousJacobian(1,1) = 0.0_WP
       firstPartialViscousJacobian(2,1) = 0.0_WP
       firstPartialViscousJacobian(3,1) = 0.0_WP
       firstPartialViscousJacobian(4,1) = 0.0_WP
       firstPartialViscousJacobian(5,1) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = 0.0_WP
       firstPartialViscousJacobian(4,2) = 0.0_WP
       firstPartialViscousJacobian(5,2) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = 0.0_WP
       firstPartialViscousJacobian(4,3) = 0.0_WP
       firstPartialViscousJacobian(5,3) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,4) = 0.0_WP
       firstPartialViscousJacobian(2,4) = 0.0_WP
       firstPartialViscousJacobian(3,4) = 0.0_WP
       firstPartialViscousJacobian(4,4) = 0.0_WP
       firstPartialViscousJacobian(5,4) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,4) = - deltaViscosity(4) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,5) = 0.0_WP
       firstPartialViscousJacobian(2,5) = 0.0_WP
       firstPartialViscousJacobian(3,5) = 0.0_WP
       firstPartialViscousJacobian(4,5) = 0.0_WP
       firstPartialViscousJacobian(5,5) = 0.0_WP
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,5) = - deltaViscosity(5) *                         &
               contravariantSpeciesFlux(k)
       end do

    end select

    ! Modify the jacobian due to flame thickening if required
    ! if (useFlameThickening)                                                                  &
    !      call flame_thickening_modify_first_partial_viscous_jacobian(gridIndex, metrics,     &
    !      firstPartialViscousJacobian)

    return
  end subroutine compute_first_partial_fluxes_jacobian


  subroutine compute_second_partial_fluxes_jacobian(dynamicViscosity, massDiffusivity,       &
       jacobian, metricsAlongFirstDir, metricsAlongSecondDir, secondPartialViscousJacobian,  &
       gridIndex)

    ! External modules
    use simulation_flags
    use solver_options

    implicit none

    ! Arguments
    real(WP), intent(in) :: dynamicViscosity, jacobian, metricsAlongFirstDir(:),             &
         metricsAlongSecondDir(:)
    integer, intent(in) :: gridIndex
    real(WP), intent(in), optional :: massDiffusivity(:)
    real(WP), intent(out) :: secondPartialViscousJacobian(:,:)

    ! Local variables
    integer :: k
    real(WP) :: temp

    ! Zero-out second partial viscous Jacobian
    secondPartialViscousJacobian = 0.0_WP

    select case (nDimensions)

    case (1)

       ! Temporary variables
       temp = metricsAlongFirstDir(1) * metricsAlongFirstDir(1)

       do k = 1, nSpecies
          secondPartialViscousJacobian(2+k,2+k) = massDiffusivity(k) * temp
       end do

    case (2)

       ! Temporary variables
       temp = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                           &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2)

       do k = 1, nSpecies
          secondPartialViscousJacobian(3+k,3+k) = massDiffusivity(k) * temp
       end do

    case (3)

       ! Temporary variables
       temp = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                           &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2) +                             &
            metricsAlongFirstDir(3) * metricsAlongSecondDir(3)

       do k = 1, nSpecies
          secondPartialViscousJacobian(4+k,4+k) = massDiffusivity(k) * temp
       end do

    end select

    ! Multiply by the Jacobian
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

    ! Modify the jacobian due to flame thickening if required
    ! if (useFlameThickening)                                                                  &
    !      call flame_thickening_modify_second_partial_viscous_jacobian(gridIndex,             &
    !      secondPartialViscousJacobian)

    return
  end subroutine compute_second_partial_fluxes_jacobian
  
end module isothermal


! ================================= !
! Setup the isothermal wall patches !
! ================================= !
subroutine isothermal_setup

  ! Internal modules
  use isothermal

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nIsothermalWalls = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_ISOTHERMAL_WALL) then
        nIsothermalWalls = nIsothermalWalls + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nIsothermalWalls .eq. 0) return

  ! Allocate the isothermal type
  allocate(isothermalWallData(nIsothermalWalls))

  ! Connect the boundary patch
  isothermalWallPatch => patches(j:j+nIsothermalWalls-1)

  ! Setup the boundary conditions
  do i = 1, nIsothermalWalls
     call setup_isothermal_patch(isothermalWallPatch(i), isothermalWallData(i))
  end do

  return
end subroutine isothermal_setup


! =================================== !
! Cleanup the isothermal wall patches !
! =================================== !
subroutine isothermal_cleanup

  ! Internal modules
  use isothermal
  ! Local variables
  integer :: i

  if (nIsothermalWalls .gt. 0) then
     do i = 1, nIsothermalWalls
        call cleanup_isothermal_patch(isothermalWallPatch(i), isothermalWallData(i))
     end do
     deallocate(isothermalWallData)
     nullify(isothermalWallPatch)
  end if

  nIsothermalWalls = 0

  return
end subroutine isothermal_cleanup


! ===================================================== !
! Add the isothermal wall source during the forward run !
! ===================================================== !
subroutine isothermal_forward(source)

  ! Internal modules
  use isothermal

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nIsothermalWalls
     call isothermal_patch_forward(isothermalWallPatch(i), isothermalWallData(i), source)
  end do

  return
end subroutine isothermal_forward


! ===================================================== !
! Add the isothermal wall source during the adjoint run !
! ===================================================== !
subroutine isothermal_adjoint(source)

  ! Internal modules
  use isothermal

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nIsothermalWalls
     call isothermal_patch_adjoint(isothermalWallPatch(i), isothermalWallData(i), source)
  end do

  return
end subroutine isothermal_adjoint
