module adiabatic

  ! External modules
  use boundary
  use impenetrable, only : t_ImpenetrableWall
  use precision

  implicit none

  integer :: nAdiabaticWalls

  type, private :: t_AdiabaticWall
     type(t_ImpenetrableWall) :: impenetrableWallData
     real(WP) :: viscousPenaltyAmounts(2)
  end type t_AdiabaticWall

  type(t_AdiabaticWall), allocatable :: adiabaticWallData(:)
  type(t_Patch), pointer :: adiabaticWallPatch(:)

contains

  subroutine setup_adiabatic_patch(patch, bc)

    ! External modules
    use parser
    use simulation_flags
    use solver_options
    use first_derivative
    use impenetrable, only : setup_impenetrable_patch 

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_AdiabaticWall), intent(inout) :: bc

    ! Local variables
    integer :: i

    ! Verify the patch type
    call verify_adiabatic_patch(patch)

    ! Adiabatic wall requires impenetrable wall data
    call setup_impenetrable_patch(patch, bc%impenetrableWallData)

    ! Viscous penalty amount
    if (useViscosity) then
       call parser_read(trim(patch%name) // ' viscous penalty amount',                       &
            bc%viscousPenaltyAmounts(1), defaultViscousPenaltyAmount)
       bc%viscousPenaltyAmounts(2) = 1.0_WP
       do i = 1, 2
          bc%viscousPenaltyAmounts(i) = bc%viscousPenaltyAmounts(i) /                        &
               firstDerivative(abs(patch%normalDirection))%normBoundary(1)
       end do
       bc%viscousPenaltyAmounts(1) = bc%viscousPenaltyAmounts(1) * reynoldsNumberInverse
       bc%viscousPenaltyAmounts(2) = sign(bc%viscousPenaltyAmounts(2),                       &
            real(patch%normalDirection, WP))
    else
       bc%viscousPenaltyAmounts = 0.0_WP
    end if

    return
  end subroutine setup_adiabatic_patch


  subroutine cleanup_adiabatic_patch(patch, bc)

    ! External modules
    use impenetrable, only : cleanup_impenetrable_patch

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_AdiabaticWall), intent(inout) :: bc

    call cleanup_impenetrable_patch(patch, bc%impenetrableWallData)

    return
  end subroutine cleanup_adiabatic_patch


  subroutine verify_adiabatic_patch(patch)

    ! External modules
    use simulation_flags
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    if (patch%patchType .ne. SAT_ADIABATIC_WALL)                                             &
         call die('verify_adiabatic_patch: patch type mismatch for SAT_ADIABATIC_WALL')

    return
  end subroutine verify_adiabatic_patch


  subroutine adiabatic_patch_forward(patch, bc, source)

    ! External modules
    use solver_options
    use geometry
    use grid
    use grid_patch
    use state
    use impenetrable, only : impenetrable_patch_forward

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_AdiabaticWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: metricsAlongNormalDirection(:), temp(:,:)

    ! Enforce the no-penetration condition
    call impenetrable_patch_forward(patch, bc%impenetrableWallData, source)

    if (.not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(metricsAlongNormalDirection(nDimensions))
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

             ! No-slip condition
             source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) -         &
                  bc%viscousPenaltyAmounts(1) * jacobian(gridIndex, 1) *                     &
                  conservedVariables(gridIndex,2:nDimensions+1)

             temp = 0.0_WP
             select case (nDimensions)

             case (1)
                temp(1:2,1) = 0.0_WP
                temp(3,1) = - heatFlux(gridIndex,1)
                do l = 1, nSpecies
                   temp(l+3,1) = - speciesFlux(gridIndex,l,1)
                end do

             case (2)
                temp(1:3,:) = 0.0_WP
                temp(4,1) = - heatFlux(gridIndex,1)
                temp(4,2) = - heatFlux(gridIndex,2)
                do l = 1, nSpecies
                   temp(l+4,1) = - speciesFlux(gridIndex,l,1)
                   temp(l+4,2) = - speciesFlux(gridIndex,l,2)
                end do

             case (3)
                temp(1:4,:) = 0.0_WP
                temp(5,1) = - heatFlux(gridIndex,1)
                temp(5,2) = - heatFlux(gridIndex,2)
                temp(5,3) = - heatFlux(gridIndex,3)
                do l = 1, nSpecies
                   temp(l+5,1) = - speciesFlux(gridIndex,l,1)
                   temp(l+5,2) = - speciesFlux(gridIndex,l,2)
                   temp(l+5,3) = - speciesFlux(gridIndex,l,3)
                end do

             end select

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             ! Zero flux condition
             source(gridIndex,:) = source(gridIndex,:) + bc%viscousPenaltyAmounts(2) *       &
                  jacobian(gridIndex, 1) * matmul(temp, metricsAlongNormalDirection)

          end do
       end do
    end do

    deallocate(metricsAlongNormalDirection)
    deallocate(temp)

    return
  end subroutine adiabatic_patch_forward


  subroutine adiabatic_patch_adjoint(patch, bc, source)

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
    type(t_AdiabaticWall), intent(inout) :: bc
    real(WP), dimension(:,:), intent(inout) :: source

    ! Local variables
    integer :: i, j, k, l, direction, gridIndex, patchIndex
    real(WP), allocatable :: metricsAlongNormalDirection(:), metricsAlongDirection2(:),      &
         localViscousFluxJacobian(:,:), localViscousFluxJacobian2(:,:), localMassFraction(:),&
         localSpeciesFlux(:,:), temp1(:,:,:), temp2(:,:)

    ! Apply no-slip conditions
    call impenetrable_patch_adjoint(patch, bc%impenetrableWallData, source)

    if (.not. useViscosity) return

    direction = abs(patch%normalDirection)

    allocate(metricsAlongNormalDirection(nDimensions))
    allocate(metricsAlongDirection2(nDimensions))
    allocate(localViscousFluxJacobian(nUnknowns, nUnknowns))
    allocate(localViscousFluxJacobian2(nUnknowns-1, nUnknowns-1))
    if (nSpecies .gt. 0) then
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

             if (useContinuousAdjoint) then
                source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) -      &
                     jacobian(gridIndex, 1) * bc%viscousPenaltyAmounts(1) *                  &
                     adjointVariables(gridIndex,2:nDimensions+1)
             else
                source(gridIndex,2:nDimensions+1) = source(gridIndex,2:nDimensions+1) +      &
                     jacobian(gridIndex, 1) * bc%viscousPenaltyAmounts(1) *                  &
                     adjointVariables(gridIndex,2:nDimensions+1)
             end if

             metricsAlongNormalDirection =                                                   &
                  metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

             if (nSpecies .gt. 0) then
                localMassFraction = massFraction(gridIndex,:)
                localSpeciesFlux = speciesFlux(gridIndex,:,:)
             end if
             call compute_first_partial_fluxes_jacobian(conservedVariables(gridIndex,:),     &
                  metricsAlongNormalDirection, specificVolume(gridIndex,1),                  &
                  velocity(gridIndex,:), temperature(gridIndex,1), heatFlux(gridIndex,:),    &
                  localViscousFluxJacobian, localMassFraction, localSpeciesFlux, gridIndex)
             source(gridIndex,:) = source(gridIndex,:) - bc%viscousPenaltyAmounts(2) *       &
                  jacobian(gridIndex, 1) * matmul(transpose(localViscousFluxJacobian),       &
                  adjointVariables(gridIndex,:))

             if (useContinuousAdjoint) cycle

             do l = 1, nDimensions
                metricsAlongDirection2 = metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)
                call compute_second_partial_fluxes_jacobian(dynamicViscosity(gridIndex,1),   &
                     thermalDiffusivity(gridIndex,1), massDiffusivity(gridIndex,:),          &
                     jacobian(gridIndex,1), metricsAlongNormalDirection,                     &
                     metricsAlongDirection2, localViscousFluxJacobian2, gridIndex)
                temp1(gridIndex,:,l) = - matmul(transpose(localViscousFluxJacobian2),        &
                     adjointVariables(gridIndex,2:nUnknowns))
             end do
          end do
       end do
    end do

    deallocate(metricsAlongNormalDirection)
    deallocate(metricsAlongDirection2)
    deallocate(localViscousFluxJacobian)
    deallocate(localViscousFluxJacobian2)
    if (allocated(localMassFraction)) deallocate(localMassFraction)
    if (allocated(localSpeciesFlux)) deallocate(localSpeciesFlux)

    if (useContinuousAdjoint) return

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
       temp2(:,i) = bc%viscousPenaltyAmounts(2) * jacobian(:, 1) * temp2(:,i)
    end do

    source(:,2:nUnknowns) = source(:,2:nUnknowns) + temp2
    source(:,1) = source(:,1) - specificVolume(:,1) * conservedVariables(:,nDimensions+2)    &
         * temp2(:,nDimensions+1) - sum(velocity * temp2(:,1:nDimensions), dim = 2)
    if (nSpecies .gt. 0) source(:,1) = source(:,1) -                                         &
         sum(massFraction * temp2(:,nDimensions+2:nUnknowns-1), dim = 2)

    deallocate(temp2)

    return
  end subroutine adiabatic_patch_adjoint


  subroutine compute_first_partial_fluxes_jacobian(conservedVariables, metrics,              &
       specificVolume, velocity, temperature, heatFlux, firstPartialViscousJacobian,         &
       massFraction, speciesFlux, gridIndex)

    ! External modules
    use simulation_flags
    use solver_options

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:), metrics(:), specificVolume, velocity(:),  &
         temperature, heatFlux(:)
    integer, intent(in) :: gridIndex
    real(WP), intent(in), optional :: massFraction(:), speciesFlux(:,:)
    real(WP), intent(out) :: firstPartialViscousJacobian(:,:)

    ! Local variables
    integer :: k
    real(WP) :: phiSquared, contravariantHeatFlux, contravariantSpeciesFlux(nSpecies),       &
         deltaViscosity(nDimensions+nSpecies+2)

    ! Zero-out first partial viscous Jacobian
    firstPartialViscousJacobian = 0.0_WP

    select case (nDimensions)

    case (1)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * velocity(1) ** 2
       contravariantHeatFlux = metrics(1) * heatFlux(1) !... not normalized
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
       firstPartialViscousJacobian(3,1) = - deltaViscosity(1) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = - deltaViscosity(2) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = - deltaViscosity(3) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(3+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

    case (2)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                               &
            (velocity(1) ** 2 + velocity(2) ** 2)
       contravariantHeatFlux = metrics(1) * heatFlux(1) +                                    &
            metrics(2) * heatFlux(2) !... not normalized
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
       firstPartialViscousJacobian(4,1) = - deltaViscosity(1) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = 0.0_WP
       firstPartialViscousJacobian(4,2) = - deltaViscosity(2) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = 0.0_WP
       firstPartialViscousJacobian(4,3) = - deltaViscosity(3) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,4) = 0.0_WP
       firstPartialViscousJacobian(2,4) = 0.0_WP
       firstPartialViscousJacobian(3,4) = 0.0_WP
       firstPartialViscousJacobian(4,4) = - deltaViscosity(4) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(4+k,4) = - deltaViscosity(4) *                         &
               contravariantSpeciesFlux(k)
       end do

    case (3)

       ! Other dependent variables
       phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                               &
            (velocity(1) ** 2 + velocity(2) ** 2 + velocity(3) ** 2)
       contravariantHeatFlux = metrics(1) * heatFlux(1) + metrics(2) * heatFlux(2) +         &
            metrics(3) * heatFlux(3) !... not normalized
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
       firstPartialViscousJacobian(5,1) = - deltaViscosity(1) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,1) = - deltaViscosity(1) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,2) = 0.0_WP
       firstPartialViscousJacobian(2,2) = 0.0_WP
       firstPartialViscousJacobian(3,2) = 0.0_WP
       firstPartialViscousJacobian(4,2) = 0.0_WP
       firstPartialViscousJacobian(5,2) = - deltaViscosity(2) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,2) = - deltaViscosity(2) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,3) = 0.0_WP
       firstPartialViscousJacobian(2,3) = 0.0_WP
       firstPartialViscousJacobian(3,3) = 0.0_WP
       firstPartialViscousJacobian(4,3) = 0.0_WP
       firstPartialViscousJacobian(5,3) = - deltaViscosity(3) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,3) = - deltaViscosity(3) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,4) = 0.0_WP
       firstPartialViscousJacobian(2,4) = 0.0_WP
       firstPartialViscousJacobian(3,4) = 0.0_WP
       firstPartialViscousJacobian(4,4) = 0.0_WP
       firstPartialViscousJacobian(5,4) = - deltaViscosity(4) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,4) = - deltaViscosity(4) *                         &
               contravariantSpeciesFlux(k)
       end do

       firstPartialViscousJacobian(1,5) = 0.0_WP
       firstPartialViscousJacobian(2,5) = 0.0_WP
       firstPartialViscousJacobian(3,5) = 0.0_WP
       firstPartialViscousJacobian(4,5) = 0.0_WP
       firstPartialViscousJacobian(5,5) = - deltaViscosity(5) * contravariantHeatFlux
       do k = 1, nSpecies
          firstPartialViscousJacobian(5+k,5) = - deltaViscosity(5) *                         &
               contravariantSpeciesFlux(k)
       end do

    end select

    return
  end subroutine compute_first_partial_fluxes_jacobian


  subroutine compute_second_partial_fluxes_jacobian(dynamicViscosity, thermalDiffusivity,    &
       massDiffusivity, jacobian, metricsAlongFirstDir, metricsAlongSecondDir,               &
       secondPartialViscousJacobian, gridIndex)

    ! External modules
    use simulation_flags
    use solver_options

    implicit none

    ! Arguments
    real(WP), intent(in) :: dynamicViscosity, thermalDiffusivity, jacobian,                  &
         metricsAlongFirstDir(:), metricsAlongSecondDir(:)
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

       secondPartialViscousJacobian(2,2) = thermalDiffusivity * temp
       do k = 1, nSpecies
          secondPartialViscousJacobian(2+k,2+k) = massDiffusivity(k) * temp
       end do

    case (2)

       ! Temporary variables
       temp = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                           &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2)

       secondPartialViscousJacobian(3,3) = thermalDiffusivity * temp
       do k = 1, nSpecies
          secondPartialViscousJacobian(3+k,3+k) = massDiffusivity(k) * temp
       end do

    case (3)

       ! Temporary variables
       temp = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                           &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2) +                             &
            metricsAlongFirstDir(3) * metricsAlongSecondDir(3)

       secondPartialViscousJacobian(4,4) = thermalDiffusivity * temp
       do k = 1, nSpecies
          secondPartialViscousJacobian(4+k,4+k) = massDiffusivity(k) * temp
       end do

    end select

    ! Multiply by the Jacobian
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

    return
  end subroutine compute_second_partial_fluxes_jacobian

end module adiabatic


! ==================================== !
! Setup the adiabatic boundary patches !
! ==================================== !
subroutine adiabatic_setup

  ! Internal modules
  use adiabatic

  ! External modules
  use parser

  implicit none

  ! Local variables
  integer :: i, j

  ! Find the number of boundaries of this type
  nAdiabaticWalls = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. SAT_ADIABATIC_WALL) then
        nAdiabaticWalls = nAdiabaticWalls + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nAdiabaticWalls .eq. 0) return

  ! Allocate the adiabatic wall type
  allocate(adiabaticWallData(nAdiabaticWalls))

  ! Connect the boundary patch
  adiabaticWallPatch => patches(j:j+nAdiabaticWalls-1)

  ! Setup the boundary conditions
  do i = 1, nAdiabaticWalls
     call setup_adiabatic_patch(adiabaticWallPatch(i), adiabaticWallData(i))
  end do

  return
end subroutine adiabatic_setup


! ====================================== !
! Cleanup the adiabatic boundary patches !
! ====================================== !
subroutine adiabatic_cleanup

  ! Internal modules
  use adiabatic

  ! Local variables
  integer :: i

  if (nAdiabaticWalls .gt. 0) then
     do i = 1, nAdiabaticWalls
        call cleanup_adiabatic_patch(adiabaticWallPatch(i), adiabaticWallData(i))
     end do
     deallocate(adiabaticWallData)
     nullify(adiabaticWallPatch)
  end if

  nAdiabaticWalls = 0

  return
end subroutine adiabatic_cleanup


! =============================================== !
! Add the adiabatic source during the forward run !
! =============================================== !
subroutine adiabatic_forward(source)

  ! Internal modules
  use adiabatic

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints,nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nAdiabaticWalls
     call adiabatic_patch_forward(adiabaticWallPatch(i), adiabaticWallData(i), source)
  end do

  return
end subroutine adiabatic_forward


! =============================================== !
! Add the adiabatic source during the adjoint run !
! =============================================== !
subroutine adiabatic_adjoint(source)

  ! Internal modules
  use adiabatic

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i

  do i = 1, nAdiabaticWalls
     call adiabatic_patch_adjoint(adiabaticWallPatch(i), adiabaticWallData(i), source)
  end do

  return
end subroutine adiabatic_adjoint
