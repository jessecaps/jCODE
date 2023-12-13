module state

  ! External modules
  use precision
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Conserved variables
  real(WP), dimension(:,:), allocatable :: conservedVariables, targetState, adjointVariables

  ! Dependent/transport variables
  real(WP), dimension(:,:), allocatable :: specificVolume, velocity, velocityGradient,       &
       pressure, temperature, massFraction, dynamicViscosity, secondCoefficientOfViscosity,  &
       thermalDiffusivity, massDiffusivity, volumeFraction, mixtureMolecularWeight

  ! Fluxes
  real(WP), dimension(:,:), allocatable :: stressTensor, heatFlux, enthalpyFlux, reynoldsStress
  real(WP), dimension(:,:,:), allocatable :: speciesFlux

  ! Subgrid-scale terms (artificial diffusion)
  real(WP), dimension(:,:), allocatable :: turbulentViscosity, artificialShearViscosity,     &
       artificialBulkViscosity, artificialThermalDiffusivity

contains

  ! Allocate the state arrays
  ! -------------------------
  subroutine allocate_state

    allocate(conservedVariables(nGridPoints, nUnknowns))
    allocate(specificVolume(nGridPoints, 1))
    allocate(velocity(nGridPoints, nDimensions))
    allocate(velocityGradient(nGridPoints, nDimensions ** 2))
    allocate(pressure(nGridPoints, 1))
    allocate(temperature(nGridPoints, 1))
    if (nSpecies .gt. 0) then
       allocate(massFraction(nGridPoints, nSpecies))
       allocate(mixtureMolecularWeight(nGridPoints, 1))
    else
       allocate(massFraction(0, 0))
       allocate(mixtureMolecularWeight(0, 0))
    end if

    if (useTargetState) allocate(targetState(nGridPoints, nUnknowns))

    if (useViscosity) then
       allocate(dynamicViscosity(nGridPoints, 1))
       allocate(secondCoefficientOfViscosity(nGridPoints, 1))
       allocate(thermalDiffusivity(nGridPoints, 1))
       allocate(massDiffusivity(nGridPoints, nSpecies))
       allocate(stressTensor(nGridPoints, nDimensions ** 2))
       allocate(heatFlux(nGridPoints, nDimensions))
       if (nSpecies .gt. 0) then
          allocate(speciesFlux(nGridPoints, nSpecies, nDimensions))
          if (equationOfState .eq. IDEAL_GAS_MIXTURE)                                        &
               allocate(enthalpyFlux(nGridPoints, nDimensions))
       else
          allocate(speciesFlux(0, 0, 0))
          allocate(enthalpyFlux(0, 0))
       end if
       if (useLES) allocate(turbulentViscosity(nGridPoints, 1))
       if (useShockCapturing) then
          allocate(artificialShearViscosity(nGridPoints, 1))
          allocate(artificialBulkViscosity(nGridPoints, 1))
          allocate(artificialThermalDiffusivity(nGridPoints, 1))
       end if
    end if

    if (.not. predictionOnly) allocate(adjointVariables(nGridPoints, nUnknowns))

    if (useParticles) allocate(volumeFraction(nGridPoints, 1))
    if (usePTKE) allocate(reynoldsStress(nGridPoints, nDimensions ** 2))

    return
  end subroutine allocate_state


  ! Update the state
  ! ----------------
  subroutine update_state(stateVector)

    ! External modules
    use grid_functions
    use equation_of_state
    use state_functions
    use grid_functions

    implicit none

    ! Arguments
    real(WP), intent(in), optional :: stateVector(nGridPoints, nUnknowns)

    ! Local variables
    integer :: i, k

    ! Start the update state timer
    call timing_start('update')

    ! Adjust for volume fraction
    if (twoWayCoupling .and. .not.present(stateVector)) then
       call compute_particle_dependent_variables
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) / volumeFraction(:,1)
       end do
    end if

    ! Compute the dependent variables via the equation of state
    if (present(stateVector)) then
       call compute_dependent_variables(stateVector, specificVolume(:,1), velocity,          &
            massFraction, mixtureMolecularWeight(:,1), pressure(:,1), temperature(:,1))
    else
       call compute_dependent_variables(conservedVariables, specificVolume(:,1), velocity,   &
            massFraction, mixtureMolecularWeight(:,1), pressure(:,1), temperature(:,1))
    end if

    ! Store the velocity gradient
    call gradient(velocity, velocityGradient)

    if (useViscosity) then

       ! Large eddy simulation
       if (useLES .and. .not.present(stateVector)) then
          call les_compute_viscosity(conservedVariables(:,1), velocityGradient,              &
               turbulentViscosity(:,1))
       end if

       ! Shock capturing model
       if (useShockCapturing .and. .not.present(stateVector)) then
          call local_artificial_diffusivity(conservedVariables(:,1), pressure(:,1),          &
               temperature(:,1), velocityGradient, artificialShearViscosity(:,1),            &
               artificialBulkViscosity(:,1), artificialThermalDiffusivity(:,1))
       end if

       ! Compute the transport variables (viscosity, diffusivity, etc.)
       if (present(stateVector)) then
          call compute_transport_variables(temperature(:,1), dynamicViscosity(:,1),          &
               secondCoefficientOfViscosity(:,1), thermalDiffusivity(:,1), massDiffusivity)
       else
          call compute_transport_variables(temperature(:,1), dynamicViscosity(:,1),          &
               secondCoefficientOfViscosity(:,1), thermalDiffusivity(:,1), massDiffusivity,  &
               turbulentViscosity(:,1), artificialShearViscosity(:,1),                       &
               artificialBulkViscosity(:,1), artificialThermalDiffusivity(:,1))
       end if

       ! Compute the stress tensor fluxes
       call compute_stress_tensor(velocityGradient, dynamicViscosity(:,1),                   &
            secondCoefficientOfViscosity(:,1), stressTensor)

       ! Compute the heat flux
       call gradient(temperature(:,1), heatFlux)
       do i = 1, nDimensions
          heatFlux(:,i) = - thermalDiffusivity(:,1) * heatFlux(:,i)
       end do

       ! Compute the species mass flux
       do k = 1, nSpecies
          call gradient(massFraction(:,k), speciesFlux(:,k,:))
          do i = 1, nDimensions
             speciesFlux(:,k,i) = - massDiffusivity(:,k) * speciesFlux(:,k,i)
          end do
       end do

       ! Compute the enthalpy flux
       if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
          enthalpyFlux = 0.0_WP
          do i = 1, nDimensions
             do k = 1, nSpecies
                enthalpyFlux(:,i) = enthalpyFlux(:,i) + speciesFlux(:,k,i) *                 &
                     temperature(:,1) / schmidtNumberInverse(k) *                            &
                     ( molecularWeightInverse(k) * schmidtNumberInverse(k) -                 &
                     molecularWeightInverse(nSpecies+1) * schmidtNumberInverse(nSpecies+1) )
             end do
          end do
       end if

    end if

    ! Multiply volume fraction back
    if (twoWayCoupling .and. .not.present(stateVector)) then
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) * volumeFraction(:,1)
       end do
    end if

    ! Compute pseudo-turbulent Reynolds stress
    if (usePTKE .and. .not.present(stateVector)) call compute_pseudo_reynolds_stress

    ! Stop the update state timer
    call timing_stop('update')

    return
  end subroutine update_state


  ! Make a quiescent flow
  ! ----------------------
  subroutine make_quiescent(stateVector)

    implicit none

    ! Arguments
    real(WP), intent(out) :: stateVector(nGridPoints, nUnknowns)

    ! Local variables
    integer :: k

    stateVector(:,1) = 1.0_WP ! ... density
    stateVector(:,2:nDimensions+1) = 0.0_WP ! ... momentum
    stateVector(:,nDimensions+2) = 1.0_WP / ratioOfSpecificHeats /                           &
         (ratioOfSpecificHeats - 1.0_WP) ! ... energy
    do k = 1, nSpecies
       stateVector(:,nDimensions+2+k) = 0.0_WP ! ... density fraction
    end do

    if (useParticles) volumeFraction = 1.0_WP

    return
  end subroutine make_quiescent


  ! Correct state after RHS has been applied
  ! ----------------------------------------
  subroutine correct_state(stateVector)

    ! External modules
    use parallel
    use grid, only : gridNorm

    implicit none

    ! Arguments
    real(WP), intent(inout) :: stateVector(nGridPoints, nUnknowns)

    ! Local variables
    integer :: i, k
    real(WP) :: rho, Yk

    ! Return if adjoint is used
    if (.not. predictionOnly) return

    ! Enforce scalar boundedness
    if (boundScalars) then
       do i = 1, nGridPoints
          rho = stateVector(i, 1)
          do k = 1, nSpecies
             Yk = min(maxSpecies, max(minSpecies, stateVector(i,nDimensions+2+k) / rho))
             stateVector(i,nDimensions+2+k) = rho * Yk
          end do
       end do
    end if

    ! Mass correction
    if (useMassCorrection) then
       rho = sum(stateVector(:,1) * gridNorm(:, 1))
       call parallel_sum(rho); rho = rho / globalGridVolume
       do i = 1, nGridPoints
          stateVector(i, 1) = stateVector(i, 1) * targetDensity / rho
       end do
    end if

    ! Strongly-imposed (non-SAT) boundary conditions
    call boundary_correct_state(stateVector)

    ! Overwrite ghost points if using IBM
    if (useIBM) call ibm_ghost_point_correct_state(stateVector)
    if (useIBM) call ibm_filter_state(stateVector)

    return
  end subroutine correct_state


  ! Add the viscous / inviscid fluxes during forward run
  ! ----------------------------------------------------
  subroutine add_fluxes_forward(source)

    ! External modules
    use grid
    use first_derivative
    use grid_functions
    use state_jacobian
    use state_functions

    implicit none

    ! Arguments
    real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

    ! Local variables
    integer :: i, j
    real(WP), dimension(nGridPoints, nUnknowns, nDimensions) :: fluxes1, fluxes2

    ! Start the fluxes timer
    call timing_start('fluxes')

    ! Adjust for volume fraction
    if (twoWayCoupling) then
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) / volumeFraction(:,1)
       end do
    end if

    ! Inviscid fluxes
    ! ---------------------------------------
    if (useUpwinding .and. (.not.useSplitConvection)) then
       ! Perform vector flux splitting if upwinding is used
       call compute_upwind_inviscid_source(conservedVariables, specificVolume(:,1),          &
            pressure(:,1), velocity, volumeFraction(:,1), source, 1.0_WP)
       fluxes1 = 0.0_WP
    else if (useUpwinding .and. useSplitConvection) then
       fluxes1 = 0.0_WP
       call compute_split_inviscid_source(conservedVariables,                                &
            specificVolume(:,1), pressure(:,1), velocity, volumeFraction(:,1), fluxes1,      &
            source)
       ! Perform vector flux splitting if upwinding is used
       call compute_upwind_inviscid_source(conservedVariables, specificVolume(:,1),          &
            pressure(:,1), velocity, volumeFraction(:,1), source, 0.5_WP)
    else
       ! Compute Cartesian form of inviscid fluxes
       call compute_cartesian_inviscid_fluxes(conservedVariables, velocity, pressure(:,1),   &
            fluxes1)

       ! Modify fluxes and add additional source terms if convective splitting is used
       if (useSplitConvection) call compute_split_inviscid_source(conservedVariables,        &
            specificVolume(:,1), pressure(:,1), velocity, volumeFraction(:,1), fluxes1,      &
            source)
    end if

    ! Viscous fluxes
    ! ------------------------------------------
    if (useViscosity) then
       ! Compute Cartesian form of viscous fluxes if viscous terms are included
       call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,               &
            fluxes2, speciesFlux, enthalpyFlux)

       ! Update Cartesian form of total fluxes (strong form)
       if (.not. useSplitViscosity) fluxes1 = fluxes1 - fluxes2

       ! Send viscous fluxes to appropriate patches
       call boundary_store_viscous_fluxes(fluxes2)
    end if

    ! Compute turbulent fluxes if LES is used
    if (usePTKE) then
       call compute_cartesian_turbulent_fluxes(velocity, reynoldsStress, fluxes2)
       fluxes1 = fluxes1 + fluxes2
    end if

    ! Multiply volume fraction back
    if (twoWayCoupling) then
       do i = 1, nUnknowns
          conservedVariables(:, i) = conservedVariables(:, i) * volumeFraction(:,1)
       end do
       do j = 1, nDimensions
          do i = 1, nUnknowns
             fluxes1(:,i,j) = fluxes1(:,i,j) * volumeFraction(:,1)
          end do
       end do
    end if

    ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian
    ! form of total fluxes... upon return, `fluxes2` has the contravariant form
    call transform_fluxes(fluxes1, fluxes2)

    ! Take derivatives of fluxes
    do i = 1, nDimensions
       call first_derivative_apply(i, fluxes2(:,:,i))
    end do
    source = source - sum(fluxes2, dim = 3)

    ! Multiply by Jacobian
    do i = 1, nUnknowns
       source(:,i) = source(:,i) * jacobian(:,1)
    end do

    ! Non-conservative (split) Laplacian-form of the viscous fluxes
    if (useSplitViscosity) call compute_split_viscous_source(velocity, velocityGradient,     &
         stressTensor, dynamicViscosity(:,1), secondCoefficientOfViscosity(:,1),             &
         temperature(:,1), thermalDiffusivity(:,1), heatFlux, massFraction, speciesFlux,     &
         massDiffusivity, volumeFraction(:,1), source)

    ! Stop the fluxes timer
    call timing_stop('fluxes')

    return
  end subroutine add_fluxes_forward


  ! Add the viscous / inviscid fluxes during adjoint run
  ! ----------------------------------------------------
  subroutine add_fluxes_adjoint(source)

    ! External modules
    use grid
    use first_derivative
    use grid_functions
    use state_jacobian
    use state_functions

    implicit none

    ! Arguments
    real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

    ! Local variables
    integer :: i, j, k
    real(WP), allocatable :: temp1(:,:,:), temp2(:,:), localFluxJacobian1(:,:),              &
         localFluxJacobian2(:,:), localConservedVariables(:), localVelocity(:),              &
         localMassFraction(:), localMetricsAlongDirection1(:),                               &
         localMetricsAlongDirection2(:), localStressTensor(:), localHeatFlux(:),             &
         localEnthalpyFlux(:), localSpeciesFlux(:,:), localAdjointDiffusion(:,:)

    ! Start the fluxes timer
    call timing_start('fluxes')

    ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates
    allocate(temp1(nGridPoints, nUnknowns, nDimensions))
    do i = 1, nDimensions
       temp1(:,:,i) = adjointVariables
       call adjoint_first_derivative_apply(i, temp1(:,:,i))
    end do

    allocate(localFluxJacobian1(nUnknowns, nUnknowns))
    allocate(localConservedVariables(nUnknowns))
    allocate(localVelocity(nDimensions))
    if (nSpecies .gt. 0) allocate(localMassFraction(nSpecies))
    allocate(localMetricsAlongDirection1(nDimensions))

    if (useViscosity) then
       allocate(localFluxJacobian2(nUnknowns, nUnknowns))
       allocate(localStressTensor(nDimensions ** 2))
       allocate(localHeatFlux(nDimensions))
       if (nSpecies .gt. 0) then
          allocate(localSpeciesFlux(nSpecies,nDimensions))
          if (equationOfState .eq. IDEAL_GAS_MIXTURE)                                     &
               allocate(localEnthalpyFlux(nDimensions))
       end if
    end if

    do j = 1, nGridPoints

       localConservedVariables = conservedVariables(j,:)
       localVelocity = velocity(j,:)
       if (nSpecies .gt. 0) localMassFraction = massFraction(j,:)
       if (useViscosity) then
          localStressTensor = stressTensor(j,:)
          localHeatFlux = heatFlux(j,:)
          if (nSpecies .gt. 0) then
             localSpeciesFlux = speciesFlux(j,:,:)
             if (equationOfState .eq. IDEAL_GAS_MIXTURE)                                  &
                  localEnthalpyFlux = enthalpyFlux(j,:)
          end if
       end if

       do i = 1, nDimensions

          localMetricsAlongDirection1 = metrics(j,1+nDimensions*(i-1):nDimensions*i)

          call compute_jacobian_of_inviscid_flux(localConservedVariables,                 &
               localMetricsAlongDirection1, localFluxJacobian1,                           &
               specificVolume = specificVolume(j,1), velocity = localVelocity,            &
               pressure = pressure(j,1), massFraction = localMassFraction)

          if (useViscosity) then
             call compute_first_partial_viscous_jacobian(localConservedVariables,         &
                  localMetricsAlongDirection1, localStressTensor, localHeatFlux,          &
                  localEnthalpyFlux, localSpeciesFlux, localFluxJacobian2,                &
                  specificVolume(j,1), localVelocity,temperature(j,1), localMassFraction)
             localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
          end if

          source(j,:) = source(j,:) + matmul(transpose(localFluxJacobian1), temp1(j,:,i))

       end do

    end do

    if (allocated(localConservedVariables)) deallocate(localConservedVariables)
    if (allocated(localFluxJacobian1)) deallocate(localFluxJacobian1)
    if (allocated(localHeatFlux)) deallocate(localHeatFlux)
    if (allocated(localEnthalpyFlux)) deallocate(localEnthalpyFlux)
    if (allocated(localStressTensor)) deallocate(localStressTensor)
    if (allocated(localMassFraction)) deallocate(localMassFraction)
    if (allocated(localSpeciesFlux)) deallocate(localSpeciesFlux)
    if (allocated(localFluxJacobian2)) deallocate(localFluxJacobian2)

    if (useViscosity) then

       allocate(temp2(nGridPoints, nUnknowns - 1)); temp2 = 0.0_WP

       allocate(localMetricsAlongDirection2(nDimensions))
       allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
       allocate(localAdjointDiffusion(nUnknowns - 1, nDimensions))

       do k = 1, nGridPoints

          localVelocity = velocity(k,:)
          localAdjointDiffusion = 0.0_WP

          do j = 1, nDimensions

             localMetricsAlongDirection2 = metrics(k,1+nDimensions*(j-1):nDimensions*j)

             do i = 1, nDimensions

                localMetricsAlongDirection1 = metrics(k,1+nDimensions*(i-1):nDimensions*i)

                call compute_second_partial_viscous_jacobian(localVelocity,               &
                     dynamicViscosity(k,1), secondCoefficientOfViscosity(k,1),            &
                     thermalDiffusivity(k,1), massDiffusivity(k,:), temperature(k,1),     &
                     jacobian(k,1), localMetricsAlongDirection1,                          &
                     localMetricsAlongDirection2, localFluxJacobian2)

                localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) +                 &
                     matmul(transpose(localFluxJacobian2), temp1(k,2:nUnknowns,i))

             end do

          end do

          do j = 1, nDimensions
             temp1(k,2:nUnknowns,j) = localAdjointDiffusion(:,j)
          end do

       end do

       do j = 1, nDimensions
          call adjoint_first_derivative_apply(j, temp1(:,2:nUnknowns,j))
       end do
       temp2 = sum(temp1(:,2:nUnknowns,:), dim = 3) !... divergence of the adjoint flux

       select case (equationOfState)

       case (IDEAL_GAS)

          temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *           &
               temp2(:,nDimensions+1)
          do i = 1, nDimensions
             temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *              &
                  temp2(:,nDimensions+1)
          end do
          do k = 1, nSpecies
             temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k)
          end do

       case (IDEAL_GAS_MIXTURE)

          temp2(:,nDimensions+1) = ratioOfSpecificHeats * specificVolume(:,1) *           &
               mixtureMolecularWeight(:,1) * temp2(:,nDimensions+1)
          do i = 1, nDimensions
             temp2(:,i) = specificVolume(:,1) * temp2(:,i) - velocity(:,i) *              &
                  temp2(:,nDimensions+1)
          end do
          do k = 1, nSpecies
             temp2(:,nDimensions+1+k) = specificVolume(:,1) * temp2(:,nDimensions+1+k) +  &
                  temperature(:,1) * (molecularWeightInverse(nSpecies+1) -                &
                  molecularWeightInverse(k)) * temp2(:,nDimensions+1) /                   &
                  ratioOfSpecificHeats
          end do

       end select

       source(:,2:nUnknowns) = source(:,2:nUnknowns) - temp2
       source(:,1) = source(:,1) + specificVolume(:,1) *                                  &
            conservedVariables(:,nDimensions+2) * temp2(:,nDimensions+1) +                &
            sum(velocity * temp2(:,1:nDimensions), dim = 2)
       if (nSpecies .gt. 0) source(:,1) = source(:,1) +                                   &
            sum(massFraction * temp2(:,ndimensions+2:nUnknowns-1), dim = 2)

       deallocate(temp2)

       deallocate(localAdjointDiffusion)
       deallocate(localFluxJacobian2)
       deallocate(localMetricsAlongDirection2)

    end if

    deallocate(localMetricsAlongDirection1)
    deallocate(localVelocity)
    deallocate(temp1)

    ! Multiply by Jacobian
    do i = 1, nUnknowns
       source(:,i) = source(:,i) * jacobian(:,1)
    end do

    ! Stop the fluxes timer
    call timing_stop('fluxes')

    return
  end subroutine add_fluxes_adjoint


  ! Add dissipation to the forward/adjoint solution
  ! -----------------------------------------------
  subroutine add_dissipation(mode, source)

    ! External modules
    use geometry
    use grid
    use dissipation
    use first_derivative
    use state_functions

    implicit none

    ! Arguments
    integer, intent(in) :: mode
    real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

    ! Local variables
    integer :: i, j
    real(WP), allocatable :: dissipationTerm(:,:)

    if (.not. useDissipation) return

    ! Start the dissipation timer
    call timing_start('dissipation')

    allocate(dissipationTerm(nGridPoints, nUnknowns))

    ! Dissipation sensor
    if (hybridDissipation) call compute_sensor(dissSensorType, pressure(:,1),                &
         specificVolume(:,1), velocity, velocityGradient, massFraction,                      &
         dissipationSensor(:,1))

    do i = 1, nDimensions

       select case (mode)
       case (FORWARD)
          if (twoWayCoupling) then
             do j = 1, nUnknowns
                dissipationTerm(:,j) = conservedVariables(:,j) / volumeFraction(:,1)
             end do
          else
             dissipationTerm = conservedVariables
          end if
       case (ADJOINT)
          dissipationTerm = adjointVariables
       end select

       call dissipation_apply(i, dissipationTerm, dissipationSensor(:,1))

       if (.not. compositeDissipation) then
          do j = 1, nUnknowns
             dissipationTerm(:,j) = -arcLengths(:,i) * dissipationTerm(:,j)
          end do
          call dissipation_transpose_apply(i, dissipationTerm)
          call first_derivative_apply_norm_inverse(i, dissipationTerm)
       end if

       ! Multiply by Jacobian
       if (twoWayCoupling) then
          do j = 1, nUnknowns
             dissipationTerm(:,j) = dissipationTerm(:,j) * jacobian(:,1) * volumeFraction(:,1)
          end do
       else
          do j = 1, nUnknowns
             dissipationTerm(:,j) = dissipationTerm(:,j) * jacobian(:,1)
          end do
       end if

       select case (mode)
       case (FORWARD)
          source = source + dissipationAmount * dissipationTerm
       case (ADJOINT)
          source = source - dissipationAmount * dissipationTerm
       end select

       ! Store for monitor
       if (all(isPeriodic(1:nDimensions))) then
          dissipationSource = dissipationAmount * dissipationTerm
       end if

    end do

    ! Cleanup
    deallocate(dissipationTerm)

    ! Stop the dissipation timer
    call timing_stop('dissipation')

    return
  end subroutine add_dissipation


  ! Filter the solution at a given timestep
  ! ---------------------------------------
  subroutine filter_solution(f, timestep)

    ! External modules
    use filter

    implicit none

    ! Arguments
    real(WP), intent(inout) :: f(:,:)
    integer, intent(in) :: timestep

    ! Local variables
    integer :: i, j
    integer, allocatable :: directions(:)

    if (.not. useFilter) return

    ! Start the filter timer
    call timing_start('filter')

    ! Alternate which direction gets filtered first each timestep
    select case (nDimensions)
    case (1)
       allocate(directions(1))
       directions(:) = (/ 1 /)
    case (2)
       allocate(directions(2))
       directions(:) = (/ 12, 21 /)
    case (3)
       allocate(directions(6))
       directions(:) = (/ 123, 231, 312, 132, 321, 213 /)
    end select

    do i = 1, nDimensions
       j = mod(directions(mod(timestep, size(directions)) + 1) / 10 ** (i - 1), 10)
       call filter_apply(j, f)
    end do

    deallocate(directions)

    ! Stop the filter timer
    call timing_stop('filter')

    return
  end subroutine filter_solution

end module state


! =============== !
! Setup the state !
! =============== !
subroutine state_setup

  ! Internal modules
  use state

  ! External modules
  use string
  use parallel

  implicit none

  ! Allocate the state arrays
  call allocate_state

  ! Make queiscent flow
  call make_quiescent(conservedVariables)

  ! LES
  call les_setup

  ! Combustion
  call combustion_setup

  ! Setup various source terms
  call sources_setup

  ! Setup the immersed boundary routine
  call ibm_setup

  return
end subroutine state_setup


! ================= !
! Cleanup the state !
! ================= !
subroutine state_cleanup

  ! Internal modules
  use state

  implicit none

  call ibm_cleanup
  call sources_cleanup
  call combustion_cleanup
  call les_cleanup

  if (allocated(conservedVariables)) deallocate(conservedVariables)
  if (allocated(targetState)) deallocate(targetState)
  if (allocated(adjointVariables)) deallocate(adjointVariables)
  if (allocated(specificVolume)) deallocate(specificVolume)
  if (allocated(velocity)) deallocate(velocity)
  if (allocated(massFraction)) deallocate(massFraction)
  if (allocated(pressure)) deallocate(pressure)
  if (allocated(temperature)) deallocate(temperature)
  if (allocated(dynamicViscosity)) deallocate(dynamicViscosity)
  if (allocated(secondCoefficientOfViscosity)) deallocate(secondCoefficientOfViscosity)
  if (allocated(turbulentViscosity)) deallocate(turbulentViscosity)
  if (allocated(artificialShearViscosity)) deallocate(artificialShearViscosity)
  if (allocated(artificialBulkViscosity)) deallocate(artificialBulkViscosity)
  if (allocated(artificialThermalDiffusivity)) deallocate(artificialThermalDiffusivity)
  if (allocated(thermalDiffusivity)) deallocate(thermalDiffusivity)
  if (allocated(massDiffusivity)) deallocate(massDiffusivity)
  if (allocated(velocityGradient)) deallocate(velocityGradient)
  if (allocated(stressTensor)) deallocate(stressTensor)
  if (allocated(heatFlux)) deallocate(heatFlux)
  if (allocated(enthalpyFlux)) deallocate(enthalpyFlux)
  if (allocated(speciesFlux)) deallocate(speciesFlux)
  if (allocated(volumeFraction)) deallocate(volumeFraction)
  if (allocated(mixtureMolecularWeight)) deallocate(mixtureMolecularWeight)
  if (allocated(reynoldsStress)) deallocate(reynoldsStress)

  return
end subroutine state_cleanup


! =========================== !
! Compute the right-hand side !
! =========================== !
subroutine state_rhs(mode, source)

  ! Internal modules
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: mode
  real(WP), dimension(nGridPoints,nUnknowns), intent(out) :: source

  ! Reset the right-hand side
  source = 0.0_WP

  ! Compute the fluxes
  if (mode .eq. FORWARD) then
     call add_fluxes_forward(source)
  else
     call add_fluxes_adjoint(source)
  end if

  ! Add forcing on actuator patches
  if (mode .eq. FORWARD) call controller_add_source(mode, source)

  ! Add adjoint forcing on cost target patches
  if (mode .eq. ADJOINT) call functional_adjoint_source(source)

  ! Combustion
  call combustion_source(mode, source)

  ! Particles
  if (mode .eq. FORWARD) then
     call particle_source_forward(source)
  else if (mode .eq. ADJOINT) then
     call particle_source_adjoint(source)
  end if

  ! Other source terms
  call add_sources(mode, source)

  ! Add boundary patch penalties
  call boundary_sources(mode, source)

  ! Immersed boundary forcing (if using weak form)
  call ibm_ghost_point_source(source)

  ! Add dissipation if required
  call add_dissipation(mode, source)

  return
end subroutine state_rhs
