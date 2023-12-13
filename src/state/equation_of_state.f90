module equation_of_state

  ! External modules
  use precision
  use simulation_flags
  use solver_options
  use geometry

  implicit none

contains

  ! Compute the dependent variables via the equation of state
  ! ---------------------------------------------------------
  subroutine compute_dependent_variables(stateVector, specificVolume, velocity,              &
       massFraction, mixtureMolecularWeight, pressure, temperature)

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:)
    real(WP), intent(out), optional :: specificVolume(:), velocity(:,:), massFraction(:,:),  &
         mixtureMolecularWeight(:), pressure(:), temperature(:)

    ! Local variables
    real(WP), allocatable :: temp(:)

    ! Local variables
    integer :: i

    ! Specific volume
    if (present(specificVolume)) then
       specificVolume = 1.0_WP / stateVector(:,1)
    end if

    ! Velocity
    if (present(velocity)) then
       if (present(specificVolume)) then
          do i = 1, nDimensions
             velocity(:,i) = specificVolume * stateVector(:,i+1)
          end do
       else
          do i = 1, nDimensions
             velocity(:,i) =stateVector(:,i+1) / stateVector(:,1)
          end do
       end if
    end if

    ! Mass fraction
    if (present(massFraction) .and. nSpecies .gt. 0) then
       if (present(specificVolume)) then
          do i = 1, nSpecies
             massFraction(:,i) = specificVolume * stateVector(:,nDimensions+2+i)
          end do
       else
          do i = 1, nSpecies
             massFraction(:,i) = stateVector(:,nDimensions+2+i) / stateVector(:,1)
          end do
       end if
    end if

    ! Mixture molecular weight
    if (present(mixtureMolecularWeight) .and. equationOfState .eq. IDEAL_GAS_MIXTURE) then
       mixtureMolecularWeight = molecularWeightInverse(nSpecies+1)
       if (present(massFraction)) then
          do i = 1, nSpecies
             mixtureMolecularWeight = mixtureMolecularWeight + massFraction(:,i) *           &
                  (molecularWeightInverse(i) - molecularWeightInverse(nSpecies+1))
          end do
       else
          if (present(specificVolume)) then
             do i = 1, nSpecies
                mixtureMolecularWeight = mixtureMolecularWeight +                            &
                     specificVolume * stateVector(:,nDimensions+2+i) *                       &
                     (molecularWeightInverse(i) - molecularWeightInverse(nSpecies+1))
             end do
          else
             do i = 1, nSpecies
                mixtureMolecularWeight = mixtureMolecularWeight +                            &
                     stateVector(:,nDimensions+2+i) / stateVector(:,1) *                     &
                     (molecularWeightInverse(i) - molecularWeightInverse(nSpecies+1))
             end do
          end if
       end if
       mixtureMolecularWeight = 1.0_WP / mixtureMolecularWeight
    end if

    ! Pressure
    if (present(pressure)) then
       if (present(velocity)) then
          pressure = (ratioOfSpecificHeats - 1.0_WP) *                                       &
               (stateVector(:,nDimensions+2) - 0.5_WP * stateVector(:,1) *                   &
               sum(velocity ** 2, dim = 2))
       else if (present(specificVolume)) then
          pressure = (ratioOfSpecificHeats - 1.0_WP) *                                       &
               (stateVector(:,nDimensions+2) - 0.5_WP *                                      &
               sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) * specificVolume)
       else
          pressure = (ratioOfSpecificHeats - 1.0_WP) * (stateVector(:,nDimensions+2) -       &
               0.5_WP * sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) /                  &
               stateVector(:,1))
       end if
       ! Account for Reynolds stresses
       if (usePTKE) pressure = pressure - (ratioOfSpecificHeats - 1.0_WP) *                  &
            stateVector(:, nDimensions + 3)
    end if

    ! Temperature from the equation of state
    if (present(temperature)) then
       select case (equationOfState)

       case (IDEAL_GAS) ! ... compute temperature using ideal gas law
          if (present(pressure)) then
             if (present(specificVolume)) then
                temperature = ratioOfSpecificHeats * pressure /                              &
                     (ratioOfSpecificHeats - 1.0_WP) * specificVolume
             else
                temperature = ratioOfSpecificHeats * pressure /                              &
                     (ratioOfSpecificHeats - 1.0_WP) / stateVector(:,1)
             end if
          else
             temperature = ratioOfSpecificHeats * (stateVector(:,nDimensions+2) -            &
                  0.5_WP * sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) /               &
                  stateVector(:,1)) / stateVector(:,1)
          end if

       case (IDEAL_GAS_MIXTURE) ! ... compute temperature using multi-component ideal gas
          if (present(mixtureMolecularWeight)) then
             if (present(pressure)) then
                if (present(specificVolume)) then
                   temperature = ratioOfSpecificHeats * pressure * mixtureMolecularWeight /  &
                        (ratioOfSpecificHeats - 1.0_WP) * specificVolume
                else
                   temperature = ratioOfSpecificHeats * pressure * mixtureMolecularWeight /  &
                        (ratioOfSpecificHeats - 1.0_WP) / stateVector(:,1)
                end if
             else
                if (present(specificVolume)) then
                   temperature = ratioOfSpecificHeats * (stateVector(:,nDimensions+2) -      &
                        0.5_WP * sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) /         &
                        stateVector(:,1)) * mixtureMolecularWeight * specificVolume
                else
                   temperature = ratioOfSpecificHeats * (stateVector(:,nDimensions+2) -      &
                        0.5_WP * sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) /         &
                        stateVector(:,1)) * mixtureMolecularWeight / stateVector(:,1)
                end if
             end if
          else
             allocate(temp(size(temperature)))
             temp = stateVector(:,1) * molecularWeightInverse(nSpecies+1)
             do i = 1, nSpecies
                temp = temp + stateVector(:,nDimensions+2+i) *                               &
                     (molecularWeightInverse(i) - molecularWeightInverse(nSpecies+1))
             end do
             temp = 1.0_WP / temp
             if (present(pressure)) then
                temperature = ratioOfSpecificHeats * pressure /                              &
                     (ratioOfSpecificHeats - 1.0_WP) * temp
             else
                temperature = ratioOfSpecificHeats * (stateVector(:,nDimensions+2) -         &
                     0.5_WP * sum(stateVector(:,2:nDimensions+1) ** 2, dim = 2) /            &
                     stateVector(:,1)) * temp
             end if
             deallocate(temp)
          end if

       end select
    end if

    return
  end subroutine compute_dependent_variables


  ! Compute the transport variables
  ! -------------------------------
  subroutine compute_transport_variables(temperature, dynamicViscosity,                      &
       secondCoefficientOfViscosity, thermalDiffusivity, massDiffusivity,                    &
       turbulentViscosity, artificialShearViscosity, artificialBulkViscosity,                &
       artificialThermalDiffusivity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: temperature(:)
    real(WP), intent(in), optional :: turbulentViscosity(:), artificialShearViscosity(:),    &
         artificialBulkViscosity(:), artificialThermalDiffusivity(:)
    real(WP), intent(out), optional :: dynamicViscosity(:), secondCoefficientOfViscosity(:), &
         thermalDiffusivity(:), massDiffusivity(:,:)

    ! Local variables
    integer :: i
    real(WP), parameter :: twoThirds = 2.0_WP / 3.0_WP

    if (powerLawExponent .le. 0.0_WP) then !... handle powerLawExponent = 0 separately

       ! Dynamic viscosity
       if (present(dynamicViscosity)) then
          dynamicViscosity = reynoldsNumberInverse
       end if

       ! Second coefficient of viscosity
       if (present(secondCoefficientOfViscosity)) then
          secondCoefficientOfViscosity =                                                     &
               (bulkViscosityRatio - twoThirds) * reynoldsNumberInverse
       end if

       ! Thermal diffusivity
       if (present(thermalDiffusivity)) then
          thermalDiffusivity = reynoldsNumberInverse * prandtlNumberInverse
       end if

       ! Mass diffusivity
       if (present(massDiffusivity) .and. nSpecies .gt. 0) then
          do i = 1, nSpecies
             massDiffusivity(:,i) = reynoldsNumberInverse * schmidtNumberInverse(i)
          end do
       end if

    else

       ! Dynamic viscosity
       if (present(dynamicViscosity)) then
          dynamicViscosity =                                                                 &
               ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** powerLawExponent *         &
               reynoldsNumberInverse
       end if

       ! Second coefficient of viscosity
       if (present(secondCoefficientOfViscosity)) then
          if (present(dynamicViscosity)) then
             secondCoefficientOfViscosity =                                                  &
                  (bulkViscosityRatio - twoThirds) * dynamicViscosity
          else
             secondCoefficientOfViscosity = (bulkViscosityRatio - twoThirds) *               &
                  ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** powerLawExponent *      &
                  reynoldsNumberInverse
          end if
       end if

       ! Thermal diffusivity
       if (present(thermalDiffusivity)) then
          if (present(dynamicViscosity)) then
             thermalDiffusivity = dynamicViscosity * prandtlNumberInverse
          else
             thermalDiffusivity =                                                            &
                  ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** powerLawExponent *      &
                  reynoldsNumberInverse * prandtlNumberInverse
          end if
       end if

       ! Mass diffusivity
       if (present(massDiffusivity) .and. nSpecies .gt. 0) then
          if (present(dynamicViscosity)) then
             do i = 1, nSpecies
                massDiffusivity(:,i) = dynamicViscosity * schmidtNumberInverse(i)
             end do
          else
             do i = 1, nSpecies
                massDiffusivity(:,i) =                                                       &
                     ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** powerLawExponent *   &
                     reynoldsNumberInverse * schmidtNumberInverse(i)
             end do
          end if
       end if

    end if

    ! Turbulent viscosity
    if (useLES .and. present(turbulentViscosity)) then
       if (present(dynamicViscosity)) then
          dynamicViscosity = dynamicViscosity + turbulentViscosity
       end if
       if (present(thermalDiffusivity)) then
          thermalDiffusivity = thermalDiffusivity + turbulentViscosity * prandtlNumberInverse
       end if
       if (present(massDiffusivity) .and. nSpecies .gt. 0) then
          do i = 1, nSpecies
             massDiffusivity(:,i) = massDiffusivity(:,i) + turbulentViscosity *              &
                  schmidtNumberInverse(i)
          end do
       end if
    end if

    ! Artificial diffusion
    if (useShockCapturing) then
       if (present(dynamicViscosity) .and. present(artificialShearViscosity))                &
            dynamicViscosity = dynamicViscosity + artificialShearViscosity
       if (present(secondCoefficientOfViscosity) .and. present(artificialBulkViscosity))     &
            secondCoefficientOfViscosity = secondCoefficientOfViscosity +                    &
            artificialBulkViscosity
       if (present(thermalDiffusivity) .and. present(artificialThermalDiffusivity))          &
            thermalDiffusivity = thermalDiffusivity + artificialThermalDiffusivity
    end if

    return
  end subroutine compute_transport_variables

end module equation_of_state
