module equation_of_state

  ! External modules
  use precision
  use simulation_flags
  use solver_options

  implicit none

contains

  subroutine compute_specific_volume(density, specificVolume)

    implicit none

    ! Arguments
    real(WP), intent(in) :: density(:)
    real(WP), intent(out) :: specificVolume(:)

    specificVolume = 1.0_WP / density

    return
  end subroutine compute_specific_volume


  subroutine compute_velocity(specificVolume, momentum, velocity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: specificVolume(:), momentum(:,:)
    real(WP), intent(out) :: velocity(:,:)

    ! Local variables
    integer :: i

    if (size(momentum) .ne. size(velocity)) call die('compute_velocity: inconsistent sizes')

    do i = 1, size(momentum, 2)
       velocity(:,i) = specificVolume * momentum(:,i)
    end do

    return
  end subroutine compute_velocity


  subroutine compute_mass_fraction(specificVolume, massDensity, massFraction)

    implicit none

    ! Arguments
    real(WP), intent(in) :: specificVolume(:), massDensity(:,:)
    real(WP), intent(out) :: massFraction(:,:)

    ! Local Variables
    integer :: i

    if (size(massDensity) .ne. size(massFraction))                                           &
         call die('compute_mass_fraction: inconsistent sizes')

    do i = 1, size(massDensity, 2)
       massFraction(:,i) = specificVolume * massDensity(:, i)
    end do

    return
  end subroutine compute_mass_fraction


  subroutine compute_mixture_molecular_weight(massFraction, mixtureMolecularWeight)

    implicit none

    ! Arguments
    real(WP), intent(in) :: massFraction(:,:)
    real(WP), intent(out) :: mixtureMolecularWeight(:)

    ! Local variables
    integer :: i

    if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
       mixtureMolecularWeight = molecularWeightInverse(nSpecies+1)
       do i = 1, size(massFraction, 2)
          mixtureMolecularWeight = mixtureMolecularWeight + massFraction(:,i) *              &
               (molecularWeightInverse(i) - molecularWeightInverse(nSpecies+1))
       end do
       mixtureMolecularWeight = 1.0_WP / mixtureMolecularWeight
    end if

    return
  end subroutine compute_mixture_molecular_weight


  subroutine compute_pressure(density, energy, velocity, pressure, turbulentKineticEnergy)

    implicit none

    ! Arguments
    real(WP), intent(in) :: density(:), energy(:), velocity(:,:)
    real(WP), intent(out) :: pressure(:)
    real(WP), intent(in), optional :: turbulentKineticEnergy

    pressure = (ratioOfSpecificHeats - 1.0_WP) *                                             &
         (energy - 0.5_WP * density * sum(velocity ** 2, 2))

    if (present(turbulentKineticEnergy)) pressure = pressure -                               &
         (ratioOfSpecificHeats - 1.0_WP) * density * turbulentKineticEnergy

    return
  end subroutine compute_pressure


  subroutine compute_temperature(specificVolume, pressure, temperature,                      &
       mixtureMolecularWeight)

    implicit none

    ! Arguments
    real(WP), intent(in) :: specificVolume(:), pressure(:)
    real(WP), intent(out) :: temperature(:)
    real(WP), intent(in), optional :: mixtureMolecularWeight(:)

    select case (equationOfState)

    case (IDEAL_GAS)
       ! Compute temperature using ideal gas law
       temperature = ratioOfSpecificHeats * pressure /                                       &
            (ratioOfSpecificHeats - 1.0_WP) * specificVolume

    case (IDEAL_GAS_MIXTURE)
       ! Compute temperature using multi-component ideal gas
       temperature = ratioOfSpecificHeats * pressure * mixtureMolecularWeight /              &
            (ratioOfSpecificHeats - 1.0_WP) * specificVolume
    end select

    return
  end subroutine compute_temperature

  subroutine compute_viscosity(temperature, dynamicViscosity, secondCoefficientOfViscosity,  &
       turbulentViscosity, bulkViscosityModel)

    implicit none

    ! Arguments
    real(WP), intent(in) :: temperature(:)
    real(WP), intent(in), optional :: bulkViscosityModel(:), turbulentViscosity(:)
    real(WP), intent(out) :: dynamicViscosity(:)
    real(WP), intent(out), optional :: secondCoefficientOfViscosity(:)

    ! Local variables
    real(WP), parameter :: twoThirds = 2.0_WP / 3.0_WP

    select case (viscosityModel)
    case (POWER_LAW)
       ! Determine viscosity from a power law
       dynamicViscosity = reynoldsNumberInverse *                                            &
            ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** powerLawExponent
       
    case (SUTHERLANDS_LAW)
       ! Determine viscosity from Sutherland's law
       dynamicViscosity = reynoldsNumberInverse *                                            &
            ((ratioOfSpecificHeats - 1.0_WP) * temperature) ** 1.5_WP *                      &
            (1.0_WP + sutherlandConstant) / ((ratioOfSpecificHeats - 1.0_WP) * temperature + &
            sutherlandConstant)

    case default
       ! Constant viscosity
       dynamicViscosity = reynoldsNumberInverse

    end select

    ! Second coefficient of viscosity
    if (present(secondCoefficientOfViscosity)) then
       secondCoefficientOfViscosity = (bulkViscosityRatio - twoThirds) * dynamicViscosity
    end if

    ! Turbulent viscosity
    if (useLES .and. present(turbulentViscosity)) dynamicViscosity = dynamicViscosity +      &
         turbulentViscosity

    ! Augment the bulk viscosity
    if (useShockCapturing .or. useLES .and. present(bulkViscosityModel))                     &
         secondCoefficientOfViscosity = secondCoefficientOfViscosity + bulkViscosityModel

    return
  end subroutine compute_viscosity


  subroutine compute_thermal_diffusivity(dynamicViscosity, thermalDiffusivity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: dynamicViscosity(:)
    real(WP), intent(out) :: thermalDiffusivity(:)

    thermalDiffusivity = dynamicViscosity * prandtlNumberInverse

    return
  end subroutine compute_thermal_diffusivity


  subroutine compute_mass_diffusivity(dynamicViscosity, massDiffusivity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: dynamicViscosity(:)
    real(WP), intent(out) :: massDiffusivity(:,:)

    ! Local variables
    integer :: i

    do i = 1, size(massDiffusivity, 2)
       massDiffusivity(:,i) = dynamicViscosity * schmidtNumberInverse
    end do

    return
  end subroutine compute_mass_diffusivity

end module equation_of_state


! ========================================================= !
! Compute the dependent variables via the equation of state !
! stateVector(:,1)                       - Density          !
! stateVector(:,2:nDimensions+1)         - Momentum         !
! stateVector(:,nDimensions+2)           - Energy           !
! stateVector(:,nDimensions+3:nUnknowns) - Mass density     !
! ========================================================= !
subroutine compute_dependent_variables(stateVector)

  ! Internal modules
  use equation_of_state

  ! External modules
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns) :: stateVector

  ! Specific volume
  call compute_specific_volume(stateVector(:,1), specificVolume(:,1))

  ! Velocity
  call compute_velocity(specificVolume(:,1), stateVector(:,2:nDimensions+1),              &
       velocity(:,1:nDimensions))

  ! Mass fraction
  call compute_mass_fraction(specificVolume(:,1),                                         &
       stateVector(:, nDimensions+2+1:nUnknowns), massFraction(:,:))

  ! Mixture molecular weight
  call compute_mixture_molecular_weight(massFraction(:,:), mixtureMolecularWeight(:,1))

  ! Pressure
  call compute_pressure(stateVector(:,1), stateVector(:,nDimensions+2),                   &
       velocity(:,1:nDimensions), pressure(:,1))

  ! Temperature
  call compute_temperature(specificVolume(:,1), pressure(:,1), temperature(:,1),          &
       mixtureMolecularWeight(:,1))

  return
end subroutine compute_dependent_variables


! =============================== !
! Compute the transport variables !
! =============================== !
subroutine compute_transport_variables

  ! Internal modules
  use equation_of_state

  ! External modules
  use state

  implicit none

  ! Viscosity
  call compute_viscosity(temperature(:,1), dynamicViscosity(:,1),                            &
       secondCoefficientOfViscosity(:,1), turbulentViscosity(:,1), bulkViscosityModel(:,1))

  ! Thermal diffusivity
  call compute_thermal_diffusivity(dynamicViscosity(:,1), thermalDiffusivity(:,1))

  ! Mass diffusivity
  call compute_mass_diffusivity(dynamicViscosity(:,1), massDiffusivity(:,:))

  return
end subroutine compute_transport_variables
