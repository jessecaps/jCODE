module boivin

  ! External modules
  use combustion
  use precision
  use thermochem
  use solver_options
  use geometry

  implicit none

  ! ============================================================= !
  ! Boivin skeletal mechanism                                     !
  ! Reduced-Kinetic Mechanisms for Hydrogen and Syngas Combustion !
  ! Including Autoignition (2011 PhD thesis)                      !
  ! ============================================================= !

  ! Species indices (8 active species + N2)
  integer :: H2, O2, H2O, H, O, OH, HO2, H2O2, N2

  ! Rate coefficients (nReactions, forward/back)
  real(WP), dimension(12,2) :: prefactor, temperatureExponent, activationTemperature

  ! Stoichiometric coefficients (nReactions, nSpecies, forward/back)
  real(WP), dimension(12,8,2) :: stoichiometricCoefficient

  ! Three body reactions
  real(WP), dimension(12,9) :: collisionEfficiency !... include contributions from N2
  real(WP), dimension(12,4) :: troeFalloffCoefficients
  logical, dimension(12) :: threeBodyReaction, falloffReaction

  ! Enthalpy of formation
  real(WP), dimension(8) :: enthalpyOfFormation

  ! Reference quantities
  real(WP) :: referenceDensity, referenceTemperature, referencePressure,                     &
       referenceMolecularWeightInverse, referenceSoundSpeed, referenceLength

contains


  ! Compute the forward & backward reaction rates
  ! ---------------------------------------------
  subroutine boivin_reaction_rates(T, concentration, reactionRate)

    implicit none

    ! Arguments
    real(WP), intent(in) :: T, concentration(nSpecies+1)
    real(WP), intent(out) :: reactionRate(nReactions,2)

    ! Local variables
    integer :: j, k
    real(WP) :: thirdBodyFactor, falloffFunction, Fc, Pr, A, T1, T2, T3

    do j = 1, nReactions
       ! Forward reaction rate (or k_0 if falloff reaction)
       reactionRate(j,1) = prefactor(j,1) * T**temperatureExponent(j,1) *                    &
            exp(-activationTemperature(j,1) / T)

       ! Backward reaction rate (or k_inf if falloff reaction)
       reactionRate(j,2) = prefactor(j,2) * T**temperatureExponent(j,2) *                    &
            exp(-activationTemperature(j,2) / T)

       ! Three body reaction
       if (threeBodyReaction(j)) then
          thirdBodyFactor = 0.0_WP
          do k = 1, nSpecies + 1
             thirdBodyFactor = thirdBodyFactor + collisionEfficiency(j,k) * concentration(k)
          end do
          if (.not. falloffReaction(j)) then
             reactionRate(j,1) = reactionRate(j,1) * thirdBodyFactor
             reactionRate(j,2) = reactionRate(j,2) * thirdBodyFactor
          end if
       end if

       ! Falloff reaction
       if (falloffReaction(j)) then
          ! Reduced pressure
          Pr = thirdBodyFactor * reactionRate(j,1) / reactionRate(j,2)

          ! Get the falloff parameter
          A  = troeFalloffCoefficients(j,1)
          T1 = troeFalloffCoefficients(j,2)
          T2 = troeFalloffCoefficients(j,3)
          T3 = troeFalloffCoefficients(j,4)
          if (T1 + T2 + T3 .gt. 0.0_WP) then
             Fc = (1.0_WP - A) * exp(-T / T3) + A * exp(-T / T1) + exp(-T2/T)
          else
             Fc = A
          end if

          ! Get Troe Falloff function
          call troe_falloff(Fc, Pr, falloffFunction)

          ! Update the reaction rates
          reactionRate(j,1) = thirdBodyFactor * reactionRate(j,1) / (1.0_WP + Pr) *          &
               falloffFunction
          reactionRate(j,2) = 0.0_WP
       end if
    end do

    return
  end subroutine boivin_reaction_rates


  ! Get Troe falloff function
  ! -------------------------
  subroutine troe_falloff(Fc, Pr, F)

    implicit none

    ! Arguments
    real(WP), intent(in) :: Fc, Pr
    real(WP), intent(out) :: F

    ! Local variables
    real(WP) :: C, N, logFc, logPr, f1

    logFc = log10(Fc)
    logPr = log10(Pr + epsilon(1.0_WP))

    C = -0.4_WP - 0.67_WP * logFc
    N = 0.75_WP - 1.27_WP * logFc
    f1 = (logPr + C) / (N - 0.14_WP * (logPr + C))

    F = Fc ** (1.0_WP / (1.0_WP + f1**2))

    return
  end subroutine troe_falloff


  ! Compute the Jacobian matrix for one-step chemistry
  ! --------------------------------------------------
  subroutine boivin_jacobian(stateVector, jacobianOfSource, specificVolume,                  &
       velocity, temperature, massFraction, molecularWeightInverse)

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:), specificVolume, velocity(:), temperature,        &
         massFraction(:), molecularWeightInverse(:)
    real(WP), intent(out) :: jacobianOfSource(:,:)

    ! Local variables

    return
  end subroutine boivin_jacobian

end module boivin


! ============================ !
! Setup the one-step chemistry !
! ============================ !
subroutine boivin_setup

  ! Internal modules
  use boivin

  ! External modules
  use parser
  use string
  use simulation_flags
  use solver_options
  use state

  implicit none

  ! Local variables
  integer :: k
  real(WP), parameter :: kJmol_K = 1.0E3_WP / universalGasConstant !... from kJ/mol to K
  character(len = str_medium) :: referenceSpecies

  if (nSpecies .eq. 0) return

  ! Ensure 8 species are used for Boivin skeletal mechanism
  if (nSpecies .ne. 8)                                                                       &
       call die('boivin:setup: 8 species must be declared for Boivin skeletal mechanism')

  ! Check the equation of state
  if (equationOfState .ne. IDEAL_GAS_MIXTURE)                                                &
       call die('boivin_setup: Boivin skeletal mechanism requires ideal gas mixture')

  ! Set species indices
  H2 = 0; O2 = 0; H2O = 0; H = 0; O = 0; OH = 0; HO2 = 0; H2O2 = 0; N2 = 0
  do k = 1, nSpecies + 1
     select case (trim(speciesName(k)))
     case ('H2')
        H2 = k
     case ('O2')
        O2 = k
     case ('H2O')
        H2O = k
     case ('H')
        H = k
     case ('O')
        O = k
     case ('OH')
        OH = k
     case ('HO2')
        HO2 = k
     case ('H2O2')
        H2O2 = k
     case ('N2')
        N2 = k
     case default
        call die("boivin_setup: unknown species: '" // trim(speciesName(k)) // "'")
     end select
  end do

  ! Make sure all of the relevant species were defined
  if (H2   .eq. 0) call die('boivin_setup: species H2 not defined')
  if (O2   .eq. 0) call die('boivin_setup: species O2 not defined')
  if (H2O  .eq. 0) call die('boivin_setup: species H2O not defined')
  if (H    .eq. 0) call die('boivin_setup: species H not defined')
  if (O    .eq. 0) call die('boivin_setup: species O not defined')
  if (OH   .eq. 0) call die('boivin_setup: species OH not defined')
  if (HO2  .eq. 0) call die('boivin_setup: species HO2 not defined')
  if (H2O2 .eq. 0) call die('boivin_setup: species H2O2 not defined')
  if (N2   .eq. 0) call die('boivin_setup: species N2 not defined')

  ! Ensure N2 is the inert species
  if (N2 .ne. nSpecies+1) call die('boivin_setup: N2 must be the inert species')

  ! Read reference quantities from the input
  call parser_read('reference length', referenceLength, 1.0_WP)
  call parser_read('reference temperature', referenceTemperature, 298.0_WP)
  call parser_read('reference pressure', referencePressure, 101325.0_WP)

  ! Look up reference quantities
  call parser_read('reference species', referenceSpecies, 'air')
  call get_density(trim(referenceSpecies), referenceDensity)
  call get_molecular_weight(trim(referenceSpecies), referenceMolecularWeightInverse)
  referenceMolecularWeightInverse = 1.0_WP / referenceMolecularWeightInverse

  ! Get dimensional speed of sound [m/s]
  referenceSoundSpeed = sqrt(ratioOfSpecificHeats * referencePressure / referenceDensity)

  ! Get enthalpies of formation [J/kg]
  do k = 1, nSpecies
     call get_enthalpy_formation(trim(speciesName(k)), enthalpyOfFormation(k))
  end do

  ! 12 reactions (6 are reversible)
  nReactions = 12

  ! Initialize the rate coefficients
  prefactor = 0.0_WP
  temperatureExponent = 0.0_WP
  activationTemperature = 0.0_WP
  stoichiometricCoefficient = 0.0_WP
  collisionEfficiency = 0.0_WP
  troeFalloffCoefficients = 0.0_WP
  threeBodyReaction = .false.
  falloffReaction = .false.

  ! Reaction 1: H + O2 <=> OH + O
  stoichiometricCoefficient(1,H,1) = 1.0_WP
  stoichiometricCoefficient(1,O2,1) = 1.0_WP
  stoichiometricCoefficient(1,OH,2) = 1.0_WP
  stoichiometricCoefficient(1,O,2) = 1.0_WP
  prefactor(1,1) = 3.52E16_WP
  prefactor(1,2) = 7.04E13_WP
  temperatureExponent(1,1) = -0.7_WP
  temperatureExponent(1,2) = -0.26_WP
  activationTemperature(1,1) = 71.42_WP * kJmol_K
  activationTemperature(1,2) = 0.6_WP * kJmol_K

  ! Reaction 2: H2 + O <=> OH + H
  stoichiometricCoefficient(2,H2,1) = 1.0_WP
  stoichiometricCoefficient(2,O,1) = 1.0_WP
  stoichiometricCoefficient(2,OH,2) = 1.0_WP
  stoichiometricCoefficient(2,H,2) = 1.0_WP
  prefactor(2,1) = 5.06E4_WP
  prefactor(2,2) = 3.03E4_WP
  temperatureExponent(2,1) = 2.67_WP
  temperatureExponent(2,2) = 2.63_WP
  activationTemperature(2,1) = 26.32_WP * kJmol_K
  activationTemperature(2,2) = 20.23_WP * kJmol_K

  ! Reaction 3: H2 + OH <=> H2O + H
  stoichiometricCoefficient(3,H2,1) = 1.0_WP
  stoichiometricCoefficient(3,OH,1) = 1.0_WP
  stoichiometricCoefficient(3,H2O,2) = 1.0_WP
  stoichiometricCoefficient(3,H,2) = 1.0_WP
  prefactor(3,1) = 1.17E9_WP
  prefactor(3,2) = 1.28E10_WP
  temperatureExponent(3,1) = 1.3_WP
  temperatureExponent(3,2) = 1.19_WP
  activationTemperature(3,1) = 15.21_WP * kJmol_K
  activationTemperature(3,2) = 78.25_WP * kJmol_K

  ! Reaction 4: H + O2 + M => HO2 + M^b
  stoichiometricCoefficient(4,H,1) = 1.0_WP
  stoichiometricCoefficient(4,O2,1) = 1.0_WP
  stoichiometricCoefficient(4,HO2,2) = 1.0_WP
  prefactor(4,1) = 5.75E19_WP
  prefactor(4,2) = 4.65E12_WP
  temperatureExponent(4,1) = -1.4_WP
  temperatureExponent(4,2) = 0.44_WP
  activationTemperature(4,1) = 0.0_WP
  activationTemperature(4,2) = 0.0_WP
  threeBodyReaction(4) = .true.
  collisionEfficiency(4,:) = 1.0_WP
  collisionEfficiency(4,H2) = 2.5_WP
  collisionEfficiency(4,H2O) = 16.0_WP
  falloffReaction(4) = .true.
  troeFalloffCoefficients(4,1) = 0.5_WP

  ! Reaction 5: HO2 + H => 2OH
  stoichiometricCoefficient(5,HO2,1) = 1.0_WP
  stoichiometricCoefficient(5,H,1) = 1.0_WP
  stoichiometricCoefficient(5,OH,2) = 2.0_WP
  prefactor(5,1) = 7.08E13_WP
  activationTemperature(5,1) = 1.23_WP * kJmol_K

  ! Reaction 6: HO2 + H <=> H2 + O2
  stoichiometricCoefficient(6,HO2,1) = 1.0_WP
  stoichiometricCoefficient(6,H,1) = 1.0_WP
  stoichiometricCoefficient(6,H2,2) = 1.0_WP
  stoichiometricCoefficient(6,O2,2) = 1.0_WP
  prefactor(6,1) = 1.66E13_WP
  prefactor(6,2) = 2.69E12_WP
  temperatureExponent(6,1) = 0.0_WP
  temperatureExponent(6,2) = 0.36_WP
  activationTemperature(6,1) = 3.44_WP * kJmol_K
  activationTemperature(6,2) = 231.86_WP * kJmol_K

  ! Reaction 7: HO2 + H => H2O + O2
  stoichiometricCoefficient(7,HO2,1) = 1.0_WP
  stoichiometricCoefficient(7,OH,1) = 1.0_WP
  stoichiometricCoefficient(7,H2O,2) = 1.0_WP
  stoichiometricCoefficient(7,O2,2) = 1.0_WP
  prefactor(7,1) = 2.89E13_WP
  activationTemperature(7,1) = -2.08 * kJmol_K

  ! Reaction 8: H + OH + M <=> H2O + M^c
  stoichiometricCoefficient(8,H,1) = 1.0_WP
  stoichiometricCoefficient(8,OH,1) = 1.0_WP
  stoichiometricCoefficient(8,H2O,2) = 1.0_WP
  prefactor(8,1) = 4.00E22_WP
  prefactor(8,2) = 1.03E23_WP
  temperatureExponent(8,1) = -2.0_WP
  temperatureExponent(8,2) = -1.75_WP
  activationTemperature(8,1) = 0.0_WP
  activationTemperature(8,2) = 496.14 * kJmol_K
  threeBodyReaction(8) = .true.
  collisionEfficiency(8,:) = 1.0_WP
  collisionEfficiency(8,H2) = 2.5_WP
  collisionEfficiency(8,H2O) = 12.0_WP

  ! Reaction 9: 2H + M <=> H2 + M^c
  stoichiometricCoefficient(9,H,1) = 2.0_WP
  stoichiometricCoefficient(9,H2,2) = 1.0_WP
  prefactor(9,1) = 1.30E18_WP
  prefactor(9,2) = 3.04E17_WP
  temperatureExponent(9,1) = -1.0_WP
  temperatureExponent(9,2) = -0.65_WP
  activationTemperature(9,1) = 0.0_WP
  activationTemperature(9,2) = 433.09_WP * kJmol_K

  ! Reaction 10: 2HO2 => H2O2 + O2
  stoichiometricCoefficient(10,HO2,1) = 2.0_WP
  stoichiometricCoefficient(10,H2O2,2) = 1.0_WP
  stoichiometricCoefficient(10,O2,2) = 1.0_WP
  prefactor(10,1) = 3.02E12_WP
  temperatureExponent(10,1) = 0.0_WP
  activationTemperature(10,1) = 5.8 * kJmol_K

  ! Reaction 11: HO2 + H2 => H2O2 + H
  stoichiometricCoefficient(11,HO2,1) = 1.0_WP
  stoichiometricCoefficient(11,H2,1) = 1.0_WP
  stoichiometricCoefficient(11,H2O2,2) = 1.0_WP
  stoichiometricCoefficient(11,H,2) = 1.0_WP
  prefactor(11,1) = 1.62E11_WP
  temperatureExponent(11,1) = 0.61_WP
  activationTemperature(11,1) = 100.14 * kJmol_K

  ! Reaction 12: H2O2 + M => 2OH + M^d
  stoichiometricCoefficient(12,H2O2,1) = 1.0_WP
  stoichiometricCoefficient(12,OH,2) = 2.0_WP
  prefactor(12,1) = 8.15E23_WP
  prefactor(12,2) = 2.62E19_WP
  temperatureExponent(12,1) = -1.9_WP
  temperatureExponent(12,2) = -1.39_WP
  activationTemperature(12,1) = 207.62_WP * kJmol_K
  activationTemperature(12,2) = 214.74_WP * kJmol_K
  threeBodyReaction(12) = .true.
  collisionEfficiency(12,:) = 1.0_WP
  collisionEfficiency(12,H2) = 2.5_WP
  collisionEfficiency(12,H2O) = 6.0_WP
  falloffReaction(12) = .true.
  troeFalloffCoefficients(12,1) = 0.735_WP
  troeFalloffCoefficients(12,2) = 1756.0_WP
  troeFalloffCoefficients(12,3) = 5182.0_WP
  troeFalloffCoefficients(12,4) = 94.0_WP

  return
end subroutine boivin_setup


! ============================== !
! Cleanup the one-step chemistry !
! ============================== !
subroutine boivin_cleanup

  ! Internal modules
  use boivin

  implicit none

  return
end subroutine boivin_cleanup


! ========================================== !
! Add one-step source during the forward run !
! ========================================== !
subroutine boivin_forward(source)

  ! Internal modules
  use boivin

  ! External modules
  use solver_options
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k
  real(WP), parameter :: kgm_gcm = 1.0E-3_WP ! ... from kg/m^3 to g/cm^3
  real(WP), parameter :: gcm_kgm = 1.0E+3_WP ! ... from kg/m^3 to g/cm^3
  real(WP) :: T, chemicalSource, dh, Qf, Qb
  real(WP), dimension(nReactions) :: progressRate
  real(WP), dimension(nReactions,2) :: reactionRate
  real(WP), dimension(nSpecies+1) :: concentration, Wi

  if (nSpecies .eq. 0 .or. nReactions .eq. 0) return

  ! Dimensionalize molecularWeightInverse
  do k = 1, nSpecies + 1
     Wi(k) = molecularWeightInverse(k) * referenceMolecularWeightInverse
  end do

  do i = 1, nGridPoints

     ! Dimensionalize the local temperature [K]
     T = temperature(i,1) * (ratioOfSpecificHeats - 1.0_WP) * referenceTemperature

     ! Compute the concentrations (dimensional)
     concentration(N2) = conservedVariables(i,1)
     do k = 1, nSpecies
        concentration(k) = referenceDensity * kgm_gcm *                                      &
             conservedVariables(i,nDimensions+2+k) * Wi(k)
        concentration(N2) = concentration(N2) - conservedVariables(i,nDimensions+2+k)
     end do
     concentration(N2) = concentration(N2) * referenceDensity * kgm_gcm * Wi(N2)

     ! Clip the concentrations
     !do k = 1, nSpecies + 1
     !   concentration(k) = max(concentration(k), 1.0E-16_WP)
     !end do

     ! Get the forward & backward reaction rates
     call boivin_reaction_rates(T, concentration, reactionRate)

     ! Compute the progress rate (law of mass action)
     do j = 1, nReactions
        Qf = reactionRate(j,1)
        Qb = reactionRate(j,2)
        do k = 1, nSpecies
           Qf = Qf * concentration(k)**stoichiometricCoefficient(j,k,1)
           Qb = Qb * concentration(k)**stoichiometricCoefficient(j,k,2)
        end do
        progressRate(j) = Qf - Qb
     end do

     ! Compute and add the chemical source terms
     do k = 1, nSpecies
        chemicalSource = 0.0_WP
        do j = 1, nReactions
           chemicalSource = chemicalSource + progressRate(j) *                               &
                (stoichiometricCoefficient(j,k,2) - stoichiometricCoefficient(j,k,1))
        end do
        chemicalSource = chemicalSource / Wi(k)

        ! Non-dimensionalize
        chemicalSource = chemicalSource * gcm_kgm  * referenceLength /                       &
             (referenceSoundSpeed * referenceDensity)
        dh = enthalpyOfFormation(k) / referenceSoundSpeed**2

        ! Update the species source term
        source(i,nDimensions+2+k) = source(i,nDimensions+2+k) + chemicalSource

        ! Heat release due to combustion
        source(i,nDimensions+2) = source(i,nDimensions+2) - dh * chemicalSource

     end do

  end do

  return
end subroutine boivin_forward


! ========================================== !
! Add one-step source during the adjoint run !
! ========================================== !
subroutine boivin_adjoint(source)

  ! Internal modules
  use boivin

  ! External modules
  use solver_options
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k
  real(WP), allocatable :: localJacobian(:,:), localConservedVariables(:),                   &
       localVelocity(:), localMassFraction(:), temp(:)
  real(WP) :: localTemperature, localMixtureMolecularWeight

  if (nSpecies .eq. 0 .or. nReactions .eq. 0) return

  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMassFraction(nSpecies))
  allocate(localJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  do j = 1, nGridPoints

     localConservedVariables = conservedVariables(j,:)
     localTemperature = temperature(j,1)
     localVelocity = velocity(j,:)
     localMassFraction = massFraction(j,:)

     call boivin_jacobian(localConservedVariables, localJacobian, specificVolume(j,1),      &
          localVelocity, localTemperature, localMassFraction, molecularWeightInverse)

     temp = matmul(transpose(localJacobian), adjointVariables(j,:))

     localMixtureMolecularWeight = mixtureMolecularWeight(j,1)

     temp(nDimensions+2) = ratioOfSpecificHeats * specificVolume(j,1) *                   &
          localMixtureMolecularWeight * temp(nDimensions+2)
     do i = 1, nDimensions
        temp(i+1) = specificVolume(j,1) * temp(i+1) - localVelocity(i) *                  &
             temp(nDimensions+2)
     end do
     do k = 1, nSpecies
        temp(nDimensions+2+k) = specificVolume(j,1) *  temp(nDimensions+2+k) +            &
             localTemperature * (molecularWeightInverse(nSpecies+1) -                     &
             molecularWeightInverse(k)) * temp(nDimensions+2) / ratioOfSpecificHeats
     end do
     temp(1) = temp(1) - specificVolume(j,1) *localConservedVariables(nDimensions+2) *    &
          temp(nDimensions+2) - sum(localVelocity * temp(2:nDimensions+1)) -              &
          sum(localMassFraction * temp(ndimensions+3:nUnknowns))

     source(j,:) = source(j,:) - temp

  end do !... j = 1, nGridPoints

  deallocate(localConservedVariables)
  deallocate(localVelocity)
  deallocate(localMassFraction)
  deallocate(localJacobian)
  deallocate(temp)

  return
end subroutine boivin_adjoint
