module onestep

  ! External modules
  use combustion
  use precision
  use solver_options
  use geometry

  implicit none

  integer :: H2, O2, N2
  real(WP) :: YFs, Damkohler, stoichiometricRatio, heatRelease, zelDovich, controlParameter
  real(WP), dimension(:), allocatable :: Y0, stoichiometricCoefficient

contains

  ! Compute the Jacobian matrix for one-step chemistry
  ! --------------------------------------------------
  subroutine onestep_jacobian(stateVector, jacobianOfSource, specificVolume,                 &
       velocity, temperature, massFraction, molecularWeightInverse)

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:), specificVolume, velocity(:), temperature,        &
         massFraction(:), molecularWeightInverse(:)
    real(WP), intent(out) :: jacobianOfSource(:,:)

    ! Local variables
    integer :: k, l
    real(WP) :: referenceTemperature, flameTemperature, activationTemperature,               &
         chemicalSource, H, temp

    ! Zero-out source Jacobian
    jacobianOfSource = 0.0_WP

    ! Dependent variables
    referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
    flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
    activationTemperature = zelDovich / heatRelease * flameTemperature
    chemicalSource = controlParameter * Damkohler * stateVector(nDimensions+2+H2) *          &
         stateVector(nDimensions+2+O2) * exp(- activationTemperature / temperature)
    H = heatRelease * flameTemperature / Yfs

    ! Jacobian of combustion source
    temp = 2.0_WP * chemicalSource * specificVolume
    jacobianOfSource(nDimensions+2,1) = temp * H
    do k = 1, nSpecies
       jacobianOfSource(nDimensions+2+k,1) = - temp * stoichiometricCoefficient(k)
    end do

    temp = activationTemperature / temperature**2
    jacobianOfSource(nDimensions+2,nDimensions+2) = H * chemicalSource * temp
    do k = 1, nSpecies
       jacobianOfSource(nDimensions+2+k,nDimensions+2) =                                     &
            - stoichiometricCoefficient(k) * chemicalSource * temp
    end do

    do k = 1, nSpecies
       if (k .eq. H2) then
          temp = controlParameter * Damkohler * stateVector(1)**2 * massFraction(O2) *       &
               exp(- activationTemperature / temperature)
       else if (k .eq. O2) then
          temp = controlParameter * Damkohler * stateVector(1)**2 * massFraction(H2) *       &
               exp(- activationTemperature / temperature)
       end if
       jacobianOfSource(nDimensions+2,nDimensions+2+k) = H * temp
       do l = 1, nSpecies
          jacobianOfSource(nDimensions+2+l,nDimensions+2+k) =                                &
               - stoichiometricCoefficient(l) * temp
       end do
    end do

    return
  end subroutine onestep_jacobian

end module onestep


! ============================ !
! Setup the one-step chemistry !
! ============================ !
subroutine onestep_setup

  ! Internal modules
  use onestep

  ! External modules
  use parser
  use simulation_flags
  use solver_options
  use state

  implicit none

  ! Local variables
  integer :: k
  real(WP) :: equivalenceRatio

  if (nSpecies .eq. 0) return

  ! Get species indices
  if (allocated(speciesName)) then
     do k = 1, nSpecies + 1
        select case (trim(speciesName(k)))
        case ('H2', 'HYDROGEN')
           H2 = k
        case ('O2', 'OXYGEN')
           O2 = k
        case ('N2', 'NITROGEN')
           N2 = k
        case default
           call die("onestep_setup: unknown species: '" // trim(speciesName(k)))
        end select
     end do
  else
     H2 = 1
     O2 = 2
     N2 = nSpecies + 1
     allocate(speciesName(nSpecies+1))
     speciesName(H2) = 'H2'
     speciesName(O2) = 'O2'
     speciesName(N2) = 'N2'
  end if

  ! One-step chemistry considers a single irreversible reaction
  nReactions = 1

  allocate(Y0(nSpecies))

  ! Read combustion parameters from input
  call parser_read('stoichiometric ratio', stoichiometricRatio)
  call parser_read('heat release', heatRelease)
  call parser_read('Zel Dovich', zelDovich)
  call parser_read('Damkohler number', Damkohler)
  call parser_read('initial fuel mass fraction', Y0(H2))
  call parser_read('initial oxidizer mass fraction', Y0(O2))

  ! Control parameter to test adjoint
  if (predictionOnly) then
     controlParameter = 1.0_WP
  else
     call parser_read('chemical source parameter', controlParameter, 1.0_WP)
  end if

  ! Stoichiometric fuel mass fraction
  equivalenceRatio = stoichiometricRatio * Y0(H2) / Y0(O2)
  Yfs = Y0(H2) / (1.0_WP + equivalenceRatio)

  ! Stoichiometric coefficients
  allocate(stoichiometricCoefficient(nSpecies))
  stoichiometricCoefficient = 0.0_WP
  stoichiometricCoefficient(H2) = 1.0_WP
  stoichiometricCoefficient(O2) = stoichiometricRatio

  return
end subroutine onestep_setup


! ============================== !
! Cleanup the one-step chemistry !
! ============================== !
subroutine onestep_cleanup

  ! Internal modules
  use onestep

  implicit none

  if (allocated(Y0)) deallocate(Y0)
  if (allocated(stoichiometricCoefficient)) deallocate(stoichiometricCoefficient)

  return
end subroutine onestep_cleanup


! ========================================== !
! Add one-step source during the forward run !
! ========================================== !
subroutine onestep_forward(source)

  ! Internal modules
  use onestep

  ! External modules
  use solver_options
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, k
  real(WP) :: referenceTemperature, flameTemperature, activationTemperature,                 &
       chemicalSource, H

  if (nSpecies .eq. 0 .or. nReactions .eq. 0) return

  ! One-step chemistry
  ! H2 + sO2 => (1+s)P

  referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
  activationTemperature = zelDovich / heatRelease * flameTemperature
  H = heatRelease * flameTemperature / Yfs

  do i = 1, nGridPoints

     chemicalSource = controlParameter * Damkohler * conservedVariables(i,nDimensions+2+H2) *&
          conservedVariables(i,nDimensions+2+O2) *                                           &
          exp(- activationTemperature / temperature(i, 1))

     ! Heat release due to combustion
     source(i,nDimensions+2) = source(i,nDimensions+2) + H * chemicalSource

     ! Species source terms
     do k = 1, nSpecies
        source(i,nDimensions+2+k) = source(i,nDimensions+2+k) -                              &
             stoichiometricCoefficient(k) * chemicalSource
     end do

  end do

  return
end subroutine onestep_forward


! ========================================== !
! Add one-step source during the adjoint run !
! ========================================== !
subroutine onestep_adjoint(source)

  ! Internal modules
  use onestep

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

     call onestep_jacobian(localConservedVariables, localJacobian, specificVolume(j,1),      &
          localVelocity, localTemperature, localMassFraction, molecularWeightInverse)

     temp = matmul(transpose(localJacobian), adjointVariables(j,:))

     select case (equationOfState)

     case (IDEAL_GAS)

        temp(nDimensions+2) = ratioOfSpecificHeats * specificVolume(j,1) *                   &
             temp(nDimensions+2)
        do i = 1, nDimensions
           temp(i+1) = specificVolume(j,1) * temp(i+1) - localVelocity(i) *                  &
                temp(nDimensions+2)
        end do
        do k = 1, nSpecies
           temp(nDimensions+2+k) = specificVolume(j,1) *  temp(nDimensions+2+k)
        end do
        temp(1) = temp(1) - specificVolume(j,1) * localConservedVariables(nDimensions+2) *   &
             temp(nDimensions+2) - sum(localVelocity * temp(2:nDimensions+1))
        if (nSpecies .gt. 0) temp(1) = temp(1) - sum(localMassFraction *                     &
             temp(ndimensions+3:nUnknowns))

     case (IDEAL_GAS_MIXTURE)

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

     end select

     source(j,:) = source(j,:) - temp

  end do !... j = 1, nGridPoints

  deallocate(localConservedVariables)
  deallocate(localVelocity)
  deallocate(localMassFraction)
  deallocate(localJacobian)
  deallocate(temp)

  return
end subroutine onestep_adjoint
