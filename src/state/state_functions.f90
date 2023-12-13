module state_functions

  ! External modules
  use precision
  use simulation_flags
  use solver_options
  use geometry

  implicit none

contains

  ! Compute the Cartesian inviscid fluxes
  ! -------------------------------------
  subroutine compute_cartesian_inviscid_fluxes(stateVector, velocity, pressure,              &
       inviscidFluxes)

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:), velocity(:,:), pressure(:)
    real(WP), intent(out) :: inviscidFluxes(:,:,:)

    ! Local variables
    integer :: i

    select case (nDimensions)

    case (1)
       inviscidFluxes(:,1,1) = stateVector(:,2)
       inviscidFluxes(:,2,1) = stateVector(:,2) * velocity(:,1) + pressure
       inviscidFluxes(:,3,1) = velocity(:,1) * (stateVector(:,3) + pressure)
       do i = 1, nSpecies
          inviscidFluxes(:,i+3,1) = stateVector(:,i+3) * velocity(:,1)
       end do

    case (2)
       inviscidFluxes(:,1,1) = stateVector(:,2)
       inviscidFluxes(:,2,1) = stateVector(:,2) * velocity(:,1) + pressure
       inviscidFluxes(:,3,1) = stateVector(:,2) * velocity(:,2)
       inviscidFluxes(:,4,1) = velocity(:,1) * (stateVector(:,4) + pressure)
       inviscidFluxes(:,1,2) = stateVector(:,3)
       inviscidFluxes(:,2,2) = inviscidFluxes(:,3,1)
       inviscidFluxes(:,3,2) = stateVector(:,3) * velocity(:,2) + pressure
       inviscidFluxes(:,4,2) = velocity(:,2) * (stateVector(:,4) + pressure)
       do i = 1, nSpecies
          inviscidFluxes(:,i+4,1) = stateVector(:,i+4) * velocity(:,1)
          inviscidFluxes(:,i+4,2) = stateVector(:,i+4) * velocity(:,2)
       end do

    case (3)
       inviscidFluxes(:,1,1) = stateVector(:,2)
       inviscidFluxes(:,2,1) = stateVector(:,2) * velocity(:,1) + pressure
       inviscidFluxes(:,3,1) = stateVector(:,2) * velocity(:,2)
       inviscidFluxes(:,4,1) = stateVector(:,2) * velocity(:,3)
       inviscidFluxes(:,5,1) = velocity(:,1) * (stateVector(:,5) + pressure)
       inviscidFluxes(:,1,2) = stateVector(:,3)
       inviscidFluxes(:,2,2) = inviscidFluxes(:,3,1)
       inviscidFluxes(:,3,2) = stateVector(:,3) * velocity(:,2) + pressure
       inviscidFluxes(:,4,2) = stateVector(:,3) * velocity(:,3)
       inviscidFluxes(:,5,2) = velocity(:,2) * (stateVector(:,5) + pressure)
       inviscidFluxes(:,1,3) = stateVector(:,4)
       inviscidFluxes(:,2,3) = inviscidFluxes(:,4,1)
       inviscidFluxes(:,3,3) = inviscidFluxes(:,4,2)
       inviscidFluxes(:,4,3) = stateVector(:,4) * velocity(:,3) + pressure
       inviscidFluxes(:,5,3) = velocity(:,3) * (stateVector(:,5) + pressure)
       do i = 1, nSpecies
          inviscidFluxes(:,i+5,1) = stateVector(:,i+5) * velocity(:,1)
          inviscidFluxes(:,i+5,2) = stateVector(:,i+5) * velocity(:,2)
          inviscidFluxes(:,i+5,3) = stateVector(:,i+5) * velocity(:,3)
       end do

    end select

    return
  end subroutine compute_cartesian_inviscid_fluxes


  ! Compute the Cartesian viscous fluxes
  ! ------------------------------------
  subroutine compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,              &
       viscousFluxes, speciesFlux, enthalpyFlux)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
    real(WP), intent(in), optional :: enthalpyFlux(:,:), speciesFlux(:,:,:)
    real(WP), intent(out) :: viscousFluxes(:,:,:)

    ! Local variables
    integer :: k

    viscousFluxes = 0.0_WP

    select case (nDimensions)

    case (1)
       viscousFluxes(:,1,1) = 0.0_WP
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = velocity(:,1) * stressTensor(:,1) - heatFlux(:,1)
       if (present(enthalpyFlux) .and. equationOfState .eq. IDEAL_GAS_MIXTURE) then
          viscousFluxes(:,3,1) = viscousFluxes(:,3,1) - enthalpyFlux(:,1)
       end if
       if (present(speciesFlux)) then
          do k = 1, nSpecies
             viscousFluxes(:,k+3,1) = - speciesFlux(:,k,1)
          end do
       end if

    case (2)
       viscousFluxes(:,1,1) = 0.0_WP
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = stressTensor(:,3)
       viscousFluxes(:,4,1) = velocity(:,1) * stressTensor(:,1) +                            &
            velocity(:,2) * stressTensor(:,3) - heatFlux(:,1)
       viscousFluxes(:,1,2) = 0.0_WP
       viscousFluxes(:,2,2) = stressTensor(:,2)
       viscousFluxes(:,3,2) = stressTensor(:,4)
       viscousFluxes(:,4,2) = velocity(:,1) * stressTensor(:,2) +                            &
            velocity(:,2) * stressTensor(:,4) - heatFlux(:,2)
       if (present(enthalpyFlux) .and. equationOfState .eq. IDEAL_GAS_MIXTURE) then
          viscousFluxes(:,4,1) = viscousFluxes(:,4,1) - enthalpyFlux(:,1)
          viscousFluxes(:,4,2) = viscousFluxes(:,4,2) - enthalpyFlux(:,2)
       end if
       if (present(speciesFlux)) then
          do k = 1, nSpecies
             viscousFluxes(:,k+4,1) = - speciesFlux(:,k,1)
             viscousFluxes(:,k+4,2) = - speciesFlux(:,k,2)
          end do
       end if

    case (3)
       viscousFluxes(:,1,1) = 0.0_WP
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = stressTensor(:,4)
       viscousFluxes(:,4,1) = stressTensor(:,7)
       viscousFluxes(:,5,1) = velocity(:,1) * stressTensor(:,1) +                            &
            velocity(:,2) * stressTensor(:,4) +                                              &
            velocity(:,3) * stressTensor(:,7) - heatFlux(:,1)
       viscousFluxes(:,1,2) = 0.0_WP
       viscousFluxes(:,2,2) = stressTensor(:,2)
       viscousFluxes(:,3,2) = stressTensor(:,5)
       viscousFluxes(:,4,2) = stressTensor(:,8)
       viscousFluxes(:,5,2) = velocity(:,1) * stressTensor(:,2) +                            &
            velocity(:,2) * stressTensor(:,5) +                                              &
            velocity(:,3) * stressTensor(:,8) - heatFlux(:,2)
       viscousFluxes(:,1,3) = 0.0_WP
       viscousFluxes(:,2,3) = stressTensor(:,3)
       viscousFluxes(:,3,3) = stressTensor(:,6)
       viscousFluxes(:,4,3) = stressTensor(:,9)
       viscousFluxes(:,5,3) = velocity(:,1) * stressTensor(:,3) +                            &
            velocity(:,2) * stressTensor(:,6) +                                              &
            velocity(:,3) * stressTensor(:,9) - heatFlux(:,3)
       if (present(enthalpyFlux) .and. equationOfState .eq. IDEAL_GAS_MIXTURE) then
          viscousFluxes(:,5,1) = viscousFluxes(:,5,1) - enthalpyFlux(:,1)
          viscousFluxes(:,5,2) = viscousFluxes(:,5,2) - enthalpyFlux(:,2)
          viscousFluxes(:,5,3) = viscousFluxes(:,5,3) - enthalpyFlux(:,3)
       end if
       if (present(speciesFlux)) then
          do k = 1, nSpecies
             viscousFluxes(:,k+5,1) = - speciesFlux(:,k,1)
             viscousFluxes(:,k+5,2) = - speciesFlux(:,k,2)
             viscousFluxes(:,k+5,3) = - speciesFlux(:,k,3)
          end do
       end if

    end select

    return
  end subroutine compute_cartesian_viscous_fluxes

  
  ! Compute the Cartesian turbulent fluxes
  ! --------------------------------------
  subroutine compute_cartesian_turbulent_fluxes(velocity, reynoldsStress, turbulentFluxes)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocity(:,:), reynoldsStress(:,:)
    real(WP), intent(out) :: turbulentFluxes(:,:,:)

    turbulentFluxes = 0.0_WP

    select case (nDimensions)

    case (1)
       
       turbulentFluxes(:,2,1) = reynoldsStress(:,1)
       turbulentFluxes(:,3,1) = velocity(:,1) * reynoldsStress(:,1)

    case (2)
       
       turbulentFluxes(:,2,1) = reynoldsStress(:,1)
       turbulentFluxes(:,3,1) = reynoldsStress(:,3)
       turbulentFluxes(:,4,1) = velocity(:,1) * reynoldsStress(:,1) +                        &
            velocity(:,2) * reynoldsStress(:,3)
       
       turbulentFluxes(:,2,2) = reynoldsStress(:,2)
       turbulentFluxes(:,3,2) = reynoldsStress(:,4)
       turbulentFluxes(:,4,2) = velocity(:,1) * reynoldsStress(:,2) +                        &
            velocity(:,2) * reynoldsStress(:,4)

    case (3)

       turbulentFluxes(:,2,1) = reynoldsStress(:,1)
       turbulentFluxes(:,3,1) = reynoldsStress(:,4)
       turbulentFluxes(:,4,1) = reynoldsStress(:,7)
       turbulentFluxes(:,5,1) = velocity(:,1) * reynoldsStress(:,1) +                        &
            velocity(:,2) * reynoldsStress(:,4) +                                            &
            velocity(:,3) * reynoldsStress(:,7)

       turbulentFluxes(:,2,2) = reynoldsStress(:,2)
       turbulentFluxes(:,3,2) = reynoldsStress(:,5)
       turbulentFluxes(:,4,2) = reynoldsStress(:,8)
       turbulentFluxes(:,5,2) = velocity(:,1) * reynoldsStress(:,2) +                        &
            velocity(:,2) * reynoldsStress(:,5) +                                            &
            velocity(:,3) * reynoldsStress(:,8)
       
       turbulentFluxes(:,2,3) = reynoldsStress(:,3)
       turbulentFluxes(:,3,3) = reynoldsStress(:,6)
       turbulentFluxes(:,4,3) = reynoldsStress(:,9)
       turbulentFluxes(:,5,3) = velocity(:,1) * reynoldsStress(:,3) +                        &
            velocity(:,2) * reynoldsStress(:,6) +                                            &
            velocity(:,3) * reynoldsStress(:,9)

    end select

    return
  end subroutine compute_cartesian_turbulent_fluxes


  ! Compute the stress tensor
  ! -------------------------
  subroutine compute_stress_tensor(velocityGradient, dynamicViscosity,                       &
       secondCoefficientOfViscosity, stressTensor)

    implicit none

    ! Arguments
    real(WP), intent(inout) :: velocityGradient(:,:)
    real(WP), intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)
    real(WP), intent(out), optional :: stressTensor(:,:)

    ! Local variables
    real(WP), allocatable :: divergenceOfVelocity(:)

    if (present(stressTensor)) then !... out-of-place computation

       select case (nDimensions)

       case (1)
          stressTensor(:,1) = (2.0_WP * dynamicViscosity +                                   &
               secondCoefficientOfViscosity) * velocityGradient(:,1)

       case (2)
          stressTensor(:,1) = 2.0_WP * dynamicViscosity * velocityGradient(:,1) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,4))
          stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) +                    &
               velocityGradient(:,3))
          stressTensor(:,3) = stressTensor(:,2)
          stressTensor(:,4) = 2.0_WP * dynamicViscosity * velocityGradient(:,4) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,4))

       case (3)
          stressTensor(:,1) = 2.0_WP * dynamicViscosity * velocityGradient(:,1) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))
          stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) +                    &
               velocityGradient(:,4))
          stressTensor(:,3) = dynamicViscosity * (velocityGradient(:,3) +                    &
               velocityGradient(:,7))
          stressTensor(:,4) = stressTensor(:,2)
          stressTensor(:,5) = 2.0_WP * dynamicViscosity * velocityGradient(:,5) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))
          stressTensor(:,6) = dynamicViscosity * (velocityGradient(:,6) +                    &
               velocityGradient(:,8))
          stressTensor(:,7) = stressTensor(:,3)
          stressTensor(:,8) = stressTensor(:,6)
          stressTensor(:,9) = 2.0_WP * dynamicViscosity * velocityGradient(:,9) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))

       end select

    else !... in-place computation

       select case (nDimensions)

       case (1)
          velocityGradient(:,1) = (2.0_WP * dynamicViscosity +                               &
               secondCoefficientOfViscosity) * velocityGradient(:,1)

       case (2)
          allocate(divergenceOfVelocity(size(velocityGradient, 1)))
          divergenceOfVelocity = secondCoefficientOfViscosity *                              &
               (velocityGradient(:,1) + velocityGradient(:,4))
          velocityGradient(:,1) = 2.0_WP * dynamicViscosity * velocityGradient(:,1) +        &
               divergenceOfVelocity
          velocityGradient(:,2) = dynamicViscosity *                                         &
               (velocityGradient(:,2) + velocityGradient(:,3))
          velocityGradient(:,3) = velocityGradient(:,2)
          velocityGradient(:,4) = 2.0_WP * dynamicViscosity * velocityGradient(:,4) +        &
               divergenceOfVelocity

       case (3)
          allocate(divergenceOfVelocity(size(velocityGradient, 1)))
          divergenceOfVelocity = secondCoefficientOfViscosity *                              &
               (velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9))
          velocityGradient(:,1) = 2.0_WP * dynamicViscosity * velocityGradient(:,1) +        &
               divergenceOfVelocity
          velocityGradient(:,2) = dynamicViscosity *                                         &
               (velocityGradient(:,2) + velocityGradient(:,4))
          velocityGradient(:,3) = dynamicViscosity *                                         &
               (velocityGradient(:,3) + velocityGradient(:,7))
          velocityGradient(:,4) = velocityGradient(:,2)
          velocityGradient(:,5) = 2.0_WP * dynamicViscosity * velocityGradient(:,5) +        &
               divergenceOfVelocity
          velocityGradient(:,6) = dynamicViscosity *                                         &
               (velocityGradient(:,6) + velocityGradient(:,8))
          velocityGradient(:,7) = velocityGradient(:,3)
          velocityGradient(:,8) = velocityGradient(:,6)
          velocityGradient(:,9) = 2.0_WP * dynamicViscosity * velocityGradient(:,9) +        &
               divergenceOfVelocity

       end select

    end if

    return
  end subroutine compute_stress_tensor


  ! Inviscid flux vector splitting for upwinding
  ! ---------------------------------------------
  subroutine compute_upwind_inviscid_source(stateVector, specificVolume, pressure,            &
       velocity, volumeFraction, source, coeffSplit)

    ! External modules
    use grid
    use first_derivative
    use grid_functions
    use state_jacobian

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:), specificVolume(:), pressure(:), velocity(:,:), &
         volumeFraction(:)
    real(WP), intent(inout) :: source(:,:)
    real(WP), intent(in), optional :: coeffSplit

    ! Local variables
    integer :: i, j, k
    real(WP), allocatable :: temp_left(:,:), temp_right(:,:),                                &
         leftJacobianOfInviscidFlux(:,:), rightJacobianOfInviscidFlux(:,:),                  &
         localMassFraction(:), localMetricsAlongDirection(:)
    real(WP) :: coeffSplit_(nUnknowns)

    ! Not yet implemented for adjoints
    if (.not.predictionOnly .or. .not.useUpwinding) return

    if (present(coeffSplit)) then
       coeffSplit_(1) = 1.0_WP
       coeffSplit_(2:nUnknowns) = coeffSplit
    else
       coeffSplit_ = 1.0_WP
    end if

    allocate(temp_left(nGridPoints, nUnknowns))
    allocate(temp_right(nGridPoints, nUnknowns))
    allocate(leftJacobianOfInviscidFlux(nUnknowns, nUnknowns))
    allocate(rightJacobianOfInviscidFlux(nUnknowns, nUnknowns))
    if (nSpecies .gt. 0) allocate(localMassFraction(nSpecies))
    allocate(localMetricsAlongDirection(nDimensions))

    do j = 1, nDimensions

       ! Compute upwind derivative in direction `j'
       if (twoWayCoupling) then
          do i = 1, nUnknowns
             temp_left(:,i) = stateVector(:,i) * volumeFraction
             temp_right(:,i) = temp_left(:,i)
          end do
       else
          temp_left(:,:) = stateVector
          temp_right(:,:) = stateVector
       end if
       call upwind_apply(j, temp_left(:,:), temp_right(:,:))

       do i = 1, nGridPoints

          !if (sensor_(i) .lt. epsilon(0.0_WP)) cycle

          do k = 1, nSpecies
             localMassFraction(k) = stateVector(i, nDimensions+2+k) * specificVolume(i)
          end do

          localMetricsAlongDirection = metrics(i,1+nDimensions*(j-1):nDimensions*j)

          call compute_Steger_Warming_jacobian_of_inviscid_flux(stateVector(i,:),            &
               localMetricsAlongDirection, leftJacobianOfInviscidFlux,                       &
               rightJacobianOfInviscidFlux, specificVolume(i), velocity(i,:),                &
               pressure(i), localMassFraction)

          ! Update source
          source(i,:) = source(i,:) -                                                        &
               ( matmul(leftJacobianOfInviscidFlux, temp_left(i,:)) +                        &
               matmul(rightJacobianOfInviscidFlux, temp_right(i,:)) ) * coeffSplit_

       end do
    end do

    ! Cleanup
    deallocate(temp_left, temp_right, localMetricsAlongDirection,                            &
         leftJacobianOfInviscidFlux, rightJacobianOfInviscidFlux)
    if (allocated(localMassFraction)) deallocate(localMassFraction)

    return
  end subroutine compute_upwind_inviscid_source


  ! Split-convective energy conserving scheme by S. Pirozzoli (2011)
  ! `Stabilized non-dissipative approximations of Euler equations in generalized
  ! curvilinear coordinates,' JCP. Apply 'FE-Split' form with total enthalpy splitting.
  ! Relies on SBP finite difference discretization when the flow is non-periodic.
  ! -----------------------------------------------------------------------------------
  subroutine compute_split_inviscid_source(stateVector, specificVolume, pressure,            &
       velocity, volumeFraction, inviscidFlux, source)

    ! External modules
    use grid
    use first_derivative

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:), specificVolume(:), pressure(:), velocity(:,:), &
         volumeFraction(:)
    real(WP), intent(inout) :: inviscidFlux(:,:,:), source(:,:)

    ! Local variables
    integer :: i, j, k
    real(WP), dimension(nGridPoints, nDimensions) :: contravariantVelocity
    real(WP), dimension(nGridPoints, nUnknowns-1) :: transportVariables
    real(WP), dimension(nGridPoints) :: alpha
    real(WP), allocatable :: temp(:,:)

    ! Not yet implemented for adjoints
    if (.not. predictionOnly) return

    ! Store the contravariant velocity
    select case (nDimensions)
    case (1)
       contravariantVelocity(:,1) = velocity(:,1) * metrics(:,1)
    case (2)
       if (isDomainCurvilinear) then
          contravariantVelocity(:,1) = sum(velocity * metrics(:,1:2), dim = 2)
          contravariantVelocity(:,2) = sum(velocity * metrics(:,3:4), dim = 2)
       else
          contravariantVelocity(:,1) = velocity(:,1) * metrics(:,1)
          contravariantVelocity(:,2) = velocity(:,2) * metrics(:,4)
       end if
    case (3)
       if (isDomainCurvilinear) then
          contravariantVelocity(:,1) = sum(velocity * metrics(:,1:3), dim = 2)
          contravariantVelocity(:,2) = sum(velocity * metrics(:,4:6), dim = 2)
          contravariantVelocity(:,3) = sum(velocity * metrics(:,7:9), dim = 2)
       else
          contravariantVelocity(:,1) = velocity(:,1) * metrics(:,1)
          contravariantVelocity(:,2) = velocity(:,2) * metrics(:,5)
          contravariantVelocity(:,3) = velocity(:,3) * metrics(:,9)
       end if
    end select

    ! Store the transport variables appearing in momentum, energy, and species
    transportVariables(:,1:nDimensions) = velocity(:,1:nDimensions)
    transportVariables(:,nDimensions+1) = (stateVector(:,nDimensions+2) + pressure) *        &
         specificVolume
    do i = 1, nSpecies
       transportVariables(:,nDimensions+1+i) = stateVector(:,nDimensions+2+i) *              &
            specificVolume
    end do

    ! Need to differentiate additional metrics if curvilinear
    if (isDomainCurvilinear) then
       allocate(temp(nGridPoints, nUnknowns+1+nDimensions))
    else
       allocate(temp(nGridPoints, nUnknowns+2))
    end if

    ! Account for volume fraction if two-way coupling is used
    if (twoWayCoupling) then
       alpha = volumeFraction
    else
       alpha = 1.0_WP
    end if

    ! Loop through each direction
    do j = 1, nDimensions

       ! Pack the array
       temp(:,1) = stateVector(:,1) * contravariantVelocity(:,j) * alpha
       temp(:,2:nUnknowns) = transportVariables
       temp(:,nUnknowns+1) = pressure * alpha
       if (isDomainCurvilinear) then
          do k = 1, nDimensions
             temp(:,nUnknowns+1+k) = metrics(:,k+nDimensions*(j-1))
          end do
       else
          temp(:,nUnknowns+2) = metrics(:,j+nDimensions*(j-1))
       end if

       ! Apply derivtive in direction `j'
       call first_derivative_apply(j, temp)       

       ! Update the original flux and add new source terms (Jacobian will be applied later)
       do i = 1, nUnknowns - 1
          inviscidFlux(:,i+1,j) = 0.5_WP * inviscidFlux(:,i+1,j)
          source(:,i+1) = source(:,i+1) - 0.5_WP * (transportVariables(:,i) * temp(:,1) +    &
               stateVector(:,1) * contravariantVelocity(:,j) * alpha * temp(:,i+1))
       end do
       
       ! Handle pressure terms separately
       if (isDomainCurvilinear) then
          do k = 1, nDimensions
             source(:,k+1) = source(:,k+1) - 0.5_WP * (                                      &
                  metrics(:,k+nDimensions*(j-1)) * temp(:,nUnknowns+1) +                     &
                  pressure * alpha * temp(:,nUnknowns+1+k))
          end do
       else
          source(:,j+1) = source(:,j+1) - 0.5_WP * (                                         &
               metrics(:,j+nDimensions*(j-1)) * temp(:,nUnknowns+1) +                        &
               pressure * alpha * temp(:,nUnknowns+2))
       end if
    end do

    ! Cleanup
    deallocate(temp)
    
    return
  end subroutine compute_split_inviscid_source


  ! Non-conservative (split) Laplacian-form of the viscous fluxes
  ! -------------------------------------------------------------
  subroutine compute_split_viscous_source(velocity, velocityGradient, stressTensor,          &
       dynamicViscosity, secondCoefficientOfViscosity, temperature, thermalDiffusivity,      &
       heatFlux, massFraction, speciesFlux, massDiffusivity, volumeFraction, source)

    ! External modules
    use grid
    use grid_functions
    use first_derivative
    use second_derivative

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocity(:,:), velocityGradient(:,:), stressTensor(:,:),         &
         dynamicViscosity(:), secondCoefficientOfViscosity(:), temperature(:),               &
         thermalDiffusivity(:), heatFlux(:,:), massFraction(:,:), speciesFlux(:,:,:),        &
         massDiffusivity(:,:), volumeFraction(:)
    real(WP), intent(inout) :: source(:,:)

    ! Local variables
    integer :: i, ii, j, jj
    real(WP), dimension(nGridPoints) :: lapScalar
    real(WP), dimension(nGridPoints, 2) :: temp
    real(WP), dimension(nGridPoints, nDimensions) :: gradScalar
    real(WP), dimension(nGridPoints, nUnknowns) :: source_
    real(WP), allocatable :: dilatation(:), strainRate(:,:)
    
    ! Not yet implemented for adjoints
    if (.not.predictionOnly .or. .not.useViscosity) return

    if (isDomainCurvilinear)                                                                 &
         call die('split viscosity not yet implemented in curvilinear coordinates')

    ! Zero-out temporary source
    source_ = 0.0_WP

    ! Compute gradient of dilatation (avoid repeated derivatives)
    ! gradScalar = d^2(u_j)/d(x_i)/d(x_j)
    do i = 1, nDimensions
       ! Diagonal terms: d^2(u_j)/d(x_i)/d(x_j) (i==j)
       ii = i + (i - 1) * nDimensions
       temp(:,1) = velocity(:,i)
       call second_derivative_apply(i, temp(:,1:1))
       gradScalar(:,i) = jacobian(:,1)**2 * metrics(:,ii)**2 * temp(:,1)
       temp(:,1) = velocity(:,i)
       temp(:,2) = jacobian(:,1) * metrics(:,ii)
       call first_derivative_apply(i, temp)
       gradScalar(:,i) = gradScalar(:,i) + jacobian(:,1) * metrics(:,ii) *                   &
            temp(:,1) * temp(:,2)

       ! Off-diagonal terms: d^2(u_j)/d(x_i)/d(x_j) (i/=j)
       temp(:,1) = 0.0_WP
       do j = 1, nDimensions
          if (i .ne. j) then
             jj = j + (j - 1) * nDimensions
             temp(:,1) = temp(:,1) + velocityGradient(:,jj)
          end if
       end do
       call first_derivative_apply(i, temp(:,1:1))
       gradScalar(:,i) = gradScalar(:,i) + jacobian(:,1) * metrics(:,ii) * temp(:,1)
      
       ! Update the source term with velocity Laplacian and dilatation gradient terms
       call laplacian(velocity(:,i), lapScalar)
       source_(:,i+1) = source_(:,i+1) + dynamicViscosity * (lapScalar + gradScalar(:,i)) +  &
            secondCoefficientOfViscosity * gradScalar(:,i)
       source_(:,nDimensions+2) = source_(:,nDimensions+2) +                                 &
            velocity(:,i) * (dynamicViscosity * (lapScalar + gradScalar(:,i)) +              &
            secondCoefficientOfViscosity * gradScalar(:,i))
    end do

    ! Add temperature Laplacian and tau:grad(u) to energy
    call laplacian(temperature, lapScalar)
    source_(:,nDimensions+2) = source_(:,nDimensions+2) + thermalDiffusivity * lapScalar +   &
         sum(stressTensor * velocityGradient, dim = 2)

    ! Mass diffusion
    do i = 1, nSpecies
       call laplacian(massFraction(:,i), lapScalar)
       source_(:,nDimensions+2+i) = source_(:,nDimensions+2+i) + massDiffusivity(:,i) *      &
            lapScalar
    end do

    ! Include gradient of diffusion coefficients if viscosity varies
    if (viscosityModel.ne.0 .or. useLES .or. useShockCapturing) then
       ! Compute dilatation and strain rate
       allocate(dilatation(nGridPoints))
       allocate(strainRate(nGridPoints, nDimensions**2))
       call compute_dilatation(velocityGradient, dilatation)
       call compute_strain_rate(velocityGradient, strainRate)

       ! Account for gradient in viscosity
       call gradient(dynamicViscosity, gradScalar)
       do i = 1, nDimensions
          source_(:,i+1) = source_(:,i+1) + 2.0_WP * sum(gradScalar *                        &
               strainRate(:,1+nDimensions*(i-1):nDimensions*i), dim = 2)
          source_(:,nDimensions+2) = source_(:,nDimensions+2) + velocity(:,i) * (2.0_WP *    &
               sum(gradScalar * strainRate(:,1+nDimensions*(i-1):nDimensions*i), dim = 2))
       end do

       ! Account for gradient in second coefficientof viscosity
       call gradient(secondCoefficientOfViscosity, gradScalar)
       do i = 1, nDimensions
          source_(:,i+1) = source_(:,i+1) + dilatation * gradScalar(:,i)
          source_(:,nDimensions+2) = source_(:,nDimensions+2) + velocity(:,i) *              &
               dilatation * gradScalar(:,i)
       end do

       ! Account for gradient in thermal diffusivity
       ! Avoid computing gradT by setting gradT = heatFlux / thermalDiffusivity
       call gradient(thermalDiffusivity, gradScalar)
       source_(:,nDimensions+2) = source_(:,nDimensions+2) - sum(gradScalar *                &
            heatFlux, dim = 2) / (thermalDiffusivity + epsilon(1.0_WP))

       ! Account for gradient in mass diffusivity
       ! Avoid computing gradY by setting gradY = speciesFlux / massDiffusivity
       do i = 1, nSpecies
          call gradient(massDiffusivity(:,i), gradScalar)
          source_(:,nDimensions+2+i) = source_(:,nDimensions+2+i) - sum(gradScalar *         &
               speciesFlux(:,i,:), dim = 2) / (massdiffusivity(:,i) + epsilon(1.0_WP))
       end do

       ! Clean up
       deallocate(dilatation, strainRate)
    end if

    ! Update the source
    if (.not. twoWayCoupling) then
       source = source + source_
    else
       ! Account for volume fraction effects
       do i = 1, nUnknowns
          source(:,i) = source(:,i) + source_(:,i) * volumeFraction
       end do
       
       call gradient(volumeFraction, gradScalar)
       
       do i = 1, nDimensions
          source(:,i+1) = source(:,i+1) + sum(gradScalar *                                   &
               stressTensor(:,1+nDimensions*(i-1):nDimensions*i), dim = 2)
          source(:,nDimensions+2) = source(:,nDimensions+2) + sum(gradScalar *               &
               stressTensor(:,1+nDimensions*(i-1):nDimensions*i), dim = 2) * velocity(:,i)
       end do
       
       source(:,nDimensions+2) = source(:,nDimensions+2) - sum(gradScalar * heatFlux, dim = 2)

       do i = 1, nSpecies
          source(:,nDimensions+2+i) = source(:,nDimensions+2+i) -                            &
               sum(gradScalar * speciesFlux(:,i,:), dim = 2)
       end do
    end if

    return
  end subroutine compute_split_viscous_source


  ! Compute bulk viscosity using model originally proposed by Cook & Cabot (2004)
  ! 'LAD-D2-0' model from: Kawai et al., Assessment of localized artificial diffusivity
  ! scheme for large-eddy simulation of compressible turbulent flows (2010), JCP
  ! Using curvilinear treatment described in: Kawai et al., Localized artificial
  ! diffusivity scheme for discontinuity capturing on curvilinear meshes (2008), JCP
  ! ---------------------------------------------------------------------------------
  subroutine local_artificial_diffusivity(density, pressure, temperature, velocityGradient,  &
       shearViscosity, bulkViscosity, thermalDiffusivity)

    ! External modules
    use math
    use grid
    use grid_functions
    use filter
    use fourth_derivative
    use parallel

    implicit none

    ! Arguments
    real(WP), intent(in) :: density(:), pressure(:), temperature(:), velocityGradient(:,:)
    real(WP), intent(out) :: shearViscosity(:), bulkViscosity(:), thermalDiffusivity(:)

    ! Local variables
    integer :: i
    real(WP), parameter :: eps = 1.0e-32_WP
    real(WP) :: heavyside, sensor, div2, vort2, cBeta, cMu, cKappa, soundSpeed
    real(WP), allocatable :: F(:), vorticity(:,:), grad(:,:), temp(:,:), disp(:,:,:),        &
         delta2(:,:)

    ! Allocate temporary arrays
    allocate(temp(nGridPoints,1))
    allocate(F(nGridPoints))
    allocate(grad(nGridPoints, nDimensions))
    allocate(disp(nGridPoints, nDimensions, nDimensions))
    allocate(delta2(nGridPoints, nDimensions))
    
    ! Constants (scale based on user-specified input)
    cBeta = shockCoefficient
    cMu = 0.0_WP!shockCoefficient * 0.002_WP
    cKappa = shockCoefficient * 0.01_WP

    ! Compute the grid displacement
    select case (nDimensions)
    case (2)
       disp(:,1,1) =  metrics(:,4)
       disp(:,1,2) = -metrics(:,2)
       disp(:,2,1) = -metrics(:,3)
       disp(:,2,2) =  metrics(:,1)
    case (3)
       disp(:,1,1) = jacobian(:,1) * (metrics(:,5)*metrics(:,9) - metrics(:,6)*metrics(:,8))
       disp(:,1,2) = jacobian(:,1) * (metrics(:,3)*metrics(:,8) - metrics(:,2)*metrics(:,9))
       disp(:,1,3) = jacobian(:,1) * (metrics(:,2)*metrics(:,6) - metrics(:,3)*metrics(:,5))
       disp(:,2,1) = jacobian(:,1) * (metrics(:,6)*metrics(:,7) - metrics(:,4)*metrics(:,9))
       disp(:,2,2) = jacobian(:,1) * (metrics(:,1)*metrics(:,9) - metrics(:,3)*metrics(:,7))
       disp(:,2,3) = jacobian(:,1) * (metrics(:,3)*metrics(:,4) - metrics(:,1)*metrics(:,6))
       disp(:,3,1) = jacobian(:,1) * (metrics(:,4)*metrics(:,8) - metrics(:,5)*metrics(:,7))
       disp(:,3,2) = jacobian(:,1) * (metrics(:,2)*metrics(:,7) - metrics(:,1)*metrics(:,8))
       disp(:,3,3) = jacobian(:,1) * (metrics(:,1)*metrics(:,5) - metrics(:,2)*metrics(:,4))
    end select

    ! Compute artificial shear viscosity
    ! ----------------------------------
    shearViscosity = 0.0_WP
    if (cMu .gt. 0.0_WP) then
       
       ! Compute magnitude of the strain rate tensor
       select case (nDimensions)
       case (1)
          F = sqrt(velocityGradient(:,1)**2)
       case (2)
          F = sqrt(velocityGradient(:,1)**2 + velocityGradient(:,4)**2 +                     &
               0.5_WP * (velocityGradient(:,2) + velocityGradient(:,3))**2 )
       case (3)
          F = sqrt(velocityGradient(:,1)**2 + velocityGradient(:,5)**2 +                     &
               velocityGradient(:,9)**2 +                                                    &
               0.5_WP * (velocityGradient(:,2) + velocityGradient(:,4))**2 +                 &
               0.5_WP * (velocityGradient(:,3) + velocityGradient(:,7))**2 +                 &
               0.5_WP * (velocityGradient(:,6) + velocityGradient(:,8))**2 )
       end select      

       ! Resolution length scale
       select case (nDimensions)
       case(1)
          delta2(:,1) = disp(:,1,1)**2
       case(2)
          delta2(:,1) = disp(:,1,1)**2 + disp(:,2,1)**2
          delta2(:,2) = disp(:,1,2)**2 + disp(:,2,2)**2
       case(3)
          delta2(:,1) = disp(:,1,1)**2 + disp(:,2,1)**2 + disp(:,3,1)**2
          delta2(:,2) = disp(:,1,2)**2 + disp(:,2,2)**2 + disp(:,3,2)**2
          delta2(:,3) = disp(:,1,3)**2 + disp(:,2,3)**2 + disp(:,3,3)**2
       end select           

       do i = 1, nDimensions
          temp(:,1) = F
          call fourth_derivative_apply(i, temp(:,1:1))
          shearViscosity = shearViscosity + temp(:,1) * delta2(:,i)
       end do

       ! Combine terms
       shearViscosity = density * abs(shearViscosity)

       ! Apply filter
       temp(:,1) = shearViscosity
       do i = 1, nDimensions
          call gaussian_filter_apply(i, temp(:,1:1))
       end do
       shearViscosity = cMu * temp(:,1)

    end if
    
    ! Compute artificial bulk viscosity
    ! ---------------------------------
    bulkViscosity = 0.0_WP
    if (cBeta .gt. 0.0_WP) then
       
       ! Compute dilatation and vorticity
       allocate(vorticity(nGridPoints, 3))
       call compute_dilatation(velocityGradient, F)
       call compute_vorticity(velocityGradient, vorticity)

       ! Compute normalized density gradient
       call gradient(density, grad)
       call normalize_vector(grad)

       ! Resolution length scale
       select case (nDimensions)
       case(1)
          delta2(:,1) = (disp(:,1,1) * grad(:,1))**2
       case(2)
          delta2(:,1) = (disp(:,1,1) * grad(:,1) + disp(:,2,1) * grad(:,2))**2
          delta2(:,2) = (disp(:,1,2) * grad(:,1) + disp(:,2,2) * grad(:,2))**2
       case(3)
          delta2(:,1) = (disp(:,1,1) * grad(:,1) + disp(:,2,1) * grad(:,2) +                 &
               disp(:,3,1) * grad(:,3))**2
          delta2(:,2) = (disp(:,1,2) * grad(:,1) + disp(:,2,2) * grad(:,2) +                 &
               disp(:,3,2) * grad(:,3))**2 
          delta2(:,3) = (disp(:,1,3) * grad(:,1) + disp(:,2,3) * grad(:,2) +                 &
               disp(:,3,3) * grad(:,3))**2 
       end select

       do i = 1, nDimensions
          temp(:,1) = F
          call fourth_derivative_apply(i, temp(:,1:1))
          bulkViscosity = bulkViscosity + temp(:,1) * delta2(:,i)
       end do

       ! Combine terms
       do i = 1, nGridPoints
          
          ! Sensor originally proposed by Ducros et al. (1999) and later improved by
          ! Hendrickson, T. R., Kartha, A., & Candler, G. V. (2018)
          heavyside = sign(0.5_WP, -F(i)) + 0.5_WP
          soundSpeed = sqrt(ratioOfSpecificHeats * pressure(i) / density(i))
          div2 = F(i)**2
          vort2 = sum(vorticity(i,:)**2)
          vort2 = max(vort2, (0.05_WP * soundSpeed / minval(gridSpacing(i,:)))**2)
          sensor = heavyside * min(4.0_WP / 3.0_WP * div2 / (div2 + vort2 + eps), 1.0_WP)

          bulkViscosity(i) = sensor * density(i) * abs(bulkViscosity(i))
       end do

       ! Apply filter
       temp(:,1) = bulkViscosity
       do i = 1, nDimensions
          call gaussian_filter_apply(i, temp(:,1:1))
       end do
       bulkViscosity = cBeta * temp(:,1)

       ! Cleanup
       deallocate(vorticity)
    end if

    ! Compute artificial diffusivity for thermal conductivity
    ! -------------------------------------------------------
    thermalDiffusivity = 0.0_WP
    if (cKappa .gt. 0.0_WP) then
       
       ! Compute internal energy and its gradient
       F = pressure / density / (ratioOfSpecificHeats - 1.0_WP)
       call gradient(F, grad)
       call normalize_vector(grad)

       ! Resolution length scale
       select case (nDimensions)
       case(1)
          delta2(:,1) = abs(disp(:,1,1) * grad(:,1))
       case(2)
          delta2(:,1) = abs(disp(:,1,1) * grad(:,1) + disp(:,2,1) * grad(:,2))
          delta2(:,2) = abs(disp(:,1,2) * grad(:,1) + disp(:,2,2) * grad(:,2))
       case(3)
          delta2(:,1) = abs(disp(:,1,1) * grad(:,1) + disp(:,2,1) * grad(:,2) +              &
               disp(:,3,1) * grad(:,3))
          delta2(:,2) = abs(disp(:,1,2) * grad(:,1) + disp(:,2,2) * grad(:,2) +              &
               disp(:,3,2) * grad(:,3)) 
          delta2(:,3) = abs(disp(:,1,3) * grad(:,1) + disp(:,2,3) * grad(:,2) +              &
               disp(:,3,3) * grad(:,3)) 
       end select

       do i = 1, nDimensions
          temp(:,1) = F
          call fourth_derivative_apply(i, temp(:,1:1))
           thermalDiffusivity =  thermalDiffusivity + temp(:,1) * delta2(:,i)
       end do

       ! Combine terms
       do i = 1, nGridPoints

          ! Sound speed
          soundSpeed = sqrt(ratioOfSpecificHeats * pressure(i) / density(i))

          thermalDiffusivity(i) = density(i) * soundSpeed / temperature(i) *                 &
               abs(thermalDiffusivity(i))
       end do

       ! Apply filter
       temp(:,1) = thermalDiffusivity
       do i = 1, nDimensions
          call gaussian_filter_apply(i, temp(:,1:1))
       end do
       thermalDiffusivity = cKappa * temp(:,1)

    end if

    ! Clean up
    deallocate(disp, delta2, F, grad, temp)

    return
  end subroutine local_artificial_diffusivity


  ! Compute a sensor that varies between 0 and 1
  ! ---------------------------------------------
  subroutine compute_sensor(sensorType, pressure, specificVolume, velocity,                  &
       velocityGradient, massFraction, sensor)

    ! External modules
    use grid

    implicit none

    ! Arguments
    integer, intent(in) :: sensorType
    real(WP), intent(in) :: pressure(:), specificVolume(:), velocity(:,:),                   &
         velocityGradient(:,:), massFraction(:,:)
    real(WP), intent(out) :: sensor(:)

    ! Local variables
    integer :: i, j
    real(WP), parameter :: eps = 1.0e-12_WP
    real(WP) :: div2, vort2, soundSpeed, minValue, maxValue
    real(WP), allocatable :: dilatation(:), vorticity(:,:)

    select case (sensorType)

    case (DUCROS_SENSOR)
       ! Sensor originally proposed by Ducros et al. (1999) and later improved by
       ! Hendrickson, T. R., Kartha, A., & Candler, G. V. (2018)

       ! Compute dilatation and vorticity
       allocate(dilatation(nGridPoints))
       allocate(vorticity(nGridPoints, 3))
       call compute_dilatation(velocityGradient, dilatation)
       call compute_vorticity(velocityGradient, vorticity)

       do i = 1, nGridPoints
          soundSpeed = sqrt(ratioOfSpecificHeats * pressure(i) * specificVolume(i))
          div2 = dilatation(i)**2
          vort2 = sum(vorticity(i,:)**2)
          vort2 = max(vort2, (0.05_WP * soundSpeed / minval(gridSpacing(i,:)))**2)
          sensor(i) = min(4.0_WP / 3.0_WP * div2 / (div2 + vort2 + eps), 1.0_WP)
       end do

       ! Clean up
       deallocate(dilatation)
       deallocate(vorticity)

    case (BOUNDED_SENSOR)
       ! Preserve scalar boundedness
       minValue = 0.0_WP
       maxValue = 1.0_WP
       sensor = 0.0_WP
       do i = 1, nGridPoints
          do j = 1, size(massFraction, 2)
             if (massFraction(i,j).lt.minValue .or. massFraction(i,j).gt.maxValue) then
                sensor(i) = 1.0_WP
                exit
             end if
          end do
       end do

    end select

    return
  end subroutine compute_sensor


  ! Compute the vorticity magnitude
  ! -------------------------------
  subroutine compute_vorticity(velocityGradient, vorticity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocityGradient(:,:)
    real(WP), intent(out) :: vorticity(:,:)

    vorticity = 0.0_WP

    select case (nDimensions)

    case (1)
       vorticity = 0.0_WP

    case (2)
       vorticity(:,3) = velocityGradient(:,3) - velocityGradient(:,2)

    case (3)
       vorticity(:,1) = velocityGradient(:,8) - velocityGradient(:,6)
       vorticity(:,2) = velocityGradient(:,3) - velocityGradient(:,7)
       vorticity(:,3) = velocityGradient(:,4) - velocityGradient(:,2)

    end select

    return
  end subroutine compute_vorticity


  ! Compute dilitation
  ! ------------------
  subroutine compute_dilatation(velocityGradient, dilatation)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocityGradient(:,:)
    real(WP), intent(out) :: dilatation(:)

    select case (nDimensions)

    case (1)
       dilatation = velocityGradient(:,1)

    case (2)
       dilatation = velocityGradient(:,1) + velocityGradient(:,4)

    case (3)
       dilatation = velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9)

    end select

    return
  end subroutine compute_dilatation


  ! Compute the Q-creiterion
  ! Originally defined for incompressible flow by Hunt et al. (1988) CTR
  ! Extended to compressible flow by Kolár et al. (2015) AIAA
  ! --------------------------------------------------------------------
  subroutine compute_Q_criterion(velocityGradient, QCriterion)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocityGradient(:,:)
    real(WP), intent(out) :: QCriterion(:)

    ! Local variables
    integer :: i
    real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
    real(WP), dimension(nDimensions ** 2) :: rateOfStrainTensor, vorticityTensor
    real(WP) :: trace

    select case (nDimensions)

    case (1)
       do i = 1, nGridPoints
          trace = velocityGradient(i,1)
          rateOfStrainTensor(1) = velocityGradient(i,1) - oneThird * trace
          
          vorticityTensor(1) = 0.0_WP
          
          QCriterion(i) = 0.5_WP * (sum(vorticityTensor**2) - sum(rateOfStrainTensor**2))
       end do

    case (2)
       do i = 1, size(QCriterion)
          trace = velocityGradient(i,1) + velocityGradient(i,4)
          rateOfStrainTensor(1) = velocityGradient(i,1) - oneThird * trace
          rateOfStrainTensor(2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,3))
          rateOfStrainTensor(3) = rateOfStrainTensor(2)
          rateOfStrainTensor(4) = velocityGradient(i,4) - oneThird * trace

          vorticityTensor(1) = 0.0_WP
          vorticityTensor(2) = 0.5_WP * (velocityGradient(i,2) - velocityGradient(i,3))
          vorticityTensor(3) = 0.5_WP * (velocityGradient(i,3) - velocityGradient(i,2))
          vorticityTensor(4) = 0.0_WP

          QCriterion(i) = 0.5_WP * (sum(vorticityTensor**2) - sum(rateOfStrainTensor**2))
       end do

    case (3)
       do i = 1, size(QCriterion)
          trace = velocityGradient(i,1) + velocityGradient(i,5) + velocityGradient(i,9)
          rateOfStrainTensor(1) = velocityGradient(i,1) - oneThird * trace
          rateOfStrainTensor(2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,4))
          rateOfStrainTensor(3) = 0.5_WP * (velocityGradient(i,3) + velocityGradient(i,7))
          rateOfStrainTensor(4) = rateOfStrainTensor(2)
          rateOfStrainTensor(5) = velocityGradient(i,5) - oneThird * trace
          rateOfStrainTensor(6) = 0.5_WP * (velocityGradient(i,6) + velocityGradient(i,8))
          rateOfStrainTensor(7) = rateOfStrainTensor(3)
          rateOfStrainTensor(8) = rateOfStrainTensor(6)
          rateOfStrainTensor(9) = velocityGradient(i,9) - oneThird * trace

          vorticityTensor(1) = 0.0_WP
          vorticityTensor(2) = 0.5_WP * (velocityGradient(i,2) - velocityGradient(i,4))
          vorticityTensor(3) = 0.5_WP * (velocityGradient(i,3) - velocityGradient(i,7))
          vorticityTensor(4) = 0.5_WP * (velocityGradient(i,4) - velocityGradient(i,2))
          vorticityTensor(5) = 0.0_WP
          vorticityTensor(6) = 0.5_WP * (velocityGradient(i,6) - velocityGradient(i,8))
          vorticityTensor(7) = 0.5_WP * (velocityGradient(i,7) - velocityGradient(i,3))
          vorticityTensor(8) = 0.5_WP * (velocityGradient(i,8) + velocityGradient(i,6))
          vorticityTensor(9) = 0.0_WP

          QCriterion(i) = 0.5_WP * (sum(vorticityTensor**2) - sum(rateOfStrainTensor**2))
       end do

    end select

    return
  end subroutine compute_Q_criterion


  ! Compute the rate-of-strain tensor
  ! ----------------------------------
  subroutine compute_strain_rate(velocityGradient, rateOfStrainTensor)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocityGradient(:,:)
    real(WP), intent(out) :: rateOfStrainTensor(:,:)

    ! Local variables
    integer :: i

    if (size(velocityGradient) .ne. size(rateOfStrainTensor))                                &
         call die('compute_strain_rate: inconsistent size of arrays')

    select case (nDimensions)

    case (1)
       do i = 1, size(velocityGradient(:,1))
          rateOfStrainTensor(i,1) = velocityGradient(i,1)
       end do

    case (2)
       do i = 1, size(velocityGradient(:,1))
          rateOfStrainTensor(i,1) = velocityGradient(i,1)
          rateOfStrainTensor(i,2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,3))
          rateOfStrainTensor(i,3) = rateOfStrainTensor(i,2)
          rateOfStrainTensor(i,4) = velocityGradient(i,4)
       end do

    case (3)
       do i = 1, size(velocityGradient(:,1))
          rateOfStrainTensor(i,1) = velocityGradient(i,1)
          rateOfStrainTensor(i,2) = 0.5_WP * (velocityGradient(i,2) + velocityGradient(i,4))
          rateOfStrainTensor(i,3) = 0.5_WP * (velocityGradient(i,3) + velocityGradient(i,7))
          rateOfStrainTensor(i,4) = rateOfStrainTensor(i,2)
          rateOfStrainTensor(i,5) = velocityGradient(i,5)
          rateOfStrainTensor(i,6) = 0.5_WP * (velocityGradient(i,6) + velocityGradient(i,8))
          rateOfStrainTensor(i,7) = rateOfStrainTensor(i,3)
          rateOfStrainTensor(i,8) = rateOfStrainTensor(i,6)
          rateOfStrainTensor(i,9) = velocityGradient(i,9)
       end do

    end select

    return
  end subroutine compute_strain_rate


  ! Compute the vorticity thickness in a given direction
  ! ----------------------------------------------------
  subroutine compute_vorticity_thickness(density, momentum, velocityDifference, direction,   &
       vorticityThickness)

    ! External modules
    use parallel
    use grid_functions

    implicit none

    ! Arguments
    real(WP), intent(in) :: density(:), momentum(:), velocityDifference
    integer, intent(in) :: direction
    real(WP), intent(out) :: vorticityThickness

    ! Local variables
    integer :: i, j, k, gridIndex
    real(WP), allocatable :: momentumGradient(:,:), meanDensity(:), dudy(:)

    ! Compute momentum gradient
    allocate(momentumGradient(nGridPoints, nDimensions))
    call gradient(momentum, momentumGradient)

    ! Favre average
    allocate(dudy(globalGridSize(direction)))
    allocate(meanDensity(globalGridSize(direction)))
    dudy = 0.0_WP
    meanDensity = 0.0_WP

    select case(direction)

    case(1)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(i) = meanDensity(i) + density(gridIndex) * gridNorm(gridIndex, 1)
                dudy(i) = dudy(i) + momentumGradient(gridIndex, direction) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (2)

       do k = iStart(3), iEnd(3)
          do i = iStart(1), iEnd(1)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(j) = meanDensity(j) + density(gridIndex) * gridNorm(gridIndex, 1)
                dudy(j) = dudy(j) + momentumGradient(gridIndex, direction) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (3)

       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             do k = iStart(3), iEnd(3)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(k) = meanDensity(k) + density(gridIndex) * gridNorm(gridIndex, 1)
                dudy(k) = dudy(k) + momentumGradient(gridIndex, direction) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    end select

    ! Sum them up
    call parallel_sum(meanDensity)
    call parallel_sum(dudy)

    ! Compute the Favre average velocity gradient
    dudy = dudy / meanDensity

    ! Compute the vorticity thickness
    vorticityThickness = velocityDifference / maxval(abs(dudy))

    ! Cleanup
    deallocate(momentumGradient, meanDensity, dudy)

    return
  end subroutine compute_vorticity_thickness


  ! Compute the momentum thickness in a given direction
  ! ----------------------------------------------------
  subroutine compute_momentum_thickness(stateVector, velocityDifference, normalDirection,    &
       flowDirection, momentumThickness)

    ! External modules
    use parallel
    use grid_functions

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:), velocityDifference
    integer, intent(in) :: normalDirection, flowDirection
    real(WP), intent(out) :: momentumThickness

    ! Local variables
    integer :: i, j, k, ijk(3), gridIndex
    real(WP), dimension(:), allocatable :: volume, meanDensity, favreVelocity
    real(WP) :: normFactor

    ! Prepare Favre average
    allocate(volume(globalGridSize(normalDirection)))
    allocate(meanDensity(globalGridSize(normalDirection)))
    allocate(favreVelocity(globalGridSize(normalDirection)))
    volume = 0.0_WP
    meanDensity = 0.0_WP
    favreVelocity = 0.0_WP

    select case(normalDirection)

    case(1)

       normFactor = real(nProcsDir(2) * nProcsDir(3), WP)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                volume(i) = volume(i) + gridNorm(gridIndex, 1)
                meanDensity(i) = meanDensity(i) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                favreVelocity(i) = favreVelocity(i) +                                        &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (2)

       normFactor = real(nProcsDir(1) * nProcsDir(3), WP)

       do k = iStart(3), iEnd(3)
          do i = iStart(1), iEnd(1)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                volume(j) = volume(j) + gridNorm(gridIndex, 1)
                meanDensity(j) = meanDensity(j) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                favreVelocity(j) = favreVelocity(j) +                                        &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (3)

       normFactor = real(nProcsDir(1) * nProcsDir(2), WP)

       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             do k = iStart(3), iEnd(3)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                volume(k) = volume(k) + gridNorm(gridIndex, 1)
                meanDensity(k) = meanDensity(k) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                favreVelocity(k) = favreVelocity(k) +                                        &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    end select

    ! Sum them up
    call parallel_sum(volume)
    call parallel_sum(meanDensity)
    call parallel_sum(favreVelocity)

    ! Normalize the Favre average velocity and mean density
    do i = 1, globalGridSize(normalDirection)
       favreVelocity(i) = favreVelocity(i) / meanDensity(i)
       meanDensity(i) = meanDensity(i) / volume(i)
    end do

    ! Compute the momentum thickness
    momentumThickness = 0.0_WP
    ijk = iStart
    do i = iStart(normalDirection), iEnd(normalDirection)
       ijk(normalDirection) = i
       momentumThickness = momentumThickness + meanDensity(i) * (0.25_WP *                   &
            velocityDifference**2 - favreVelocity(i)**2) *                                   &
            gridSpacing(grid_index(ijk(1), ijk(2), ijk(3)), normalDirection)
    end do
    call parallel_sum(momentumThickness)

    momentumThickness = momentumThickness / velocityDifference**2 / normFactor

    ! Cleanup
    deallocate(volume, meanDensity, favreVelocity)

    return
  end subroutine compute_momentum_thickness

  
  ! Compute the time rate of change of the momentum thickness in a given direction
  ! From 'Compressible mixing layer growth rate and turbulence characteristics' by
  ! Vreman (1996)
  ! ------------------------------------------------------------------------------
  subroutine compute_shear_growth_rate(stateVector, velocityDifference,                      &
       normalDirection, flowDirection, growthRate)

    ! External modules
    use parallel
    use grid_functions

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:,:), velocityDifference
    integer, intent(in) :: normalDirection, flowDirection
    real(WP), intent(out) :: growthRate

    ! Local variables
    integer :: i, j, k, ijk(3), gridIndex
    real(WP), dimension(:), allocatable :: meanDensity, reynoldsStress, dudy
    real(WP), dimension(:,:), allocatable :: favreVelocity, momentumGradient
    real(WP) :: normFactor, u1, u2

    ! Compute momentum gradient
    allocate(momentumGradient(nGridPoints, nDimensions))
    call gradient(stateVector(:,flowDirection+1), momentumGradient)

    ! Prepare Favre average
    allocate(meanDensity(globalGridSize(normalDirection)))
    allocate(dudy(globalGridSize(normalDirection)))
    allocate(favreVelocity(globalGridSize(normalDirection),2))
    allocate(reynoldsStress(globalGridSize(normalDirection)))
    meanDensity = 0.0_WP
    dudy = 0.0_WP
    favreVelocity = 0.0_WP

    select case(normalDirection)

    case(1)

       normFactor = real(nProcsDir(2) * nProcsDir(3), WP)

       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(i) = meanDensity(i) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                dudy(i) = dudy(i) + momentumGradient(gridIndex, normalDirection) *           &
                     gridNorm(gridIndex, 1)
                favreVelocity(i,1) = favreVelocity(i,1) +                                    &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
                favreVelocity(i,2) = favreVelocity(i,2) +                                    &
                     stateVector(gridIndex, normalDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (2)

       normFactor = real(nProcsDir(1) * nProcsDir(3), WP)

       do k = iStart(3), iEnd(3)
          do i = iStart(1), iEnd(1)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(j) = meanDensity(j) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                dudy(j) = dudy(j) + momentumGradient(gridIndex, normalDirection) *           &
                     gridNorm(gridIndex, 1)
                favreVelocity(j,1) = favreVelocity(j,1) +                                    &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
                favreVelocity(j,2) = favreVelocity(j,2) +                                    &
                     stateVector(gridIndex, normalDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case (3)

       normFactor = real(nProcsDir(1) * nProcsDir(2), WP)

       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             do k = iStart(3), iEnd(3)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                meanDensity(k) = meanDensity(k) + stateVector(gridIndex, 1) *                &
                     gridNorm(gridIndex, 1)
                dudy(k) = dudy(k) + momentumGradient(gridIndex, normalDirection) *           &
                     gridNorm(gridIndex, 1)
                favreVelocity(k,1) = favreVelocity(k,1) +                                    &
                     stateVector(gridIndex, flowDirection+1) * gridNorm(gridIndex, 1)
                favreVelocity(k,2) = favreVelocity(k,2) +                                    &
                     stateVector(gridIndex, normalDirection+1) * gridNorm(gridIndex, 1)
             end do
          end do
       end do

    end select

    ! Sum them up
    call parallel_sum(meanDensity)
    call parallel_sum(dudy)
    call parallel_sum(favreVelocity(:,1))
    call parallel_sum(favreVelocity(:,2))

    ! Normalize the Favre average velocity and mean density
    do i = 1, globalGridSize(normalDirection)
       dudy(i) = dudy(i) / meanDensity(i)
       favreVelocity(i,1) = favreVelocity(i,1) / meanDensity(i)
       favreVelocity(i,2) = favreVelocity(i,2) / meanDensity(i)
    end do

    ! Compute the Reynolds Stress
    reynoldsStress = 0.0_WP
    select case (normalDirection)

    case(1)
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                u1 = stateVector(gridIndex, flowDirection+1) / stateVector(gridIndex, 1)
                u2 = stateVector(gridIndex, normalDirection+1) / stateVector(gridIndex, 1)
                reynoldsStress(i) = reynoldsStress(i) + stateVector(gridIndex, 1) *          &
                     (u1 - favreVelocity(i,1)) * (u2 - favreVelocity(i,2)) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case(2)
       do k = iStart(3), iEnd(3)
          do i = iStart(1), iEnd(1)
             do j = iStart(2), iEnd(2)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                u1 = stateVector(gridIndex, flowDirection+1) / stateVector(gridIndex, 1)
                u2 = stateVector(gridIndex, normalDirection+1) / stateVector(gridIndex, 1)
                reynoldsStress(j) = reynoldsStress(j) + stateVector(gridIndex, 1) *          &
                     (u1 - favreVelocity(j,1)) * (u2 - favreVelocity(j,2)) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    case(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             do k = iStart(3), iEnd(3)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                u1 = stateVector(gridIndex, flowDirection+1) / stateVector(gridIndex, 1)
                u2 = stateVector(gridIndex, normalDirection+1) / stateVector(gridIndex, 1)
                reynoldsStress(k) = reynoldsStress(k) + stateVector(gridIndex, 1) *          &
                     (u1 - favreVelocity(k,1)) * (u2 - favreVelocity(k,2)) *                 &
                     gridNorm(gridIndex, 1)
             end do
          end do
       end do

    end select
    do i = 1, globalGridSize(normalDirection)
       reynoldsStress(i) = reynoldsStress(i) / meanDensity(i)
    end do

    ! Integrate
    growthRate = 0.0_WP
    ijk = iStart
    do i = iStart(normalDirection), iEnd(normalDirection)
       ijk(normalDirection) = i
       growthRate = growthRate + reynoldsStress(i) * dudy(i) *                               &
            gridSpacing(grid_index(ijk(1), ijk(2), ijk(3)), normalDirection)
    end do
    call parallel_sum(growthRate)

    growthRate = -2.0_WP * growthRate / velocityDifference**2 / normFactor

    ! Cleanup
    deallocate(meanDensity, favreVelocity, reynoldsStress, dudy, momentumGradient)

    return
  end subroutine compute_shear_growth_rate


  subroutine compute_schlieren(density, direction, component)

    ! External modules
    use parallel
    use grid_functions
    use grid_levelset

    implicit none

    ! Arguments
    real(WP), intent(inout) :: density(:)
    integer, intent(in) :: direction, component

    ! Local variables
    integer :: i, j, k, ii, jj, kk, gridIndex
    real(WP), dimension(nGridPoints, 3) :: densityGradient
    real(WP), dimension(:,:), allocatable :: schlieren

    ! Compute density gradient
    densityGradient = 0.0_WP
    call gradient(density, densityGradient(:, 1:nDimensions))

    ! Select component of density gradient
    if (component .gt. 0) then
       density = densityGradient(:, component)
    else
       select case (direction)
       case (0)
          density = sqrt(sum(densityGradient(:,1:nDimensions)**2, dim = 2))
       case (1)
          density = sqrt(densityGradient(:, 2)**2 + densityGradient(:, 3)**2)
       case (2)
          density = sqrt(densityGradient(:, 1)**2 + densityGradient(:, 3)**2)
       case (3)
          density = sqrt(densityGradient(:, 1)**2 + densityGradient(:, 2)**2)
       end select
    end if

    ! Nothing left to do if direction = 0
    if (direction .eq. 0) return

    ! Integrate along the line of site
    if (globalGridSize(direction) .gt. 1) then

       ! Not yet implemented for curvilinear coordinates
       if (isDomainCurvilinear)                                                              &
            call die('compute_schlieren not implemented for curvilinear coordinates!')

       ! Integrate along line-of-sight
       select case (direction)
       case (1)
          allocate(schlieren(localGridSize(2), localGridSize(3)))
       case (2)
          allocate(schlieren(localGridSize(1), localGridSize(3)))
       case (3)
          allocate(schlieren(localGridSize(1), localGridSize(2)))
       end select
       schlieren = 0.0_WP
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                if (allocated(levelset)) then
                   if (levelset(gridIndex,1) .le. 0.0_WP) cycle
                end if
                ii = i - gridOffset(1); jj = j - gridOffset(2); kk = k - gridOffset(3)
                select case (direction)
                case (1)
                   schlieren(jj,kk) = schlieren(jj,kk) + density(gridIndex) *                &
                        gridSpacing(gridIndex,1)
                case (2)
                   schlieren(ii,kk) = schlieren(ii,kk) + density(gridIndex) *                &
                        gridSpacing(gridIndex,2)
                case (3)
                   schlieren(ii,jj) = schlieren(ii,jj) + density(gridIndex) *                &
                        gridSpacing(gridIndex,3)
                end select
             end do
          end do
       end do

       ! Sum up over processors along line-of-sight
       call parallel_sum_dir(schlieren, direction)

       ! Transfer back
       do k = iStart(3), iEnd(3)
          do j = iStart(2), iEnd(2)
             do i = iStart(1), iEnd(1)
                gridIndex = i - gridOffset(1) + localGridSize(1) *                           &
                     (j - 1 - gridOffset(2) + localGridSize(2) *                             &
                     (k - 1 - gridOffset(3)))
                ii = i - gridOffset(1); jj = j - gridOffset(2); kk = k - gridOffset(3)
                select case (direction)
                case (1)
                   density(gridIndex) = schlieren(jj,kk)
                case (2)
                   density(gridIndex) = schlieren(ii,kk)
                case (3)
                   density(gridIndex) = schlieren(ii,jj)
                end select
             end do
          end do
       end do

       ! Clean up
       deallocate(schlieren)

    end if

    ! Scale the value for proper shading
    if (component .eq. 0) call scale_schlieren(density)

  contains

    ! Non-linear scaling
    ! Quirk JJ (1994) A contribution to the great Riemann solver debate
    ! Int J Number Methods in Fluids 18: 555-574
    subroutine scale_schlieren(shade)
      implicit none

      ! Arguments
      real(WP), dimension(:), intent(inout) :: shade

      ! Local avariables
      integer :: i
      real(WP), parameter :: c0 = -0.001_WP, c1 = 0.05_WP, ck = 5.0_WP
      real(WP) :: rmax

      rmax = maxval(shade)
      call parallel_max(rmax)

      do i = 1, nGridPoints
         shade(i) = exp(-ck * (shade(i) - c0 * rmax) / (c1 * rmax - c0 * rmax))
      end do

      return
    end subroutine scale_schlieren

  end subroutine compute_schlieren

end module state_functions
