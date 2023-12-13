module state_jacobian

  ! External modules
  use precision
  use simulation_flags
  use solver_options
  use geometry

  implicit none

contains

  pure subroutine compute_delta_variables(conservedVariables, dynamicViscosity,              &
       deltaConservedVariables, deltaSpecificVolume, deltaVelocity, deltaPressure,           &
       deltaTemperature, deltaMassFraction, deltaViscosity)

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:)
    real(WP), intent(in), optional :: dynamicViscosity
    real(WP), intent(out), optional :: deltaConservedVariables(:,:),                         &
         deltaSpecificVolume(:), deltaVelocity(:,:), deltaPressure(:), deltaTemperature(:),  &
         deltaMassFraction(:,:), deltaViscosity(:)

    ! Local variables
    integer :: i, k, nUnknowns
    real(WP) :: temp

    nUnknowns = nDimensions + nSpecies + 2

    ! Variation of conservedVariables
    if (present(deltaConservedVariables)) then
       deltaConservedVariables = 0.0_WP
       do i = 1, nUnknowns
          deltaConservedVariables(i,i) = 1.0_WP
       end do
    end if

    ! Variation of specific volume
    if (present(deltaSpecificVolume)) then
       deltaSpecificVolume = 0.0_WP
       deltaSpecificVolume(1) = -1.0_WP / conservedVariables(1) ** 2
    end if

    ! Variation of velocity
    if (present(deltaVelocity)) then
       deltaVelocity = 0.0_WP
       do i = 1, nDimensions
          deltaVelocity(i,1) = - conservedVariables(1+i) / conservedVariables(1) ** 2
          deltaVelocity(i,i+1) =  1.0_WP / conservedVariables(1)
       end do
    end if

    ! Variation of pressure
    if (present(deltaPressure)) then
       deltaPressure = 0.0_WP
       deltaPressure(1) = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                         &
            sum(conservedVariables(2:nDimensions+1) **2) / conservedVariables(1) ** 2
       do i = 1, nDimensions
          deltaPressure(i+1) =  - (ratioOfSpecificHeats - 1.0_WP) *                          &
               conservedVariables(i+1) / conservedVariables(1)
       end do
       deltaPressure(nDimensions+2) = (ratioOfSpecificHeats - 1.0_WP)
    end if

    ! Variation of mass fraction
    if (present(deltaMassFraction)) then
       deltaMassFraction = 0.0_WP
       do k = 1, nSpecies
          deltaMassFraction(k,1) = - conservedVariables(nDimensions+2+k) /                   &
               conservedVariables(1) ** 2
          deltaMassFraction(nDimensions+2+k,nDimensions+2+k) = 1.0_WP / conservedVariables(1)
       end do
    end if

    ! Variation of temperature
    if (present(deltaTemperature)) then
       deltaTemperature = 0.0_WP

       if (equationOfState .eq. IDEAL_GAS) then ! ... single-component gas

          deltaTemperature(1) = ratioOfSpecificHeats / conservedVariables(1)**2 * (          &
               sum(conservedVariables(2:nDimensions+1) **2) / conservedVariables(1) -        &
               conservedVariables(nDimensions+2) )
          do i = 1, nDimensions
             deltaTemperature(i+1) = - ratioOfSpecificHeats * conservedVariables(i+1) /      &
                  conservedVariables(1)**2
          end do
          deltaTemperature(nDimensions+2) = ratioOfSpecificHeats / conservedVariables(1)

       else ! ... ideal gas mixture
          temp = conservedVariables(1) * molecularWeightInverse(nSpecies+1)
          do k = 1, nSpecies
             temp = temp + conservedVariables(nDimensions+2+k) *                             &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
          end do
          temp = 1.0_WP / temp
          deltaTemperature(1) = 0.5_WP * ratioOfSpecificHeats *                              &
               sum(conservedVariables(2:nDimensions+1)**2) / conservedVariables(1)**2 *      &
               temp - ratioOfSpecificHeats * molecularWeightInverse(nSpecies+1) * temp**2 *  &
               (conservedVariables(nDimensions+2) - 0.5_WP *                                 &
               sum(conservedVariables(2:nDimensions+1)**2) / conservedVariables(1))
          do i = 1, nDimensions
             deltaTemperature(i+1) = - ratioOfSpecificHeats * conservedVariables(i+1) /      &
                  conservedVariables(1) * temp
          end do
          deltaTemperature(nDimensions+2) = ratioOfSpecificHeats * temp
          do k = 1, nSpecies
             deltaTemperature(nDimensions+2+k) = - ratioOfSpecificHeats * temp**2 *          &
                  (conservedVariables(nDimensions+2) - 0.5_WP *                              &
                  sum(conservedVariables(2:nDimensions+1) **2) / conservedVariables(1)) *    &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
          end do

       end if
    end if

    return
  end subroutine compute_delta_variables

  pure subroutine compute_jacobian_of_inviscid_flux(conservedVariables, metrics,             &
       jacobianOfInviscidFlux, deltaConservedVariables, specificVolume, velocity, pressure,  &
       massFraction, deltaJacobianOfInviscidFlux)

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:), metrics(:)
    real(WP), intent(out) :: jacobianOfInviscidFlux(:,:)
    real(WP), intent(in), optional :: deltaConservedVariables(:,:), specificVolume,         &
         velocity(:), pressure, massFraction(:)
    real(WP), intent(out), optional :: deltaJacobianOfInviscidFlux(:,:,:)

    ! Local variables
    integer :: k
    real(WP) :: specificVolume_, velocity_(nDimensions), pressure_,                          &
         massFraction_(nSpecies), contravariantVelocity, phiSquared, enthalpy,               &
         deltaConservedVariables_(nDimensions+nSpecies+2, nDimensions+nSpecies+2),           &
         deltaSpecificVolume(nDimensions+nSpecies+2),                                        &
         deltaVelocity(nDimensions, nDimensions+nSpecies+2),                                 &
         deltaPressure(nDimensions+nSpecies+2),                                              &
         deltaMassFraction(nSpecies, nDimensions+nSpecies+2),                                &
         deltaContravariantVelocity(nDimensions+nSpecies+2),                                 &
         deltaPhiSquared(nDimensions+nSpecies+2), deltaEnthalpy(nDimensions+nSpecies+2)

    ! Compute specific volume if it was not specified
    if (present(specificVolume)) then
       specificVolume_ = specificVolume
    else
       specificVolume_ = 1.0_WP / conservedVariables(1)
    end if

    ! Compute velocity if it was not specified
    if (present(velocity)) then
       velocity_ = velocity
    else
       do k = 1, nDimensions
          velocity_(k) = specificVolume_ * conservedVariables(k+1)
       end do
    end if

    ! Compute pressure if it was not specified
    if (present(pressure)) then
       pressure_ = pressure
    else
       pressure_ = (ratioOfSpecificHeats - 1.0_WP) *                                         &
            (conservedVariables(nDimensions+2) - 0.5_WP * conservedVariables(1) *            &
            sum(velocity_ ** 2))
    end if

    ! Compute mass fraction if it was not specified
    if (present(massFraction) .and. nSpecies .gt. 0) then
       massFraction_ = massFraction
    else
       do k = 1, nSpecies
          massFraction_(k) = conservedVariables(nDimensions+2+k) * specificVolume_
       end do
    end if

    ! Other dependent variables
    contravariantVelocity = sum(metrics * velocity_) !... not normalized
    phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * sum(velocity_**2)
    enthalpy = specificVolume_ * (conservedVariables(nDimensions+2) + pressure_)

    ! Zero-out Jacobian of inviscid flux
    jacobianOfInviscidFlux = 0.0_WP

    ! Compute Jacobian of inviscid flux
    select case (nDimensions)

    case (1)

       jacobianOfInviscidFlux(1,1) = 0.0_WP
       jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                               &
            contravariantVelocity * velocity_(1)
       jacobianOfInviscidFlux(3,1) = contravariantVelocity * (phiSquared - enthalpy)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(3+k,1) = -massFraction_(k) * contravariantVelocity
       end do

       jacobianOfInviscidFlux(1,2) = metrics(1)
       jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(1) * metrics(1)
       jacobianOfInviscidFlux(3,2) = enthalpy * metrics(1) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(1)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(3+k,2) = massFraction_(k) * metrics(1)
       end do

       jacobianOfInviscidFlux(1,3) = 0.0_WP
       jacobianOfInviscidFlux(2,3) = (ratioOfSpecificHeats - 1.0_WP) * metrics(1)
       jacobianOfInviscidFlux(3,3) = ratioOfSpecificHeats * contravariantVelocity

       do k = 1, nSpecies
          jacobianOfInviscidFlux(3+k,3+k) = contravariantVelocity
       end do

    case (2)

       jacobianOfInviscidFlux(1,1) = 0.0_WP
       jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                               &
            contravariantVelocity * velocity_(1)
       jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                               &
            contravariantVelocity * velocity_(2)
       jacobianOfInviscidFlux(4,1) = contravariantVelocity * (phiSquared - enthalpy)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(4+k,1) = -massFraction_(k) * contravariantVelocity
       end do

       jacobianOfInviscidFlux(1,2) = metrics(1)
       jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(1) * metrics(1)
       jacobianOfInviscidFlux(3,2) = velocity_(2) * metrics(1) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(1) * metrics(2)
       jacobianOfInviscidFlux(4,2) = enthalpy * metrics(1) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(1)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(4+k,2) = massFraction_(k) * metrics(1)
       end do

       jacobianOfInviscidFlux(1,3) = metrics(2)
       jacobianOfInviscidFlux(2,3) = velocity_(1) * metrics(2) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(2) * metrics(1)
       jacobianOfInviscidFlux(3,3) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(2) * metrics(2)
       jacobianOfInviscidFlux(4,3) = enthalpy * metrics(2) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(2)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(4+k,3) = massFraction_(k) * metrics(2)
       end do

       jacobianOfInviscidFlux(1,4) = 0.0_WP
       jacobianOfInviscidFlux(2,4) = (ratioOfSpecificHeats - 1.0_WP) * metrics(1)
       jacobianOfInviscidFlux(3,4) = (ratioOfSpecificHeats - 1.0_WP) * metrics(2)
       jacobianOfInviscidFlux(4,4) = ratioOfSpecificHeats * contravariantVelocity

       do k = 1, nSpecies
          jacobianOfInviscidFlux(4+k,4+k) = contravariantVelocity
       end do

    case (3)

       jacobianOfInviscidFlux(1,1) = 0.0_WP
       jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                               &
            contravariantVelocity * velocity_(1)
       jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                               &
            contravariantVelocity * velocity_(2)
       jacobianOfInviscidFlux(4,1) = phiSquared * metrics(3) -                               &
            contravariantVelocity * velocity_(3)
       jacobianOfInviscidFlux(5,1) = contravariantVelocity * (phiSquared - enthalpy)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(5+k,1) = -massFraction_(k) * contravariantVelocity
       end do

       jacobianOfInviscidFlux(1,2) = metrics(1)
       jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(1) * metrics(1)
       jacobianOfInviscidFlux(3,2) = velocity_(2) * metrics(1) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(1) * metrics(2)
       jacobianOfInviscidFlux(4,2) = velocity_(3) * metrics(1) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(1) * metrics(3)
       jacobianOfInviscidFlux(5,2) = enthalpy * metrics(1) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(1)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(5+k,2) = massFraction_(k) * metrics(1)
       end do

       jacobianOfInviscidFlux(1,3) = metrics(2)
       jacobianOfInviscidFlux(2,3) = velocity_(1) * metrics(2) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(2) * metrics(1)
       jacobianOfInviscidFlux(3,3) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(2) * metrics(2)
       jacobianOfInviscidFlux(4,3) = velocity_(3) * metrics(2) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(2) * metrics(3)
       jacobianOfInviscidFlux(5,3) = enthalpy * metrics(2) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(2)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(5+k,3) = massFraction_(k) * metrics(2)
       end do

       jacobianOfInviscidFlux(1,4) = metrics(3)
       jacobianOfInviscidFlux(2,4) = velocity_(1) * metrics(3) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(3) * metrics(1)
       jacobianOfInviscidFlux(3,4) = velocity_(2) * metrics(3) -                             &
            (ratioOfSpecificHeats - 1.0_WP) * velocity_(3) * metrics(2)
       jacobianOfInviscidFlux(4,4) = contravariantVelocity -                                 &
            (ratioOfSpecificHeats - 2.0_WP) * velocity_(3) * metrics(3)
       jacobianOfInviscidFlux(5,4) = enthalpy * metrics(3) -                                 &
            (ratioOfSpecificHeats - 1.0_WP) * contravariantVelocity * velocity_(3)
       do k = 1, nSpecies
          jacobianOfInviscidFlux(5+k,4) = massFraction_(k) * metrics(3)
       end do

       jacobianOfInviscidFlux(1,5) = 0.0_WP
       jacobianOfInviscidFlux(2,5) = (ratioOfSpecificHeats - 1.0_WP) * metrics(1)
       jacobianOfInviscidFlux(3,5) = (ratioOfSpecificHeats - 1.0_WP) * metrics(2)
       jacobianOfInviscidFlux(4,5) = (ratioOfSpecificHeats - 1.0_WP) * metrics(3)
       jacobianOfInviscidFlux(5,5) = ratioOfSpecificHeats * contravariantVelocity

       do k = 1, nSpecies
          jacobianOfInviscidFlux(5+k,5+k) = contravariantVelocity
       end do

    end select

    ! Compute variations in Jacobian of inviscid flux
    if (present(deltaJacobianOfInviscidFlux)) then

       ! Zero-out the variations of the Jacobian of inviscid flux
       deltaJacobianOfInviscidFlux = 0.0_WP

       select case (nDimensions)

       case (1)

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(3+k,3+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = -1.0_WP / conservedVariables(1) ** 2 *                       &
               deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables_(2,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(3,:) - &
               velocity_(1) * deltaVelocity(1,:)) + phiSquared * deltaConservedVariables_(1,:)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(3+k) +        &
                  specificVolume_ * deltaConservedVariables_(3+k,:)
          end do

          ! Compute variations of other dependent variables:
          deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(3,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(3) + pressure_)

          deltaJacobianOfInviscidFlux(1,1,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                &
               deltaContravariantVelocity * velocity_(1) - contravariantVelocity *           &
               deltaVelocity(1,:)
          deltaJacobianOfInviscidFlux(3,1,:) = deltaContravariantVelocity * (phiSquared -    &
               enthalpy) + contravariantVelocity * (deltaPhiSquared - deltaEnthalpy)
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(3+k,1,:) = -massFraction_(k) *                      &
                  deltacontravariantVelocity - deltaMassFraction(k,:) *                      &
                  contravariantVelocity
          end do

          deltaJacobianOfInviscidFlux(1,2,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(1,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,2,:) = deltaEnthalpy * metrics(1) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(1) +&
               contravariantVelocity * deltaVelocity(1,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(3+k,2,:) = deltaMassFraction(k,:) * metrics(1)
          end do

          deltaJacobianOfInviscidFlux(1,3,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,3,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(3,3,:) = ratioOfSpecificHeats *                        &
               deltaContravariantVelocity
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(3+k,3,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(3+k,3+k,:) = deltaContravariantVelocity
          end do

       case (2)

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             deltaConservedVariables_(4,4) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(4+k,4+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = -1.0_WP / conservedVariables(1) ** 2 *                       &
               deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables(2,:)
          deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                 &
               specificVolume_ * deltaConservedVariables(3,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(4,:) - &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))) +    &
               phiSquared * deltaConservedVariables_(1,:)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(4+k) +        &
                  specificVolume_ * deltaConservedVariables_(4+k,:)
          end do

          ! Compute variations of other dependent variables:
          deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:) + metrics(2) *        &
               deltaVelocity(2,:)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(4,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(4) + pressure_)

          deltaJacobianOfInviscidFlux(1,1,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                &
               deltaContravariantVelocity * velocity_(1) -                                   &
               contravariantVelocity * deltaVelocity(1,:)
          deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                &
               deltaContravariantVelocity * velocity_(2) -                                   &
               contravariantVelocity * deltaVelocity(2,:)
          deltaJacobianOfInviscidFlux(4,1,:) = deltaContravariantVelocity * (phiSquared -    &
               enthalpy) + contravariantVelocity * (deltaPhiSquared - deltaEnthalpy)
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(4+k,1,:) = -massFraction_(k) *                      &
                  deltacontravariantVelocity - deltaMassFraction(k,:) *                      &
                  contravariantVelocity
          end do

          deltaJacobianOfInviscidFlux(1,2,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(1,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(1,:) * metrics(2)
          deltaJacobianOfInviscidFlux(4,2,:) = deltaEnthalpy * metrics(1) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(1) +&
               contravariantVelocity * deltaVelocity(1,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(4+k,2,:) = deltaMassFraction(k,:) * metrics(1)
          end do

          deltaJacobianOfInviscidFlux(1,3,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(2,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(2,:) * metrics(2)
          deltaJacobianOfInviscidFlux(4,3,:) = deltaEnthalpy * metrics(2) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(2) +&
               contravariantVelocity * deltaVelocity(2,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(4+k,3,:) = deltaMassFraction(k,:) * metrics(2)
          end do

          deltaJacobianOfInviscidFlux(1,4,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,4,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(3,4,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(4,4,:) = ratioOfSpecificHeats *                        &
               deltaContravariantVelocity
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(4+k,4,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(4+k,4+k,:) = deltaContravariantVelocity
          end do

       case (3)

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             deltaConservedVariables_(4,4) = 1.0_WP
             deltaConservedVariables_(5,5) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(5+k,5+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = -1.0_WP / conservedVariables(1) ** 2 *                       &
               deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables(2,:)
          deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                 &
               specificVolume_ * deltaConservedVariables(3,:)
          deltaVelocity(3,:) = deltaSpecificVolume * conservedVariables(4) +                 &
               specificVolume_ * deltaConservedVariables(4,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(5,:) - &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +      &
               velocity_(3) * deltaVelocity(3,:))) + phiSquared * deltaConservedVariables_(1,:)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(5+k) +        &
                  specificVolume_ * deltaConservedVariables_(5+k,:)
          end do

          ! Compute variations of other dependent variables:
          deltaContravariantVelocity =                                                       &
               metrics(1) * deltaVelocity(1,:) +                                             &
               metrics(2) * deltaVelocity(2,:) +                                             &
               metrics(3) * deltaVelocity(3,:)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +      &
               velocity_(3) * deltaVelocity(3,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(5,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(5) + pressure_)

          deltaJacobianOfInviscidFlux(1,1,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                &
               deltaContravariantVelocity * velocity_(1) - contravariantVelocity *           &
               deltaVelocity(1,:)
          deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                &
               deltaContravariantVelocity * velocity_(2) - contravariantVelocity *           &
               deltaVelocity(2,:)
          deltaJacobianOfInviscidFlux(4,1,:) = deltaPhiSquared * metrics(3) -                &
               deltaContravariantVelocity * velocity_(3) - contravariantVelocity *           &
               deltaVelocity(3,:)
          deltaJacobianOfInviscidFlux(5,1,:) = deltaContravariantVelocity * (phiSquared -    &
               enthalpy) + contravariantVelocity * (deltaPhiSquared - deltaEnthalpy)
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,1,:) = -massFraction_(k) *                      &
                  deltacontravariantVelocity - deltaMassFraction(k,:) *                      &
                  contravariantVelocity
          end do

          deltaJacobianOfInviscidFlux(1,2,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(1,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(1,:) * metrics(2)
          deltaJacobianOfInviscidFlux(4,2,:) = deltaVelocity(3,:) * metrics(1) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(1,:) * metrics(3)
          deltaJacobianOfInviscidFlux(5,2,:) = deltaEnthalpy * metrics(1) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(1) +&
               contravariantVelocity * deltaVelocity(1,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,2,:) = deltaMassFraction(k,:) * metrics(1)
          end do

          deltaJacobianOfInviscidFlux(1,3,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(2,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(2,:) * metrics(2)
          deltaJacobianOfInviscidFlux(4,3,:) = deltaVelocity(3,:) * metrics(2) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(2,:) * metrics(3)
          deltaJacobianOfInviscidFlux(5,3,:) = deltaEnthalpy * metrics(2) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(2) +&
               contravariantVelocity * deltaVelocity(2,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,3,:) = deltaMassFraction(k,:) * metrics(2)
          end do

          deltaJacobianOfInviscidFlux(1,4,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,4,:) = deltaVelocity(1,:) * metrics(3) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(3,:) * metrics(1)
          deltaJacobianOfInviscidFlux(3,4,:) = deltaVelocity(2,:) * metrics(3) -             &
               (ratioOfSpecificHeats - 1.0_WP) * deltaVelocity(3,:) * metrics(2)
          deltaJacobianOfInviscidFlux(4,4,:) = deltaContravariantVelocity -                  &
               (ratioOfSpecificHeats - 2.0_WP) * deltaVelocity(3,:) * metrics(3)
          deltaJacobianOfInviscidFlux(5,4,:) = deltaEnthalpy * metrics(3) -                  &
               (ratioOfSpecificHeats - 1.0_WP) * (deltaContravariantVelocity * velocity_(3) +&
               contravariantVelocity * deltaVelocity(3,:))
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,4,:) = deltaMassFraction(k,:) * metrics(3)
          end do

          deltaJacobianOfInviscidFlux(1,5,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(2,5,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(3,5,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(4,5,:) = 0.0_WP
          deltaJacobianOfInviscidFlux(5,5,:) = ratioOfSpecificHeats *                        &
               deltaContravariantVelocity
          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,5,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaJacobianOfInviscidFlux(5+k,5+k,:) = deltaContravariantVelocity
          end do

       end select

    end if

    return
  end subroutine compute_jacobian_of_inviscid_flux

  pure subroutine compute_incoming_jacobian_of_inviscid_flux(conservedVariables, metrics,    &
       incomingDirection, incomingJacobianOfInviscidFlux,                                    &
       deltaIncomingJacobianOfInviscidFlux, deltaConservedVariables,                         &
       specificVolume, velocity, pressure, massFraction)

    implicit none

    ! Arguments
    integer, intent(in) :: incomingDirection
    real(WP), intent(in) :: conservedVariables(:), metrics(:)
    real(WP), intent(out) :: incomingJacobianOfInviscidFlux(:,:)
    real(WP), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(:,:,:)
    real(WP), intent(in), optional :: deltaConservedVariables(:,:), specificVolume,          &
         velocity(:), pressure, massFraction(:)

    ! Local variables
    integer :: i, j, k
    real(WP) :: arcLength, normalizedMetrics(nDimensions), specificVolume_,                  &
         velocity_(nDimensions), pressure_, massFraction_(nSpecies), temperature_,           &
         contravariantVelocity, speedOfSound, phiSquared, enthalpy,                          &
         rightEigenvectors(nDimensions+nSpecies+2, nDimensions+nSpecies+2),                  &
         eigenvalues(nDimensions+nSpecies+2),                                                &
         leftEigenvectors(nDimensions+nSpecies+2, nDimensions+nSpecies+2),                   &
         deltaConservedVariables_(nDimensions+nSpecies+2, nDimensions+nSpecies+2),           &
         deltaSpecificVolume(nDimensions+nSpecies+2),                                        &
         deltaVelocity(nDimensions, nDimensions+nSpecies+2),                                 &
         deltaPressure(nDimensions+nSpecies+2), deltaTemperature(nDimensions+nSpecies+2),    &
         deltaMassFraction(nSpecies, nDimensions+nSpecies+2),                                &
         deltaContravariantVelocity(nDimensions+nSpecies+2),                                 &
         deltaSpeedOfSound(nDimensions+nSpecies+2),                                          &
         deltaPhiSquared(nDimensions+nSpecies+2), deltaEnthalpy(nDimensions+nSpecies+2),     &
         deltaRightEigenvectors                                                              &
         (nDimensions+nSpecies+2, nDimensions+nSpecies+2,nDimensions+nSpecies+2),            &
         deltaEigenvalues(nDimensions+nSpecies+2, nDimensions+nSpecies+2),                   &
         deltaLeftEigenvectors                                                               &
         (nDimensions+nSpecies+2, nDimensions+nSpecies+2, nDimensions+nSpecies+2),           &
         temp(nDimensions+nSpecies+2)

    ! Normalize the metrics
    arcLength = sqrt(sum(metrics ** 2))
    normalizedMetrics = metrics / arcLength

    ! Compute specific volume if it was not specified
    if (present(specificVolume)) then
       specificVolume_ = specificVolume
    else
       specificVolume_ = 1.0_WP / conservedVariables(1)
    end if

    ! Compute velocity if it was not specified
    if (present(velocity)) then
       velocity_ = velocity
    else
       do i = 1, nDimensions
          velocity_(i) = specificVolume_ * conservedVariables(i + 1)
       end do
    end if

    ! Compute pressure if it was not specified
    if (present(pressure)) then
       pressure_ = pressure
    else
       pressure_ = (ratioOfSpecificHeats - 1.0_WP) *                                         &
            (conservedVariables(nDimensions+2) - 0.5_WP * conservedVariables(1) *            &
            sum(velocity_ ** 2))
    end if

    ! Compute mass fraction if it was not specified
    if (present(massFraction) .and. nSpecies .gt. 0) then
       massFraction_ = massFraction
    else
       do k = 1, nSpecies
          massFraction_(k) = specificVolume_ * conservedVariables(nDimensions+2+k)
       end do
    end if

    ! Temperature (this is not the actual temperature, assumes ideal gas law)
    temperature_ = ratioOfSpecificHeats * pressure_ * specificVolume_ /                      &
         (ratioOfSpecificHeats - 1.0_WP)

    ! Other dependent variables
    contravariantVelocity = sum(normalizedMetrics * velocity_)
    speedOfSound = sqrt(ratioOfSpecificHeats * pressure_ * specificVolume_)
    phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * sum(velocity_ ** 2)
    enthalpy = specificVolume_ * (conservedVariables(nDimensions+2) + pressure_)

    select case (nDimensions)

    case (1)

       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity + speedOfSound
       eigenvalues(3) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(3+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(3+k,3+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables_(2,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(3,:) - &
               velocity_(1) * deltaVelocity(1,:)) + phiSquared * deltaConservedVariables_(1,:)
          deltaTemperature = ratioOfSpecificHeats * deltaPressure * specificVolume_ /        &
               (ratioOfSpecificHeats - 1.0_WP) + ratioOfSpecificHeats * pressure_ *          &
               deltaSpecificVolume / (ratioOfSpecificHeats - 1.0_WP)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(3+k) +        &
                  specificVolume_ * deltaConservedVariables_(3+k,:)
          end do

          ! Compute variations of other dependent variables
          deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:)
          deltaSpeedOfSound = 0.5_WP * ratioOfSpecificHeats / speedOfSound *                 &
               (deltaPressure * specificVolume_ + deltaSpecificVolume * pressure_)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(3,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(3) + pressure_)

          ! Variation of matrix containing eigenvalues
          deltaEigenvalues(1,:) = deltaContravariantVelocity
          deltaEigenvalues(2,:) = deltaContravariantVelocity + deltaSpeedOfSound
          deltaEigenvalues(3,:) = deltaContravariantVelocity - deltaSpeedOfSound
          do k = 1, nSpecies
             deltaEigenvalues(3+k,:) = deltaContravariantVelocity
          end do
          deltaEigenvalues = arcLength * deltaEigenvalues

       end if

       ! Zero-out the eigenvalues corresponding to outgoing characteristics and
       ! corresponding variations
       do i = 1, 3 + nSpecies
          if (real(incomingDirection, WP) * real(eigenvalues(i), WP) .lt. 0.0_WP) then
             eigenvalues(i) = 0.0_WP
             if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_WP
          end if
       end do

       ! Matrix whose columns are the right eigenvectors:

       rightEigenvectors = 0.0_WP

       rightEigenvectors(1,1) = 1.0_WP
       rightEigenvectors(2,1) = velocity_(1)
       rightEigenvectors(3,1) = phiSquared / (ratioOfSpecificHeats - 1.0_WP)
       do k = 1, nSpecies
          rightEigenvectors(3+k,1) = massFraction_(k)
       end do

       rightEigenvectors(1,2) = 1.0_WP
       rightEigenvectors(2,2) = velocity_(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,2) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(3+k,2) = massFraction_(k)
       end do

       rightEigenvectors(1,3) = 1.0_WP
       rightEigenvectors(2,3) = velocity_(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,3) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(3+k,3) = massFraction_(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(3+k,3+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:

       leftEigenvectors = 0.0_WP

       leftEigenvectors(1,1) = 1.0_WP - phiSquared / speedOfSound ** 2
       leftEigenvectors(2,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(3,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(3+k,1) = -massFraction(k)
       end do

       leftEigenvectors(1,2) = velocity_(1) / temperature_
       leftEigenvectors(2,2) = - 0.5_WP * (velocity_(1) / temperature_ -                     &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(3,2) = - 0.5_WP * (velocity_(1) / temperature_ +                     &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(3+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = - 1.0_WP / temperature_
       leftEigenvectors(2,3) = 0.5_WP / temperature_
       leftEigenvectors(3,3) = 0.5_WP / temperature_
       do k = 1, nSpecies
          leftEigenvectors(3+k,3) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(3+k,3+k) = 1.0_WP
       end do

       ! ``Incoming'' part
       do j = 1, 3 + nSpecies
          do i = 1, 3 + nSpecies
             incomingJacobianOfInviscidFlux(i,j) =                                           &
                  rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +          &
                  rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +          &
                  rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j)
             do k = 1, nSpecies
                incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +  &
                     rightEigenvectors(i,3+k) * eigenvalues(3+k) * leftEigenvectors(3+k,j)
             end do
          end do
       end do

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! Variation of the matrix whose columns are the right eigenvectors:

          deltaRightEigenvectors = 0.0_WP

          deltaRightEigenvectors(1,1,:) = 0.0_WP
          deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
          deltaRightEigenvectors(3,1,:) = deltaPhiSquared /                                  &
               (ratioOfSpecificHeats - 1.0_WP)
          do k = 1, nSpecies
             deltaRightEigenvectors(3+k,1,:) = deltaMassFraction(k,:)
          end do

          deltaRightEigenvectors(1,2,:) = 0.0_WP
          deltaRightEigenvectors(2,2,:) = deltaVelocity(1,:) +                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,2,:) = deltaEnthalpy + deltaSpeedOfSound *                &
               contravariantVelocity + speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(3+k,2,:) = deltaMassFraction(k,:)
          end do

          deltaRightEigenvectors(1,3,:) = 0.0_WP
          deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) -                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,3,:) = deltaEnthalpy - deltaSpeedOfSound *                &
               contravariantVelocity - speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(3+k,3,:) = deltaMassFraction(k,:)
          end do

          do k = 1, nSpecies
             deltaRightEigenvectors(3+k,3+k,:) = 0.0_WP
          end do

          ! Variation of the matrix whose rows are the left eigenvectors:

          deltaLeftEigenvectors = 0.0_WP

          temp = deltaPhiSquared / speedOfSound ** 2 -                                       &
               2.0_WP * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
          deltaLeftEigenvectors(1,1,:) = -temp
          deltaLeftEigenvectors(2,1,:) = 0.5_WP * (temp -                                    &
               deltaContravariantVelocity / speedOfSound +                                   &
               contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(3,1,:) = 0.5_WP * (temp +                                    &
               deltaContravariantVelocity / speedOfSound -                                   &
               contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(3+k,1,:) = -deltaMassFraction(k,:)
          end do

          temp = deltaVelocity(1,:) / temperature_ -                                         &
               velocity_(1) / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,2,:) = temp
          deltaLeftEigenvectors(2,2,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(3,2,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(3+k,2,:) = 0.0_WP
          end do

          temp =  -1.0_WP / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,3,:) = -temp
          deltaLeftEigenvectors(2,3,:) = 0.5_WP * temp
          deltaLeftEigenvectors(3,3,:) = 0.5_WP * temp
          do k = 1, nSpecies
             deltaLeftEigenvectors(3+k,3,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaLeftEigenvectors(3+k,3+k,:) = 0.0_WP
          end do

          ! Variation of the ``incoming'' part
          do j = 1, 3 + nSpecies
             do i = 1, 3 + nSpecies
                deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                 &
                     deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +&
                     deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +&
                     deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +&
                     rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +&
                     rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +&
                     rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +&
                     rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +&
                     rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +&
                     rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:)
                do k = 1, nSpecies
                   deltaIncomingJacobianOfInviscidFlux(i,j,:) =                              &
                        deltaIncomingJacobianOfInviscidFlux(i,j,:) +                         &
                        deltaRightEigenvectors(i,3+k,:) * eigenvalues(3+k) *                 &
                        leftEigenvectors(3+k,j) + rightEigenvectors(i,3+k) *                 &
                        deltaEigenvalues(3+k,:) * leftEigenvectors(3+k,j) +                  &
                        rightEigenvectors(i,3+k) * eigenvalues(3+k) *                        &
                        deltaLeftEigenvectors(3+k,j,:)
                end do
             end do
          end do

       end if

    case (2)

       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity
       eigenvalues(3) = contravariantVelocity + speedOfSound
       eigenvalues(4) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(4+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             deltaConservedVariables_(4,4) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(4+k,4+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables_(2,:)
          deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                 &
               specificVolume_ * deltaConservedVariables_(3,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(4,:) - &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))) +    &
               phiSquared * deltaConservedVariables_(1,:)
          deltaTemperature = ratioOfSpecificHeats * deltaPressure * specificVolume_ /        &
               (ratioOfSpecificHeats - 1.0_WP) + ratioOfSpecificHeats * pressure_ *          &
               deltaSpecificVolume / (ratioOfSpecificHeats - 1.0_WP)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(4+k) +        &
                  specificVolume_ * deltaConservedVariables_(4+k,:)
          end do

          ! Compute variations of other dependent variables
          deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +           &
               normalizedMetrics(2) * deltaVelocity(2,:)
          deltaSpeedOfSound = 0.5_WP * ratioOfSpecificHeats / speedOfSound *                 &
               (deltaPressure * specificVolume_ + deltaSpecificVolume * pressure_)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(4,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(4) + pressure_)

          ! Variation of matrix containing eigenvalues
          deltaEigenvalues(1,:) = deltaContravariantVelocity
          deltaEigenvalues(2,:) = deltaContravariantVelocity
          deltaEigenvalues(3,:) = deltaContravariantVelocity + deltaSpeedOfSound
          deltaEigenvalues(4,:) = deltaContravariantVelocity - deltaSpeedOfSound
          do k = 1, nSpecies
             deltaEigenvalues(4+k,:) = deltaContravariantVelocity
          end do
          deltaEigenvalues = arcLength * deltaEigenvalues

       end if

       ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
       ! variations
       do i = 1, 4 + nSpecies
          if (real(incomingDirection, WP) * real(eigenvalues(i), WP) .lt. 0.0_WP) then
             eigenvalues(i) = 0.0_WP
             if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_WP
          end if
       end do

       ! Matrix whose columns are the right eigenvectors:

       rightEigenvectors = 0.0_WP

       rightEigenvectors(1,1) = 1.0_WP
       rightEigenvectors(2,1) = velocity_(1)
       rightEigenvectors(3,1) = velocity_(2)
       rightEigenvectors(4,1) = phiSquared / (ratioOfSpecificHeats - 1.0_WP)
       do k = 1, nSpecies
          rightEigenvectors(4+k,1) = massFraction_(k)
       end do

       rightEigenvectors(1,2) = 0.0_WP
       rightEigenvectors(2,2) = normalizedMetrics(2) * conservedVariables(1)
       rightEigenvectors(3,2) = - normalizedMetrics(1) * conservedVariables(1)
       rightEigenvectors(4,2) = conservedVariables(1) * (normalizedMetrics(2) *              &
            velocity_(1) - normalizedMetrics(1) * velocity_(2))
       do k = 1, nSpecies
          rightEigenvectors(4+k,2) = 0.0_WP
       end do

       rightEigenvectors(1,3) = 1.0_WP
       rightEigenvectors(2,3) = velocity_(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,3) = velocity_(2) + normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,3) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(4+k,3) = massFraction_(k)
       end do

       rightEigenvectors(1,4) = 1.0_WP
       rightEigenvectors(2,4) = velocity_(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,4) = velocity_(2) - normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,4) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(4+k,4) = massFraction_(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(4+k,4+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:

       leftEigenvectors = 0.0_WP

       leftEigenvectors(1,1) = 1.0_WP - phiSquared / speedOfSound ** 2
       leftEigenvectors(2,1) = - specificVolume_ * (normalizedMetrics(2) * velocity_(1) -    &
            normalizedMetrics(1) * velocity_(2))
       leftEigenvectors(3,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(4,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,1) = -massFraction_(k)
       end do

       leftEigenvectors(1,2) = velocity_(1) / temperature_
       leftEigenvectors(2,2) = specificVolume_ * normalizedMetrics(2)
       leftEigenvectors(3,2) = - 0.5_WP * (velocity_(1) / temperature_ -                     &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(4,2) = - 0.5_WP * (velocity_(1) / temperature_ +                     &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = velocity_(2) / temperature_
       leftEigenvectors(2,3) = - specificVolume_ * normalizedMetrics(1)
       leftEigenvectors(3,3) = - 0.5_WP * (velocity_(2) / temperature_ -                     &
            normalizedMetrics(2) / speedOfSound)
       leftEigenvectors(4,3) = - 0.5_WP * (velocity_(2) / temperature_ +                     &
            normalizedMetrics(2) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,3) = 0.0_WP
       end do

       leftEigenvectors(1,4) = - 1.0_WP / temperature_
       leftEigenvectors(2,4) = 0.0_WP
       leftEigenvectors(3,4) = 0.5_WP / temperature_
       leftEigenvectors(4,4) = 0.5_WP / temperature_
       do k = 1, nSpecies
          leftEigenvectors(4+k,4) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(4+k,4+k) = 1.0_WP
       end do

       ! ``Incoming'' part
       do j = 1, 4 + nSpecies
          do i = 1, 4 + nSpecies
             incomingJacobianOfInviscidFlux(i,j) =                                           &
                  rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +          &
                  rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +          &
                  rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +          &
                  rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j)
             do k = 1, nSpecies
                incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +  &
                     rightEigenvectors(i,4+k) * eigenvalues(4+k) * leftEigenvectors(4+k,j)
             end do
          end do
       end do

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! Variation of the matrix whose columns are the right eigenvectors:

          deltaRightEigenvectors = 0.0_WP

          deltaRightEigenvectors(1,1,:) = 0.0_WP
          deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
          deltaRightEigenvectors(3,1,:) = deltaVelocity(2,:)
          deltaRightEigenvectors(4,1,:) = deltaPhiSquared /                                  &
               (ratioOfSpecificHeats - 1.0_WP)
          do k = 1, nSpecies
             deltaRightEigenvectors(4+k,1,:) = deltaMassFraction(k,:)
          end do

          deltaRightEigenvectors(1,2,:) = 0.0_WP
          deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) *                             &
               deltaConservedVariables_(1,:)
          deltaRightEigenvectors(3,2,:) = -normalizedMetrics(1) *                            &
               deltaConservedVariables_(1,:)
          deltaRightEigenvectors(4,2,:) = deltaConservedVariables_(1,:) *                    &
               (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) + &
               conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -          &
               normalizedMetrics(1) * deltaVelocity(2,:))
          do k = 1, nSpecies
             deltaRightEigenvectors(4+k,2,:) = 0.0_WP
          end do

          deltaRightEigenvectors(1,3,:) = 0.0_WP
          deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) +                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,3,:) = deltaVelocity(2,:) +                               &
               normalizedMetrics(2) * deltaSpeedOfSound
          deltaRightEigenvectors(4,3,:) = deltaEnthalpy + deltaSpeedOfSound *                &
               contravariantVelocity + speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(4+k,3,:) = deltaMassFraction(k,:)
          end do

          deltaRightEigenvectors(1,4,:) = 0.0_WP
          deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) -                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) -                               &
               normalizedMetrics(2) * deltaSpeedOfSound
          deltaRightEigenvectors(4,4,:) = deltaEnthalpy - deltaSpeedOfSound *                &
               contravariantVelocity - speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(4+k,4,:) = deltaMassFraction(k,:)
          end do

          do k = 1, nSpecies
             deltaRightEigenvectors(4+k,4+k,:) = 0.0_WP
          end do

          ! Variation of the matrix whose rows are the left eigenvectors:

          deltaLeftEigenvectors = 0.0_WP

          temp = deltaPhiSquared / speedOfSound ** 2 -                                       &
               2.0_WP * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
          deltaLeftEigenvectors(1,1,:) = -temp
          deltaLeftEigenvectors(2,1,:) = -deltaSpecificVolume *                              &
               (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) - &
               specificVolume_ * (normalizedMetrics(2) * deltaVelocity(1,:) -                &
               normalizedMetrics(1) * deltaVelocity(2,:))
          deltaLeftEigenvectors(3,1,:) = 0.5_WP * (temp -                                    &
               deltaContravariantVelocity / speedOfSound +                                   &
               contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(4,1,:) = 0.5_WP * (temp +                                    &
               deltaContravariantVelocity / speedOfSound -                                   &
               contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(4+k,1,:) = -deltaMassFraction(k,:)
          end do

          temp = deltaVelocity(1,:) / temperature_ -                                         &
               velocity_(1) / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,2,:) = temp
          deltaLeftEigenvectors(2,2,:) = deltaSpecificVolume * normalizedMetrics(2)
          deltaLeftEigenvectors(3,2,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(4,2,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)

          temp = deltaVelocity(2,:) / temperature_ -                                         &
               velocity_(2) / temperature_ ** 2 * deltaTemperature
          do k = 1, nSpecies
             deltaLeftEigenvectors(4+k,2,:) = 0.0_WP
          end do

          deltaLeftEigenvectors(1,3,:) = temp
          deltaLeftEigenvectors(2,3,:) = -deltaSpecificVolume * normalizedMetrics(1)
          deltaLeftEigenvectors(3,3,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(4,3,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(4+k,3,:) = 0.0_WP
          end do

          temp =  -1.0_WP / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,4,:) = -temp
          deltaLeftEigenvectors(2,4,:) = 0.0_WP
          deltaLeftEigenvectors(3,4,:) = 0.5_WP * temp
          deltaLeftEigenvectors(4,4,:) = 0.5_WP * temp
          do k = 1, nSpecies
             deltaLeftEigenvectors(4+k,4,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaLeftEigenvectors(4+k,4+k,:) = 0.0_WP
          end do

          ! Variation of the ``incoming'' part
          do j = 1, 4 + nSpecies
             do i = 1, 4 + nSpecies
                deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                 &
                     deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +&
                     deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +&
                     deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +&
                     deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +&
                     rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +&
                     rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +&
                     rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +&
                     rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +&
                     rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +&
                     rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +&
                     rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +&
                     rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:)
                do k = 1, nSpecies
                   deltaIncomingJacobianOfInviscidFlux(i,j,:) =                              &
                        deltaIncomingJacobianOfInviscidFlux(i,j,:) +                         &
                        deltaRightEigenvectors(i,4+k,:) * eigenvalues(4+k) *                 &
                        leftEigenvectors(4+k,j) + rightEigenvectors(i,4+k) *                 &
                        deltaEigenvalues(4+k,:) * leftEigenvectors(4+k,j) +                  &
                        rightEigenvectors(i,4+k) * eigenvalues(4+k) *                        &
                        deltaLeftEigenvectors(4+k,j,:)
                end do
             end do
          end do

       end if

    case (3)

       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity
       eigenvalues(3) = contravariantVelocity
       eigenvalues(4) = contravariantVelocity + speedOfSound
       eigenvalues(5) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(5+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! If not specified, use identity matrix for the variation of conservedVariables
          if (present(deltaConservedVariables)) then
             deltaConservedVariables_ = deltaConservedVariables
          else
             deltaConservedVariables_ = 0.0_WP
             deltaConservedVariables_(1,1) = 1.0_WP
             deltaConservedVariables_(2,2) = 1.0_WP
             deltaConservedVariables_(3,3) = 1.0_WP
             deltaConservedVariables_(4,4) = 1.0_WP
             deltaConservedVariables_(5,5) = 1.0_WP
             do k = 1, nSpecies
                deltaConservedVariables_(5+k,5+k) = 1.0_WP
             end do
          end if

          ! Compute variations of specific volume, velocity, pressure and mass fraction
          deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
          deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                 &
               specificVolume_ * deltaConservedVariables_(2,:)
          deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                 &
               specificVolume_ * deltaConservedVariables_(3,:)
          deltaVelocity(3,:) = deltaSpecificVolume * conservedVariables(4) +                 &
               specificVolume_ * deltaConservedVariables_(4,:)
          deltaPressure = (ratioOfSpecificHeats - 1.0_WP) * (deltaConservedVariables_(5,:) - &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +      &
               velocity_(3) * deltaVelocity(3,:))) + phiSquared * deltaConservedVariables_(1,:)
          deltaTemperature = ratioOfSpecificHeats * deltaPressure * specificVolume_ /        &
               (ratioOfSpecificHeats - 1.0_WP) + ratioOfSpecificHeats * pressure_ *          &
               deltaSpecificVolume / (ratioOfSpecificHeats - 1.0_WP)
          do k = 1, nSpecies
             deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(5+k) +        &
                  specificVolume_ * deltaConservedVariables_(5+k,:)
          end do

          ! Compute variations of other dependent variables
          deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +           &
               normalizedMetrics(2) * deltaVelocity(2,:) +                                   &
               normalizedMetrics(3) * deltaVelocity(3,:)
          deltaSpeedOfSound = 0.5_WP * ratioOfSpecificHeats / speedOfSound *                 &
               (deltaPressure * specificVolume_ + deltaSpecificVolume * pressure_)
          deltaPhiSquared = (ratioOfSpecificHeats - 1.0_WP) *                                &
               (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +      &
               velocity_(3) * deltaVelocity(3,:))
          deltaEnthalpy = specificVolume_ * (deltaConservedVariables_(5,:) + deltaPressure) +&
               deltaSpecificVolume * (conservedVariables(5) + pressure_)

          ! Variation of matrix containing eigenvalues
          deltaEigenvalues(1,:) = deltaContravariantVelocity
          deltaEigenvalues(2,:) = deltaContravariantVelocity
          deltaEigenvalues(3,:) = deltaContravariantVelocity
          deltaEigenvalues(4,:) = deltaContravariantVelocity + deltaSpeedOfSound
          deltaEigenvalues(5,:) = deltaContravariantVelocity - deltaSpeedOfSound
          do k = 1, nSpecies
             deltaEigenvalues(5+k,:) = deltaContravariantVelocity
          end do
          deltaEigenvalues = arcLength * deltaEigenvalues

       end if

       ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
       ! variations
       do i = 1, 5 + nSpecies
          if (real(incomingDirection, WP) * real(eigenvalues(i), WP) .lt. 0.0_WP) then
             eigenvalues(i) = 0.0_WP
             if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_WP
          end if
       end do

       ! Matrix whose columns are the right eigenvectors:

       rightEigenvectors = 0.0_WP

       rightEigenvectors(1,1) = normalizedMetrics(1)
       rightEigenvectors(2,1) = normalizedMetrics(1) * velocity_(1)
       rightEigenvectors(3,1) = normalizedMetrics(1) * velocity_(2) +                        &
            conservedVariables(1) * normalizedMetrics(3)
       rightEigenvectors(4,1) = normalizedMetrics(1) * velocity_(3) -                        &
            conservedVariables(1) * normalizedMetrics(2)
       rightEigenvectors(5,1) = conservedVariables(1) * (normalizedMetrics(3) *              &
            velocity_(2) - normalizedMetrics(2) * velocity_(3)) + phiSquared /               &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(1)
       do k = 1, nSpecies
          rightEigenvectors(5+k,1) = massFraction_(k) * normalizedMetrics(1)
       end do

       rightEigenvectors(1,2) = normalizedMetrics(2)
       rightEigenvectors(2,2) = normalizedMetrics(2) * velocity_(1) -                        &
            conservedVariables(1) * normalizedMetrics(3)
       rightEigenvectors(3,2) = normalizedMetrics(2) * velocity_(2)
       rightEigenvectors(4,2) = normalizedMetrics(2) * velocity_(3) +                        &
            conservedVariables(1) * normalizedMetrics(1)
       rightEigenvectors(5,2) = conservedVariables(1) * (normalizedMetrics(1) *              &
            velocity_(3) - normalizedMetrics(3) * velocity_(1)) + phiSquared /               &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(2)
       do k = 1, nSpecies
          rightEigenvectors(5+k,2) = massFraction_(k) * normalizedMetrics(2)
       end do

       rightEigenvectors(1,3) = normalizedMetrics(3)
       rightEigenvectors(2,3) = normalizedMetrics(3) * velocity_(1) +                        &
            conservedVariables(1) * normalizedMetrics(2)
       rightEigenvectors(3,3) = normalizedMetrics(3) * velocity_(2) -                        &
            conservedVariables(1) * normalizedMetrics(1)
       rightEigenvectors(4,3) = normalizedMetrics(3) * velocity_(3)
       rightEigenvectors(5,3) = conservedVariables(1) * (normalizedMetrics(2) *              &
            velocity_(1) - normalizedMetrics(1) * velocity_(2)) + phiSquared /               &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(3)
       do k = 1, nSpecies
          rightEigenvectors(5+k,3) = massFraction_(k) * normalizedMetrics(3)
       end do

       rightEigenvectors(1,4) = 1.0_WP
       rightEigenvectors(2,4) = velocity_(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,4) = velocity_(2) + normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,4) = velocity_(3) + normalizedMetrics(3) * speedOfSound
       rightEigenvectors(5,4) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(5+k,4) = massFraction_(k)
       end do

       rightEigenvectors(1,5) = 1.0_WP
       rightEigenvectors(2,5) = velocity_(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,5) = velocity_(2) - normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,5) = velocity_(3) - normalizedMetrics(3) * speedOfSound
       rightEigenvectors(5,5) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(5+k,5) = massFraction_(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(5+k,5+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:

       leftEigenvectors = 0.0_WP

       leftEigenvectors(1,1) = normalizedMetrics(1) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume_ *                    &
            (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3))
       leftEigenvectors(2,1) = normalizedMetrics(2) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume_ *                    &
            (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1))
       leftEigenvectors(3,1) = normalizedMetrics(3) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume_ *                    &
            (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2))
       leftEigenvectors(4,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(5,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,1) = -massFraction_(k)
       end do

       leftEigenvectors(1,2) = normalizedMetrics(1) * velocity_(1) / temperature_
       leftEigenvectors(2,2) = normalizedMetrics(2) * velocity_(1) / temperature_ -          &
            specificVolume_ * normalizedMetrics(3)
       leftEigenvectors(3,2) = normalizedMetrics(3) * velocity_(1) / temperature_ +          &
            specificVolume_ * normalizedMetrics(2)
       leftEigenvectors(4,2) = - 0.5_WP * (velocity_(1) / temperature_ -                     &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(5,2) = - 0.5_WP * (velocity_(1) / temperature_ +                     &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = normalizedMetrics(1) * velocity_(2) / temperature_ +          &
            specificVolume_ * normalizedMetrics(3)
       leftEigenvectors(2,3) = normalizedMetrics(2) * velocity_(2) / temperature_
       leftEigenvectors(3,3) = normalizedMetrics(3) * velocity_(2) / temperature_ -          &
            specificVolume_ * normalizedMetrics(1)
       leftEigenvectors(4,3) = - 0.5_WP * (velocity_(2) / temperature_ -                     &
            normalizedMetrics(2) / speedOfSound)
       leftEigenvectors(5,3) = - 0.5_WP * (velocity_(2) / temperature_ +                     &
            normalizedMetrics(2) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,3) = 0.0_WP
       end do

       leftEigenvectors(1,4) = normalizedMetrics(1) * velocity_(3) / temperature_ -          &
            specificVolume_ * normalizedMetrics(2)
       leftEigenvectors(2,4) = normalizedMetrics(2) * velocity_(3) / temperature_ +          &
            specificVolume_ * normalizedMetrics(1)
       leftEigenvectors(3,4) = normalizedMetrics(3) * velocity_(3) / temperature_
       leftEigenvectors(4,4) = - 0.5_WP * (velocity_(3) / temperature_ -                     &
            normalizedMetrics(3) / speedOfSound)
       leftEigenvectors(5,4) = - 0.5_WP * (velocity_(3) / temperature_ +                     &
            normalizedMetrics(3) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,4) = 0.0_WP
       end do

       leftEigenvectors(1,5) = - normalizedMetrics(1) / temperature_
       leftEigenvectors(2,5) = - normalizedMetrics(2) / temperature_
       leftEigenvectors(3,5) = - normalizedMetrics(3) / temperature_
       leftEigenvectors(4,5) = 0.5_WP / temperature_
       leftEigenvectors(5,5) = 0.5_WP / temperature_
       do k = 1, nSpecies
          leftEigenvectors(5+k,5) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(5+k,5+k) = 1.0_WP
       end do

       ! ``Incoming'' part
       do j = 1, 5 + nSpecies
          do i = 1, 5 + nSpecies
             incomingJacobianOfInviscidFlux(i,j) =                                           &
                  rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +          &
                  rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +          &
                  rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +          &
                  rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j) +          &
                  rightEigenvectors(i,5) * eigenvalues(5) * leftEigenvectors(5,j)
             do k = 1, nSpecies
                incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +  &
                     rightEigenvectors(i,5+k) * eigenvalues(5+k) * leftEigenvectors(5+k,j)
             end do
          end do
       end do

       if (present(deltaIncomingJacobianOfInviscidFlux)) then

          ! Variation of the matrix whose columns are the right eigenvectors:

          deltaRightEigenvectors = 0.0_WP

          deltaRightEigenvectors(1,1,:) = 0.0_WP
          deltaRightEigenvectors(2,1,:) = normalizedMetrics(1) * deltaVelocity(1,:)
          deltaRightEigenvectors(3,1,:) = normalizedMetrics(1) * deltaVelocity(2,:) +        &
               deltaConservedVariables_(1,:) * normalizedMetrics(3)
          deltaRightEigenvectors(4,1,:) = normalizedMetrics(1) * deltaVelocity(3,:) -        &
               deltaConservedVariables_(1,:) * normalizedMetrics(2)
          deltaRightEigenvectors(5,1,:) = deltaConservedVariables_(1,:) *                    &
               (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3)) + &
               conservedVariables(1) * (normalizedMetrics(3) * deltaVelocity(2,:) -          &
               normalizedMetrics(2) * deltaVelocity(3,:)) + deltaPhiSquared /                &
               (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(1)
          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,1,:) = 0.0_WP
          end do

          deltaRightEigenvectors(1,2,:) = 0.0_WP
          deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) * deltaVelocity(1,:) -        &
               deltaConservedVariables_(1,:) * normalizedMetrics(3)
          deltaRightEigenvectors(3,2,:) = normalizedMetrics(2) * deltaVelocity(2,:)
          deltaRightEigenvectors(4,2,:) = normalizedMetrics(2) * deltaVelocity(3,:) +        &
               deltaConservedVariables_(1,:) * normalizedMetrics(1)
          deltaRightEigenvectors(5,2,:) = deltaConservedVariables_(1,:) *                    &
               (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1)) + &
               conservedVariables(1) * (normalizedMetrics(1) * deltaVelocity(3,:) -          &
               normalizedMetrics(3) * deltaVelocity(1,:)) + deltaPhiSquared /                &
               (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(2)
          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,2,:) = 0.0_WP
          end do

          deltaRightEigenvectors(1,3,:) = 0.0_WP
          deltaRightEigenvectors(2,3,:) = normalizedMetrics(3) * deltaVelocity(1,:) +        &
               deltaConservedVariables_(1,:) * normalizedMetrics(2)
          deltaRightEigenvectors(3,3,:) = normalizedMetrics(3) * deltaVelocity(2,:) -        &
               deltaConservedVariables_(1,:) * normalizedMetrics(1)
          deltaRightEigenvectors(4,3,:) = normalizedMetrics(3) * deltaVelocity(3,:)
          deltaRightEigenvectors(5,3,:) = deltaConservedVariables_(1,:) *                    &
               (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) + &
               conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -          &
               normalizedMetrics(1) * deltaVelocity(2,:)) + deltaPhiSquared /                &
               (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(3)
          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,3,:) = 0.0_WP
          end do

          deltaRightEigenvectors(1,4,:) = 0.0_WP
          deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) +                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) +                               &
               normalizedMetrics(2) * deltaSpeedOfSound
          deltaRightEigenvectors(4,4,:) = deltaVelocity(3,:) +                               &
               normalizedMetrics(3) * deltaSpeedOfSound
          deltaRightEigenvectors(5,4,:) = deltaEnthalpy + deltaSpeedOfSound *                &
               contravariantVelocity + speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,4,:) = deltaMassFraction(k,:)
          end do

          deltaRightEigenvectors(1,5,:) = 0.0_WP
          deltaRightEigenvectors(2,5,:) = deltaVelocity(1,:) -                               &
               normalizedMetrics(1) * deltaSpeedOfSound
          deltaRightEigenvectors(3,5,:) = deltaVelocity(2,:) -                               &
               normalizedMetrics(2) * deltaSpeedOfSound
          deltaRightEigenvectors(4,5,:) = deltaVelocity(3,:) -                               &
               normalizedMetrics(3) * deltaSpeedOfSound
          deltaRightEigenvectors(5,5,:) = deltaEnthalpy - deltaSpeedOfSound *                &
               contravariantVelocity - speedOfSound * deltaContravariantVelocity
          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,5,:) = deltaMassFraction(k,:)
          end do

          do k = 1, nSpecies
             deltaRightEigenvectors(5+k,5+k,:) = 0.0_WP
          end do

          ! Variation of the matrix whose rows are the left eigenvectors:

          deltaLeftEigenvectors = 0.0_WP

          temp = deltaPhiSquared / speedOfSound ** 2 -                                       &
               2.0_WP * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
          deltaLeftEigenvectors(1,1,:) = -normalizedMetrics(1) * temp - deltaSpecificVolume *&
               (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3)) - &
               specificVolume_ * (normalizedMetrics(3) * deltaVelocity(2,:) -                &
               normalizedMetrics(2) * deltaVelocity(3,:))
          deltaLeftEigenvectors(2,1,:) = -normalizedMetrics(2) * temp - deltaSpecificVolume *&
               (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1)) - &
               specificVolume_ * (normalizedMetrics(1) * deltaVelocity(3,:) -                &
               normalizedMetrics(3) * deltaVelocity(1,:))
          deltaLeftEigenvectors(3,1,:) = -normalizedMetrics(3) * temp - deltaSpecificVolume *&
               (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) - &
               specificVolume_ * (normalizedMetrics(2) * deltaVelocity(1,:) -                &
               normalizedMetrics(1) * deltaVelocity(2,:))
          deltaLeftEigenvectors(4,1,:) = 0.5_WP * (temp - deltaContravariantVelocity /       &
               speedOfSound + contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(5,1,:) = 0.5_WP * (temp + deltaContravariantVelocity /       &
               speedOfSound - contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,1,:) = -deltaMassFraction(k,:)
          end do

          temp = deltaVelocity(1,:) / temperature_ -                                         &
               velocity_(1) / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,2,:) = normalizedMetrics(1) * temp
          deltaLeftEigenvectors(2,2,:) = normalizedMetrics(2) * temp -                       &
               deltaSpecificVolume * normalizedMetrics(3)
          deltaLeftEigenvectors(3,2,:) = normalizedMetrics(3) * temp +                       &
               deltaSpecificVolume * normalizedMetrics(2)
          deltaLeftEigenvectors(4,2,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(5,2,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,2,:) = 0.0_WP
          end do

          temp = deltaVelocity(2,:) / temperature_ -                                         &
               velocity_(2) / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,3,:) = normalizedMetrics(1) * temp +                       &
               deltaSpecificVolume * normalizedMetrics(3)
          deltaLeftEigenvectors(2,3,:) = normalizedMetrics(2) * temp
          deltaLeftEigenvectors(3,3,:) = normalizedMetrics(3) * temp -                       &
               deltaSpecificVolume * normalizedMetrics(1)
          deltaLeftEigenvectors(4,3,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(5,3,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,3,:) = 0.0_WP
          end do

          temp = deltaVelocity(3,:) / temperature_ -                                         &
               velocity_(3) / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,4,:) = normalizedMetrics(1) * temp -                       &
               deltaSpecificVolume * normalizedMetrics(2)
          deltaLeftEigenvectors(2,4,:) = normalizedMetrics(2) * temp +                       &
               deltaSpecificVolume * normalizedMetrics(1)
          deltaLeftEigenvectors(3,4,:) = normalizedMetrics(3) * temp
          deltaLeftEigenvectors(4,4,:) = -0.5_WP * (temp +                                   &
               normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)
          deltaLeftEigenvectors(5,4,:) = -0.5_WP * (temp -                                   &
               normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)
          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,4,:) = 0.0_WP
          end do

          temp =  -1.0_WP / temperature_ ** 2 * deltaTemperature
          deltaLeftEigenvectors(1,5,:) = -normalizedMetrics(1) * temp
          deltaLeftEigenvectors(2,5,:) = -normalizedMetrics(2) * temp
          deltaLeftEigenvectors(3,5,:) = -normalizedMetrics(3) * temp
          deltaLeftEigenvectors(4,5,:) = 0.5_WP * temp
          deltaLeftEigenvectors(5,5,:) = 0.5_WP * temp
          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,5,:) = 0.0_WP
          end do

          do k = 1, nSpecies
             deltaLeftEigenvectors(5+k,5+k,:) = 0.0_WP
          end do

          ! Variation of the ``incoming'' part
          do j = 1, 5 + nSpecies
             do i = 1, 5 + nSpecies
                deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                 &
                     deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +&
                     deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +&
                     deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +&
                     deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +&
                     deltaRightEigenvectors(i,5,:) * eigenvalues(5) * leftEigenvectors(5,j) +&
                     rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +&
                     rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +&
                     rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +&
                     rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +&
                     rightEigenvectors(i,5) * deltaEigenvalues(5,:) * leftEigenvectors(5,j) +&
                     rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +&
                     rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +&
                     rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +&
                     rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:) +&
                     rightEigenvectors(i,5) * eigenvalues(5) * deltaLeftEigenvectors(5,j,:)
                do k = 1, nSpecies
                   deltaIncomingJacobianOfInviscidFlux(i,j,:) =                              &
                        deltaIncomingJacobianOfInviscidFlux(i,j,:) +                         &
                        deltaRightEigenvectors(i,5+k,:) * eigenvalues(5+k) *                 &
                        leftEigenvectors(5+k,j) + rightEigenvectors(i,5+k) *                 &
                        deltaEigenvalues(5+k,:) * leftEigenvectors(5+k,j) +                  &
                        rightEigenvectors(i,5+k) * eigenvalues(5+k) *                        &
                        deltaLeftEigenvectors(5+k,j,:)
                end do
             end do
          end do

       end if

    end select

    return
  end subroutine compute_incoming_jacobian_of_inviscid_flux

  pure subroutine compute_first_partial_viscous_jacobian(conservedVariables, metrics,        &
       stressTensor, heatFlux, enthalpyFlux, speciesFlux, firstPartialViscousJacobian,       &
       specificVolume, velocity, temperature, massFraction)

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:), metrics(:), stressTensor(:), heatFlux(:)
    real(WP), intent(out) :: firstPartialViscousJacobian(:,:)
    real(WP), intent(in), optional :: specificVolume, velocity(:), temperature,              &
         massFraction(:), enthalpyFlux(:), speciesFlux(:,:)

    ! Local variables
    integer :: i, k
    real(WP) :: specificVolume_, velocity_(nDimensions), temperature_, phiSquared,           &
         contravariantStressTensor(nDimensions), contravariantHeatFlux, temp,                &
         massFraction_(nSpecies), contravariantSpeciesFlux(nSpecies),                        &
         contravariantEnthalpyFlux, deltaViscosity(nDimensions+nSpecies+2),                  &
         deltaTemperature(nDimensions+nSpecies+2)

    ! Compute specific volume if it was not specified
    if (present(specificVolume)) then
       specificVolume_ = specificVolume
    else
       specificVolume_ = 1.0_WP / conservedVariables(1)
    end if

    ! Compute velocity if it was not specified
    if (present(velocity)) then
       velocity_ = velocity
    else
       do i = 1, nDimensions
          velocity_(i) = specificVolume_ * conservedVariables(i + 1)
       end do
    end if

    ! Compute mass fraction if it was not specified
    if (present(massFraction) .and. nSpecies .gt. 0) then
       massFraction_ = massFraction
    else
       do k = 1, nSpecies
          massFraction_(k) = specificVolume_ * conservedVariables(nDimensions+2+k)
       end do
    end if

    ! Compute temperature if it was not specified
    if (present(temperature)) then
       temperature_ = temperature
    else
       if (equationOfState .eq. IDEAL_GAS) then ! ... single-component gas
          temperature_ = ratioOfSpecificHeats * (specificVolume_ *                           &
               conservedVariables(nDimensions+2) - 0.5_WP * sum(velocity_ ** 2))
       else ! ... ideal gas mixture
          temp = conservedVariables(1) * molecularWeightInverse(nSpecies+1)
          do k = 1, nSpecies
             temp = temp + conservedVariables(nDimensions+2+k) * (molecularWeightInverse(k) -&
                  molecularWeightInverse(nSpecies+1))
          end do
          temperature_ = ratioOfSpecificHeats * (conservedVariables(nDimensions+2) - 0.5_WP *&
               sum(conservedVariables(2:nDimensions+1)**2) / conservedVariables(1)) / temp
       end if
    end if

    ! Zero-out first partial viscous Jacobian
    firstPartialViscousJacobian = 0.0_WP

    select case (nDimensions)

    case (1)

       if (equationOfState .eq. IDEAL_GAS) then ! ... single-component gas

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * velocity_(1) ** 2
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) !... not normalized
          end do
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature_ / &
               ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(1)
          deltaViscosity(3) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_
          do k = 1, nSpecies
             deltaViscosity(3+k) = 0.0_WP
          end do
          temp = velocity(1) * contravariantStressTensor(1) - contravariantHeatFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * temp - specificVolume_ *    &
               (velocity(1) * contravariantStressTensor(1))
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * temp +                      &
               specificVolume_ * contravariantStressTensor(1)
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * temp
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

       else ! ... ideal gas mixture

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * velocity_(1) ** 2
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) !... not normalized
          end do
          contravariantEnthalpyFlux = metrics(1) * enthalpyFlux(1) !... not normalized
          temp = molecularWeightInverse(nSpecies+1)
          do k = 1, nSpecies
             temp = temp + massFraction_(k) * (molecularWeightInverse(k) -                   &
                  molecularWeightInverse(nSpecies+1))
          end do
          temp = 1.0_WP / temp
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaTemperature(1) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(1) * temp
          deltaTemperature(2) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(1) * temp
          deltaViscosity(3) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp
          deltaTemperature(3) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp
          do k = 1, nSpecies
             deltaViscosity(3+k) = - powerLawExponent * specificVolume_ * temp *             &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
             deltaTemperature(3+k) = - specificVolume_ * temp *                              &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
          end do
          temp = velocity(1) * contravariantStressTensor(1) - contravariantHeatFlux -        &
               contravariantEnthalpyFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * temp -                      &
               specificVolume_ * velocity(1) * contravariantStressTensor(1) -                &
               deltaTemperature(1) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * temp +                      &
               specificVolume_ * contravariantStressTensor(1) -                              &
               deltaTemperature(2) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * temp -                      &
               deltaTemperature(3) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(3+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

          do k = 1, nSpecies
             firstPartialViscousJacobian(1,3+K) = 0.0_WP
             firstPartialViscousJacobian(2,3+k) = deltaViscosity(3+k) *                      &
                  contravariantStressTensor(1)
             firstPartialViscousJacobian(3,3+k) = deltaViscosity(3+k) * temp -               &
               deltaTemperature(3+k) * contravariantEnthalpyFlux
             do i = 1, nSpecies
                firstPartialViscousJacobian(3+i,3+k) = - deltaViscosity(3+k) *               &
                     contravariantSpeciesFlux(i)
             end do
          end do

       end if

    case (2)

       if (equationOfState .eq. IDEAL_GAS) then ! ... single-component gas

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                            &
               (velocity_(1) ** 2 + velocity_(2) ** 2)
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                      &
               metrics(2) * stressTensor(2) !... not normalized
          contravariantStressTensor(2) = metrics(1) * stressTensor(3) +                      &
               metrics(2) * stressTensor(4) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) +                                 &
               metrics(2) * heatFlux(2) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                   &
                  metrics(2) * speciesFlux(k,2) !... not normalized
          end do
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature_ / &
               ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats *                    &
               specificVolume_ / temperature_ * velocity(1)
          deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats *                    &
               specificVolume_ / temperature_ * velocity(2)
          deltaViscosity(4) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_
          do k = 1, nSpecies
             deltaViscosity(4+k) = 0.0_WP
          end do
          temp = velocity(1) * contravariantStressTensor(1) +                                &
               velocity(2) * contravariantStressTensor(2) - contravariantHeatFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,1) = deltaViscosity(1) * temp - specificVolume_ *    &
               (velocity(1) * contravariantStressTensor(1) +                                 &
               velocity(2) * contravariantStressTensor(2))
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,2) = deltaViscosity(2) * temp +                      &
               specificVolume_ * contravariantStressTensor(1)
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,3) = deltaViscosity(3) * temp +                      &
               specificVolume_ * contravariantStressTensor(2)
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,4) = 0.0_WP
          firstPartialViscousJacobian(2,4) = deltaViscosity(4) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,4) = deltaViscosity(4) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,4) = deltaViscosity(4) * temp
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,4) = - deltaViscosity(4) *                      &
                  contravariantSpeciesFlux(k)
          end do

       else ! ... ideal gas mixture

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                            &
               (velocity_(1) ** 2 + velocity_(2) ** 2)
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                      &
               metrics(2) * stressTensor(2) !... not normalized
          contravariantStressTensor(2) = metrics(1) * stressTensor(3) +                      &
               metrics(2) * stressTensor(4) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) +                                 &
               metrics(2) * heatFlux(2) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                   &
                  metrics(2) * speciesFlux(k,2) !... not normalized
          end do
          contravariantEnthalpyFlux = metrics(1) * enthalpyFlux(1) +                         &
               metrics(2) * enthalpyFlux(2) !... not normalized
          temp = molecularWeightInverse(nSpecies+1)
          do k = 1, nSpecies
             temp = temp + massFraction_(k) * (molecularWeightInverse(k) -                   &
                  molecularWeightInverse(nSpecies+1))
          end do
          temp = 1.0_WP / temp
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaTemperature(1) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(1) * temp
          deltaTemperature(2) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(1) * temp
          deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(2) * temp
          deltaTemperature(3) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(2) * temp
          deltaViscosity(4) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp
          deltaTemperature(4) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp
          do k = 1, nSpecies
             deltaViscosity(4+k) = - powerLawExponent * specificVolume_ * temp *             &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
             deltaTemperature(4+k) = - specificVolume_ * temp *                              &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
          end do
          temp = velocity(1) * contravariantStressTensor(1) + velocity(2) *                  &
               contravariantStressTensor(2) - contravariantHeatFlux -                        &
               contravariantEnthalpyFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,1) = deltaViscosity(1) * temp - specificVolume_ *    &
               (velocity(1) * contravariantStressTensor(1) + velocity(2) *                   &
               contravariantStressTensor(2)) - deltaTemperature(1) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,2) = deltaViscosity(2) * temp + specificVolume_ *    &
               contravariantStressTensor(1) - deltaTemperature(2) * contravariantEnthalpyFlux

          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,3) = deltaViscosity(3) * temp + specificVolume_ *    &
               contravariantStressTensor(2) - deltaTemperature(3) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,4) = 0.0_WP
          firstPartialViscousJacobian(2,4) = deltaViscosity(4) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,4) = deltaViscosity(4) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,4) = deltaViscosity(4) * temp -                      &
               deltaTemperature(4) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(4+k,4) = - deltaViscosity(4) *                      &
                  contravariantSpeciesFlux(k)
          end do

          do k = 1, nSpecies
             firstPartialViscousJacobian(1,4+k) = 0.0_WP
             firstPartialViscousJacobian(2,4+k) = deltaViscosity(4+k) *                      &
                  contravariantStressTensor(1)
             firstPartialViscousJacobian(3,4+k) = deltaViscosity(4+k) *                      &
                  contravariantStressTensor(2)
             firstPartialViscousJacobian(4,4+k) = deltaViscosity(4+k) * temp -               &
                  deltaTemperature(4+k) * contravariantEnthalpyFlux
             do i = 1, nSpecies
                firstPartialViscousJacobian(4+i,4+k) = - deltaViscosity(4+k) *               &
                     contravariantSpeciesFlux(i)
             end do
          end do

       end if

    case (3)

       if (equationOfState .eq. IDEAL_GAS) then ! ... single-component gas

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                            &
               (velocity_(1) ** 2 + velocity_(2) ** 2 + velocity_(3) ** 2)
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                      &
               metrics(2) * stressTensor(2) + metrics(3) * stressTensor(3) !... not normalized
          contravariantStressTensor(2) = metrics(1) * stressTensor(4) +                      &
               metrics(2) * stressTensor(5) + metrics(3) * stressTensor(6) !... not normalized
          contravariantStressTensor(3) = metrics(1) * stressTensor(7) +                      &
               metrics(2) * stressTensor(8) + metrics(3) * stressTensor(9) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) + metrics(2) * heatFlux(2) +      &
               metrics(3) * heatFlux(3) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                   &
                  metrics(2) * speciesFlux(k,2) +                                            &
                  metrics(3) * speciesFlux(k,3) !... not normalized
          end do
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) - temperature_ / &
               ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats *                    &
               specificVolume_ / temperature_ * velocity(1)
          deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats *                    &
               specificVolume_ / temperature_ * velocity(2)
          deltaViscosity(4) = - powerLawExponent * ratioOfSpecificHeats *                    &
               specificVolume_ / temperature_ * velocity(3)
          deltaViscosity(5) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_
          do k = 1, nSpecies
             deltaViscosity(5+k) = 0.0_WP
          end do
          temp = velocity(1) * contravariantStressTensor(1) +                                &
               velocity(2) * contravariantStressTensor(2) +                                  &
               velocity(3) * contravariantStressTensor(3) - contravariantHeatFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,1) = deltaViscosity(1) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,1) = deltaViscosity(1) * temp - specificVolume_ *    &
               (velocity(1) * contravariantStressTensor(1) +                                 &
               velocity(2) * contravariantStressTensor(2) +                                  &
               velocity(3) * contravariantStressTensor(3))
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,2) = deltaViscosity(2) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,2) = deltaViscosity(2) * temp +                      &
               specificVolume_ * contravariantStressTensor(1)
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,3) = deltaViscosity(3) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,3) = deltaViscosity(3) * temp +                      &
               specificVolume_ * contravariantStressTensor(2)
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,4) = 0.0_WP
          firstPartialViscousJacobian(2,4) = deltaViscosity(4) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,4) = deltaViscosity(4) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,4) = deltaViscosity(4) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,4) = deltaViscosity(4) * temp +                      &
               specificVolume_ * contravariantStressTensor(3)
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,4) = - deltaViscosity(4) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,5) = 0.0_WP
          firstPartialViscousJacobian(2,5) = deltaViscosity(5) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,5) = deltaViscosity(5) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,5) = deltaViscosity(5) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,5) = deltaViscosity(5) * temp
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,5) = - deltaViscosity(5) *                      &
                  contravariantSpeciesFlux(k)
          end do

       else ! ... ideal gas mixture

          ! Other dependent variables
          phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                            &
               (velocity_(1) ** 2 + velocity_(2) ** 2 + velocity_(3) ** 2)
          contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                      &
               metrics(2) * stressTensor(2) + metrics(3) * stressTensor(3) !... not normalized
          contravariantStressTensor(2) = metrics(1) * stressTensor(4) +                      &
               metrics(2) * stressTensor(5) + metrics(3) * stressTensor(6) !... not normalized
          contravariantStressTensor(3) = metrics(1) * stressTensor(7) +                      &
               metrics(2) * stressTensor(8) + metrics(3) * stressTensor(9) !... not normalized
          contravariantHeatFlux = metrics(1) * heatFlux(1) + metrics(2) * heatFlux(2) +      &
               metrics(3) * heatFlux(3) !... not normalized
          do k = 1, nSpecies
             contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                   &
                  metrics(2) * speciesFlux(k,2) +                                            &
                  metrics(3) * speciesFlux(k,3) !... not normalized
          end do
          contravariantEnthalpyFlux = metrics(1) * enthalpyFlux(1) +                         &
               metrics(2) * enthalpyFlux(2) + metrics(3) * enthalpyFlux(3) !... not normalized
          temp = molecularWeightInverse(nSpecies+1)
          do k = 1, nSpecies
             temp = temp + massFraction_(k) * (molecularWeightInverse(k) -                   &
                  molecularWeightInverse(nSpecies+1))
          end do
          temp = 1.0_WP / temp
          deltaViscosity(1) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaTemperature(1) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp * (phiSquared / (ratioOfSpecificHeats - 1.0_WP) -         &
               temperature_ * molecularWeightInverse(nSpecies+1) / ratioOfSpecificHeats)
          deltaViscosity(2) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(1) * temp
          deltaTemperature(2) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(1) * temp
          deltaViscosity(3) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(2) * temp
          deltaTemperature(3) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(2) * temp
          deltaViscosity(4) = - powerLawExponent * ratioOfSpecificHeats * specificVolume_ /  &
               temperature_ * velocity(3) * temp
          deltaTemperature(4) = - ratioOfSpecificHeats * specificVolume_ /                   &
               temperature_ * velocity(3) * temp
          deltaViscosity(5) = powerLawExponent * ratioOfSpecificHeats * specificVolume_ /    &
               temperature_ * temp
          deltaTemperature(5) = ratioOfSpecificHeats * specificVolume_ /                     &
               temperature_ * temp
          do k = 1, nSpecies
             deltaViscosity(5+k) = - powerLawExponent * specificVolume_ * temp *             &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
             deltaTemperature(5+k) = - specificVolume_ * temp *                              &
                  (molecularWeightInverse(k) - molecularWeightInverse(nSpecies+1))
          end do
          temp = velocity(1) * contravariantStressTensor(1) +                                &
               velocity(2) * contravariantStressTensor(2) +                                  &
               velocity(3) * contravariantStressTensor(3) - contravariantHeatFlux -          &
               contravariantEnthalpyFlux

          firstPartialViscousJacobian(1,1) = 0.0_WP
          firstPartialViscousJacobian(2,1) = deltaViscosity(1) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,1) = deltaViscosity(1) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,1) = deltaViscosity(1) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,1) = deltaViscosity(1) * temp - specificVolume_ *    &
               (velocity(1) * contravariantStressTensor(1) + velocity(2) *                   &
               contravariantStressTensor(2) + velocity(3) * contravariantStressTensor(3)) -  &
               deltaTemperature(1) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,1) = - deltaViscosity(1) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,2) = 0.0_WP
          firstPartialViscousJacobian(2,2) = deltaViscosity(2) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,2) = deltaViscosity(2) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,2) = deltaViscosity(2) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,2) = deltaViscosity(2) * temp + specificVolume_ *    &
               contravariantStressTensor(1) - deltaTemperature(2) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,2) = - deltaViscosity(2) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,3) = 0.0_WP
          firstPartialViscousJacobian(2,3) = deltaViscosity(3) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,3) = deltaViscosity(3) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,3) = deltaViscosity(3) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,3) = deltaViscosity(3) * temp + specificVolume_ *    &
               contravariantStressTensor(2) - deltaTemperature(3) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,3) = - deltaViscosity(3) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,4) = 0.0_WP
          firstPartialViscousJacobian(2,4) = deltaViscosity(4) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,4) = deltaViscosity(4) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,4) = deltaViscosity(4) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,4) = deltaViscosity(4) * temp + specificVolume_ *    &
               contravariantStressTensor(3) - deltaTemperature(4) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,4) = - deltaViscosity(4) *                      &
                  contravariantSpeciesFlux(k)
          end do

          firstPartialViscousJacobian(1,5) = 0.0_WP
          firstPartialViscousJacobian(2,5) = deltaViscosity(5) * contravariantStressTensor(1)
          firstPartialViscousJacobian(3,5) = deltaViscosity(5) * contravariantStressTensor(2)
          firstPartialViscousJacobian(4,5) = deltaViscosity(5) * contravariantStressTensor(3)
          firstPartialViscousJacobian(5,5) = deltaViscosity(5) * temp -                      &
               deltaTemperature(5) * contravariantEnthalpyFlux
          do k = 1, nSpecies
             firstPartialViscousJacobian(5+k,5) = - deltaViscosity(5) *                      &
                  contravariantSpeciesFlux(k)
          end do

          do k = 1, nSpecies
             firstPartialViscousJacobian(1,5+k) = 0.0_WP
             firstPartialViscousJacobian(2,5+k) = deltaViscosity(5+k) *                      &
                  contravariantStressTensor(1)
             firstPartialViscousJacobian(3,5+k) = deltaViscosity(5+k) *                      &
                  contravariantStressTensor(2)
             firstPartialViscousJacobian(4,5+k) = deltaViscosity(5+k) *                      &
                  contravariantStressTensor(3)
             firstPartialViscousJacobian(5,5+k) = deltaViscosity(5+k) * temp -               &
                  deltaTemperature(5+k) * contravariantEnthalpyFlux
             do i = 1, nSpecies
                firstPartialViscousJacobian(5+i,5+k) = - deltaViscosity(5+k) *               &
                     contravariantSpeciesFlux(i)
             end do
          end do

       end if

    end select

    return
  end subroutine compute_first_partial_viscous_jacobian

  pure subroutine compute_second_partial_viscous_jacobian(velocity, dynamicViscosity,        &
       secondCoefficientOfViscosity, thermalDiffusivity, massDiffusivity, temperature,       &
       jacobian, metricsAlongFirstDir, metricsAlongSecondDir, secondPartialViscousJacobian)

    implicit none

    ! Arguments
    real(WP), intent(in) :: velocity(:), dynamicViscosity, secondCoefficientOfViscosity,     &
         thermalDiffusivity, jacobian, metricsAlongFirstDir(:)
    real(WP), intent(in), optional :: massDiffusivity(:), temperature,                       &
         metricsAlongSecondDir(:)
    real(WP), intent(out) :: secondPartialViscousJacobian(:,:)

    ! Local variables
    integer :: k
    real(WP) :: temp1, temp2, temp3

    ! Zero-out second partial viscous Jacobian
    secondPartialViscousJacobian = 0.0_WP

    select case (nDimensions)

    case (1)

       ! Temporary variables
       temp1 = metricsAlongFirstDir(1) * metricsAlongFirstDir(1)
       temp2 = dynamicViscosity * metricsAlongFirstDir(1) * velocity(1)
       temp3 = secondCoefficientOfViscosity * metricsAlongFirstDir(1) * velocity(1)

       secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(1) * metricsAlongFirstDir(1)
       secondPartialViscousJacobian(2,1) = dynamicViscosity * temp1 * velocity(1) +          &
            metricsAlongFirstDir(1) * temp2 + metricsAlongFirstDir(1) * temp3

       secondPartialViscousJacobian(1,2) = 0.0_WP
       secondPartialViscousJacobian(2,2) = thermalDiffusivity * temp1

       if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
          do k = 1, nSpecies
             secondPartialViscousJacobian(2,2+k) = massDiffusivity(k) * temperature /        &
                  schmidtNumberInverse(k) * (molecularWeightInverse(k) *                     &
                  schmidtNumberInverse(k) - molecularWeightInverse(nSpecies+1) *             &
                  schmidtNumberInverse(nSpecies+1)) * temp1
          end do
       end if

       do k = 1, nSpecies
          secondPartialViscousJacobian(2+k,2+k) = massDiffusivity(k) * temp1
       end do

    case (2)

       ! Temporary variables
       temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                          &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
       temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                  &
            metricsAlongSecondDir(2) * velocity(2))
       temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +       &
            metricsAlongFirstDir(2) * velocity(2))

       secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
       secondPartialViscousJacobian(2,1) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
       secondPartialViscousJacobian(3,1) = dynamicViscosity * temp1 * velocity(1) +          &
            metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

       secondPartialViscousJacobian(1,2) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
       secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
       secondPartialViscousJacobian(3,2) = dynamicViscosity * temp1 * velocity(2) +          &
            metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

       secondPartialViscousJacobian(1,3) = 0.0_WP
       secondPartialViscousJacobian(2,3) = 0.0_WP
       secondPartialViscousJacobian(3,3) = thermalDiffusivity * temp1

       if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
          do k = 1, nSpecies
             secondPartialViscousJacobian(3,3+k) = massDiffusivity(k) * temperature /        &
                  schmidtNumberInverse(k) * (molecularWeightInverse(k) *                     &
                  schmidtNumberInverse(k) - molecularWeightInverse(nSpecies+1) *             &
                  schmidtNumberInverse(nSpecies+1)) * temp1
          end do
       end if

       do k = 1, nSpecies
          secondPartialViscousJacobian(3+k,3+k) = massDiffusivity(k) * temp1
       end do

    case (3)

       ! Temporary variables
       temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                          &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2) +                             &
            metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
       temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                  &
            metricsAlongSecondDir(2) * velocity(2) + metricsAlongSecondDir(3) * velocity(3))
       temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +       &
            metricsAlongFirstDir(2) * velocity(2) + metricsAlongFirstDir(3) * velocity(3))

       secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
       secondPartialViscousJacobian(2,1) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
       secondPartialViscousJacobian(3,1) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1)
       secondPartialViscousJacobian(4,1) = dynamicViscosity * temp1 * velocity(1) +          &
            metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

       secondPartialViscousJacobian(1,2) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
       secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
       secondPartialViscousJacobian(3,2) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2)
       secondPartialViscousJacobian(4,2) = dynamicViscosity * temp1 * velocity(2) +          &
            metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

       secondPartialViscousJacobian(1,3) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3)
       secondPartialViscousJacobian(2,3) =                                                   &
            dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2) +          &
            secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3)
       secondPartialViscousJacobian(3,3) = dynamicViscosity * temp1 +                        &
            (dynamicViscosity + secondCoefficientOfViscosity) *                              &
            metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
       secondPartialViscousJacobian(4,3) = dynamicViscosity * temp1 * velocity(3) +          &
            metricsAlongFirstDir(3) * temp2 + metricsAlongSecondDir(3) * temp3

       secondPartialViscousJacobian(1,4) = 0.0_WP
       secondPartialViscousJacobian(2,4) = 0.0_WP
       secondPartialViscousJacobian(3,4) = 0.0_WP
       secondPartialViscousJacobian(4,4) = thermalDiffusivity * temp1

       if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
          do k = 1, nSpecies
             secondPartialViscousJacobian(4,4+k) = massDiffusivity(k) * temperature /        &
                  schmidtNumberInverse(k) * (molecularWeightInverse(k) *                     &
                  schmidtNumberInverse(k) - molecularWeightInverse(nSpecies+1) *             &
                  schmidtNumberInverse(nSpecies+1)) * temp1
          end do
       end if

       do k = 1, nSpecies
          secondPartialViscousJacobian(4+k,4+k) = massDiffusivity(k) * temp1
       end do

    end select

    ! Multiply by the Jacobian
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

    return
  end subroutine compute_second_partial_viscous_jacobian


  subroutine compute_Steger_Warming_jacobian_of_inviscid_flux(conservedVariables,            &
       metrics, leftJacobianOfInviscidFlux, rightJacobianOfInviscidFlux,                     &
       specificVolume, velocity, pressure, massFraction)

    implicit none

    ! Arguments
    real(WP), intent(in) :: conservedVariables(:), metrics(:), specificVolume, velocity(:),  &
         pressure, massFraction(:)
    real(WP), intent(out) :: leftJacobianOfInviscidFlux(:,:), rightJacobianOfInviscidFlux(:,:)

    ! Local variables
    integer :: i, j, k
    real(WP) :: arcLength, normalizedMetrics(nDimensions), contravariantVelocity,            &
         speedOfSound, phiSquared, enthalpy, temperature,                                    &
         rightEigenvectors(nDimensions+nSpecies+2, nDimensions+nSpecies+2),                  &
         leftEigenvectors(nDimensions+nSpecies+2, nDimensions+nSpecies+2),                   &
         eigenvalues(nDimensions+nSpecies+2), posEigenvalues(nDimensions+nSpecies+2),        &
         negEigenvalues(nDimensions+nSpecies+2)

    ! Normalize the metrics
    arcLength = sqrt(sum(metrics ** 2))
    normalizedMetrics = metrics / arcLength

    ! Temperature (this is not the actual temperature, assumes ideal gas law)
    temperature = ratioOfSpecificHeats * pressure * specificVolume /                         &
         (ratioOfSpecificHeats - 1.0_WP)

    ! Dependent variables
    contravariantVelocity = sum(normalizedMetrics * velocity)
    speedOfSound = sqrt(ratioOfSpecificHeats * pressure * specificVolume)
    phiSquared = 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) * sum(velocity ** 2)
    enthalpy = specificVolume * (conservedVariables(nDimensions+2) + pressure)

    select case (nDimensions)

    case (1)

       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity + speedOfSound
       eigenvalues(3) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(3+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       ! Select positive and negative eigenvalues
       posEigenvalues = 0.5_WP * (eigenvalues + abs(eigenvalues))
       negEigenvalues = 0.5_WP * (eigenvalues - abs(eigenvalues))

       ! Matrix whose columns are the right eigenvectors:
       rightEigenvectors(1,1) = 1.0_WP
       rightEigenvectors(2,1) = velocity(1)
       rightEigenvectors(3,1) = phiSquared / (ratioOfSpecificHeats - 1.0_WP)
       do k = 1, nSpecies
          rightEigenvectors(3+k,1) = massFraction(k)
       end do

       rightEigenvectors(1,2) = 1.0_WP
       rightEigenvectors(2,2) = velocity(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,2) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(3+k,2) = massFraction(k)
       end do

       rightEigenvectors(1,3) = 1.0_WP
       rightEigenvectors(2,3) = velocity(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,3) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(3+k,3) = massFraction(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(3+k,3+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:
       leftEigenvectors(1,1) = 1.0_WP - phiSquared / speedOfSound ** 2
       leftEigenvectors(2,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(3,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(3+k,1) = -massFraction(k)
       end do

       leftEigenvectors(1,2) = velocity(1) / temperature
       leftEigenvectors(2,2) = - 0.5_WP * (velocity(1) / temperature -                       &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(3,2) = - 0.5_WP * (velocity(1) / temperature +                       &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(3+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = - 1.0_WP / temperature
       leftEigenvectors(2,3) = 0.5_WP / temperature
       leftEigenvectors(3,3) = 0.5_WP / temperature
       do k = 1, nSpecies
          leftEigenvectors(3+k,3) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(3+k,3+k) = 1.0_WP
       end do

       ! Perform splitting
       leftJacobianOfInviscidFlux  = 0.0_WP
       rightJacobianOfInviscidFlux = 0.0_WP
       do j = 1, 3 + nSpecies
          do i = 1, 3 + nSpecies
             do k = 1, 3 + nSpecies
                leftJacobianOfInviscidFlux(i,j) = leftJacobianOfInviscidFlux(i,j) +          &
                     rightEigenvectors(i,k) * posEigenvalues(k) * leftEigenvectors(k,j)

                rightJacobianOfInviscidFlux(i,j) = rightJacobianOfInviscidFlux(i,j) +        &
                     rightEigenvectors(i,k) * negEigenvalues(k) * leftEigenvectors(k,j)
             end do
          end do
       end do

    case (2)
       
       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity
       eigenvalues(3) = contravariantVelocity + speedOfSound
       eigenvalues(4) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(4+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       ! Select positive and negative eigenvalues
       posEigenvalues = 0.5_WP * (eigenvalues + abs(eigenvalues))
       negEigenvalues = 0.5_WP * (eigenvalues - abs(eigenvalues))
       
       ! Matrix whose columns are the right eigenvectors:
       rightEigenvectors(1,1) = 1.0_WP
       rightEigenvectors(2,1) = velocity(1)
       rightEigenvectors(3,1) = velocity(2)
       rightEigenvectors(4,1) = phiSquared / (ratioOfSpecificHeats - 1.0_WP)
       do k = 1, nSpecies
          rightEigenvectors(4+k,1) = massFraction(k)
       end do

       rightEigenvectors(1,2) = 0.0_WP
       rightEigenvectors(2,2) = normalizedMetrics(2) * conservedVariables(1)
       rightEigenvectors(3,2) = - normalizedMetrics(1) * conservedVariables(1)
       rightEigenvectors(4,2) = conservedVariables(1) * (normalizedMetrics(2) *              &
            velocity(1) - normalizedMetrics(1) * velocity(2))
       do k = 1, nSpecies
          rightEigenvectors(4+k,2) = 0.0_WP
       end do

       rightEigenvectors(1,3) = 1.0_WP
       rightEigenvectors(2,3) = velocity(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,3) = velocity(2) + normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,3) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(4+k,3) = massFraction(k)
       end do

       rightEigenvectors(1,4) = 1.0_WP
       rightEigenvectors(2,4) = velocity(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,4) = velocity(2) - normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,4) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(4+k,4) = massFraction(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(4+k,4+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:
       leftEigenvectors(1,1) = 1.0_WP - phiSquared / speedOfSound ** 2
       leftEigenvectors(2,1) = - specificVolume * (normalizedMetrics(2) * velocity(1) -      &
            normalizedMetrics(1) * velocity(2))
       leftEigenvectors(3,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(4,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,1) = -massFraction(k)
       end do

       leftEigenvectors(1,2) = velocity(1) / temperature
       leftEigenvectors(2,2) = specificVolume * normalizedMetrics(2)
       leftEigenvectors(3,2) = - 0.5_WP * (velocity(1) / temperature -                       &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(4,2) = - 0.5_WP * (velocity(1) / temperature +                       &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = velocity(2) / temperature
       leftEigenvectors(2,3) = - specificVolume * normalizedMetrics(1)
       leftEigenvectors(3,3) = - 0.5_WP * (velocity(2) / temperature -                       &
            normalizedMetrics(2) / speedOfSound)
       leftEigenvectors(4,3) = - 0.5_WP * (velocity(2) / temperature +                       &
            normalizedMetrics(2) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(4+k,3) = 0.0_WP
       end do

       leftEigenvectors(1,4) = - 1.0_WP / temperature
       leftEigenvectors(2,4) = 0.0_WP
       leftEigenvectors(3,4) = 0.5_WP / temperature
       leftEigenvectors(4,4) = 0.5_WP / temperature
       do k = 1, nSpecies
          leftEigenvectors(4+k,4) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(4+k,4+k) = 1.0_WP
       end do

       ! Perform splitting
       leftJacobianOfInviscidFlux  = 0.0_WP
       rightJacobianOfInviscidFlux = 0.0_WP
       do j = 1, 4 + nSpecies
          do i = 1, 4 + nSpecies
             do k = 1, 4 + nSpecies
                leftJacobianOfInviscidFlux(i,j) = leftJacobianOfInviscidFlux(i,j) +          &
                     rightEigenvectors(i,k) * posEigenvalues(k) * leftEigenvectors(k,j)

                rightJacobianOfInviscidFlux(i,j) = rightJacobianOfInviscidFlux(i,j) +        &
                     rightEigenvectors(i,k) * negEigenvalues(k) * leftEigenvectors(k,j)
             end do
          end do
       end do

    case (3)

       ! Eigenvalues
       eigenvalues(1) = contravariantVelocity
       eigenvalues(2) = contravariantVelocity
       eigenvalues(3) = contravariantVelocity
       eigenvalues(4) = contravariantVelocity + speedOfSound
       eigenvalues(5) = contravariantVelocity - speedOfSound
       do k = 1, nSpecies
          eigenvalues(5+k) = contravariantVelocity
       end do
       eigenvalues = arcLength * eigenvalues

       ! Select positive and negative eigenvalues
       posEigenvalues = 0.5_WP * (eigenvalues + abs(eigenvalues))
       negEigenvalues = 0.5_WP * (eigenvalues - abs(eigenvalues))

       ! Matrix whose columns are the right eigenvectors:
       rightEigenvectors(1,1) = normalizedMetrics(1)
       rightEigenvectors(2,1) = normalizedMetrics(1) * velocity(1)
       rightEigenvectors(3,1) = normalizedMetrics(1) * velocity(2) +                         &
            conservedVariables(1) * normalizedMetrics(3)
       rightEigenvectors(4,1) = normalizedMetrics(1) * velocity(3) -                         &
            conservedVariables(1) * normalizedMetrics(2)
       rightEigenvectors(5,1) = conservedVariables(1) * (normalizedMetrics(3) *              &
            velocity(2) - normalizedMetrics(2) * velocity(3)) + phiSquared /                 &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(1)
       do k = 1, nSpecies
          rightEigenvectors(5+k,1) = massFraction(k) * normalizedMetrics(1)
       end do

       rightEigenvectors(1,2) = normalizedMetrics(2)
       rightEigenvectors(2,2) = normalizedMetrics(2) * velocity(1) -                         &
            conservedVariables(1) * normalizedMetrics(3)
       rightEigenvectors(3,2) = normalizedMetrics(2) * velocity(2)
       rightEigenvectors(4,2) = normalizedMetrics(2) * velocity(3) +                         &
            conservedVariables(1) * normalizedMetrics(1)
       rightEigenvectors(5,2) = conservedVariables(1) * (normalizedMetrics(1) *              &
            velocity(3) - normalizedMetrics(3) * velocity(1)) + phiSquared /                 &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(2)
       do k = 1, nSpecies
          rightEigenvectors(5+k,2) = massFraction(k) * normalizedMetrics(2)
       end do

       rightEigenvectors(1,3) = normalizedMetrics(3)
       rightEigenvectors(2,3) = normalizedMetrics(3) * velocity(1) +                         &
            conservedVariables(1) * normalizedMetrics(2)
       rightEigenvectors(3,3) = normalizedMetrics(3) * velocity(2) -                         &
            conservedVariables(1) * normalizedMetrics(1)
       rightEigenvectors(4,3) = normalizedMetrics(3) * velocity(3)
       rightEigenvectors(5,3) = conservedVariables(1) * (normalizedMetrics(2) *              &
            velocity(1) - normalizedMetrics(1) * velocity(2)) + phiSquared /                 &
            (ratioOfSpecificHeats - 1.0_WP) * normalizedMetrics(3)
       do k = 1, nSpecies
          rightEigenvectors(5+k,3) = massFraction(k) * normalizedMetrics(3)
       end do

       rightEigenvectors(1,4) = 1.0_WP
       rightEigenvectors(2,4) = velocity(1) + normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,4) = velocity(2) + normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,4) = velocity(3) + normalizedMetrics(3) * speedOfSound
       rightEigenvectors(5,4) = enthalpy + speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(5+k,4) = massFraction(k)
       end do

       rightEigenvectors(1,5) = 1.0_WP
       rightEigenvectors(2,5) = velocity(1) - normalizedMetrics(1) * speedOfSound
       rightEigenvectors(3,5) = velocity(2) - normalizedMetrics(2) * speedOfSound
       rightEigenvectors(4,5) = velocity(3) - normalizedMetrics(3) * speedOfSound
       rightEigenvectors(5,5) = enthalpy - speedOfSound * contravariantVelocity
       do k = 1, nSpecies
          rightEigenvectors(5+k,5) = massFraction(k)
       end do

       do k = 1, nSpecies
          rightEigenvectors(5+k,5+k) = 1.0_WP
       end do

       ! Matrix whose rows are the left eigenvectors:
       leftEigenvectors(1,1) = normalizedMetrics(1) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume *                     &
            (normalizedMetrics(3) * velocity(2) - normalizedMetrics(2) * velocity(3))
       leftEigenvectors(2,1) = normalizedMetrics(2) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume *                     &
            (normalizedMetrics(1) * velocity(3) - normalizedMetrics(3) * velocity(1))
       leftEigenvectors(3,1) = normalizedMetrics(3) *                                        &
            (1.0_WP - phiSquared / speedOfSound ** 2) - specificVolume *                     &
            (normalizedMetrics(2) * velocity(1) - normalizedMetrics(1) * velocity(2))
       leftEigenvectors(4,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 -                    &
            contravariantVelocity / speedOfSound)
       leftEigenvectors(5,1) = 0.5_WP * (phiSquared / speedOfSound ** 2 +                    &
            contravariantVelocity / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,1) = -massFraction(k)
       end do

       leftEigenvectors(1,2) = normalizedMetrics(1) * velocity(1) / temperature
       leftEigenvectors(2,2) = normalizedMetrics(2) * velocity(1) / temperature -            &
            specificVolume * normalizedMetrics(3)
       leftEigenvectors(3,2) = normalizedMetrics(3) * velocity(1) / temperature +            &
            specificVolume * normalizedMetrics(2)
       leftEigenvectors(4,2) = - 0.5_WP * (velocity(1) / temperature -                       &
            normalizedMetrics(1) / speedOfSound)
       leftEigenvectors(5,2) = - 0.5_WP * (velocity(1) / temperature +                       &
            normalizedMetrics(1) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,2) = 0.0_WP
       end do

       leftEigenvectors(1,3) = normalizedMetrics(1) * velocity(2) / temperature +            &
            specificVolume * normalizedMetrics(3)
       leftEigenvectors(2,3) = normalizedMetrics(2) * velocity(2) / temperature
       leftEigenvectors(3,3) = normalizedMetrics(3) * velocity(2) / temperature -            &
            specificVolume * normalizedMetrics(1)
       leftEigenvectors(4,3) = - 0.5_WP * (velocity(2) / temperature -                       &
            normalizedMetrics(2) / speedOfSound)
       leftEigenvectors(5,3) = - 0.5_WP * (velocity(2) / temperature +                       &
            normalizedMetrics(2) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,3) = 0.0_WP
       end do

       leftEigenvectors(1,4) = normalizedMetrics(1) * velocity(3) / temperature -            &
            specificVolume * normalizedMetrics(2)
       leftEigenvectors(2,4) = normalizedMetrics(2) * velocity(3) / temperature +            &
            specificVolume * normalizedMetrics(1)
       leftEigenvectors(3,4) = normalizedMetrics(3) * velocity(3) / temperature
       leftEigenvectors(4,4) = - 0.5_WP * (velocity(3) / temperature -                       &
            normalizedMetrics(3) / speedOfSound)
       leftEigenvectors(5,4) = - 0.5_WP * (velocity(3) / temperature +                       &
            normalizedMetrics(3) / speedOfSound)
       do k = 1, nSpecies
          leftEigenvectors(5+k,4) = 0.0_WP
       end do

       leftEigenvectors(1,5) = - normalizedMetrics(1) / temperature
       leftEigenvectors(2,5) = - normalizedMetrics(2) / temperature
       leftEigenvectors(3,5) = - normalizedMetrics(3) / temperature
       leftEigenvectors(4,5) = 0.5_WP / temperature
       leftEigenvectors(5,5) = 0.5_WP / temperature
       do k = 1, nSpecies
          leftEigenvectors(5+k,5) = 0.0_WP
       end do

       do k = 1, nSpecies
          leftEigenvectors(5+k,5+k) = 1.0_WP
       end do

       ! Perform splitting
       leftJacobianOfInviscidFlux  = 0.0_WP
       rightJacobianOfInviscidFlux = 0.0_WP
       do j = 1, 5 + nSpecies
          do i = 1, 5 + nSpecies
             do k = 1, 5 + nSpecies
                leftJacobianOfInviscidFlux(i,j) = leftJacobianOfInviscidFlux(i,j) +          &
                     rightEigenvectors(i,k) * posEigenvalues(k) * leftEigenvectors(k,j)

                rightJacobianOfInviscidFlux(i,j) = rightJacobianOfInviscidFlux(i,j) +        &
                     rightEigenvectors(i,k) * negEigenvalues(k) * leftEigenvectors(k,j)
             end do
          end do
       end do

    end select

    return
  end subroutine compute_Steger_Warming_jacobian_of_inviscid_flux

end module state_jacobian
