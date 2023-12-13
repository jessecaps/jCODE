module well_stirred

  ! External modules
  use combustion
  use precision
  use solver_options
  use geometry

  implicit none

  real(WP) :: Tin, residenceTime
  real(WP), dimension(:), allocatable :: Yin
  
contains

  ! Compute the Jacobian matrix for the WSR
  ! ---------------------------------------
  subroutine well_stirred_jacobian(stateVector, jacobianOfSource, specificVolume,            &
       velocity, temperature, massFraction, molecularWeightInverse, Wmix)

    implicit none

    ! Arguments
    real(WP), intent(in) :: stateVector(:), specificVolume, velocity(:), temperature,        &
       massFraction(:), molecularWeightInverse(:), Wmix
    real(WP), intent(out) :: jacobianOfSource(:,:)

    ! Local variables
    integer :: k
    real(WP) :: temp

    ! Zero-out source Jacobian
    jacobianOfSource = 0.0_WP

    temp = -residenceTime

    select case (equationOfState)

    case(IDEAL_GAS)

       jacobianOfSource(nDimensions+2,1) = jacobianOfSource(nDimensions+2,1) + temp *        &
            (specificVolume * stateVector(nDimensions+2) - Tin / ratioOfSpecificHeats)
       do k = 1, nDimensions
          jacobianOfSource(nDimensions+2,k+1) = jacobianOfSource(nDimensions+2,k+1) + temp * &
               0.5_WP * stateVector(k+1)
       end do
       jacobianOfSource(nDimensions+2,nDimensions+2) =                                       &
            jacobianOfSource(nDimensions+2,nDimensions+2) + temp * stateVector(1) /          &
            ratioOfSpecificHeats

    case (IDEAL_GAS_MIXTURE)

       jacobianOfSource(nDimensions+2,1) = jacobianOfSource(nDimensions+2,1) + temp *        &
            (specificVolume * stateVector(nDimensions+2) - Tin / ratioOfSpecificHeats)
       do k = 1, nDimensions
          jacobianOfSource(nDimensions+2,k+1) = jacobianOfSource(nDimensions+2,k+1) + temp * &
               0.5_WP * stateVector(k+1)
       end do
       jacobianOfSource(nDimensions+2,nDimensions+2) =                                       &
            jacobianOfSource(nDimensions+2,nDimensions+2) + temp * stateVector(1) /          &
            Wmix / ratioOfSpecificHeats
       do k = 1, nSpecies
          jacobianOfSource(nDimensions+2,nDimensions+2+k) =                                  &
               jacobianOfSource(nDimensions+2,nDimensions+2+k) + temp *                      &
               stateVector(1) * (temperature - Tin) * (molecularWeightInverse(k) -           &
               molecularWeightInverse(nSpecies+1)) / ratioOfSpecificHeats
       end do

    end select

    do k = 1, nSpecies
       jacobianOfSource(nDimensions+2+k,1) = jacobianOfSource(nDimensions+2+k,1) + temp *    &
            (massFraction(k) - Yin(k))
       jacobianOfSource(nDimensions+2+k,nDimensions+2+k) =                                   &
            jacobianOfSource(nDimensions+2+k,nDimensions+2+k) + temp * stateVector(1)
    end do

    return
  end subroutine well_stirred_jacobian

end module well_stirred


! ============================== !
! Setup the well-stirred reactor !
! ============================== !
subroutine well_stirred_setup

  ! Internal modules
  use well_stirred

  ! External modules
  use string
  use parser
  use state

  implicit none

  ! Local variables
  integer :: nInput

  if (.not. wellStirredReactor) return

  if (any(.not. isPeriodic(1:nDimensions)))                                                  &
       call die ('well_stirred_setup: WSR requires fully periodic boundaries')

  allocate(Yin(nSpecies))

  call parser_read('residence time', residenceTime)
  residenceTime = 1.0_WP / residenceTime
  call parser_read('inlet temperature', Tin)
  call parser_getsize('inlet mass fraction', nInput)
  if (nInput .ne. nSpecies) call die("Error reading value of 'inlet mass fraction'.&
                & Please input nSpecies values")
  call parser_read('inlet mass fraction', Yin)

  return
end subroutine well_stirred_setup


! ================================ !
! Cleanup the well-stirred reactor !
! ================================ !
subroutine well_stirred_cleanup

  ! Internal modules
  use well_stirred

  implicit none

  if (allocated(Yin)) deallocate(Yin)

  return
end subroutine well_stirred_cleanup


! ========================================= !
! Add the WSR source during the forward run !
! ========================================= !
subroutine well_stirred_forward(source)

  ! Internal modules
  use well_stirred

  ! External modules
  use solver_options
  use geometry
  use state

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, k
  real(WP) :: density

  if (.not. wellStirredReactor) return

  do i = 1, nGridPoints
     density = conservedVariables(i, 1)
     select case (equationOfState)
     case (IDEAL_GAS)
        source(i,nDimensions+2) = source(i,nDimensions+2) + residenceTime *                  &
             (density * Tin / ratioOfSpecificHeats - conservedVariables(i,nDimensions+2))
     case (IDEAL_GAS_MIXTURE)
        source(i,nDimensions+2) = source(i,nDimensions+2) + residenceTime *                  &
             (density * Tin / ratioOfSpecificHeats / mixtureMolecularWeight(i,1) -           &
             conservedVariables(i,nDimensions+2))
     end select
     do k = 1, nSpecies
        source(i,nDimensions+2+k) = source(i,nDimensions+2+k) + residenceTime *              &
             (density * Yin(k) - conservedVariables(i,nDimensions+2+k))
     end do
  end do

  return
end subroutine well_stirred_forward


! ========================================= !
! Add the WSR source during the adjoint run !
! ========================================= !
subroutine well_stirred_adjoint(source)

  ! Internal modules
  use well_stirred

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

  if (.not. wellStirredReactor) return

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

     call well_stirred_jacobian(localConservedVariables, localJacobian,                      &
          specificVolume(j,1), localVelocity, localTemperature, localMassFraction,           &
          molecularWeightInverse, mixtureMolecularWeight(j,1))

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
        temp(1) = temp(1) - specificVolume(j,1) * localConservedVariables(nDimensions+2) *   &
             temp(nDimensions+2) - sum(localVelocity * temp(2:nDimensions+1)) -              &
             sum(localMassFraction * temp(ndimensions+3:nUnknowns))

     end select

     source(j,:) = source(j,:) - temp

  end do

  deallocate(localConservedVariables)
  deallocate(localVelocity)
  deallocate(localMassFraction)
  deallocate(localJacobian)
  deallocate(temp)

  return
end subroutine well_stirred_adjoint
