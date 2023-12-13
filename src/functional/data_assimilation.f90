module data_assimilation

  ! External modules
  use functional

  implicit none

  ! Assimilation types
  integer, parameter ::                                                                      &
       ASSIMILATE_VELOCITY     = 1,                                                          &
       ASSIMILATE_TEMPERATURE  = 2,                                                          &
       ASSIMILATE_HEAT_RELEASE = 3

  integer :: assimilationType, power, velocityComponent
  real(WP), dimension(:,:), allocatable :: data

contains

  subroutine verify_data_assimilation_patch(patch)

    ! External modules
    use solver_options

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions

    do i = 1, nDimensions
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    n = nDimensions
    if (n .ne. nDimensions) then
       write(message, '(2(A,I0.0),A)') 'verify_data_assimilation_patch: Expected a ',         &
            nDimensions, 'D patch, but extent represents a ', n, 'D patch!'
       call die(trim(message))
    end if

    return
  end subroutine verify_data_assimilation_patch

end module data_assimilation


! =========================================== !
! Setup the data assimilation cost functional !
! =========================================== !
subroutine data_assimilation_setup

  ! Internal modules
  use data_assimilation

  ! External modules
  use parser
  use solver_options
  use simulation_flags
  use combustion

  implicit none

  ! Local variables
  character(len = str_medium) :: filename, type

  ! Verify the patch type
  call verify_data_assimilation_patch(functionalPatch)

  ! Get the norm power
  call parser_read('norm power', power, 2)
  if (power .le. 1)                                                                          &
       call die('data_assimilation_setup: norm power must be > 1')

  ! Get the assimilation type
  call parser_read('data assimilation type', type)
  select case (trim(type))

  case ('velocity')

     assimilationType = ASSIMILATE_VELOCITY

     call parser_read('data assimilation velocity component', velocityComponent)

  case ('temperature')

     assimilationType = ASSIMILATE_TEMPERATURE

  case ('heat release')

     assimilationType = ASSIMILATE_HEAT_RELEASE

     if (chemistryModel .ne. ONE_STEP)                                                          &
          call die('data_assimilation_setup: data assimilation requires one-step chemistry!')

  case default

     call die("data_assimilation_setup: unknown assimilation type '" // trim(type) // "'")

  end select

  ! Allocate the data array
  allocate(data(nGridPoints, nUnknowns))

  ! Get the data
  call parser_read('data assimilation file', filename)
  if (useSerialIO) then
     call state_read_serial(data, trim(filename))
  else
     call state_read_parallel(data, trim(filename))
  end if

  return
end subroutine data_assimilation_setup


! ============================================= !
! Cleanup the data assimilation cost functional !
! ============================================= !
subroutine data_assimilation_cleanup

  ! Internal modules
  use data_assimilation

  implicit none

  ! Cleanup
  if (allocated(data)) deallocate(data)

  return
end subroutine data_assimilation_cleanup


! ============================================= !
! Compute the data assimilation cost functional !
! ============================================= !
subroutine data_assimilation_compute

  ! Internal modules
  use data_assimilation

  ! External modules
  use grid_functions, only : inner_product
  use state
  use equation_of_state
  use onestep

  implicit none

  ! Local variables
  integer :: i
  real(WP), allocatable :: F(:), dataTemperature(:)
  real(WP) :: chemicalSource, chemicalSourceData, referenceTemperature, flameTemperature,    &
       activationTemperature

  allocate(F(nGridPoints))

  select case (assimilationType)

  case (ASSIMILATE_VELOCITY)

     F = (velocity(:, velocityComponent) - data(:,velocityComponent+1)/data(:,1))**power

  case (ASSIMILATE_TEMPERATURE)

     ! Get temperature from the data file
     allocate(dataTemperature(nGridPoints))
     call compute_dependent_variables(data, temperature = dataTemperature)

     F = (temperature(:, 1) - dataTemperature(:))**power

  case (ASSIMILATE_HEAT_RELEASE)

     ! One-step chemistry
     ! H2 + sO2 => (1+s)P

     ! Get temperature from the data file
     allocate(dataTemperature(nGridPoints))
     call compute_dependent_variables(data, temperature = dataTemperature)

     referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
     flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
     activationTemperature = zelDovich / heatRelease * flameTemperature

     do i = 1, nGridPoints

        ! Chemical source from the simulation
        chemicalSource = -controlParameter * Damkohler * stoichiometricCoefficient(H2) *     &
             conservedVariables(i,nDimensions+2+H2) *                                        &
             conservedVariables(i,nDimensions+2+O2) *                                        &
             exp(- activationTemperature / temperature(i, 1))

        ! Chemical source from the data
        chemicalSourceData = -Damkohler * stoichiometricCoefficient(H2) *                    &
             data(i,nDimensions+2+H2) * data(i,nDimensions+2+O2) *                           &
             exp(- activationTemperature / dataTemperature(i))

        F(i) = (chemicalSource - chemicalSourceData)**power

     end do

     deallocate(dataTemperature)

  end select

  instantaneousCostFunctional = inner_product(F, targetMollifier(:,1))

  deallocate(F)

  return
end subroutine data_assimilation_compute


! ============================================= !
! Compute the data assimilation adjoint forcing !
! ============================================= !
subroutine data_assimilation_adjoint_source(source)

  ! Internal modules
  use data_assimilation

  ! External modules
  use simulation_flags
  use solver_options
  use state_jacobian
  use geometry
  use grid
  use state
  use equation_of_state
  use onestep

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: forcingFactor, F, chemicalSource, chemicalSourceData, referenceTemperature,    &
       flameTemperature, activationTemperature, dataVelocity, temp
  real(WP), allocatable :: dataTemperature(:), deltaTemperature(:)

  if (useContinuousAdjoint) then
     forcingFactor = 1.0_WP
  else
     forcingFactor = adjointForcingFactor
  end if

  select case (assimilationType)

  case (ASSIMILATE_VELOCITY)

     do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
        do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
           do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              dataVelocity = data(gridIndex, velocityComponent+1) / data(gridIndex, 1)

              F = - forcingFactor * targetMollifier(gridIndex, 1) *                          &
                   specificVolume(gridIndex, 1)

              source(gridIndex, 1) = source(gridIndex, 1) + F * (-real(power,WP)) *          &
                   velocity(gridIndex, velocityComponent) *                                  &   
                   (velocity(gridIndex, velocityComponent) - dataVelocity)**(power-1)
              source(gridIndex, velocityComponent+1) =                                       &
                   source(gridIndex, velocityComponent+1) + F * real(power,WP) *             &
                   (velocity(gridIndex, velocityComponent) - dataVelocity)**(power-1)

           end do
        end do
     end do

  case (ASSIMILATE_TEMPERATURE)

     ! Get temperature from the data file
     allocate(dataTemperature(nGridPoints))
     call compute_dependent_variables(data, temperature = dataTemperature)

     allocate(deltaTemperature(nDimensions + nSpecies + 2))

     do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
        do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
           do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              ! Compute temperature variations
              call compute_delta_variables(conservedVariables(gridIndex,:),                  &
                   deltaTemperature = deltaTemperature)

              source(gridIndex,:) = source(gridIndex,:) - forcingFactor *                    &
                   targetMollifier(gridIndex, 1) * real(power, WP) *                         &
                   (temperature(gridIndex, 1) - dataTemperature(gridIndex))**(power-1) *     &
                   deltatemperature

           end do
        end do
     end do

     ! Cleanup
     deallocate(dataTemperature)
     deallocate(deltaTemperature)

  case (ASSIMILATE_HEAT_RELEASE)

     ! Get temperature from the data file
     allocate(dataTemperature(nGridPoints))
     call compute_dependent_variables(data, temperature = dataTemperature)

     allocate(deltaTemperature(nDimensions + nSpecies + 2))

     referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
     flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
     activationTemperature = zelDovich / heatRelease * flameTemperature

     do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
        do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
           do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              ! Chemical source from the simulation
              temp = -controlParameter * Damkohler * stoichiometricCoefficient(H2) *         &
                   exp(- activationTemperature / temperature(gridIndex, 1))
              chemicalSource = temp * conservedVariables(gridIndex,nDimensions+2+H2) *       &
                   conservedVariables(gridIndex,nDimensions+2+O2)

              ! Chemical source from the data
              chemicalSourceData = -Damkohler * stoichiometricCoefficient(H2) *              &
                   data(gridIndex,nDimensions+2+H2) * data(gridIndex,nDimensions+2+O2) *     &
                   exp(- activationTemperature / dataTemperature(gridIndex))

              ! Compute temperature variations
              call compute_delta_variables(conservedVariables(gridIndex,:),                  &
                   deltaTemperature = deltaTemperature)

              F = - forcingFactor * targetMollifier(gridIndex, 1) * real(power, WP) *        &
                   (chemicalSource - chemicalSourceData)**(power-1)

              source(gridIndex,:) = source(gridIndex,:) + F * chemicalSource *               &
                   activationTemperature / temperature(gridIndex, 1)**2 * deltaTemperature

              source(gridIndex,nDimensions+2+H2) = source(gridIndex,nDimensions+2+H2) +      &
                   F * temp * conservedVariables(gridIndex,nDimensions+2+O2)

              source(gridIndex,nDimensions+2+O2) = source(gridIndex,nDimensions+2+O2) +      &
                   F * temp * conservedVariables(gridIndex,nDimensions+2+H2)

           end do
        end do
     end do

     ! Cleanup
     deallocate(dataTemperature)
     deallocate(deltaTemperature)

  end select

  return
end subroutine data_assimilation_adjoint_source


! ===================================================================== !
! Compute sensitivity if data assimilation is a function of the control !
! ===================================================================== !
subroutine data_assimilation_compute_sensitivity

  ! Internal modules
  use data_assimilation

  ! External modules
  use simulation_flags
  use solver_options
  use state_jacobian
  use geometry
  use grid
  use grid_functions, only : inner_product
  use state
  use equation_of_state
  use onestep
  use controller

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: chemicalSource, chemicalSourceData, referenceTemperature,       &
       flameTemperature, activationTemperature 
  real(WP), allocatable :: dataTemperature(:), F(:)

  instantaneousSensitivity = 0.0_WP
  
  if (assimilationType .ne. ASSIMILATE_HEAT_RELEASE .or.                                     &
       controllerType .ne. CHEMICAL_ACTUATOR) return

  allocate(F(nGridPoints)); F = 0.0_WP

  ! Get temperature from the data file
  allocate(dataTemperature(nGridPoints))
  call compute_dependent_variables(data, temperature = dataTemperature)

  referenceTemperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  flameTemperature = referenceTemperature / (1.0_WP - heatRelease)
  activationTemperature = zelDovich / heatRelease * flameTemperature

  do k = functionalPatch%iStart(3), functionalPatch%iEnd(3)
     do j = functionalPatch%iStart(2), functionalPatch%iEnd(2)
        do i = functionalPatch%iStart(1), functionalPatch%iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))

           ! Chemical source from the simulation
           chemicalSource = -Damkohler * stoichiometricCoefficient(H2) *                     &
                conservedVariables(gridIndex,nDimensions+2+H2) *                             &
                conservedVariables(gridIndex,nDimensions+2+O2) *                             &
                exp(- activationTemperature / temperature(gridIndex, 1))

           ! Chemical source from the data
           chemicalSourceData = -Damkohler * stoichiometricCoefficient(H2) *                 &
                data(gridIndex,nDimensions+2+H2) * data(gridIndex,nDimensions+2+O2) *        &
                exp(- activationTemperature / dataTemperature(gridIndex))

           F(gridIndex) = real(power, WP) * (controlParameter * chemicalSource -             &
                chemicalSourceData)**(power-1) * chemicalSource

        end do
     end do
  end do

  instantaneousSensitivity = instantaneousSensitivity + inner_product(F, targetMollifier(:,1))

  ! Cleanup
  deallocate(dataTemperature)
  deallocate(F)

  return
end subroutine data_assimilation_compute_sensitivity
