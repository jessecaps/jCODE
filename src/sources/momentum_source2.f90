module momentum_source

  ! External modules
  use precision

  implicit none

  integer :: nMomentumSources

  type :: t_MomentumSource
     real(WP) :: initialMomentum
  end type t_MomentumSource

  type(t_MomentumSource), allocatable, dimension(:) :: momentumSources
  logical :: forceMomentum(3), forceDensity, forceTemperature, useTimeRamp
  real(WP) :: initialDensity, initialTemperature, timeRampDuration

end module momentum_source


! ========================== !
! Setup the momentum sources !
! ========================== !
subroutine momentum_source_setup

  ! Internal modules
  use momentum_source

  ! External modules
  use parser
  use string
  use geometry

  implicit none

  ! Local variables
  integer :: i, n
  character(len = str_short) :: forceDir

  ! Find the direction to force momentum
  forceMomentum = .false.
  call parser_read('force momentum', forceDir, 'none')
  if (index(forceDir,'x') .ne. 0) then
     forceMomentum(1) = .true.
     if (.not. isPeriodic(1)) call die('momentum_source_setup: &
          & Cannot force momentum in x, direction not periodic')
  end if
  if (index(forceDir,'y') .ne. 0) then
     forceMomentum(2) = .true.
     if (.not. isPeriodic(2)) call die('momentum_source_setup: &
          & Cannot force momentum in y, direction not periodic')
  end if
  if (index(forceDir,'z') .ne. 0) then
     forceMomentum(3) = .true.
     if (.not. isPeriodic(3)) call die('momentum_source_setup: &
          & Cannot force momentum in z, direction not periodic')
  end if
  
  ! Allocate the source term
  n = 0
  do i = 1, nDimensions
     if (forceMomentum(i)) n = n + 1
  end do
  nMomentumSources = n
  if (n .le. 0) return
  allocate(momentumSources(n))

  ! Determine if mass/heat should be corrected
  call parser_read('force density', forceDensity, .false.)
  call parser_read('force temperature', forceTemperature, .false.)

  ! Should we ramp up the source in time?
  call parser_read('use time ramp', useTimeRamp, .false.)
  if (useTimeRamp) call parser_read('time ramp duration', timeRampDuration)

  return
end subroutine momentum_source_setup


! ======================== !
! Compute initial momentum !
! ======================== !
subroutine momentum_source_init

  ! Internal modules
  use momentum_source

  ! External modules
  use parallel
  use parser
  use geometry
  use grid
  use state
  use equation_of_state

  implicit none

  ! Local variables
  integer :: i, n
  real(WP) :: targetMomentum(nMomentumSources)
  logical :: useTarget

  if (.not. allocated(momentumSources)) return

  ! Determine if a target momentum is specified
  call parser_is_defined('target momentum', useTarget)

  if (useTarget) then

     call parser_read('target momentum', targetMomentum)

     n = 1
     do i = 1, nDimensions
        if (forceMomentum(i)) then
           momentumSources(n)%initialMomentum = targetMomentum(i)
           n = n + 1
        end if
     end do

  else

     ! Compute domain-averaged momentum
     n = 1
     do i = 1, nDimensions
        if (forceMomentum(i)) then
           momentumSources(n)%initialMomentum = sum(conservedVariables(:,1+i) * gridNorm(:,1))
           call parallel_sum(momentumSources(n)%initialMomentum)
           momentumSources(n)%initialMomentum = momentumSources(n)%initialMomentum /         &
                globalGridVolume
           n = n + 1
        end if
     end do

  end if

  ! Get initial density and temperature
  if (forceDensity) then
     initialDensity = sum(conservedVariables(:,1) * gridNorm(:,1))
     call parallel_sum(initialDensity); initialDensity = initialDensity / globalGridVolume
  end if
  !if (forceTemperature) then
     call compute_dependent_variables(conservedVariables, temperature = temperature(:,1))
     initialTemperature = sum(temperature(:,1) * gridNorm(:,1))
     call parallel_sum(initialTemperature); initialTemperature = initialTemperature /        &
          globalGridVolume
 ! end if
  
  return
end subroutine momentum_source_init


! ============================ !
! Cleanup the momentum sources !
! ============================ !
subroutine momentum_source_cleanup

  ! Internal modules
  use momentum_source

  implicit none

  if (allocated(momentumSources)) deallocate(momentumSources)

  return
end subroutine momentum_source_cleanup


! =================================================== !
! Add the the momentum sources during the forward run !
! =================================================== !
subroutine momentum_source_forward(source)

  ! Internal modules
  use momentum_source

  ! External modules
  use parallel
  use geometry
  use grid
  use state
  use time_info

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j, n
  real(WP) :: currentMomentum, currentDensity, currentTemperature, source_, timeRamp,        &
       meanKineticEnergy, buf, src(nGridPoints, nUnknowns), meanSrc(2)

  if (.not. allocated(momentumSources)) return

  ! Time ramp
  timeRamp = 1.0_WP
  if (useTimeRamp .and. time.lt.timeRampDuration) then
     timeRamp = tanh(5.0_WP * time / timeRampDuration)
  end if

  ! Add the source term
  n = 1
  meanKineticEnergy = 0.0_WP
  do i = 1, nDimensions
     if (forceMomentum(i)) then
        currentMomentum = sum(conservedVariables(:,1+i) * gridNorm(:,1))
        call parallel_sum(currentMomentum)
        currentMomentum = currentMomentum / globalGridVolume

        source_ = (momentumSources(n)%initialMomentum - currentMomentum) / timeStepSize *    &
             timeRamp

        do j = 1, nGridPoints
           
           ! Momentum treatment
           source(j, i + 1) = source(j, i + 1) + source_

           ! Energy treatment
           source(j, nDimensions+2) = source(j, nDimensions+2) + source_ * velocity(j, i)
           meanKineticEnergy = meanKineticEnergy + source_ * velocity(j, i) * gridNorm(j,1)
        end do

        n = n + 1
     end if
     
  end do

  ! Remove the mean kinetic energy contribution
  call parallel_sum(meanKineticEnergy)
  meanKineticEnergy = meanKineticEnergy / globalGridVolume
  source(:, nDimensions+2) = source(:, nDimensions+2) - meanKineticEnergy

!!$  ! Add mass correction
!!$  if (forceDensity) then
!!$     currentDensity = sum(conservedVariables(:,1) * gridNorm(:,1))
!!$     call parallel_sum(currentDensity)
!!$     currentDensity = currentDensity / globalGridVolume
!!$     source_ = (initialDensity - currentDensity) / timeStepSize * timeRamp
!!$     do j = 1, nGridPoints
!!$        source(j,1) = source(j,1) + source_
!!$        source(j,nDimensions+2) = source(j,nDimensions+2) + source_ * (temperature(j,1) /    &
!!$             ratioOfSpecificHeats - 0.5_WP * sum(velocity(j,:)**2))
!!$       ! buf = buf + source_ * (temperature(j,1) /                                            &
!!$       !      ratioOfSpecificHeats - 0.5_WP * sum(velocity(j,:)**2)) * gridNorm(j,1)
!!$     end do
!!$  end if
!!$
!!$  ! Add temperature correction
!!$  if (forceTemperature) then
!!$     currentTemperature = sum(temperature(:,1) * gridNorm(:,1))
!!$     call parallel_sum(currentTemperature)
!!$     currentTemperature = currentTemperature / globalGridVolume
!!$     do j = 1, nGridPoints
!!$        source_ = conservedVariables(j,1) * (initialTemperature - currentTemperature) /      &
!!$             ratioOfSpecificHeats / timeStepSize * timeRamp
!!$        source(j,nDimensions+2) = source(j,nDimensions+2) + source_
!!$      !  buf = buf + source_ * gridNorm(j,1)
!!$     end do
!!$  end if


  ! Add temperature correction
  if (forceTemperature) then

     ! Recompute sources
     src = 0.0_WP
     call add_fluxes_forward(src)
     call ibm_source(src)
     call add_dissipation(FORWARD, src)

     ! Remove from density and energy
     meanSrc = 0.0_WP
     do j = 1, nGridPoints
        meanSrc(1) = meanSrc(1) + src(j,1) * gridNorm(j,1)
        meanSrc(2) = meanSrc(2) + src(j,nDimensions+2) * gridNorm(j,1)
     end do
     call parallel_sum(meanSrc(1)); meanSrc(1) = meanSrc(1) / globalGridVolume
     call parallel_sum(meanSrc(2)); meanSrc(2) = meanSrc(2) / globalGridVolume
     do j = 1, nGridPoints
        source(j,1) = source(j,1) - meanSrc(1)
        !source(j,nDimensions+2) = source(j,nDimensions+2) - meanSrc(1) * (temperature(j,1) /    &
        !     ratioOfSpecificHeats - 0.5_WP * sum(velocity(j,:)**2))
        source(j,nDimensions+2) = source(j,nDimensions+2) - meanSrc(2)
     end do

     currentTemperature = sum(temperature(:,1) * gridNorm(:,1))
     call parallel_sum(currentTemperature)
     currentTemperature = currentTemperature / globalGridVolume
     do j = 1, nGridPoints
        source_ = conservedVariables(j,1) * (initialTemperature - currentTemperature) /      &
             ratioOfSpecificHeats / timeStepSize * timeRamp
        source(j,nDimensions+2) = source(j,nDimensions+2) + source_
     end do
  end if

  return
end subroutine momentum_source_forward


! =================================================== !
! Add the the momentum sources during the adjoint run !
! =================================================== !
subroutine momentum_source_adjoint(source)

  ! Internal modules
  use momentum_source

  ! External modules
  use solver_options
  use geometry

  implicit none

  ! Arguments
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  if (.not. allocated(momentumSources)) return

  call die('momentum_source_adjoint not yet implemented!')

  return
end subroutine momentum_source_adjoint
