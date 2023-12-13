module time_info

  ! External modules
  use precision

  implicit none

  integer :: timestep, timeStage
  real(WP) :: timeStepSize
  real(WP) :: time, adjointCoefficientTime
  real(WP) :: cfl, acousticCfl, viscousCfl, particleCfl, ibmCfl

contains

  ! Compute the CFL based on the input timestep
  ! -------------------------------------------
  subroutine get_cfl

    ! External modules
    use parallel
    use simulation_flags
    use grid
    use grid_levelset
    use state

    implicit none

    ! Local variables
    integer :: i, j
    real(WP) :: dt, localSpeedOfSound, localWaveSpeed, maxViscosity

    ! Get the time step
    if (useConstantCfl) then
       dt = timeStepSize
    else
       dt = inputTimeStepSize
    end if

    ! Reset various CFLs
    acousticCfl = 0.0_WP
    viscousCfl = 0.0_WP
    particleCfl = 0.0_WP
    ibmCfl = 0.0_WP

    ! Advection
    do i = 1, nGridPoints
       if (allocated(levelset)) then
          if (levelset(i,1) .lt. 0.0_WP) cycle ! ... skip IB points
       end if
       localSpeedOfSound = sqrt(ratioOfSpecificHeats * pressure(i,1) * specificVolume(i,1))
       localWaveSpeed = 0.0_WP
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed +                                                  &
               sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
       end do
       localWaveSpeed = localSpeedOfSound * sqrt(localWaveSpeed)
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed + abs(dot_product(velocity(i,:),                   &
               metrics(i,1+nDimensions*(j-1):nDimensions*j)))
       end do
       localWaveSpeed = jacobian(i,1) * localWaveSpeed
       acousticCfl = max(acousticCfl, localWaveSpeed * dt)
    end do
    
    ! Diffusion
    if (useViscosity) then
       do i = 1, nGridPoints
          if (allocated(levelset)) then
             if (levelset(i,1) .lt. 0.0_WP) cycle ! ... skip IB points
          end if
          localWaveSpeed = 0.0_WP
          maxViscosity = max(2.0_WP * dynamicViscosity(i,1),                                 &
               max(secondCoefficientOfViscosity(i,1), thermalDiffusivity(i,1)))
          if (nSpecies .gt. 0) maxViscosity = max(maxViscosity, massDiffusivity(i,1))
          localWaveSpeed = jacobian(i,1) ** 2 * sum(metrics(i,:) ** 2) * maxViscosity
          viscousCfl = max(viscousCfl, localWaveSpeed * dt)
       end do
    end if

    ! Particles
    if (useParticles) call get_particle_cfl(dt)

    ! IBM
    if (useIBM) call get_ibm_cfl(dt)

    ! Get global max
    call parallel_max(acousticCfl)
    call parallel_max(viscousCfl)
    call parallel_max(particleCfl)
    call parallel_max(ibmCfl)

    ! Set the CFL
    if (useConstantCfl) then
       cfl = inputCFL
    else
       cfl = max(acousticCfl, viscousCfl, particleCfl, ibmCfl)
    end if

    return
  end subroutine get_cfl


  ! Compute the timestep size based on the input CFL
  ! ------------------------------------------------
  subroutine get_timestep_size(timeStepSizeOld)

    ! External modules
    use parallel
    use simulation_flags
    use grid
    use grid_levelset
    use state

    implicit none

    ! Arguments
    real(WP), intent(in), optional :: timeStepSizeOld

    ! Local variables
    integer :: i, j
    real(WP) :: localSpeedOfSound, localWaveSpeed, maxViscosity
    real(WP), parameter :: alpha = 0.7_WP

    if (useConstantCfl) then

       timeStepSize = huge(0.0_WP)

       ! Advection
       do i = 1, nGridPoints
          if (allocated(levelset)) then
             if (levelset(i,1) .lt. 0.0_WP) cycle ! ... skip IB points
          end if
          localSpeedOfSound = sqrt(ratioOfSpecificHeats * pressure(i,1) * specificVolume(i,1))
          localWaveSpeed = 0.0_WP
          do j = 1, nDimensions
             localWaveSpeed = localWaveSpeed +                                               &
                  sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
          end do
          localWaveSpeed = localSpeedOfSound * sqrt(localWaveSpeed)
          do j = 1, nDimensions
             localWaveSpeed = localWaveSpeed + abs(dot_product(velocity(i,:),                &
                  metrics(i,1+nDimensions*(j-1):nDimensions*j)))
          end do
          localWaveSpeed = jacobian(i,1) * localWaveSpeed
          timeStepSize = min(timeStepSize, inputCFL / localWaveSpeed)
       end do

       ! Diffusion
       if (useViscosity) then
          do i = 1, nGridPoints
             if (allocated(levelset)) then
                if (levelset(i,1) .lt. 0.0_WP) cycle ! ... skip IB points
             end if
             localWaveSpeed = 0.0_WP
             maxViscosity = max(2.0_WP * dynamicViscosity(i,1),                              &
                  max(secondCoefficientOfViscosity(i,1), thermalDiffusivity(i,1)))
             if (nSpecies .gt. 0) maxViscosity = max(maxViscosity, massDiffusivity(i,1))
             localWaveSpeed = jacobian(i,1) ** 2 * sum(metrics(i,:) ** 2) * maxViscosity
             timeStepSize = min(timeStepSize, inputCFL / localWaveSpeed)
          end do
       end if

       ! Particles
       if (useParticles) call get_particle_timestep_size

       ! IBM
       if (useIBM) call get_ibm_timestep_size

       ! Communicate
       call parallel_min(timeStepSize)

    else

       timeStepSize = inputTimeStepSize

    end if

    if (present(timeStepSizeOld)) then
       if (timeStepSize .gt. timeStepSizeOld) timeStepSize = alpha * timeStepSize +          &
            (1.0_WP - alpha) * timeStepSizeOld
    end if

    return
  end subroutine get_timestep_size

end module time_info
