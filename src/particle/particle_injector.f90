module particle_injector

  ! External modules
  use particle
  use particle_exchange

  implicit none

  ! Injection type
  integer :: injectionType
  integer, parameter ::                                                                      &
       NO_INJECT   = 0,                                                                      &
       BULK_INJECT = 1

  ! Global variables
  integer :: nParticlesMax
  real(WP) :: dMean, dMin, dMax, dShift, dStd
  real(WP) :: particleMassFlowRate, injectionVelocity(3), injectionPosition(3,2),            &
       injectionTemperature, injectionRadius, injectionVelocityStd(3)
  logical :: useInjectionVelocity, useRadialInjection

contains

  ! Compute particle diameter
  function get_diameter() result(dp)

    ! External modules
    use math
    use random

    implicit none

    ! Arguments
    real(WP) :: dp

    dp=random_lognormal(m = dMean - dShift, sd = dStd) + dShift
    do while (dp.gt.dMax+epsilon(1.0_WP).or.dp.lt.dMin-epsilon(1.0_WP))
       dp=random_lognormal(m = dMean - dShift, sd = dStd) + dShift
    end do

  end function get_diameter

end module particle_injector


! =========================== !
! Setup the particle injector !
! =========================== !
subroutine particle_injector_setup

  ! Internal modules
  use particle_injector

  ! External modules
  use parser
  use string
  use parallel
  use solver_options

  implicit none

  ! Local variables
  character(len=str_medium) :: val

  ! Maximum allowable number of particles
  call parser_read('maximum number of particles', nParticlesMax, huge(nParticlesMax))

  ! Set number of particles in (for monitoring)
  nParticlesIn = 0

  call parser_read('particle injection', val, 'none')
  select case (trim(val))
  case ('none', 'NONE', 'None')

     ! Set type
     injectionType = NO_INJECT

  case ('bulk', 'Bulk', 'BULK')

     ! Set type
     injectionType = BULK_INJECT

     ! Mass flow rate to impose
     call parser_read('particle mass flow rate', particleMassFlowRate)

     ! Injection velocity
     call parser_is_defined('particle u velocity', useInjectionVelocity)
     if (useInjectionVelocity) then
        call parser_read('particle u velocity', injectionVelocity(1))
        call parser_read('particle v velocity', injectionVelocity(2), 0.0_WP)
        call parser_read('particle w velocity', injectionVelocity(3), 0.0_WP)
        call parser_read('particle u std', injectionVelocityStd(1), 0.0_WP)
        call parser_read('particle v std', injectionVelocityStd(2), injectionVelocityStd(1))
        call parser_read('particle w std', injectionVelocityStd(3), injectionVelocityStd(1))
     end if

     ! Injection temperature
     call parser_read('particle temperature', injectionTemperature, 2.5_WP)

     ! Read particle diameter from input file
     call parser_read('particle mean diameter', dMean)
     call parser_read('particle std diameter', dStd, 0.0_WP)
     call parser_read('particle min diameter', dMin, 0.0_WP)
     call parser_read('particle max diameter', dMax, 0.0_WP)
     call parser_read('particle diameter shift', dShift, 0.0_WP)
     if (dStd.le.epsilon(1.0_WP)) then
        dMin = dMean
        dMax = dMean
     end if
     
     ! Read min/max particle position to inject
     call parser_read('particle xmin', injectionPosition(1,1), domainExtent(1,1))
     call parser_read('particle xmax', injectionPosition(1,2), domainExtent(1,2))
     call parser_is_defined('injection radius', useRadialInjection)
     if (useRadialInjection) then
        call parser_read('injection radius', injectionRadius)
        call parser_read('injection y0', injectionPosition(2,1), 0.0_WP)
        call parser_read('injection z0', injectionPosition(3,1), 0.0_WP)
     else
        call parser_read('particle ymin', injectionPosition(2,1), domainExtent(2,1))
        call parser_read('particle ymax', injectionPosition(2,2), domainExtent(2,2))
        call parser_read('particle zmin', injectionPosition(3,1), domainExtent(3,1))
        call parser_read('particle zmax', injectionPosition(3,2), domainExtent(3,2))
     end if

  case default

     call die('particle_injector_setup: unknown particle injection type.')

  end select

  return
end subroutine particle_injector_setup


! ============================= !
! Cleanup the particle injector !
! ============================= !
subroutine particle_injector_cleanup

  ! Internal modules
  use particle_injector

  implicit none

  return
end subroutine particle_injector_cleanup


! ==================== !
! Inject new particles !
! ==================== !
subroutine particle_injector_inject

  ! Internal modules
  use particle_injector

  implicit none

  ! Return if exceeded particle count
  if (nParticlesGlobal .ge. nParticlesMax) return

  ! Choose injection type
  select case (injectionType)

  case (NO_INJECT)

     ! No injection
     return

  case (BULK_INJECT)

     ! Uniform particle injection
     call particle_injector_bulk

  end select

  return
end subroutine particle_injector_inject


! =========================== !
! Add new particles uniformly !
! =========================== !
subroutine particle_injector_bulk

  ! Internal modules
  use particle_injector

  ! External modules
  use math
  use random
  use time_info
  use parallel
  use state
  use particle_comm
  use particle_io, only : communicate_globally

  implicit none

  ! Local variables
  integer :: count, n, nParticles0, nParticlesTmp, ijk(3), jp, ierror
  integer(kind=8) :: localMaxid, maxid
  real(WP) :: Mgoal, Madded, Mtmp, buf
  real(WP), save :: previousError = 0.0_WP

  ! Initial number of particles
  nParticles0 = nParticles
  nParticlesIn = 0

  ! Get the particle mass that should be added to the system
  Mgoal  = particleMassFlowRate * timeStepSize + previousError
  Madded = 0.0_WP

  ! Determine id to assign to particle
  localMaxid = 0
  do n = 1, nParticles
     localMaxid = max(localMaxid, particles(n)%id)
  end do
  call MPI_ALLREDUCE(localMaxid, maxid, 1, MPI_INTEGER8, MPI_MAX, comm, ierror)

  ! Add new particles until desired mass is achieved
  do while (Madded .lt. Mgoal)

     if (iRank .eq. iRoot) then
        ! Initialize parameters
        Mtmp = 0.0_WP
        nParticlesTmp = 0

        ! Loop while the added volume is not sufficient
        new_part: do while (Mtmp .lt. Mgoal - Madded)

           ! Increment counter
           nParticlesTmp = nParticlesTmp + 1
           count = nParticles0 + nParticlesTmp

           ! Create space for new particle
           call resize_particles(particles, count)

           ! Generate a diameter
           particles(count)%diameter = get_diameter()

           ! Set various parameters for the particle
           particles(count)%id              = maxid + int(nParticlesTmp, 8)
           particles(count)%temperature     = injectionTemperature
           particles(count)%collision       = 0.0_WP
           particles(count)%torque          = 0.0_WP
           particles(count)%velocity        = 0.0_WP
           particles(count)%angularVelocity = 0.0_WP
           
           ! Make it an "official" particle
           particles(count)%stop = 0

           ! Give a position at the injector to the particle
           call get_position(particles(count)%position, 0.6_WP * particles(count)%diameter)

           ! Check for collisions
           if (collisionsOn) then
              ! Loop over previously injected particles
              do n = 1, nParticlesTmp - 1
                 jp = n + nParticles0
                 buf = 0.5_WP * (particles(count)%diameter + particles(jp)%diameter) -       &
                      sqrt(sum((particles(count)%position - particles(jp)%position)**2))
                 if (buf .gt. 0.0_WP) then
                    ! Reset, try again
                    nParticlesTmp = nParticlesTmp - 1
                    cycle new_part
                 end if
              end do
           end if

           ! Localize the particle
           ijk = iStart
           do n = 1, nDimensions
              call bisection(particles(count)%position(n), ijk(n),                           &
                   cartesianCoordinates(1:globalGridSize(n), n), 1, globalGridSize(n))
           end do
           particles(count)%gridIndex = ijk

           ! Give it a velocity
           if (useInjectionVelocity) then
              do n = 1, nDimensions
                 particles(count)%velocity(n) = random_normal(injectionVelocity(n),          &
                      injectionVelocityStd(n))
              end do
           end if

           ! Update the added mass for the timestep
           Mtmp = Mtmp + particleDensity * Pi / 6.0_WP * particles(count)%diameter**3
        end do new_part
     end if

     ! Communicate particles
     call communicate_globally

     ! Loop through newly created particles
     buf = 0.0_WP
     do n = nParticles0 + 1, nParticles
        if (particles(n)%stop.eq.0) then
           ! Interpolate fluid velocity to particle position
           if (.not. useInjectionVelocity) then
              call interpolate_fluid_to_particle(particles(n)%gridIndex,                     &
                   particles(n)%position, velocity = particles(n)%velocity)
           end if
           
           ! Increment counter
           nParticlesIn = nParticlesIn + 1
           
           ! Update the added mass for the timestep
           buf = buf + particleDensity * Pi / 6.0_WP * particles(n)%diameter**3
           
           ! Update the max particle id
           maxid = max(maxid, particles(n)%id)
        end if
     end do

     ! Total mass added
     call parallel_sum(buf, Mtmp); Madded = Madded + Mtmp

     ! Clean up particles
     call recycle_particles(particles)

     ! Update initial npart
     nParticles0 = nParticles

     ! Maximum particle id
     call MPI_ALLREDUCE(maxid, localMaxid, 1, MPI_INTEGER8, MPI_MAX, comm, ierror)
     maxid = localMaxid
  end do

  ! Remember the error
  previousError = Mgoal - Madded

  ! Sum up injected particles
  call parallel_sum(nParticlesIn)

contains

  ! Position for bulk injection of particles
  subroutine get_position(position, mind)

    ! Internal modules
    use particle_injector

    ! External modules
    use random
    use math

    implicit none

    ! Arguments
    real(WP), intent(out) :: position(3)
    real(WP), intent(in) :: mind

    ! Local variables
    integer :: i
    real(WP) :: r, theta, rand

    ! Random position within a specified region
    if (useRadialInjection) then
       call random_number(rand)
       r = sqrt(rand) * injectionRadius ! sqrt(rand) avoids accumulation near the center
       call random_number(rand)
       theta = rand * twoPi
       position(1) = injectionPosition(1,1) + mind + rand *                                  &
            (injectionPosition(1,2) - injectionPosition(1,1) - 2.0_WP * mind)
       position(2) = r * cos(theta) + injectionPosition(2,1)
       if (nDimensions.gt.2) position(3) = r * sin(theta) + injectionPosition(3,1)
    else
       do i = 1, nDimensions
          call random_number(rand)
          position(i) = injectionPosition(i,1) + mind + rand *                               &
               (injectionPosition(i,2) - injectionPosition(i,1) - 2.0_WP * mind)
       end do
    end if

    return
  end subroutine get_position

end subroutine particle_injector_bulk
