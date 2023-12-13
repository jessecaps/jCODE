module particle_integrator

  ! External modules
  use particle
  use particle_comm

  implicit none

  ! RK4 parameters
  type(t_Particle), allocatable :: particleBuffer1(:), particleBuffer2(:)

end module particle_integrator


! ============================= !
! Setup the particle integrator !
! ============================= !
subroutine particle_integrator_setup

  ! Internal modules
  use particle_integrator

  ! External modules
  use simulation_flags

  implicit none

  if (.not. useParticles) return

  return
end subroutine particle_integrator_setup


! =============================== !
! Cleanup the particle integrator !
! =============================== !
subroutine particle_integrator_cleanup

  ! Internal modules
  use particle_integrator

  implicit none

  if (allocated(particleBuffer1)) deallocate(particleBuffer1)
  if (allocated(particleBuffer2)) deallocate(particleBuffer2)

  return
end subroutine particle_integrator_cleanup


! ========================================== !
! Firat-order Euler scheme for the particles !
! ========================================== !
subroutine particle_substep_euler

  ! Internal modules
  use particle_integrator

  ! External modules
  use math
  use parallel
  use simulation_flags
  use time_info
  use particle_solver, only : minDragTime, minRep, maxRep, minMap, maxMap, minKnp, maxKnp

  implicit none

  ! Local variables
  integer :: i, p
  real(WP), dimension(3) :: dxdt, dudt, dwdt
  real(WP) :: mass, dTdt, dmdt
  real(WP), parameter :: oneHalf = 1.0_WP / 2.0_WP
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP

  if (.not. useParticles) return

  ! Start the particle timer
  call timing_start('particles')

  ! Reset min/max quantities
  minDragTime = huge(1.0_WP)
  minRep = huge(1.0_WP); maxRep = 0.0_WP
  minMap = huge(1.0_WP); maxMap = 0.0_WP
  minKnp = huge(1.0_WP); maxKnp = 0.0_WP

  ! Prepare the particle solver
  call prepare_cartesian_arrays
  call particle_injector_inject
  call compute_particle_collisions

  do p = 1, nParticles

     ! Compute particle right-hand side
     call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

     ! Only update particles if id > 0
     if (particles(p)%id.le.0) cycle

     do i = 1, nDimensions

        ! Update particle position
        particles(p)%position(i) = particles(p)%position(i) + timeStepSize * dxdt(i)

        ! Update particle velocity
        particles(p)%velocity(i) = particles(p)%velocity(i) + timeStepSize * dudt(i)

     end do

     ! Update particle temperature
     particles(p)%temperature = particles(p)%temperature + timeStepSize * dTdt

     ! Update the particle angular velocity
     if (useFriction) then
        particles(p)%angularVelocity = particles(p)%angularVelocity + timeStepSize * dwdt
     end if

     ! Update particle diameter
     if (usePhaseChange) then
        mass = particleDensity*pi/6.0_WP*particles(p)%diameter**3
        mass = mass + timeStepSize * dmdt
        particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
     end if

     ! Remove particles that left the domain and update local grid index
     call check_particle_bounds(particles(p))
     call localize_particle(particles(p))

  end do

  ! Communicate the particles
  call communicate_particles

  ! Get global min/max values
  if (twoWayCoupling) call parallel_min(minDragTime)
  if (nParticlesGlobal .gt. 0) then
     call parallel_min(minRep)
     call parallel_min(minMap)
     call parallel_min(minKnp)
     call parallel_max(maxRep)
     call parallel_max(maxMap)
     call parallel_max(maxKnp)
  else
     minRep = 0.0_WP; maxRep = 0.0_WP
     minMap = 0.0_WP; maxMap = 0.0_WP
     minKnp = 0.0_WP; maxKnp = 0.0_WP
  end if

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_substep_euler


! ============================================== !
! 2nd-order Runge-Kutta scheme for the particles !
! ============================================== !
subroutine particle_substep_rk2(stage)

  ! Internal modules
  use particle_integrator

  ! External modules
  use math
  use parallel
  use simulation_flags
  use time_info
  use particle_solver, only : minDragTime, minRep, maxRep, minMap, maxMap, minKnp, maxKnp

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p
  real(WP), dimension(3) :: dxdt, dudt, dwdt
  real(WP) :: mass, dTdt, dmdt

  if (.not. useParticles) return

  ! Start the particle timer
  call timing_start('particles')

  ! Reset min/max quantities
  minDragTime = huge(1.0_WP)
  minRep = huge(1.0_WP); maxRep = 0.0_WP
  minMap = huge(1.0_WP); maxMap = 0.0_WP
  minKnp = huge(1.0_WP); maxKnp = 0.0_WP

  select case (stage)

  case (1)

     ! Prepare the particle solver
     call prepare_cartesian_arrays
     call particle_injector_inject
     call compute_particle_collisions

     ! Prepare temporary particle vectors
     if (allocated(particleBuffer1)) deallocate(particleBuffer1)
     allocate(particleBuffer1(nParticles))
     particleBuffer1 = particles

     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = particleBuffer1(p)%position(i) +                       &
                0.5_WP * timeStepSize * dxdt(i)

           ! Update particle velocity
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                0.5_WP * timeStepSize * dudt(i)

        end do

        ! Update particle temperature
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             0.5_WP * timeStepSize * dTdt

        ! Update the particle angular velocity
        if (useFriction) then
           particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +               &
                0.5_WP * timeStepSize * dwdt
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + 0.5_WP * timeStepSize * dmdt
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (2)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = particleBuffer1(p)%position(i) +                       &
                timeStepSize * dxdt(i)

           ! Update particle velocity
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                timeStepSize * dudt(i)

        end do

        ! Update particle temperature
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             timeStepSize * dTdt

        ! Update the particle angular velocity
        if (useFriction) then
              particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +            &
                   timeStepSize * dwdt
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + timeStepSize * dmdt
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

     ! Communicate the particles
     call communicate_particles

  end select

  ! Get global min/max values
  if (twoWayCoupling) call parallel_min(minDragTime)
  if (nParticlesGlobal .gt. 0) then
     call parallel_min(minRep)
     call parallel_min(minMap)
     call parallel_min(minKnp)
     call parallel_max(maxRep)
     call parallel_max(maxMap)
     call parallel_max(maxKnp)
  else
     minRep = 0.0_WP; maxRep = 0.0_WP
     minMap = 0.0_WP; maxMap = 0.0_WP
     minKnp = 0.0_WP; maxKnp = 0.0_WP
  end if

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_substep_rk2


! ========================================================== !
! Optimal 3rd-order TVD Runge-Kutta scheme for the particles !
! ========================================================== !
subroutine particle_substep_rk3_tvd(stage)

  ! Internal modules
  use particle_integrator

  ! External modules
  use math
  use parallel
  use simulation_flags
  use time_info
  use particle_solver, only : minDragTime, minRep, maxRep, minMap, maxMap, minKnp, maxKnp

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p
  real(WP), dimension(3) :: dxdt, dudt, dwdt
  real(WP) :: mass, dTdt, dmdt

  if (.not. useParticles) return

  ! Start the particle timer
  call timing_start('particles')

  ! Reset min/max quantities
  minDragTime = huge(1.0_WP)
  minRep = huge(1.0_WP); maxRep = 0.0_WP
  minMap = huge(1.0_WP); maxMap = 0.0_WP
  minKnp = huge(1.0_WP); maxKnp = 0.0_WP

  select case (stage)

  case (1)

     ! Prepare the particle solver
     call prepare_cartesian_arrays
     call particle_injector_inject
     call compute_particle_collisions

     ! Prepare temporary particle vectors
     if (allocated(particleBuffer1)) deallocate(particleBuffer1)
     allocate(particleBuffer1(nParticles))
     particleBuffer1 = particles

     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = particleBuffer1(p)%position(i) +                       &
                timeStepSize * dxdt(i)

           ! Update particle velocity
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                timeStepSize * dudt(i)

        end do

        ! Update particle temperature
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             timeStepSize * dTdt

        ! Update the particle angular velocity
        if (useFriction) then
           particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +               &
                timeStepSize * dwdt
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + timeStepSize * dmdt
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (2)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = 0.75_WP * particleBuffer1(p)%position(i) +             &
                0.25_WP * particles(p)%position(i) + 0.25_WP *  timeStepSize * dxdt(i)

           ! Update particle velocity
           particles(p)%velocity(i) = 0.75_WP * particleBuffer1(p)%velocity(i) +             &
                0.25_WP * particles(p)%velocity(i) + 0.25_WP * timeStepSize * dudt(i)

        end do

        ! Update particle temperature
        particles(p)%temperature = 0.75_WP * particleBuffer1(p)%temperature +                &
             0.25_WP * particles(p)%temperature + 0.25_WP *  timeStepSize * dTdt

        ! Update the particle angular velocity
        if (useFriction) then
           particles(p)%angularVelocity = 0.75_WP * particleBuffer1(p)%angularVelocity +     &
                0.25_WP * particles(p)%angularVelocity + 0.25_WP * timeStepSize * dwdt
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = 0.75_WP * particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3 +       &
                0.25_WP * particleDensity*pi/6.0_WP*particles(p)%diameter**3 +               &
                0.25_WP * timeStepSize * dmdt
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (3)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = 1.0_WP / 3.0_WP * particleBuffer1(p)%position(i) +     &
                2.0_WP / 3.0_WP * (particles(p)%position(i) + timeStepSize * dxdt(i))

           ! Update particle velocity
           particles(p)%velocity(i) = 1.0_WP / 3.0_WP * particleBuffer1(p)%velocity(i) +     &
                2.0_WP / 3.0_WP * (particles(p)%velocity(i) + timeStepSize * dudt(i))

        end do

        ! Update particle temperature
        particles(p)%temperature = 1.0_WP / 3.0_WP * particleBuffer1(p)%temperature +        &
             2.0_WP / 3.0_WP * (particles(p)%temperature + timeStepSize * dTdt)

        ! Update the particle angular velocity
        if (useFriction) then
           particles(p)%angularVelocity = 1.0_WP / 3.0_WP *                                  &
                particleBuffer1(p)%angularVelocity + 2.0_WP / 3.0_WP *                       &
                (particles(p)%angularVelocity + timeStepSize * dwdt)
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = 1.0_WP / 3.0_WP * particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3+&
                2.0_WP / 3.0_WP * (particleDensity*pi/6.0_WP*particles(p)%diameter**3 +      &
                timeStepSize * dmdt)
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

     ! Communicate the particles
     call communicate_particles

  end select

  ! Get global min/max values
  if (twoWayCoupling) call parallel_min(minDragTime)
  if (nParticlesGlobal .gt. 0) then
     call parallel_min(minRep)
     call parallel_min(minMap)
     call parallel_min(minKnp)
     call parallel_max(maxRep)
     call parallel_max(maxMap)
     call parallel_max(maxKnp)
  else
     minRep = 0.0_WP; maxRep = 0.0_WP
     minMap = 0.0_WP; maxMap = 0.0_WP
     minKnp = 0.0_WP; maxKnp = 0.0_WP
  end if

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_substep_rk3_tvd


! ============================================== !
! 4th-order Runge-Kutta scheme for the particles !
! ============================================== !
subroutine particle_substep_rk4(stage)

  ! Internal modules
  use particle_integrator

  ! External modules
  use math
  use parallel
  use simulation_flags
  use time_info
  use particle_solver, only : minDragTime, minRep, maxRep, minMap, maxMap, minKnp, maxKnp

  implicit none

  ! Arguments
  integer, intent(in) :: stage

  ! Local variables
  integer :: i, p
  real(WP), dimension(3) :: dxdt, dudt, dwdt
  real(WP) :: mass, dTdt, dmdt
  real(WP), parameter :: oneHalf = 1.0_WP / 2.0_WP
  real(WP), parameter :: oneThird = 1.0_WP / 3.0_WP
  real(WP), parameter :: oneSixth = 1.0_WP / 6.0_WP

  if (.not. useParticles) return

  ! Start the particle timer
  call timing_start('particles')

  ! Reset min/max quantities
  minDragTime = huge(1.0_WP)
  minRep = huge(1.0_WP); maxRep = 0.0_WP
  minMap = huge(1.0_WP); maxMap = 0.0_WP
  minKnp = huge(1.0_WP); maxKnp = 0.0_WP

  select case (stage)

  case (1)

     ! Prepare the particle solver
     call prepare_cartesian_arrays
     call particle_injector_inject
     call compute_particle_collisions

     ! Prepare temporary particle vectors
     if (allocated(particleBuffer1)) deallocate(particleBuffer1)
     if (allocated(particleBuffer2)) deallocate(particleBuffer2)
     allocate(particleBuffer1(nParticles))
     allocate(particleBuffer2(nParticles))
     particleBuffer1 = particles

     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particleBuffer2(p)%position(i) = particles(p)%position(i) +                       &
                timeStepSize * dxdt(i) * oneSixth
           particles(p)%position(i) = particleBuffer1(p)%position(i) +                       &
                timeStepSize * dxdt(i) * oneHalf

           ! Update particle velocity
           particleBuffer2(p)%velocity(i) = particles(p)%velocity(i) +                       &
                timeStepSize * dudt(i) * oneSixth
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                timeStepSize * dudt(i) * oneHalf

        end do

        ! Update particle temperature
        particleBuffer2(p)%temperature = particles(p)%temperature +                          &
             timeStepSize * dTdt * oneSixth
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             timeStepSize * dTdt * oneHalf

        ! Update the particle angular velocity
        if (useFriction) then
           particleBuffer2(p)%angularVelocity = particles(p)%angularVelocity +               &
                timeStepSize * dwdt * oneSixth
           particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +               &
                timeStepSize * dwdt * oneHalf
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particles(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneSixth
           particleBuffer2(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneHalf
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (2)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particleBuffer2(p)%position(i) = particleBuffer2(p)%position(i) +                 &
                timeStepSize * dxdt(i) * oneThird
           particles(p)%position(i) = particleBuffer1(p)%position(i) +                       &
                timeStepSize * dxdt(i) * oneHalf

           ! Update particle velocity
           particleBuffer2(p)%velocity(i) = particleBuffer2(p)%velocity(i) +                 &
                timeStepSize * dudt(i) * oneThird
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                timeStepSize * dudt(i) * oneHalf

        end do

        ! Update particle temperature
        particleBuffer2(p)%temperature = particleBuffer2(p)%temperature +                    &
             timeStepSize * dTdt * oneThird
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             timeStepSize * dTdt * oneHalf

        ! Update the particle angular velocity
        if (useFriction) then
           particleBuffer2(p)%angularVelocity = particleBuffer2(p)%angularVelocity +         &
                timeStepSize * dwdt * oneThird
           particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +               &
                timeStepSize * dwdt * oneHalf
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer2(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneThird
           particleBuffer2(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneHalf
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (3)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particleBuffer2(p)%position(i) = particleBuffer2(p)%position(i) +                 &
                timeStepSize * dxdt(i) * oneThird
           particles(p)%position(i) = particleBuffer1(p)%position(i) + timeStepSize * dxdt(i)

           ! Update particle velocity
           particleBuffer2(p)%velocity(i) = particleBuffer2(p)%velocity(i) +                 &
                timeStepSize * dudt(i) * oneThird
           particles(p)%velocity(i) = particleBuffer1(p)%velocity(i) +                       &
                timeStepSize * dudt(i)

        end do

        ! Update particle temperature
        particleBuffer2(p)%temperature = particleBuffer2(p)%temperature +                    &
             timeStepSize * dTdt * oneThird
        particles(p)%temperature = particleBuffer1(p)%temperature +                          &
             timeStepSize * dTdt

        ! Update the particle angular velocity
        if (useFriction) then
           particleBuffer2(p)%angularVelocity = particleBuffer2(p)%angularVelocity +         &
                timeStepSize * dwdt * oneThird
           particles(p)%angularVelocity = particleBuffer1(p)%angularVelocity +               &
                timeStepSize * dwdt
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer2(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneThird
           particleBuffer2(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
           mass = particleDensity*pi/6.0_WP*particleBuffer1(p)%diameter**3
           mass = mass + timeStepSize * dmdt
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

  case (4)

     call prepare_cartesian_arrays
     do p = 1, nParticles

        ! Compute particle right-hand side
        call particle_rhs(particles(p), dxdt, dudt, dwdt, dTdt, dmdt)

        ! Only update particles if id > 0
        if (particles(p)%id.le.0) cycle

        do i = 1, nDimensions

           ! Update particle position
           particles(p)%position(i) = particleBuffer2(p)%position(i) +                       &
                timeStepSize * dxdt(i) * oneSixth

           ! Update particle velocity
           particles(p)%velocity(i) = particleBuffer2(p)%velocity(i) +                       &
                timeStepSize * dudt(i) * oneSixth

        end do

        ! Update particle temperature
        particles(p)%temperature = particleBuffer2(p)%temperature +                          &
             timeStepSize * dTdt * oneSixth

        ! Update the particle angular velocity
        if (useFriction) then
              particles(p)%angularVelocity = particleBuffer2(p)%angularVelocity +            &
                   timeStepSize * dwdt * oneSixth
        end if

        ! Update particle diameter
        if (usePhaseChange) then
           mass = particleDensity*pi/6.0_WP*particleBuffer2(p)%diameter**3
           mass = mass + timeStepSize * dmdt * oneSixth
           particles(p)%diameter = (6.0_WP/(pi*particleDensity)*mass)**(1.0_WP/3.0_WP)
        end if

        ! Remove particles that left the domain and update local grid index
        call check_particle_bounds(particles(p))
        call localize_particle(particles(p))

     end do

     ! Communicate the particles
     call communicate_particles

  end select

  ! Get global min/max values
  if (twoWayCoupling) call parallel_min(minDragTime)
  if (nParticlesGlobal .gt. 0) then
     call parallel_min(minRep)
     call parallel_min(minMap)
     call parallel_min(minKnp)
     call parallel_max(maxRep)
     call parallel_max(maxMap)
     call parallel_max(maxKnp)
  else
     minRep = 0.0_WP; maxRep = 0.0_WP
     minMap = 0.0_WP; maxMap = 0.0_WP
     minKnp = 0.0_WP; maxKnp = 0.0_WP
  end if

  ! Stop the particle timer
  call timing_stop('particles')

  return
end subroutine particle_substep_rk4


! ============================================== !
! Compute the CFL based on the particle velocity !
! ============================================== !
subroutine get_particle_cfl(dt)

  ! Internal modules
  use particle_integrator

  ! External modules
  use grid
  use state
  use time_info
  use particle_solver

  implicit none

  ! Arguments
  real(WP), intent(in) :: dt

  ! Local variables
  integer :: i, j, p
  real(WP) :: localWaveSpeed

  if (.not. useParticles) return

  particleCfl = 0.0_WP

  ! Advection
  do p = 1, nParticles
     i = particles(p)%gridIndex(1) - gridOffset(1) + localGridSize(1) *                      &
          (particles(p)%gridIndex(2) - 1 - gridOffset(2) + localGridSize(2) *                &
          (particles(p)%gridIndex(3) - 1 - gridOffset(3)))
     localWaveSpeed = 0.0_WP
     do j = 1, nDimensions
        localWaveSpeed = localWaveSpeed +                                                    &
             abs(dot_product(particles(p)%velocity(1:nDimensions),                           &
             metrics(i,1+nDimensions*(j-1):nDimensions*j)))
     end do
     localWaveSpeed = jacobian(i,1) * localWaveSpeed
     particleCfl = max(particleCfl, localWaveSpeed * dt)
  end do

  ! Collisions
  if (collisionsOn) then
     do p = 1, nParticles
        particleCfl = max(particleCfl,                                                       &
             10.0_WP * sqrt(sum(particles(p)%velocity(1:nDimensions) **2)) *                 &
             dt / particles(p)%diameter)
     end do
  end if

  ! Drag
  if (twoWayCoupling) then
     particleCfl = max(particleCfl, 0.5_WP * dt / minDragTime)
  end if

  return
end subroutine get_particle_cfl


! ======================================================== !
! Compute the timestep size based on the particle velocity !
! ======================================================== !
subroutine get_particle_timestep_size

  ! Internal modules
  use particle_integrator

  ! External modules
  use grid
  use state
  use time_info
  use particle_solver

  implicit none

  ! Local variables
  integer :: i, j, p
  real(WP) :: localWaveSpeed

  if (.not. useConstantCfl .or. .not. useParticles) return

  ! Advection
  do p = 1, nParticles
     i = particles(p)%gridIndex(1) - gridOffset(1) + localGridSize(1) *                      &
          (particles(p)%gridIndex(2) - 1 - gridOffset(2) + localGridSize(2) *                &
          (particles(p)%gridIndex(3) - 1 - gridOffset(3)))
     localWaveSpeed = 0.0_WP
     do j = 1, nDimensions
        localWaveSpeed = localWaveSpeed +                                                    &
             abs(dot_product(particles(p)%velocity(1:nDimensions),                           &
             metrics(i,1+nDimensions*(j-1):nDimensions*j)))
     end do
     localWaveSpeed = jacobian(i,1) * localWaveSpeed
     timeStepSize = min(timeStepSize, inputCFL / localWaveSpeed)
  end do

  ! Collisions
  if (collisionsOn) then
     do p = 1, nParticles
        timeStepSize =                                                                       &
             min(timeStepSize, 0.1_WP * inputCFL * particles(p)%diameter /                   &
             (sqrt(sum(particles(p)%velocity(1:nDimensions) **2)) + epsilon(1.0_WP)))
     end do
  end if

  ! Drag
  if (twoWayCoupling) timeStepSize = min(timeStepSize, 2.0_WP * minDragTime)

  return
end subroutine get_particle_timestep_size
