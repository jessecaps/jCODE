module shock_tube

  ! External modules
  use precision
  use string
  use parser
  use random
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Global variables
  integer :: nx, ny, nz
  real(WP) :: Lx, Ly, Lz, dx, dy, dz

end module shock_tube


subroutine shock_tube_grid

  ! Internal modules
  use shock_tube

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Ly', Ly, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Compute the grid spacing
  if (periodicityType(1) .eq. PLANE) then
     dx = Lx / real(nx, WP)
  else
     dx = Lx / real(nx-1, WP)
  end if
  if (periodicityType(2) .eq. PLANE) then
     dy = Ly / real(ny, WP)
  else
     dy = Ly / real(ny-1, WP)
  end if
  if (periodicityType(3) .eq. PLANE) then
     dz = Lz / real(nz, WP)
  else
     dz = Lz / real(nz-1, WP)
  end if

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create x
           coordinates(grid_index(i,j,k), 1) = Lx * real(i-1, WP) / real(nx-1, WP)

           ! Uniform and periodic in y
           coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) /                &
                real(ny - 1, WP) - 0.5_WP * Ly

           ! Uniform and periodic in z
           if (nz .gt. 1) then
              coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /             &
                   real(nz - 1, WP) - 0.5_WP * Lz
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine shock_tube_grid


subroutine shock_tube_data

  ! Internal modules
  use shock_tube

  ! External modules
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, x, x0, p1, p2, rho1, rho2, u1, u2, density, pressure, velocity, Yv, sigma
  logical :: smoothShock, useVapor

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  call parser_is_defined('vapor mass fraction', useVapor)
  if (useVapor) nUnknowns = nUnknowns + 1
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read in shock parameters
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('shock position', x0)
  call parser_read('post-shock pressure', p1)
  call parser_read('pre-shock pressure', p2)
  call parser_read('post-shock density', rho1)
  call parser_read('pre-shock density', rho2)
  call parser_read('post-shock velocity', u1, 0.0_WP)
  call parser_read('pre-shock velocity', u2, 0.0_WP)
  call parser_read('smooth shock', smoothShock, .false.)

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Get local coordinate
           x = coordinates(grid_index(i,j,k), 1)

           ! Set the density & pressure
           if (smoothShock) then
              ! Smooth transition
              sigma=2.0_WP/dx
              density = rho1 + 0.5_WP * (rho2 - rho1) * (1.0_WP + tanh(sigma * (x-x0)))
              pressure = p1 + 0.5_WP * (p2 - p1) * (1.0_WP + tanh(sigma * (x-x0)))
              velocity = u1 + 0.5_WP * (u2 - u1) * (1.0_WP + tanh(sigma * (x-x0)))
           else
              ! Sharp (discontinuous) transition
              if (x.le.x0) then
                 density = rho1
                 pressure = p1
                 velocity = u1
              else
                 density = rho2
                 pressure = p2
                 velocity = u2
              end if
           end if

           ! Set the conserved variables
           conservedVariables(grid_index(i,j,k), 1) = density
           conservedVariables(grid_index(i,j,k), 2) = density * velocity
           conservedVariables(grid_index(i,j,k), 3:nDimensions+1) = 0.0_WP
           conservedVariables(grid_index(i,j,k), nDimensions+2) = pressure /                 &
                (gamma - 1.0_WP) + 0.5_WP * density * velocity ** 2
        end do
     end do
  end do

  ! Set vapor mass fraction if used
  if (useVapor) then
     call parser_read('vapor mass fraction', Yv)
     conservedVariables(:, nDimensions+3) = conservedVariables(:, 1) * Yv
  end if

  return
end subroutine shock_tube_data

subroutine shock_tube_particles

  ! Internal modules
  use shock_tube

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, j, k, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(WP) :: volume, volumeFraction, particleVolume, particleThickness, std,                &
       volumeFactor, sumParticleVolume, Lp, Lpy, Lpz, x0, rand, r, distance(3),              &
       dMean, dStd, dMin, dMax, dShift, Tp
  character(len = str_medium) :: filename, particleDistribution
  logical :: success, preventOverlap

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return 

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('particle layer thickness', particleThickness, Ly) 
  call parser_read('particle layer position', x0, 0.0_WP)
  call parser_read('particle temperature', Tp, 2.5_WP)
  call parser_read('particle layer std', std, 0.0_WP)

  ! Read particle diameter
  call parser_read('particle mean diameter', dMean)
  call parser_read('particle std diameter', dStd, 0.0_WP)
  call parser_read('particle min diameter', dMin, 0.0_WP)
  call parser_read('particle max diameter', dMax, 0.0_WP)
  call parser_read('particle diameter shift', dShift, 0.0_WP)

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')
  call parser_read('prevent particle overlap', preventOverlap, .true.)
  if (preventOverlap .and. nProcs.gt.1)                                                      &
       call die('Prevent particle overlap only implemented in serial!')

  ! Initial Volume that particles occupy
  if (nz.gt.1) then
     volume = particleThickness * Ly * Lz
  else
     volume = particleThickness * Ly
  end if

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * dMean ** nDimensions

  select case (trim(particleDistribution))

  case ('random')

     ! Get the number of particles based on the volume fraction and distribute to processors
     nParticlesGlobal = int(volumeFraction * volume / particleVolume)

  case ('uniform')

     ! Mean interparticle distance
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     npartx = int(particleThickness / Lp)
     Lp = particleThickness / real(npartx, WP)
     if (ny .gt. 1) nparty = int(Ly / Lp)
     Lpy = Ly / real(nparty, WP)
     if (nz .gt. 1) npartz = int(Lz / Lp)
     Lpz = Lz / real(npartz, WP)
     nParticlesGlobal = npartx * nparty * npartz

  case default

     call die("Unknown particle distribution '" // trim(particleDistribution) // "'")

  end select

  ! Distribute particles to processors
  call pigeon_hole(nParticlesGlobal, nProcs, iRank, particleOffset, nParticles)

  ! Allocate the particle vector
  if (allocated(particles)) deallocate(particles)
  allocate(particles(nParticles))

  ! Initialize the particle parameters
  sumParticleVolume = 0.0_WP
  do i = 1, nParticles

     ! Unique id
     particles(i)%id = int(particleOffset + i, kind = 8)

     ! Set particle size distribution (compact support lognormal)
     particles(i)%diameter = random_lognormal(m=dMean,sd=dStd) + dShift
     if (dStd.gt.0.0_WP) then
        do while (particles(i)%diameter.gt.dMax+epsilon(1.0_WP) .or.                         &
             particles(i)%diameter.lt.dMin-epsilon(1.0_WP))
           particles(i)%diameter=random_lognormal(m=dMean,sd=dStd) + dShift
        end do
     else
        particles(i)%diameter = dMean
     end if

     ! Position (this will be reset later)
     particles(i)%position = 0.0_WP

     ! Particle velocity
     particles(i)%velocity = 0.0_WP

     ! Particle temperature
     particles(i)%temperature = Tp

     ! Other parameters
     particles(i)%angularVelocity = 0.0_WP
     particles(i)%collision = 0.0_WP
     particles(i)%torque = 0.0_WP

     ! The particle exists
     particles(i)%stop = 0

     ! Sum particle volume
     sumParticleVolume = sumParticleVolume + pi * volumeFactor *                             &
          particles(i)%diameter ** nDimensions
  end do

  ! Compute effective volume fraction
  call parallel_sum(sumParticleVolume)
  volumeFraction = sumParticleVolume / volume

  ! Distribute the particles
  select case (trim(particleDistribution))
  case ('random')
     i=1
     do while (i.le.nParticles)

        ! Give the particle a random position
        if (nx .gt. 1) then
           if (std .gt. 0.0_WP) then
              particles(i)%position(1) = random_normal(x0 + 0.5_WP * particleThickness, std)
           else
              call random_number(rand)
              particles(i)%position(1) = 0.5_WP * particles(i)%diameter + x0 +               &
                   (particleThickness - 0.5_WP * particles(i)%diameter) * rand
           end if
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = Ly * (1.0_WP * rand - 0.5_WP) 
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = Lz * (1.0_WP * rand - 0.5_WP) 
        end if

        ! Prevent particle overlap
        success = .true.
        if (preventOverlap) then
           part2: do j = 1, i-1
              ! Compute separation distance (account for periodicity)
              distance = 0.0_WP
              do k = 1, nDimensions
                 distance(k) = abs(particles(i)%position(k) - particles(j)%position(k))
                 if (isPeriodic(k)) then
                    distance(k) = min(distance(k), periodicLength(k) -                       &
                         abs(particles(i)%position(k) - particles(j)%position(k)))
                 end if
              end do
              r = sqrt(sum(distance**2))
              if (r.le.0.55_WP*(particles(i)%diameter+particles(j)%diameter)) then
                 i=i-1
                 success = .false.
                 exit part2
              end if
           end do part2
        end if
        i=i+1

        if (success .and. modulo(i,1000).eq.0)                                               &
             print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'%'

     end do

  case ('uniform')
     do i = 1, nParticles
        ix = (particleOffset + i - 1) / (nparty * npartz)
        iy = (particleOffset + i - 1 - nparty * npartz * ix) / npartz
        iz = particleOffset + i - 1 - nparty * npartz * ix - npartz * iy
        particles(i)%position(1) = (ix + 0.5_WP) * Lp + x0
        particles(i)%position(2) = (iy + 0.5_WP) * Lpy - 0.5_WP * Ly
        particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz
     end do

  end select

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  return
end subroutine shock_tube_particles
