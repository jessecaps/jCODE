module impulsive_plate

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

end module impulsive_plate


subroutine impulsive_plate_grid

  ! Internal modules
  use impulsive_plate

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP):: ytilde, r
  logical :: stretchGrid

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

  ! Should we stretch the mesh?
  call parser_read('grid stretching', stretchGrid, .false.)
  if (stretchGrid) r = 2.0_WP

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Uniform and periodic in x
           coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /             &
                   real(nx - 1, WP) - 0.5_WP * Lx

           ! Create Y
           if (stretchGrid) then
              ytilde = real(ny-j, WP) / real(ny-1, WP)
              coordinates(grid_index(i,j,k), 2) = Ly * (1.0_WP - tanh(r * ytilde) / tanh(r))
           else
              coordinates(grid_index(i,j,k), 2) = Ly * real(j-1, WP) / real(ny-1, WP)
           end if

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

end subroutine impulsive_plate_grid


subroutine impulsive_plate_data

  ! Internal modules
  use impulsive_plate

  ! External modules
  use parallel
  use parser
  use random
  use state, only : conservedVariables, targetState

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, U0, p0, T0, rho0, velocity(3)

  ! Set the number of conserved variables and allocate the arrays
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))
  allocate(targetState(nGridPoints, nUnknowns))

  ! Set flow quantities
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('free stream velocity', U0)
  p0 = 1.0_WP / gamma
  T0 = 1.0_WP / (gamma - 1.0_WP)
  rho0 = 1.0_WP
  velocity = 0.0_WP

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Set the conserved variables
           velocity(1) = U0
           conservedVariables(grid_index(i,j,k), 1) = rho0
           conservedVariables(grid_index(i,j,k), 2:nDimensions+1) =                          &
                rho0 * velocity(1:nDimensions)
           conservedVariables(grid_index(i,j,k), nDimensions+2) = p0 /                       &
                (gamma - 1.0_WP) + 0.5_WP * rho0 * sum(velocity ** 2)

           ! Set the target state
           if (j .eq. 1) then
              velocity(1) = 0.0_WP
           else
              velocity(1) = U0
           end if
           targetState(grid_index(i,j,k), 1) = rho0
           targetState(grid_index(i,j,k), 2:nDimensions+1) =                                 &
                rho0 * velocity(1:nDimensions)
           targetState(grid_index(i,j,k), nDimensions+2) = p0 /                              &
                (gamma - 1.0_WP) + 0.5_WP * rho0 * sum(velocity ** 2)

        end do
     end do
  end do

  return
end subroutine impulsive_plate_data

subroutine impulsive_plate_particles

  ! Internal modules
  use impulsive_plate

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, j, k, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(wp) :: diameter, volume, volumeFraction, particleVolume, particleThickness,           &
       particleVelocity, volumeFactor, sumParticleVolume, velocityFluctuation,               &
       Lp, Lpy, Lpz, y0, rand, Tp
  character(len = str_medium) :: filename, particleDistribution

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle diameter', diameter)
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('particle velocity', particleVelocity)
  call parser_read('particle velocity fluctuations', velocityFluctuation, 0.0_WP)
  call parser_read('particle layer thickness', particleThickness, Ly) 
  call parser_read('particle layer position', y0, 0.0_WP)
  call parser_read('particle temperature', Tp, 2.5_WP)

  ! Initial Volume that particles occupy
  volume = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           if (coordinates(grid_index(i,j,k), 2).ge.y0 .and.                                 &
                coordinates(grid_index(i,j,k), 2).le.y0+particleThickness)                   &
                volume = volume + gridNorm(grid_index(i,j,k), 1)
        end do
     end do
  end do
  call parallel_sum(volume)

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * diameter ** nDimensions

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')

  select case (trim(particleDistribution))

  case ('random')

     ! Get the number of particles based on the volume fraction and distribute to processors
     nParticlesGlobal = int(volumeFraction * volume / particleVolume)

  case ('uniform')

     ! Mean interparticle distance
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     npartx = int(Lx / Lp)
     Lp = Lx / real(npartx, WP)
     if (ny .gt. 1) nparty = int(particleThickness / Lp)
     Lpy = particleThickness / real(nparty, WP)
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

  ! Initialize the particles
  sumParticleVolume = 0.0_WP
  do i = 1, nParticles

     ! Unique id
     particles(i)%id = int(particleOffset + i, kind = 8)

     ! Particle diameter
     particles(i)%diameter = diameter

     ! Particle temperature
     particles(i)%temperature = Tp

     ! Other parameters
     particles(i)%angularVelocity = 0.0_WP
     particles(i)%collision = 0.0_WP
     particles(i)%torque = 0.0_WP

     ! Distribute the particles
     particles(i)%position = 0.0_WP
     select case (trim(particleDistribution))
     case ('random')
        if (nx .gt. 1) then
           call random_number(rand)
           particles(i)%position(1) = Lx * (1.0_WP * rand - 0.5_WP)  
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = 0.5_WP * particles(i)%diameter + y0 +                  &
                (particleThickness - 0.5_WP * particles(i)%diameter) * rand
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = Lz * (1.0_WP * rand - 0.5_WP) 
        end if

     case ('uniform')
        ix = (particleOffset + i - 1) / (nparty * npartz)
        iy = (particleOffset + i - 1 - nparty * npartz * ix) / npartz
        iz = particleOffset + i - 1 - nparty * npartz * ix - npartz * iy
        particles(i)%position(1) = (ix + 0.5_WP) * Lp - 0.5_WP * Lx
        particles(i)%position(2) = (iy + 0.5_WP) * Lpy + y0
        particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz

     end select

     ! Particle velocity
     particles(i)%velocity = 0.0_WP
     particles(i)%velocity(1) = particleVelocity
     do j = 1, nDimensions
        if (globalGridSize(j) .gt. 1 .and. velocityFluctuation .gt. 0.0_WP) then
           call random_number(rand)
           rand = 2.0_WP * rand - 1.0_WP
           particles(i)%velocity(j) = particles(i)%velocity(j) + velocityFluctuation * rand
        end if
     end do

     ! The particle exists
     particles(i)%stop = 0

     ! Sum particle volume
     sumParticleVolume = sumParticleVolume + pi * volumeFactor *                             &
          particles(i)%diameter ** nDimensions
  end do

  ! Compute effective volume fraction
  call parallel_sum(sumParticleVolume)
  volumeFraction = sumParticleVolume / volume

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  return
end subroutine impulsive_plate_particles
