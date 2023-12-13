module particle_box

  ! External modules
  use precision
  use string
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Global variables
  integer :: nx, ny, nz
  real(WP) :: Lx, Ly, Lz, dx, dy, dz

end module particle_box

subroutine particle_box_grid

  ! Internal modules
  use particle_box

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

           ! Create X
           if (nx .gt. 1) then
              if (periodicityType(1) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /          &
                      real(nx - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /                  &
                      real(nx - 1, WP)
              end if
           end if

           ! Create Y
           if (ny .gt. 1) then
              if (periodicityType(2) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) /          &
                      real(ny - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                  &
                      real(ny - 1, WP)
              end if
           end if

           ! Create Z
           if (nz .gt. 1) then
              if (periodicityType(3) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /          &
                      real(nz - 1, WP)
              else
                 coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /                  &
                      real(nz - 1, WP)
              end if
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

  return
end subroutine particle_box_grid

subroutine particle_box_data

  ! Internal modules
  use particle_box

  ! External modules
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: gamma, P0, rho0, meanVelocity(3)

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  if (usePTKE) nUnknowns = nUnknowns + 1
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Initialize fluid velocity
  call parser_read('Mean U Velocity', meanVelocity(1), 0.0_WP)
  call parser_read('Mean V Velocity', meanVelocity(2), 0.0_WP)
  call parser_read('Mean W Velocity', meanVelocity(3), 0.0_WP)

  ! Ambient pressure
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  P0 = 1.0_WP / gamma

  ! Constant density
  rho0 = 1.0_WP

  ! Set the density
  conservedVariables(:, 1) = rho0

  ! Set the momentum
  do i = 1, nDimensions
     conservedVariables(:, i+1) = rho0 * meanVelocity(i)
  end do

  ! Set the energy
  conservedVariables(:, nDimensions+2) = P0 / (gamma - 1.0_WP) + 0.5_WP * rho0 *             &
       sum(meanVelocity(1:nDimensions)**2)

  ! Set the pseudo-TKE
  if (usePTKE) conservedVariables(:, nDimensions+3) = 0.0_WP

  return
end subroutine particle_box_data

subroutine particle_box_particles

  ! Internal modules
  use particle_box

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, j, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(WP) :: domainVolume, volumeFraction, particleVolume, volumeFactor,          &
       sumParticleVolume, velocityFluctuation, Lp, Lpy, Lpz, rand, dMean, dStd,    &
       dMin, dMax, dShift, Tp,  meanVelocity(3)
  character(len = str_medium) :: filename, particleDistribution

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return 

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('particle velocity fluctuations', velocityFluctuation, 0.0_WP)

  ! Read particle diameter
  call parser_read('particle mean diameter', dMean)
  call parser_read('particle std diameter', dStd, 0.0_WP)
  call parser_read('particle min diameter', dMin, 0.0_WP)
  call parser_read('particle max diameter', dMax, 0.0_WP)
  call parser_read('particle diameter shift', dShift, 0.0_WP)
  call parser_read('particle temperature', Tp, 2.5_WP)

  ! Compute domain volume
  domainVolume = sum(gridNorm(:,1))
  call parallel_sum(domainVolume)

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * dMean ** nDimensions

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')

  select case (trim(particleDistribution))

  case ('random')

     ! Get the number of particles based on the volume fraction and distribute to processors
     nParticlesGlobal = int(volumeFraction * domainVolume / particleVolume)

  case ('uniform')

     ! Mean interparticle distance
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     npartx = int(Lx / Lp)
     Lp = Lx / real(npartx, WP)
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

  ! Initialize the particles
  sumParticleVolume = 0.0_WP
  meanVelocity = 0.0_WP
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
           particles(i)%position(1) = Lx * rand
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = Ly * rand
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = Lz * rand
        end if

     case ('uniform')
        ix = (particleOffset + i - 1) / (nparty * npartz)
        iy = (particleOffset + i - 1 - nparty * npartz * ix) / npartz
        iz = particleOffset + i - 1 - nparty * npartz * ix - npartz * iy
        particles(i)%position(1) = (ix + 0.5_WP) * Lp
        particles(i)%position(2) = (iy + 0.5_WP) * Lpy
        particles(i)%position(3) = (iz + 0.5_WP) * Lpz

     end select

     ! Particle velocity
     particles(i)%velocity = 0.0_WP
     do j = 1, nDimensions
        if (globalGridSize(j) .gt. 1 .and. velocityFluctuation .gt. 0.0_WP) then
           call random_number(rand)
           rand = 2.0_WP * rand - 1.0_WP
           particles(i)%velocity(j) = particles(i)%velocity(j) + velocityFluctuation * rand
           meanVelocity(j) = meanVelocity(j) + velocityFluctuation * rand
        end if
     end do

     ! The particle exists
     particles(i)%stop = 0

     ! Sum particle volume
     sumParticleVolume = sumParticleVolume + pi * volumeFactor *                             &
          particles(i)%diameter ** nDimensions
  end do

  ! Remove mean fluctuations
  do j=1, nDimensions
     call parallel_sum(meanVelocity(j))
     meanVelocity(j) = meanVelocity(j) / real(nParticlesGlobal, WP)
  end do
  do i = 1, nParticles
     particles(i)%velocity(1:nDimensions) = particles(i)%velocity(1:nDimensions) -           &
          meanVelocity(1:nDimensions)
  end do

  ! Compute effective volume fraction
  call parallel_sum(sumParticleVolume)
  volumeFraction = sumParticleVolume / domainVolume

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  return
end subroutine particle_box_particles
