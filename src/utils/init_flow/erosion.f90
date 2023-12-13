module erosion

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
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, Hb

end module erosion

subroutine erosion_grid

  ! Internal modules
  use erosion

  implicit none

  ! Local variables
  integer :: i, j, k, ny_bed, offset
  real(WP) :: dyMin, dyMax, r, y0
  logical :: stretch

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
  call parser_read('stretch grid above bed', stretch, .false.)

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
           if (ny.gt.1 .and. .not.stretch) then
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

  ! Apply stretching above from the bed
  ! Applies uniform grid spacing with dy=dyMax inside the bed and dy=dyMin above the bed
  if (stretch) then
     call parser_read('bed height', Hb)
     call parser_read('min dy', dyMin)
     call parser_read('max dy', dyMax)
     call parser_read('stretching coefficient',r, 2.0_WP)
     call parser_read('stretching offset',offset, 2)
     ny_bed = floor(Hb / dyMax)

     ! Grid spacing follows tanh profile
     j = 1
     y0 = r * log(cosh((real(ny_bed + offset - j, WP)) / r)) *                               &
          0.5_WP * real(dyMin - dyMax, WP) -                                                 &
          0.5_WP * real(j, WP) * real(dyMin - dyMax, WP) +                                   &
          real(dyMin, WP) * real(j, WP)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              coordinates(grid_index(i,j,k), 2) = r *                                        &
                   log(cosh((real(ny_bed + offset - j, WP)) / r)) *                          &
                   0.5_WP * real(dyMin - dyMax, WP) -                                        &
                   0.5_WP * real(j, WP) * real(dyMin - dyMax, WP) +                          &
                   real(dyMin, WP) * real(j, WP) - y0
           end do
        end do
     end do
  end if

  ! Setup the grid metrics
  call grid_metrics_setup

  return
end subroutine erosion_grid

subroutine erosion_data

  ! Internal modules
  use erosion

  ! External modules
  use math
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, U0, P0, rho0, x, y, z, amp, u, v, w, coeff, rnd, spongeThickness
  character(len = str_medium) :: val

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read from input
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('Mean U Velocity', U0, 0.0_WP)
  call parser_read('profile', val, 'none')
  call parser_read('fluctuation amplitude', amp, 0.0_WP)
  call parser_read('bed height', Hb)
  call parser_read('sponge thickness', SpongeThickness, 0.0_WP)

  ! Get ambient conditions
  P0 = 1.0_WP / gamma
  rho0 = 1.0_WP

  select case (trim(val))
  case ('none', 'None', 'NONE')

     do i = 1, nGridPoints
        conservedVariables(i, 1) = rho0
        conservedVariables(i, 2) = rho0 * U0
        conservedVariables(i, 3:nDimensions+1) = 0.0_WP
        conservedVariables(i, nDimensions+2) = P0 / (gamma - 1.0_WP) + 0.5_WP * rho0 * U0**2
     end do
     
  case ('couette', 'Couette', 'COUETTE')

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)

              ! Local coordinates
              x  = coordinates(grid_index(i,j,k), 1)
              y  = coordinates(grid_index(i,j,k), 2)
              if (nz.gt.1) z = coordinates(grid_index(i,j,k), 3)

              if (y.le.Hb) then
                 ! Zero out the velocity in the bed
                 u = 0.0_WP
                 v = 0.0_WP
                 w = 0.0_WP
              else if (y.gt.Hb .and. y.le.Ly-spongeThickness) then
                 ! Linear profile above the bed
                 u = u0 * (y - Hb) / (Ly - spongeThickness - Hb)

                 ! Add fluctuations (avoid walls)
                 if (j.gt.1 .and. j.lt.ny) then
                    
                    ! For faster transition
                    if (nz.gt.1) then
                       ! Fluctuations in x for w
                       w = w + amp * abs(u0) * cos(8.0_WP * twoPi * x / Lx)
                       ! Fluctuations in z for u
                       u = u + amp * abs(u0) * cos(8.0_WP * twoPi * z / Lz)
                    end if

                    ! Random values
                    call random_number(rnd)
                    u = u + (rnd - 0.5_WP) * amp * u
                    call random_number(rnd)
                    v = v + (rnd - 0.5_WP) * amp * v
                    call random_number(rnd)
                    w = w + (rnd - 0.5_WP) * amp * w
                 end if
              
              else
                 ! Constant velocity in the sponge
                 u = u0
              end if

              ! Set the density
              conservedVariables(grid_index(i,j,k), 1) = rho0

              ! Set the momentum
              conservedVariables(grid_index(i,j,k), 2) = rho0 * u
              conservedVariables(grid_index(i,j,k), 3) = rho0 * v
              if (nz.gt.1) conservedVariables(grid_index(i,j,k), 4) = rho0 * w

              ! Set the energy
              conservedVariables(grid_index(i,j,k), nDimensions+2) =                         &
                   p0 / (gamma - 1.0_wp) + 0.5_WP *                                          &
                   sum(conservedVariables(grid_index(i,j,k), 2:nDimensions+1)**2) /          &
                   conservedVariables(grid_index(i,j,k), 1)

           end do
        end do
     end do
     
  case ('poiseuille', 'Poiseuille', 'POISEUILLE')

     ! Set the conserved variables
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)

              ! Local coordinates
              x  = coordinates(grid_index(i,j,k), 1)
              y  = coordinates(grid_index(i,j,k), 2)
              if (nz.gt.1) z = coordinates(grid_index(i,j,k), 3)

              if (y.le.Hb) then
                 ! Zero out the velocity in the bed
                 u = 0.0_WP
                 v = 0.0_WP
                 w = 0.0_WP
              else
                 ! Laminar profile
                 coeff = 3.0_WP / 2.0_WP
                 u = u0 * coeff * (y - Hb) * (2.0_WP * (Ly - Hb) - (y - Hb)) / (Ly - Hb)**2
                 w = 0.0_WP
                 
                 ! For faster transition
                 if (nz.gt.1) then
                    ! Fluctuations in x for w
                    w = w + amp * abs(u0) * cos(8.0_WP * twoPi * x / Lx)
                    ! Fluctuations in z for u
                    u = u + amp * abs(u0) * cos(8.0_WP * twoPi * z / Lz)
                 end if
              
                 ! Random values
                 call random_number(rnd)
                 u = u + (rnd - 0.5_WP) * amp * u
                 call random_number(rnd)
                 v = v + (rnd - 0.5_WP) * amp * v
                 call random_number(rnd)
                 w = w + (rnd - 0.5_WP) * amp * w
              end if

              ! Set the density
              conservedVariables(grid_index(i,j,k), 1) = rho0

              ! Set the momentum
              conservedVariables(grid_index(i,j,k), 2) = rho0 * u
              conservedVariables(grid_index(i,j,k), 3) = rho0 * v
              if (nz.gt.1) conservedVariables(grid_index(i,j,k), 4) = rho0 * w

              ! Set the energy
              conservedVariables(grid_index(i,j,k), nDimensions+2) =                         &
                   p0 / (gamma - 1.0_wp) + 0.5_WP *                                          &
                   sum(conservedVariables(grid_index(i,j,k), 2:nDimensions+1)**2) /          &
                   conservedVariables(grid_index(i,j,k), 1)

           end do
        end do
     end do

  case default
     
     call die("Unknown profile '" // trim(val) // "'")
     
  end select
  
  return
end subroutine erosion_data


subroutine erosion_particles

  ! Internal modules
  use erosion

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, j, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(WP) :: diameter, domainVolume, volumeFraction, particleVolume, volumeFactor,          &
       sumParticleVolume, velocityFluctuation, Lp, Lpy, Lpz, rand, hbed, Tp
  character(len = str_medium) :: filename, particleDistribution

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return 

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle diameter', diameter)
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('particle velocity fluctuations', velocityFluctuation, 0.0_WP)
  call parser_read('initial bed height', hbed)
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
  particleVolume = pi * volumeFactor * diameter ** nDimensions

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
     if (ny .gt. 1) nparty = int(hbed / Lp)
     Lpy = hbed / real(nparty, WP)
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
           particles(i)%position(1) = Lx * rand
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = (hbed - 0.5_WP*diameter) * rand + 0.5_WP*diameter
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
        particles(i)%position(2) = (iy + 0.5_WP) * Lpy + 0.5_WP*diameter
        particles(i)%position(3) = (iz + 0.5_WP) * Lpz

     end select

     ! Particle velocity
     particles(i)%velocity = 0.0_WP
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
  volumeFraction = sumParticleVolume / domainVolume

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  return
end subroutine erosion_particles
