module el_ibm

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

  ! Global variables (grid)
  integer :: nx, ny, nz
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, maxDxDyDz

  ! Global variables (particles)
  real(WP) :: curtainThickness, curtainPosition

end module el_ibm


subroutine el_ibm_grid

  ! Internal modules
  use el_ibm

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

  maxDxDyDz = dx
  if (nDimensions .gt. 1) maxDxDyDz = max(maxDxDyDz, dy)
  if (nDimensions .gt. 2) maxDxDyDz = max(maxDxDyDz, dz)

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
                      real(ny - 1, WP) - 0.5_WP * Ly
              else
                 coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                  &
                      real(ny - 1, WP) - 0.5_WP * Ly
              end if
           end if

           ! Create Z
           if (nz .gt. 1) then
              if (periodicityType(3) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /          &
                      real(nz - 1, WP) - 0.5_WP * Lz
              else
                 coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /                  &
                      real(nz - 1, WP) - 0.5_WP * Lz
              end if
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine el_ibm_grid


subroutine el_ibm_data

  ! Internal modules
  use el_ibm

  ! External modules
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, x, x0, p1, p2, rho1, rho2, u1, u2, density, pressure, velocity

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Specific heat ratio
  call parser_read('ratio of specific heats', gamma, 1.4_WP)

  ! Read in shock parameters
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('shock position', x0, 0.0_WP)
  call parser_read('post-shock pressure', p1, 1.0_WP / gamma)
  call parser_read('pre-shock pressure', p2, 1.0_WP / gamma)
  call parser_read('post-shock density', rho1, 1.0_WP)
  call parser_read('pre-shock density', rho2, 1.0_WP)
  call parser_read('post-shock velocity', u1, 0.0_WP)
  call parser_read('pre-shock velocity', u2, 0.0_WP)

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Get local coordinate
           x = coordinates(grid_index(i,j,k), 1)

           ! Set the density & pressure
           if (x.le.x0) then
              density = rho1
              pressure = p1
              velocity = u1
           else
              density = rho2
              pressure = p2
              velocity = u2
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

  return
end subroutine el_ibm_data


subroutine el_ibm_objects

  ! Internal modules
  use el_ibm

  ! External modules
  use random
  use math
  use parallel
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, npartx, nparty, npartz, ix, iy, iz, ierr
  integer :: nParticles
  real(WP) :: volume, volumeFraction, particleVolume, volumeFactor, dp,                      &
       sumParticleVolume, rand, r, distance(3), particleSpacing,                             &
       coordCenter(nDimensions), inBound(1:nDimensions),                                     &
       maxDistCenter(1:nDimensions), lShift, minDistances
  character(len = str_medium) :: filename
  logical :: success

  ! Particle type
  type :: t_Particle
     real(WP) :: diameter
     real(WP), dimension(3) :: position
     real(WP), dimension(3) :: velocity
  end type t_Particle
  type(t_Particle), dimension(:), allocatable :: particles_ibm

  ! Return if not writing an ibm file
  call parser_read('init ibm file', filename, '')

  if (len_trim(filename).eq.0) return

  ! Read in the particle parameters
  call parser_read('curtain thickness', curtainThickness, Lx) 
  call parser_read('curtain position', curtainPosition, 0.0_WP)

  ! Factor for computing volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if

  ! Read particle diameter
  call parser_read('big particle diameter', dp)
  particleVolume = pi * volumeFactor * dp ** nDimensions

  ! Initialize the random number generator
  call random_setup

  ! Initial Volume that particles occupy
  if (nz.gt.1) then
     volume = curtainThickness * Ly * Lz
  else
     volume = curtainThickness * Ly
  end if

  ! Get the number of particles based on the volume fraction and distribute to processors
  call parser_read('big particle volume fraction', volumeFraction)
  nParticles = int(volumeFraction * volume / particleVolume)

  ! Allocate the object vector
  allocate(particles_ibm(nParticles))

  ! Only root process generates IBM particle data
  if (iRank .eq. iRoot) then

     ! Initialize the particle parameters
     sumParticleVolume = 0.0_WP
     do i = 1, nParticles

        ! Store diameter
        particles_ibm(i)%diameter = dp

        ! Position (this will be reset later)
        particles_ibm(i)%position = 0.0_WP

        ! Particle velocity
        particles_ibm(i)%velocity = 0.0_WP

        ! Sum particle volume
        sumParticleVolume = sumParticleVolume + pi * volumeFactor *                          &
             particles_ibm(i)%diameter ** nDimensions
     end do

     ! Compute effective volume fraction
     volumeFraction = sumParticleVolume / volume

     ! Distribute the particles
     i=1
     do while (i.le.nParticles)

        ! Give the particle a random position
        if (nx .gt. 1) then
           call random_number(rand)
           particles_ibm(i)%position(1) = curtainThickness * rand + curtainPosition
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles_ibm(i)%position(2) = Ly * (1.0_WP * rand - 0.5_WP) 
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles_ibm(i)%position(3) = Lz * (1.0_WP * rand - 0.5_WP) 
        end if

        ! Prevent particle overlap
        success = .true.
        part2: do j = 1, i-1
           ! Compute separation distance (account for periodicity)
           distance = 0.0_WP
           do k = 1, nDimensions
              distance(k) = abs(particles_ibm(i)%position(k) - particles_ibm(j)%position(k))
              if (isPeriodic(k)) then
                 distance(k) = min(distance(k), periodicLength(k) -                          &
                      abs(particles_ibm(i)%position(k) - particles_ibm(j)%position(k)))
              end if
           end do
           r = sqrt(sum(distance**2))
           particleSpacing = 0.5_WP*(particles_ibm(i)%diameter+particles_ibm(j)%diameter)    &
                + 3.0_WP * maxDxDyDz
           if (r.le.particleSpacing) then
              i=i-1
              success = .false.
              exit part2
           end if
        end do part2
        i=i+1

        if (success .and. modulo(i,1000).eq.0)                                               &
             print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'% (big)'
     end do

     ! Output stuff to the screen
     print *
     print *, 'Number of big particles: ', nParticles
     print *, 'Big particle diameter: ', dp
     print *, 'Volume fraction (big): ', volumeFraction
     print *
  end if

  ! Set IBM object data for i/o
  nObjects = nParticles
  allocate(object(nObjects))
  if (iRank .eq. iRoot) then
     do i = 1, nObjects
        select case (nDimensions)
        case (2)
           object(i)%volume = 0.25_WP * pi * particles_ibm(i)%diameter**2
        case (3)
           object(i)%volume = pi * particles_ibm(i)%diameter**3 / 6.0_WP
        end select
        object(i)%position = particles_ibm(i)%position
        object(i)%velocity = particles_ibm(i)%velocity
        object(i)%angularVelocity = 0.0_WP
        object(i)%remove = .false.
     end do
  end if
  call MPI_BCAST(object, nObjects, MPI_OBJECT, iRoot, MPI_COMM_WORLD, ierr)

  return
end subroutine el_ibm_objects

subroutine el_ibm_particles

  ! Internal modules
  use el_ibm

  ! External modules
  use random
  use math
  use parallel
  use particle
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, particleOffset, ierror
  real(WP) :: volume, volumeFraction, particleVolume, volumeFactor, sumParticleVolume, rand, &
       r, distance(3), dp, dpBig, Tp, myLength, xStart, lengthProc(nprocs)
  character(len = str_medium) :: filename, particleDistribution
  logical :: success

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return 

  ! Read in the (small) particle parameters
  call parser_read('small particle volume fraction', volumeFraction)
  call parser_read('small particle diameter', dp)
  call parser_read('small particle temperature', Tp, 2.5_WP)

  ! Initial Volume that particles occupy
  if (nz.gt.1) then
     volume = curtainThickness * Ly * Lz
  else
     volume = curtainThickness * Ly
  end if

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * dp ** nDimensions

  ! Get the number of particles based on the volume fraction and distribute to processors
  nParticlesGlobal = int(volumeFraction * volume / particleVolume)

  ! Distribute particles to processors
  call pigeon_hole(nParticlesGlobal, nProcs, iRank, particleOffset, nParticles)

  ! Determine local extent to distribute particles
    myLength = curtainThickness * real(nParticles, WP) / real(nParticlesGlobal, WP)
  call MPI_Allgather(myLength, 1, MPI_REAL_WP, lengthProc, 1, MPI_REAL_WP, comm, ierror)
  xStart = curtainPosition
  do i = 1, iRank
     xStart = xStart + lengthProc(i)
  end do

  ! Allocate the particle vector
  if (allocated(particles)) deallocate(particles)
  allocate(particles(nParticles))

  ! Initialize the particle parameters
  sumParticleVolume = 0.0_WP
  do i = 1, nParticles

     ! Unique id
     particles(i)%id = int(particleOffset + i, kind = 8)

     ! Set particle size distribution (compact support lognormal)
     particles(i)%diameter = dp

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
     sumParticleVolume = sumParticleVolume + pi * volumeFactor * dp ** nDimensions
  end do

  ! Compute effective volume fraction
  call parallel_sum(sumParticleVolume)
  volumeFraction = sumParticleVolume / volume
  
  ! Distribute the particles
  i=1
  do while (i.le.nParticles)

     ! Give the particle a random position
     if (nx .gt. 1) then
        call random_number(rand)
        particles(i)%position(1) = xStart + 0.5_WP * dp + (myLength - dp) * rand
     end if
     if (ny .gt. 1) then
        call random_number(rand)
        particles(i)%position(2) = Ly * (1.0_WP * rand - 0.5_WP) 
     end if
     if (nz .gt. 1) then
        call random_number(rand)
        particles(i)%position(3) = Lz * (1.0_WP * rand - 0.5_WP) 
     end if

     ! Prevent overlap with big particles
     success = .true.
     part2: do j = 1, nObjects
        ! Compute separation distance (account for periodicity)
        dpBig = (2.0_WP * real(nDimensions, WP) * object(j)%volume / pi)                    &
             ** (1.0_WP / real(nDimensions, WP))
        distance = 0.0_WP
        do k = 1, nDimensions
           distance(k) = abs(particles(i)%position(k) - object(j)%position(k))
           if (isPeriodic(k)) then
              distance(k) = min(distance(k), periodicLength(k) -                             &
                   abs(particles(i)%position(k) - object(j)%position(k)))
           end if
        end do
        r = sqrt(sum(distance**2))
        if (r.le.0.55_WP*(dp+dpBig)) then
           i=i-1
           success = .false.
           exit part2
        end if
     end do part2

     ! Prevent inter-particle overlap
     if (success) then
        part3: do j=1, i-1
           distance=0.0_WP
           ! Compute separation distance 
           do k = 1, nDimensions
              distance(k) = abs(particles(i)%position(k) - particles(j)%position(k))
              if (isPeriodic(k)) then
                 distance(k) = min(distance(k), periodicLength(k) -                          &
                      abs(particles(i)%position(k) - particles(j)%position(k)))
              end if
           end do
           r = sqrt(sum(distance**2))
           if (r.le.0.55_WP*2.0_WP*dp) then
              i=i-1
              success = .false.
              exit part3
           end if
        end do part3
     end if

     ! Update counter
     i=i+1

     if (iRank .eq. iRoot .and. success .and. modulo(i,1000).eq.0)                            &
          print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'% (small)'

  end do

     ! Output stuff to the screen.
     if (iRank .eq. iRoot) then
     print *
     print *, 'Number of small particles: ', nParticlesGlobal
     print *, 'Small particle diameter: ', dp
     print *, 'Volume fraction (small): ', volumeFraction
     print *
  end if

  return
end subroutine el_ibm_particles
