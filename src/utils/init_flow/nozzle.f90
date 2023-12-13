module nozzle

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
  real(WP) :: L, ID1, ID2,OD1, OD2, steepness, L0, wallThickness

  ! Nozzle types
  integer :: nozzleType
  integer, parameter ::                                                                      &
       LINEAR_NOZZLE = 1,                                                                    &
       TANH_NOZZLE   = 2

end module nozzle

subroutine nozzle_grid

  ! Internal modules
  use nozzle

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: sigma, b, c
  real(WP), allocatable, dimension(:) :: s, g
  logical :: stretchGrid

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx)
  call parser_read('Ly', Ly)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Compute the grid spacing
  dx = Lx / real(nx-1, WP)
  dy = Ly / real(ny-1, WP)
  if (nz.gt.1) dz = Lz / real(nz-1, WP)

  ! Should we stretch the mesh?
  call parser_read('grid stretching', stretchGrid, .false.)

  ! Create the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) / real(nx - 1, WP)
           coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP)                          &
                / real(ny - 1, WP) -  0.5_WP * Ly
           if (nz.gt.1) coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP)             &
                / real(nz - 1, WP) -  0.5_WP * Lz
        end do
     end do
  end do

  ! Grid stretching
  if (stretchGrid) then
     ! Grid stretching parameters
     sigma = 0.21_WP
     b = 12.0_WP
     c = 0.6_WP

     ! Stretch in y
     !--------------------------------------------------
     if (ny.gt.1) then
        ! Create uniform spacing
        allocate(s(ny))
        do j = 1, ny
           s(j) = real(j - 1, WP) / real(ny - 1, WP)
        end do

        ! Compute mapping g(s)
        allocate(g(ny))
        call mapping_function(s, b, c, sigma, g)

        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              do i = iStart(1), iEnd(1)
                 ! Create y
                 coordinates(grid_index(i,j,k), 2) = 0.5_WP * Ly * (1.0_WP + g(j)) - 0.5_WP * Ly
              end do
           end do
        end do
        deallocate(s)
        deallocate(g)
     end if

     ! Stretch in z
     !--------------------------------------------------
     if (nz.gt.1) then
        ! Create uniform spacing
        allocate(s(nz))
        do k = 1, nz
           s(k) = real(k - 1, WP) / real(nz - 1, WP)
        end do

        ! Compute mapping g(s)
        allocate(g(nz))
        call mapping_function(s, b, c, sigma, g)

        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              do i = iStart(1), iEnd(1)
                 ! Create z
                 coordinates(grid_index(i,j,k), 3) = 0.5_WP * Lz * (1.0_WP + g(k)) - 0.5_WP * Lz
              end do
           end do
        end do
        deallocate(s)
        deallocate(g)
     end if
  end if

  ! Setup the grid metrics
  call grid_metrics_setup

contains

  ! Mapping functional for grid stretching
  ! --------------------------------------
  subroutine mapping_function(s, b, c, sigma, g)

    ! External modules
    use math, only : pi

    implicit none

    real(WP), intent(in) :: s(:), b, c, sigma
    real(WP), intent(out) :: g(size(s))

    g = ((s - 0.5_WP) * (1.0_WP + 2.0_WP * b) - b * sigma *                                  &
         (exp(- ((s - 0.5_WP + c) / sigma) ** 2) / sqrt(pi) +                                &
         ((s - 0.5_WP + c) / sigma) * erf((s - 0.5_WP + c) / sigma) -                        &
         exp(- ((s - 0.5_WP - c) / sigma) ** 2) / sqrt(pi) -                                 &
         ((s - 0.5_WP - c) / sigma) * erf((s - 0.5_WP - c) / sigma))) /                      &
         (0.5_WP + b - b * sigma * (exp(- ((0.5_WP + c) / sigma) ** 2) /                     &
         sqrt(pi) + ((0.5_WP + c) / sigma) * erf((0.5_WP + c) / sigma) -                     &
         exp(- ((0.5_WP - c) / sigma) ** 2) / sqrt(pi) - ((0.5_WP - c) / sigma) *            &
         erf((0.5_WP - c) / sigma)))

    return
  end subroutine mapping_function

end subroutine nozzle_grid


subroutine nozzle_data

  ! Internal modules
  use nozzle

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, density, pressure, velocity, x, r, x0, U0, p0, T0, rho0, D, buf
  character(len = str_medium) :: inputString

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Select the nozzle type and read in geometric parameters
  call parser_read('nozzle type', inputString, 'linear')
  select case (trim(inputString))
  case ('linear', 'Linear', 'LINEAR')
     nozzleType = LINEAR_NOZZLE

     call parser_read('inner diameter 1', ID1)
     call parser_read('inner diameter 2', ID2)
     call parser_read('length', L)
     call parser_read('wall thickness', wallThickness)
     call parser_read('straight length', L0)

  case ('tanh', 'TANH')
     nozzleType = TANH_NOZZLE

     call parser_read('inner diameter 1', ID1)
     call parser_read('inner diameter 2', ID2)
     call parser_read('outer diameter 1', OD1)
     call parser_read('outer diameter 2', OD2)
     call parser_read('length', L)
     call parser_read('steepness', steepness, 10.0_WP)
     call parser_read('convergence location', L0, 0.5_WP * L)

  case default
     call die("nozzle_data: invalid nozzle type '" // trim(inputString) // "'!")

  end select

  ! Read in fluid parameters
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('nozzle temperature', T0, 1.0_WP / (gamma-1.0_WP))
  call parser_read('nozzle pressure', p0)
  call parser_read('nozzle velocity', U0)

  ! Apply ideal gas law
  rho0 = gamma / (gamma - 1.0_WP) * p0 / T0

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Get position
           x = coordinates(grid_index(i,j,k), 1)
           r = sqrt(sum(coordinates(grid_index(i,j,k), 2:nDimensions)**2))

           ! Set state variables based on nozzle geometry
           select case (nozzleType)
           case (LINEAR_NOZZLE)
              x0 = L0
              if ((r.le.0.5_WP*ID1 .and. x.le.x0)) then
                 buf = 1.0_WP
              else
                 buf = 0.0_WP
              end if

           case (TANH_NOZZLE)
              x0 = 1.0_WP*L
              D = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (x - L0)))
              if (r.le.0.5_WP * D .and. x.le.L) then
                 buf = 1.0_WP
              else
                 buf = 0.0_WP
              end if
           end select
           
           velocity = (U0 - 0.0_WP) * buf + 0.0_WP
           density = (rho0 - 1.0_WP) * buf + 1.0_WP
           pressure = (p0 - 1.0_WP / gamma) * buf + 1.0_WP/gamma

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
end subroutine nozzle_data


subroutine nozzle_levelset

  ! Internal modules
  use nozzle
  
  ! External modules
  use math
  use grid_levelset
  use parallel

  implicit none

  ! General
  integer :: gridIndex, i, j, k, ip
  real(WP), parameter :: eps = 1.0e-9_WP
  real(WP) :: r, x, y, z, rp, xp, rin, rout, ID1e, ID2e, xComp, rComp, straightHeightLevel
  integer, dimension(nGridPoints) :: normIndex
  real(WP), dimension(nGridPoints) :: thetaGrid, normIndicator
  real(WP), dimension(2) :: tempVec

  ! Levelset
  integer :: count
  real(WP) :: mydist

  ! Progress monitoring
  integer :: togo, prog
  integer :: iratio_new, iratio_old

  ! Levelset file
  real(WP) :: tmp

  ! Allocate distance array
  allocate(levelset(nGridPoints,1)); levelset=huge(1.0_WP)

  select case (nozzleType)

  case (LINEAR_NOZZLE)

     ! Define the nozzle levelset using straight lines
     ! -----------------------------------------------

     select case (nDimensions)   

     case(2)
        xComp = sin(atan((ID1*0.5_WP - ID2*0.5_WP)/L))
        rComp = sqrt(1.0_WP**2 - xComp**2)

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = 0.0_WP
           if (y .le. -1e-15_WP) thetaGrid(i) = 1.5_WP*pi
           if (y .ge. 1e-15) thetaGrid(i) = 0.5_WP*pi

           if (abs(y) .lt. 1e-15_WP) then
              thetaGrid(i) = 0.0_WP
           end if

           levelset(i,1) = y**2 + z**2
        end do

        levelset(:,1) = sqrt(levelset(:,1))

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)

           if (x .le. L0) then
              rin = 0.5_WP * ID1
           else
              rin = (ID2*0.5_WP - ID1*0.5_WP)/ L * x + ID1 * 0.5_WP                          &
                   - (ID2*0.5_WP - ID1*0.5_WP) / L * L0
           end if
           rout = rin + wallThickness

           normIndicator(i) = sign(1.0_WP,levelset(i,1) - (rin + rout)*0.5_WP)

           if(abs(levelset(i,1) - (rin + rout)*0.5_WP) .lt. 1e-15_WP) then 
              normIndicator(i) = 1.0_WP
           end if

           tempVec(1) = abs(levelset(i,1) - (rin + rout)*0.5_WP) - (rout - rin)*0.5_WP
           tempVec(2) = x - (L + L0) 
           straightHeightLevel = x - L0

           normIndex(i) = maxloc(tempVec(:), DIM = 1)
           levelset(i,1) = tempVec(normIndex(i))
           if(abs(tempVec(1) - tempVec(2)) .lt. 1e-15_WP) then
              levelset(i,1) = tempVec(2)
              normIndex(i) = 2    
           end if

        end do

     case(3)
        xComp = sin(atan((ID1*0.5_WP - ID2*0.5_WP)/L))
        rComp = sqrt(1.0_WP**2 - xComp**2)   

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = coordinates(i,3)

           thetaGrid(i) = pi - 0.5_WP*pi*(1.0_WP + sign(1.0_WP,z))                           &
                * (1.0_WP - sign(1.0_WP,y**2)) - 0.25_WP * pi * (2.0_WP + sign(1.0_WP, z))   &
                * (sign(1.0_WP,y)) - sign(1.0_WP, y*z)                                       &
                * atan((abs(z) - abs(y))/(abs(z) + abs(y)))
           if (abs(y) .lt. 1e-14_WP .and. abs(z) .lt. 1e-14_WP) then
              thetaGrid(i) = 0.0_WP
           end if

           levelset(i,1) = y**2 + z**2
        end do

        levelset(:,1) = sqrt(levelset(:,1))

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = coordinates(i,3)

           if (x .le. L0) then
              rin = ID1*0.5_WP
              rout = rin + wallThickness
           else
              rin = (ID2*0.5_WP - ID1*0.5_WP) / L * x + ID1*0.5_WP                           &
                   - (ID2*0.5_WP - ID1*0.5_WP)/L * L0
              rout = rin + wallThickness
           end if

           normIndicator(i) = sign(1.0_WP,levelset(i,1) - (rin + rout)*0.5_WP)

           if(abs(levelset(i,1) - (rin + rout)*0.5_WP) .lt. 4e-14_WP) then 
              normIndicator(i) = 1.0_WP
           end if

           tempVec(1) = abs(levelset(i,1) - (rin + rout)*0.5_WP) - (rout - rin)*0.5_WP
           tempVec(2) = x - (L + L0)
           straightHeightLevel = x - L0

           normIndex(i) = maxloc(tempVec(:), DIM = 1)
           levelset(i,1) = tempVec(normIndex(i))
           if(abs(tempVec(1) - tempVec(2)) .lt. 4e-14_WP) then
              levelset(i,1) = tempVec(2)
              normIndex(i) = 2    
           end if

        end do
     end select

  case (TANH_NOZZLE)

     ! Define the nozzle levelset as a tanh function
     ! ---------------------------------------------

     ! Prepare counter
     prog = 0
     togo = nGridPoints
     iratio_old = 0

     ! Store effective inner diameters
     ID1e = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (0.0_WP - L0)))
     ID2e = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (L - L0)))

     ! Compute projection
     if (iRank .eq. iRoot) print *,'Computing levelset analytically...'
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              ! Prepare projections
              count = 0
              mydist = huge(1.0_WP)

              ! Get position
              x = coordinates(gridIndex, 1)
              r = sqrt(sum(coordinates(gridIndex, 2:nDimensions)**2))

              ! Determine radius to inner and outer nozzle geometry
              rin = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                                   &
                   (1.0_WP - tanh(steepness / L * (x - L0))))
              rout = 0.5_WP * OD1 + 0.5_WP * (OD2 - OD1) * x / L

              ! Calculate distance
              if (r.lt.rin .and. x.le.L) then
                 ! Inner nozzle
                 do ip = 1,1000
                    xp = L * real(ip - 1,WP) / real(1000 - 1, WP)
                    rp = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                              &
                         (1.0_WP - tanh(steepness / L * (xp - L0))))
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp .lt. mydist) mydist = tmp
                 end do
              elseif (r .gt. rout .and. x.le.L) then
                 ! Cone
                 mydist = (r-rout) * cos(atan(0.5_WP*(OD2 - ID2e) / L))
              elseif (r.ge.rin .and. r.le.rout .and. x.le.L) then
                 ! In the solid
                 do ip = 1,1000
                    xp = L * real(ip - 1,WP) / real(1000 - 1, WP)
                    rp = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                              &
                         (1.0_WP - tanh(steepness / L * (xp - L0))))
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp.lt.mydist) mydist = tmp
                 end do
                 mydist = min(mydist, abs((rout-r) * cos(atan(0.5_WP * (OD2 - ID2e) / L))))
                 mydist = min(mydist, abs(x-L))
                 mydist = -mydist
              elseif (r.gt.0.5_WP * OD2 .and. x.gt.L) then
                 mydist = sqrt((x - L)**2 + (r - 0.5_WP * OD2)**2)
              elseif (r.le.0.5_WP * OD2 .and. r.ge.0.5_WP * ID2e .and. x.gt.L) then
                 mydist = x - L
              else
                 mydist = sqrt((x - L)**2 + (r - 0.5_WP * ID2e)**2)
              end if

              ! Store the minimum distance
              levelset(gridIndex,1) = mydist

              ! Add point to counter
              if (irank.eq.iroot) then
                 prog=prog+1
                 iratio_new=int(real(prog,WP)/real(togo,WP)*100.0_WP)
                 if (iratio_new.ge.iratio_old+20) then
                    iratio_old=iratio_new
                    write(*,'(i3,x,a1)') iratio_new,'%'
                 end if
              end if

           end do
        end do
     end do

  end select

  return
end subroutine nozzle_levelset


subroutine nozzle_particles

  ! Internal modules
  use nozzle

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(WP) :: diameter, domainVolume, volumeFraction, particleVolume, volumeFactor,          &
       sumParticleVolume, Lp, Lpx, Lpz, rand, hbed, Tp
  character(len = str_medium) :: filename, particleDistribution

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).eq.0) return 

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle diameter', diameter)
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('bed height', hbed)
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
     nparty = int(Ly / Lp)
     Lp = Ly / real(nparty, WP)
     npartx = int(hbed / Lp)
     Lpx = hbed / real(npartx, WP)
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
     particles(i)%velocity = 0.0_WP
     particles(i)%angularVelocity = 0.0_WP
     particles(i)%collision = 0.0_WP
     particles(i)%torque = 0.0_WP

     ! Distribute the particles
     particles(i)%position = 0.0_WP
     select case (trim(particleDistribution))
     case ('random')
        if (nx .gt. 1) then
           call random_number(rand)
           particles(i)%position(1) = Lx - hbed * rand
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = 0.5_WP * Ly * (2.0_WP * rand - 1.0_WP)
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = 0.5_WP * Lz * (2.0_WP * rand - 1.0_WP)
        end if

     case ('uniform')
        iy = (particleOffset + i - 1) / (npartx * npartz)
        ix = (particleOffset + i - 1 - npartx * npartz * iy) / npartz
        iz = particleOffset + i - 1 - npartx * npartz * iy - npartz * ix
        particles(i)%position(2) = (iy + 0.5_WP) * Lp - 0.5_WP * Ly
        particles(i)%position(1) = Lx - hbed - 0.5_WP * diameter + (ix + 0.5_WP) * Lpx
        particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz
     end select

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
  end if

  return
end subroutine nozzle_particles


subroutine nozzle_ibm_particles

  ! Internal modules
  use nozzle

  ! External modules
  use random
  use math
  use parallel
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, npartx, nparty, npartz, ix, iy, iz, ierr
  real(WP) :: dMean, Up
  character(len = str_medium) :: filename, particleDistribution

  ! Particle type
  type :: t_Particle
     real(WP) :: diameter
     real(WP), dimension(3) :: position
     real(WP), dimension(3) :: velocity
  end type t_Particle
  type(t_Particle), dimension(:), allocatable :: particles
  integer :: nParticles, nParticles1, nParticles2

  ! Return if not writing an ibm file
  call parser_read('init ibm file', filename, '')
  if (len_trim(filename).eq.0) return

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')
  call parser_read('particle mean diameter', dMean)
  call parser_read('particle velocity', Up)

  select case (trim(particleDistribution))

  case ('single')
     nParticles = 1

  case ('double')
     nParticles = 2

  case default

     call die("Unknown particle distribution '" // trim(particleDistribution) // "'")

  end select

  ! Allocate the object vector
  allocate(particles(nParticles))
  
  ! Only root process generates IBM particle data
  if (iRank .eq. iRoot) then

     ! Initialize the particle parameters
     do i = 1, nParticles
        particles(i)%diameter = dMean
        particles(i)%velocity = 0.0_WP
        particles(i)%velocity(1) = Up
        particles(i)%position = 0.0_WP
     end do

     ! Distribute the particles
     select case (trim(particleDistribution))

     case ('single')
        call parser_read('particle center', particles(1)%position(1:nDimensions))

     case ('double')
        call parser_read('particle 1 center', particles(1)%position(1:nDimensions))
        call parser_read('particle 2 center', particles(2)%position(1:nDimensions))

     end select

     ! Output stuff to the screen
     print *
     print *, 'Number of particles: ', nParticles
  end if

  ! Set IBM object data for i/o
  nObjects = nParticles
  allocate(object(nObjects))
  if (iRank .eq. iRoot) then
     do i = 1, nObjects
        select case (nDimensions)
        case (2)
           object(i)%volume = 0.25_WP * pi * particles(i)%diameter**2
        case (3)
           object(i)%volume = pi * particles(i)%diameter**3 / 6.0_WP
        end select
        object(i)%position = particles(i)%position
        object(i)%velocity = particles(i)%velocity
        object(i)%angularVelocity = 0.0_WP
     end do
  end if
  call MPI_BCAST(object, nObjects, MPI_OBJECT, iRoot, MPI_COMM_WORLD, ierr)

  return
end subroutine nozzle_ibm_particles
