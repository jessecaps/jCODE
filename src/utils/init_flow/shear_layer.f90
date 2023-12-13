module shear_layer

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
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, initialThickness

end module shear_layer

subroutine shear_layer_grid

  ! Internal modules
  use shear_layer

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: sigma, b, c, minMeshsize, maxMeshsize, y1,  y2
  real(WP), allocatable, dimension(:) :: s, g
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
           if (ny .gt. 1 .and. .not. stretchGrid) then
              coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                     &
                   real(ny - 1, WP) - 0.5_WP * Ly
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

  ! Grid stretching
  if (stretchGrid) then
     ! Grid stretching parameters
     sigma = 0.21_WP
     b = 12.0_WP
     c = 0.6_WP

     ! Create uniform spacing
     allocate(s(ny))
     do j = 1, ny
        s(j) = real(j - 1, WP) / real(ny - 1, WP)
     end do

     ! Compute mapping g(s)
     allocate(g(ny))
     call mapping_function(s, b, c, sigma, g)

     ! Find min/max spacing
     minMeshsize =  huge(1.0_WP)
     maxMeshsize = -huge(1.0_WP)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Create y
              coordinates(grid_index(i,j,k), 2) = 0.5_WP * Ly * (1.0_WP + g(j)) - 0.5_WP * Ly

              ! Find min/max spacing
              if (j .gt. iStart(2) + 1) then
                 y1 = coordinates(grid_index(i,j  ,k), 2)
                 y2 = coordinates(grid_index(i,j-1,k), 2)
                 minMeshsize = min(minMeshsize, abs(y2 - y1))
                 maxMeshsize = max(maxMeshsize, abs(y2 - y1))
              end if
           end do
        end do
     end do
     call parallel_max(maxMeshsize)
     call parallel_min(minMeshsize)
     if (iRank .eq. iRoot) then
        print *
        print *, 'min/max y-spacing:', minMeshsize, maxMeshsize
        print *
     end if

     deallocate(s)
     deallocate(g)
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

end subroutine shear_layer_grid

subroutine shear_layer_data

  ! Internal modules
  use shear_layer

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, H2, O2, N2, VAP, AIR
  integer, parameter :: nkx = 10
  integer, parameter :: nkz = 5
  integer:: midplaneCount, zint, xint, numBlocks = 5
  real(WP) :: density, temperature, pressure, velocity, velocityDifference,                  &
       lowerTemperature, upperTemperature, T0, P0, x, y, z, Lpx, Lpz,                        &
       scalar, fuel, oxidizer, inert, mixtureFraction, Z0, YF0, Yo0, rand,                   &
       shear_layer_initial_rms_fraction, midplaneRMS, minDensity, maxDensity, Yv
  real(WP), dimension(:,:), allocatable :: gradA, phaseJitter, velocityFluctuations
  real(WP), dimension(:), allocatable :: A, Wi
  real(WP), dimension(nkx, -nkz:nkz) :: phases
  logical :: useSolenoidalPerturbations, useH2combustion, useVapor

  ! Get the number of species
  call parser_read('number of species', nSpecies, 0)

  call parser_is_defined('vapor mass fraction', useVapor)
  if (useVapor) nSpecies = nSpecies + 1

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2 + nSpecies

  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Setup solver_options to get species data
  call solver_options_setup

  ! Get species indices
  H2 = 0
  O2 = 0
  N2 = 0
  VAP= 0
  if (nSpecies .gt. 0) then
     if (allocated(speciesName)) then
        do k = 1, nSpecies + 1
           select case (trim(speciesName(k)))
           case ('H2', 'HYDROGEN')
              H2 = k
           case ('O2', 'OXYGEN')
              O2 = k
           case ('N2', 'NITROGEN')
              N2 = k
           case ('water', 'WATER', 'H2O')
              VAP = k
           case ('air', 'AIR')
              AIR = k
           case default
              call die("shear_layer_data: unknown species: '" //                             &
                   trim(speciesName(k)) // "!")

           end select
        end do
     else
        ! Hard-code H2-O2 mechanism
        H2 = 1
        O2 = 2
        N2 = nSpecies + 1
     end if
  end if

  ! Get molecular weights
  if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
     allocate(Wi(nSpecies+1))
     Wi = molecularWeightInverse
  end if

  ! Initialize random number generator
  call random_setup

  ! Zero-out the conserved variables
  conservedVariables = 0.0_WP

  ! Reference temperature and pressure
  T0 = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  P0 = 1.0_WP / ratioOfSpecificHeats

  ! Read in the shear layer parameters
  call parser_read('shear layer thickness', initialThickness)
  call parser_read('shear layer velocity difference', velocityDifference)
  call parser_read('shear layer lower temperature', lowerTemperature, T0)
  call parser_read('shear layer upper temperature', upperTemperature, T0)
  call parser_read('shear layer pressure', pressure, P0)

  ! Mixture properties
  useH2combustion = .false.
  if (H2.gt.0 .and. O2.gt.0) then
     useH2combustion = .true.
     call parser_read('initial mixture fraction', Z0, 1.0_WP)
     call parser_read('initial fuel mass fraction', Yf0)
     call parser_read('initial oxidizer mass fraction', Yo0)
  end if

  ! Set vapor mass fraction if used
  if (useVapor) call parser_read('vapor mass fraction', Yv)

  ! Mean profiles
  minDensity = huge(1.0_WP); maxDensity = -huge(1.0_WP)
  do i = 1, nGridPoints

     ! Vertical coordinate
     y = coordinates(i, 2)

     ! Velocity profile
     velocity = 0.5_WP * velocityDifference * tanh(0.5_WP * y / initialThickness)

     ! Temperature profile
     temperature = lowerTemperature + 0.5_WP * (upperTemperature - lowerTemperature) *       &
          (1.0_WP + tanh(0.5_WP * y / initialThickness))

     ! Hydrogen combustion
     if (useH2combustion) then
        ! Mixture fraction
        mixtureFraction = 0.5_WP * Z0 * ( 1.0_WP + tanh(-0.5_WP * y / initialThickness))

        ! Mass fractions
        fuel = YF0 * mixtureFraction
        oxidizer = YO0 * (1.0_WP - mixtureFraction)

        ! Correct species mass fractions
        inert = 1.0_WP - fuel - oxidizer
        if (inert .lt. 0.0_WP) oxidizer = oxidizer + inert
     end if

     ! Assign the density
     if (equationOfState .eq. IDEAL_GAS_MIXTURE) then
        ! Account for mixture of moleculure weights
        density = ratioOfSpecificHeats * pressure / (ratioOfSpecificHeats - 1.0_WP) /        &
             temperature
        if (useH2combustion) then
           density = density /                                                               &
                (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2))
        elseif (useVapor) then
           density = density / (Yv * (Wi(VAP) - Wi(AIR)) + Wi(AIR))
        end if
     else
        ! Crocco-Busemann (Sandham 1990 & Vaghefi 2014)
        density = 1.0_WP + 0.5_WP * (ratioOfSpecificHeats - 1.0_WP) *                        &
             (0.25_WP * velocityDifference**2 - velocity**2)
        density = 1.0_WP / density
     end if

     ! Keep track of min/max density
     minDensity = min(minDensity, density)
     maxDensity = max(maxDensity, density)

     ! Assign the density
     conservedVariables(i,1) = density

     ! Assign the mean velocity
     conservedVariables(i,2) = velocity

     ! Assign the species
     if (H2 .gt. 0) conservedVariables(i,nDimensions+2+H2) = fuel * conservedVariables(i,1)
     if (O2 .gt. 0) conservedVariables(i,nDimensions+2+O2) = oxidizer * conservedVariables(i,1)
     if (VAP .gt. 0) conservedVariables(i,nDimensions+2+VAP) = Yv * conservedVariables(i,1)
  end do
  call parallel_min(minDensity)
  call parallel_max(maxDensity)
  if (iRank .eq. iRoot) then
     print *
     print *, 'min/max density:', minDensity, maxDensity
     print *
  end if

  ! Should we add solenoidal perturbations?
  call parser_read('include solenoidal perturbations', useSolenoidalPerturbations, .false.)

  if (useSolenoidalPerturbations) then

     ! Storage for vector potential
     allocate(A(nGridPoints))
     allocate(gradA(nGridPoints, nDimensions))

     ! Vector potential
     Lpx = Lx
     Lpz = 0.0_WP
     if (nDimensions .eq. 3) Lpz = Lz

     phases = 0.0_WP
     do j = 1, nkx
        do i = -nkz, nkz
           call random_number(rand)
           rand = 2.0_WP * rand - 1.0_WP
           phases(j,i) = rand * pi
        end do
     end do

     allocate(phaseJitter(0:numBlocks, 0:numBlocks))
     phaseJitter = 0.0_WP
     do i = 0, numBlocks
        do j = 0, numBlocks
           call random_number(rand)
           phaseJitter(i,j) = rand
        end do
     end do

     z = 0.0_WP
     zint = 1
     do i = 1, nGridPoints
        ! Get the local coordinates
        x = coordinates(i, 1)
        y = coordinates(i, 2)
        z = 0.0_WP
        if (nDimensions .eq. 3) z = coordinates(i, 3) + 0.5_WP * Lz

        xint = int(x / (Lpx / real(numBlocks, WP))) 
        if (nDimensions .eq. 3) zint = int(z / (Lpz / real(numBlocks, WP)))

        ! Broadbanded scalar with Phase Jittering similar to Kim thesis and Lui thesis,
        ! except theirs were phase jitted in time.
        if (phaseJitter(xint, zint) .lt. 0.33_WP) then
           call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases + pi, nkx, nkz,   &
                scalar)
        else if (phaseJitter(xint, zint) .lt. 0.667_WP .and.                                 &
             phaseJitter(xint, zint) .ge. 0.33_WP) then
           call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases - pi, nkx, nkz,   &
                scalar)
        else
           call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases, nkx, nkz, scalar)
        end if
        A(i) = scalar
     end do

     deallocate(phaseJitter) ! ... not needed anymore

     ! Vector cross product to get a solenoidal VELOCITY field
     call gradient(A, gradA)

     allocate(velocityFluctuations(nGridPoints, 3))

     if (nDimensions .eq. 3) then
        velocityFluctuations(:, 1) = (gradA(:, 2) - gradA(:, 3))
        velocityFluctuations(:, 2) =-(gradA(:, 1) - gradA(:, 3))
        velocityFluctuations(:, 3) = (gradA(:, 1) - gradA(:, 2))
     else if (nDimensions .eq. 2) then
        velocityFluctuations(:, 1) = gradA(:, 2)
        velocityFluctuations(:, 2) =-gradA(:, 1)
     end if

     ! Rescale the velocity perturbations to have rms of x% of delta U
     midplaneRMS = 0.0_WP; midplaneCount = 0
     do i = 1, nGridPoints
        y = coordinates(i, 2)
        if (y .ge. -initialThickness .and. y .le. initialThickness) then   
           midplaneRMS = midplaneRMS + sum(velocityFluctuations(i,:)**2)
           midplaneCount = midplaneCount + 1
        end if
     end do
     call parallel_sum(midplaneRMS)
     call parallel_sum(midplaneCount)
     midplaneRMS = midplaneRMS / real(midplaneCount, WP)
     midplaneRMS = midplaneRMS**0.5

     ! Rescale the velocity
     call parser_read('shear layer rms fraction', shear_layer_initial_rms_fraction)

     velocityFluctuations = velocityFluctuations * shear_layer_initial_rms_fraction /        &
          midplaneRMS * velocityDifference

     ! Mollify the perturbations after making it solenoidal
     do i = 1, nGridPoints
        y = coordinates(i, 2)
        scalar = 0.5_WP * (tanh(5.0_WP * (y + 1.0_WP)) - (tanh(5.0_WP * (y - 1.0_WP))))
        velocityFluctuations(i,:) = velocityFluctuations(i,:) * scalar
     end do

     ! Add perturbations to the mean specified from above.
     conservedVariables(:, 2:nDimensions+1) = conservedVariables(:, 2:nDimensions+1) +       &
          velocityFluctuations(:, 1:nDimensions)

     ! Clean up
     deallocate(A)
     deallocate(gradA)
     deallocate(velocityFluctuations)
     if (allocated(Wi)) deallocate(Wi)

  end if

  do i = 1, nGridPoints

     ! Set the momentum
     conservedVariables(i,2:nDimensions+1) = conservedVariables(i,1) *                       &
          conservedVariables(i,2:nDimensions+1)

     ! Assign the energy
     conservedVariables(i, nDimensions+2) = pressure / (ratioOfSpecificHeats - 1.0_WP) +     &
          0.5_WP * conservedVariables(i,2)**2 / conservedVariables(i,1)
     if (nDimensions .gt. 1) conservedVariables(i, nDimensions+2) =                          &
          conservedVariables(i, nDimensions+2) + 0.5_WP * conservedVariables(i,3)**2 /       &
          conservedVariables(i,1)
     if (nDimensions .eq. 3) conservedVariables(i, nDimensions+2) =                          &
          conservedVariables(i, nDimensions+2) + 0.5_WP * conservedVariables(i,4)**2 /       &
          conservedVariables(i,1)

  end do

contains

  ! Routine for generating realistic velocity perturbations
  ! -------------------------------------------------------
  subroutine broadband_real(x, y, z, Lpx, Lpz, lengthScale, phases,num_kx, num_kz, val)

    ! External modules
    use math, only : pi

    implicit none

    ! Arguments
    real(WP), intent(in) :: x, y, z, Lpx, Lpz, lengthScale
    real(WP), intent(out) :: val
    integer,intent(in):: num_kx, num_kz
    real(WP), dimension(num_kx,-num_kz:num_kz),intent(in) :: phases

    ! Local variables
    real(WP) :: kxo, kzo, kyo
    integer:: i, k

    kxo = 0.0_WP; kzo = 0.0_WP
    if (Lpx .gt. 0.0_WP) kxo = 2.0_WP / (Lpx)
    if (Lpz .gt. 0.0_WP) kzo = 2.0_WP / (Lpz)
    kyo = 1.0_WP / (10.0_WP * lengthScale)
    val = 0.0_WP
    do i = 1, num_kx
       do k = -num_kz, num_kz
          val = val + cos(2.0_WP * pi * real(i, WP) * kxo * x +                              &
               2.0_WP * pi * real(k, WP) * kzo * z +                                         &
               2.0_WP * pi * kyo * y + phases(i,k)) 
       end do
    end do

    return
  end subroutine broadband_real

end subroutine shear_layer_data

subroutine shear_layer_particles

  ! Internal modules
  use shear_layer

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, j, k, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(wp) :: diameter, volume, volumeFraction, particleVolume, volumeFactor,                &
       sumParticleVolume, velocityFluctuation, Lp, Lpy, Lpz, rand, particleThickness, Tp
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
  call parser_read('particle thickness', particleThickness, initialThickness * 10)
  call parser_read('particle temperature', Tp, 2.5_WP)

  ! Initial Volume that particles occupy
  volume = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           if (abs(coordinates(grid_index(i,j,k), 2)) .le. 0.5_WP * particleThickness)       &
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
           particles(i)%position(1) = Lx * rand
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = particleThickness * (rand - 0.5_WP)
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = Lz * (rand - 0.5_WP)
        end if

     case ('uniform')
        ix = (particleOffset + i - 1) / (nparty * npartz)
        iy = (particleOffset + i - 1 - nparty * npartz * ix) / npartz
        iz = particleOffset + i - 1 - nparty * npartz * ix - npartz * iy
        particles(i)%position(1) = (ix + 0.5_WP) * Lp
        particles(i)%position(2) = (iy + 0.5_WP) * Lpy - 0.5_WP * particleThickness
        if (nz .gt. 1) particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz

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
  volumeFraction = sumParticleVolume / volume

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
     print *
  end if

  return
end subroutine shear_layer_particles
