module ibm_part

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
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, maxDxDyDz

end module ibm_part


subroutine ibm_part_grid

  ! Internal modules
  use ibm_part

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max
  real(WP) :: ymin, ymax, ymin_uni, ymax_uni, dy_min, dy_max
  real(WP) :: zmin, zmax, zmin_uni, zmax_uni, dz_min, dz_max
  real(WP) :: cornerMin(3), cornerMax(3), cornerMinU(3), cornerMaxU(3), hMin(3), hMax(3)
  real(WP), dimension(:), allocatable :: x1D, y1D, z1D
  logical :: stretch


  call parser_read('stretch grid', stretch, .false.)

  if (stretch) then

     ! Should not be used for periodic cases
     do i = 1, nDimensions
        if (periodicityType(i) .ne. NONE)                                                    &
             call die('Stretched grid only for non-periodic domains')
     end do

     call parser_read('corner min', cornerMin)
     call parser_read('corner max', cornerMax)
     call parser_read('corner uniform min', cornerMinU)
     call parser_read('corner uniform max', cornerMaxU)
     call parser_read('h min', hMax)
     call parser_read('h max', hMin)

     ! Create x
     nx = 1
     if (globalGridSize(1) .gt. 1) then
        xmin = min(cornerMin(1), cornerMax(1))
        xmax = max(cornerMin(1), cornerMax(1))
        xmin_uni = min(cornerMinU(1), cornerMaxU(1))
        xmax_uni = max(cornerMinU(1), cornerMaxU(1))
        dx_min = min(hMin(1), hMax(1))
        dx_max = max(hMin(1), hMax(1))
        call grid_stretch(xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max, x1D)
        nx = size(x1d)
        Lx = xmax - xmin
     end if

     ! Create y
     ny = 1
     if (globalGridSize(2) .gt. 1) then
        ymin = min(cornerMin(2), cornerMax(2))
        ymax = max(cornerMin(2), cornerMax(2))
        ymin_uni = min(cornerMinU(2), cornerMaxU(2))
        ymax_uni = max(cornerMinU(2), cornerMaxU(2))
        dy_min = min(hMin(2), hMax(2))
        dy_max = max(hMin(2), hMax(2))
        call grid_stretch(ymin, ymax, ymin_uni, ymax_uni, dy_min, dy_max, y1D)
        ny = size(y1d)
        Ly = ymax - ymin
     end if

     ! Create z
     nz = 1
     if (globalGridSize(3) .gt. 1) then
        zmin = min(cornerMin(3), cornerMax(3))
        zmax = max(cornerMin(3), cornerMax(3))
        zmin_uni = min(cornerMinU(3), cornerMaxU(3))
        zmax_uni = max(cornerMinU(3), cornerMaxU(3))
        dz_min = min(hMin(3), hMax(3))
        dz_max = max(hMin(3), hMax(3))
        call grid_stretch(zmin, zmax, zmin_uni, zmax_uni, dz_min, dz_max, z1D)
        nz = size(z1d)
        Lz = zmax - zmin
     end if

     ! Output new grid size
     if (nx.gt.1) print *, 'xmin: ', xmin, ', xmax: ', xmax
     if (ny.gt.1) print *, 'ymin: ', ymin, ', ymax: ', ymax
     if (nz.gt.1) print *, 'zmin: ', zmin, ', zmax: ', zmax
     print *, 'nx: ', nx, ', ny: ', ny,', nz: ', nz
     if (nx.gt.1) print *, 'x_start: ',x1D(1), ' x_end: ', x1D(nx)
     if (ny.gt.1) print *, 'y_start: ',y1D(1), ' y_end: ', y1D(ny)
     if (nz.gt.1) print *, 'z_start: ',z1D(1), ' z_end: ', z1D(nz)

     ! Update grid metrics
     globalGridSize(1) = nx
     globalGridSize(2) = ny
     globalGridSize(3) = nz
     call grid_cleanup
     call operator_cleanup
     call geometry_setup
     call grid_setup
     call operator_setup

     ! Generate the grid
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              if (nx.gt.1) coordinates(grid_index(i,j,k), 1) = x1D(i)
              if (ny.gt.1) coordinates(grid_index(i,j,k), 2) = y1D(j)
              if (nz.gt.1) coordinates(grid_index(i,j,k), 3) = z1D(k)
           end do
        end do
     end do

     dx = dx_min
     if (nDimensions .gt. 1) dx = min(dx_min, dy_min)
     if (nDimensions .gt. 2) dx = min(dx_min, dy_min, dz_min)

  else

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

  end if

  ! Setup the grid metrics
  call grid_metrics_setup

contains

  ! --------------------------------------------------------- !
  ! Generate a stretched grid (in 1D)                         !
  !                                                           !
  !   dx_max                 dx_min                   dx_max  !
  !  |-----|---|--|-|-|-|-|-|-|-|-|-|-|-|-|-|-|--|---|-----|  !
  ! xmin        xmin_uni                   xmax_uni      xmax !
  ! --------------------------------------------------------- !

  subroutine grid_stretch(xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max, x)

    implicit none

    ! Arguments
    real(WP), intent(in) :: xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max
    real(WP), dimension(:), allocatable, intent(out) :: x

    ! Local variables
    integer :: nu, ns1, ns2, nx
    real(WP) :: ratio1, ratio2, ratio, sc

    ! sanity check
    if (dx_min .gt. dx_max) call die('dx_min is larger than dx_max')
    if (xmin_uni .gt. xmax_uni) call die('xmin_uni is larger than xmax_uni')

    sc = min(dx_min, abs(xmin_uni-xmin), abs(xmax-xmax_uni), abs(xmax_uni-xmin_uni))
    if (abs(sc-dx_min) .gt. epsilon(1.0_WP)) call die('dx_min too large')

    sc = min(dx_max, abs(xmin_uni-xmin), abs(xmax-xmax_uni), abs(xmax_uni-xmin_uni))
    if (abs(sc-dx_max) .gt. epsilon(1.0_WP)) call die('dx_max too large')

    ! Determine number of grid points and stretching rate
    nu = ceiling((xmax_uni-xmin_uni) / dx_min) ! Number of uniform grid
    if (abs(dx_min-dx_max) .lt. epsilon(1.0_WP)) then
       ratio1 = 1.0_WP; ratio2 = 1.0_WP
       ns1 = ceiling((xmin_uni-xmin) / dx_min) + 1
       ns2 = ceiling((xmax-xmax_uni) / dx_min) + 1
    else
       ! Stretching ratio on the left
       ratio1 = (xmin_uni-xmin)/(xmin_uni-xmin-dx_max+dx_min)
       ! Stretching ratio on the right
       ratio2 = (xmax-xmax_uni)/(xmax-xmax_uni-dx_max+dx_min)

       ! Preserve the same ratio around the uniform zone
       ratio = min(ratio1,ratio2)

       !nb of points from xmin to xmin_uni (xmin_uni included)
       ns1 = ceiling(1  + log(1+(ratio-1)*(xmin_uni-xmin)/(dx_min*ratio))/log(ratio))
       !nb of points from xmax to xmax_uni included (xmax_uni included)
       ns2 = ceiling(1  + log(1+(ratio-1)*(xmax-xmax_uni)/(dx_min*ratio))/log(ratio))

       ! Update (instead of using a newton algorithm to get the exact ratio)
       ratio = (dx_max/dx_min)**(1/(real(max(ns1,ns2)-1,WP)))
       ns1 = ceiling(1  + log(1+(ratio-1)*(xmin_uni-xmin)/(dx_min*ratio))/log(ratio))
       ns2 = ceiling(1  + log(1+(ratio-1)*(xmax-xmax_uni)/(dx_min*ratio))/log(ratio))

    end if

    ! Total number of grid points
    nx = nu + ns1 + ns2 - 1

    ! Generate the grid
    allocate(x(nx))

    ! left
    x(ns1) = xmin_uni
    do i = 1,ns1-1
       x(ns1-i) = x(ns1-i+1) - dx_min * ratio**real(i,WP)
    end do

    ! uniform
    do i=ns1+1,ns1+nu-1
       x(i) = x(i-1) + dx_min
    end do
    x(ns1+nu) = xmax_uni

    !right
    do i = ns1+nu+1,nx
       x(i) = x(i-1) + dx_min * ratio**real(i-(ns1+nu),WP)
    end do

    !adjust the boundaries
    x(1) = xmin
    x(nx) = xmax

    !prevent the first and last interval to be too small or larger than dx_max
    if (x(2)-xmin < x(3)-x(2)) then
       x(2) = (xmin + x(3)) / 2.0_WP
    end if

    if (xmax-x(nx-1) < x(nx-1)-x(nx-2)) then
       x(nx-1) = (xmax + x(nx-2)) / 2.0_WP
    end if

    return
  end subroutine grid_stretch

end subroutine ibm_part_grid


subroutine ibm_part_data

  ! Internal modules
  use ibm_part

  ! External modules
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, x, x0, p1, p2, rho1, rho2, u1, u2, density, pressure, velocity
  logical :: useShock

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Uniform flow or shock?
  call parser_read('include shock', useShock, .false.)

  ! Specific heat ratio
  call parser_read('ratio of specific heats', gamma, 1.4_WP)

  if (useShock) then
     ! Read in shock parameters
     call parser_read('ratio of specific heats', gamma, 1.4_WP)
     call parser_read('shock position', x0, 0.0_WP)
     call parser_read('post-shock pressure', p1, 1.0_WP / gamma)
     call parser_read('pre-shock pressure', p2, 1.0_WP / gamma)
     call parser_read('post-shock density', rho1, 1.0_WP)
     call parser_read('pre-shock density', rho2, 1.0_WP)
     call parser_read('post-shock velocity', u1, 0.0_WP)
     call parser_read('pre-shock velocity', u2, 0.0_WP)
  else
     call parser_read('crossflow mach number', velocity)
     call parser_read('density', density, 1.0_WP)
     pressure = density / gamma
  end if

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           if (useShock) then
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
end subroutine ibm_part_data


subroutine ibm_part_objects

  ! Internal modules
  use ibm_part

  ! External modules
  use random
  use math
  use parallel
  use ibm

  implicit none

  ! Local variables
  integer :: i, j, k, npartx, nparty, npartz, ix, iy, iz, ierr
  integer :: iter, iterMax
  real(WP) :: volume,volumeFraction, volumeFraction1, volumeFraction2, particleVolume,       &
       particleVolume1, particleVolume2, particleThickness, volumeFactor, sumParticleVolume, &
       Lp, Lpy, Lpz, x0, rand, r, distance(3), dMean, dStd, dMin, dMax, dShift, d1, d2,      &
       particleSpacing, coordCenter(nDimensions), inBound(1:nDimensions),                    &
       maxDistCenter(1:nDimensions), lShift, minDistances
  character(len = str_medium) :: filename, particleDistribution
  logical :: success, preventOverlap, setSeed
  integer :: mySeed(33)=123456789

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

  ! Read in the particle parameters
  call parser_read('particle layer thickness', particleThickness, Lx) 
  call parser_read('particle layer position', x0, 0.0_WP)

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')
  call parser_read('prevent particle overlap', preventOverlap, .true.)

  ! Factor for computing volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if

  ! Read particle diameter
  if (trim(adjustl(particleDistribution)) .eq. 'bidisperse') then
     call parser_read('particle diameter 1', d1)
     call parser_read('particle diameter 2', d2)
     particleVolume1 = pi * volumeFactor * d1 ** nDimensions
     particleVolume2 = pi * volumeFactor * d2 ** nDimensions
  else
     call parser_read('particle mean diameter', dMean)
     call parser_read('particle std diameter', dStd, 0.0_WP)
     call parser_read('particle min diameter', dMin, 0.0_WP)
     call parser_read('particle max diameter', dMax, 0.0_WP)
     call parser_read('particle diameter shift', dShift, 0.0_WP)
     particleVolume = pi * volumeFactor * dMean ** nDimensions
  end if

  ! Initialize the random number generator
  call parser_read('set seed', setSeed, .false.)
  if (.not. setSeed) then
     call random_setup
  else
     call random_seed(put=mySeed)
  end if

  ! Initial Volume that particles occupy
  if (nz.gt.1) then
     volume = particleThickness * Ly * Lz
  else
     volume = particleThickness * Ly
  end if

  select case (trim(particleDistribution))

  case ('single')
     nParticles = 1

  case ('double')
     nParticles = 2

  case ('random', 'random neighbor exploration')
     ! Get the number of particles based on the volume fraction and distribute to processors
     call parser_read('particle volume fraction', volumeFraction)
     nParticles = int(volumeFraction * volume / particleVolume)

  case ('uniform', 'uniform perturbed')
     ! Mean interparticle distance
     call parser_read('particle volume fraction', volumeFraction)
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     npartx = int(particleThickness / Lp)
     Lp = particleThickness / real(npartx, WP)
     if (ny .gt. 1) nparty = int(Ly / Lp)
     Lpy = Ly / real(nparty, WP)
     if (nz .gt. 1) npartz = int(Lz / Lp)
     Lpz = Lz / real(npartz, WP)
     nParticles = npartx * nparty * npartz

  case ('bidisperse')
     call parser_read('particle volume fraction 1', volumeFraction1)
     call parser_read('particle volume fraction 2', volumeFraction2)
     nParticles1 = int(volumeFraction1 * volume / particleVolume1) 
     nParticles2 = int(volumeFraction2 * volume / particleVolume2)
     nParticles = nParticles1 + nParticles2

  case default

     call die("Unknown particle distribution '" // trim(particleDistribution) // "'")

  end select

  ! Allocate the object vector
  allocate(particles(nParticles))
  
  ! Only root process generates IBM particle data
  if (iRank .eq. iRoot) then

     ! Initialize the particle parameters
     sumParticleVolume = 0.0_WP
     do i = 1, nParticles

        if(trim(particleDistribution) .eq. 'bidisperse') then
           if (i .lt. nParticles1) then
              particles(i)%diameter = d1
           else
              particles(i)%diameter = d2
           end if
        else
           ! Set particle size distribution (compact support lognormal)
           particles(i)%diameter = random_lognormal(m=dMean,sd=dStd) + dShift
           if (dStd.gt.0.0_WP) then
              do while (particles(i)%diameter.gt.dMax+epsilon(1.0_WP) .or.                   &
                   particles(i)%diameter.lt.dMin-epsilon(1.0_WP))
                 particles(i)%diameter=random_lognormal(m=dMean,sd=dStd) + dShift
              end do
           else
              particles(i)%diameter = dMean
           end if
        end if

        ! Position (this will be reset later)
        particles(i)%position = 0.0_WP

        ! Particle velocity
        particles(i)%velocity = 0.0_WP

        ! Sum particle volume
        sumParticleVolume = sumParticleVolume + pi * volumeFactor *                          &
             particles(i)%diameter ** nDimensions
     end do

     ! Compute effective volume fraction
     volumeFraction = sumParticleVolume / volume

     ! Distribute the particles
     select case (trim(particleDistribution))

     case ('single')
        call parser_read('particle center', particles(1)%position(1:nDimensions))
        call parser_read('particle velocity', particles(1)%velocity(1:nDimensions))

     case ('double')
        call parser_read('particle 1 center', particles(1)%position(1:nDimensions))
        call parser_read('particle 2 center', particles(2)%position(1:nDimensions))
        call parser_read('particle 1 velocity', particles(1)%velocity(1:nDimensions))
        call parser_read('particle 2 velocity', particles(2)%velocity(1:nDimensions))

     case ('random', 'bidisperse')
        i=1
        do while (i.le.nParticles)

           ! Give the particle a random position
           if (nx .gt. 1) then
              call random_number(rand)
              particles(i)%position(1) = particleThickness * rand + x0
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
                       distance(k) = min(distance(k), periodicLength(k) -                    &
                            abs(particles(i)%position(k) - particles(j)%position(k)))
                    end if
                 end do
                 r = sqrt(sum(distance**2))
                 particleSpacing = 0.5_WP*(particles(i)%diameter+particles(j)%diameter)      &
                      + 3.0_WP * maxDxDyDz
                 if (setSeed) then
                    particleSpacing = 0.5375_WP*(particles(i)%diameter+particles(j)%diameter)
                 end if
                 if (r.le.particleSpacing) then
                    i=i-1
                    success = .false.
                    exit part2
                 end if
              end do part2
           end if
           i=i+1

           if (success .and. modulo(i,1000).eq.0)                                            &
                print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'%'

        end do

     case ('random neighbor exploration')
        coordCenter(1)   = Lx/2.0_WP
        maxDistCenter(1) = Lx/2.0_WP
        if (ny .gt. 1) then
           coordCenter(2)   = 0.0_WP
           maxDistCenter(2) = Ly/2.0_WP
        end if
        if (nz .gt. 1) then
           coordCenter(3)   = 0.0_WP
           maxDistCenter(3) = Lz/2.0_WP
        end if

        ! parameters requiring tuning - need a thorough study
        ! by default, no exploration of the neighborhood
        iterMax = 0
        ! 2D
        if ( (nDimensions .eq. 2) .and. (volumeFraction .gt. 0.44_WP )) then
           iterMax = 16
        end if
        ! 3D
        if (nDimensions .eq. 3) then
           if (volumeFraction .gt. 0.28_WP) iterMax = 16
           if (volumeFraction .gt. 0.33_WP) iterMax = 512
        end if

        i=1
        do while (i.le.nParticles)

           ! Give the particle a random position
           if (nx .gt. 1) then
              call random_number(rand)
              particles(i)%position(1) = Lx * rand
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
           if (preventOverlap) then
              j = 1
              iter = 0
              part3: do while (j.le.i-1)

                 ! Compute separation distance (account for periodicity)
                 distance = particles(i)%position(1:nDimensions) -                           &
                      particles(j)%position(1:nDimensions)
                 do k=1, nDimensions
                    if (abs(distance(k))>periodicLength(k)-abs(distance(k))) then
                       distance(k) = particles(i)%position(k)-particles(j)%position(k) -     &
                            sign(periodicLength(k), particles(i)%position(k) -               &
                            particles(j)%position(k))
                    end if
                 end do

                 ! compute distance
                 r = sqrt(sum(distance**2))

                 ! compute interparticle spacing
                 particleSpacing = 0.5_WP*(particles(i)%diameter+particles(j)%diameter)      &
                      + 3.0_WP * maxDxDyDz
                 if (setSeed) then
                    ! maybe we need to set the need the constant in the input file? 
                    particleSpacing = 0.5375_WP*(particles(i)%diameter+particles(j)%diameter)
                 end if

                 if ((r.lt.particleSpacing) .and. (r.gt.epsilon(1.0_WP))) then
                    ! compute minimum shifting length possible
                    ! eps is required because distance/norm(distance) may not be exactly of norm 1
                    ! eps need to be adjusted for a different computer/architecture?
                    ! real(iter/iterMax,WP) is used to accelerate the discovery of a new available spot
                    !lShift = (1+0.5*real(iter/iterMax,WP))*(epsilon(1.0_WP) + particleSpacing - r)
                    lShift = (1+0.5*real(iter/iterMax,WP))*(1.0E-10_WP + particleSpacing - r)

                    ! shift the particle
                    particles(i)%position(1:nDimensions) =                                   &
                         particles(i)%position(1:nDimensions) + lShift * distance / r

                    ! check if the shifted particle is inside the domain
                    inBound = particles(i)%position(1:nDimensions) - coordCenter(1:nDimensions)
                    do k=1, nDimensions
                       if (abs(inBound(k))>maxDistCenter(k)) then
                          particles(i)%position(k) = particles(i)%position(k)                &
                               - sign(periodicLength(k), inBound(k))
                       end if
                    end do

                    ! reset j to test the new location against all the other particles
                    j = 1

                    ! increment number of shifting
                    iter = iter + 1

                    ! check maximum number of shifting for this particle (for pathological case)
                    if(iter>iterMax) then
                       j = nParticles
                       i = i-1
                    end if
                 else if (r.le.epsilon(1.0_WP)) then
                    j = nParticles
                    i = i-1
                 else
                    j = j + 1
                 end if
              end do part3
           end if
           i = i + 1

           if (modulo(i,10).eq.0) print *, real(i,SP)/real(nParticles,SP)*100.0_SP,'%'

        end do

        ! Check if the distribution is fine (minimal distance should be larger than particle interspacing)
        minDistances = huge(1.0_WP)
        do i=1, nParticles
           do j=1,i-1
              do k=1, nDimensions
                 distance(k) = abs(particles(i)%position(k)-particles(j)%position(k))
                 distance(k) = min(distance(k), periodicLength(k)-distance(k))
              end do
              r = sqrt(sum(distance**2))
              minDistances = min(minDistances, r)
           end do
        end do
        print*, "The minimal distance is : ", minDistances

     case ('uniform')
        do i = 1, nParticles
           ix = (i - 1) / (nparty * npartz)
           iy = (i - 1 - nparty * npartz * ix) / npartz
           iz = i - 1 - nparty * npartz * ix - npartz * iy
           particles(i)%position(1) = (ix + 0.5_WP) * Lp + x0
           particles(i)%position(2) = (iy + 0.5_WP) * Lpy - 0.5_WP * Ly
           particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz
        end do

     case ('uniform perturbed')
        do i = 1, nParticles
           ix = (i - 1) / (nparty * npartz)
           iy = (i - 1 - nparty * npartz * ix) / npartz
           iz = i - 1 - nparty * npartz * ix - npartz * iy
           particles(i)%position(1) = (ix + 0.5_WP) * Lp + x0
           particles(i)%position(2) = (iy + 0.5_WP) * Lpy - 0.5_WP * Ly
           particles(i)%position(3) = (iz + 0.5_WP) * Lpz - 0.5_WP * Lz
           if (nx .gt. 1) then
              call random_number(rand)
              rand = 2.0_WP * rand - 1.0_WP
              particles(i)%position(1) = particles(i)%position(1) + rand *                   &
                   (0.5_WP * Lp - 0.5_WP * particles(i)%diameter - 2.0_WP * dx)
           end if
           if (ny .gt. 1) then
              call random_number(rand)
              rand = 2.0_WP * rand - 1.0_WP
              particles(i)%position(2) = particles(i)%position(2) + rand *                   &
                   (0.5_WP * Lpy - 0.5_WP * particles(i)%diameter - 2.0_WP * dy)
           end if
           if (nz .gt. 1) then
              call random_number(rand)
              rand = 2.0_WP * rand - 1.0_WP
              particles(i)%position(3) = particles(i)%position(3) + rand *                   &
                   (0.5_WP * Lpz - 0.5_WP * particles(i)%diameter - 2.0_WP * dz)
           end if
        end do

     end select

     ! Output stuff to the screen
     print *
     print *, 'Number of particles: ', nParticles
     print *, 'Volume fraction: ', volumeFraction
     print *
     !print *, 'Particles coordinates: '
     !do i=1, nParticlesGlobal
     !   print *,i, particles(i)%position(1), particles(i)%position(2), particles(i)%position(3)
     !end do
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
        object(i)%remove = .false.
     end do
  end if
  call MPI_BCAST(object, nObjects, MPI_OBJECT, iRoot, MPI_COMM_WORLD, ierr)

  return
end subroutine ibm_part_objects
