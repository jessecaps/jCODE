module jet_crossflow

  ! External modules
  use precision
  use string
  use parser
  use random
  use simulation_flags
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Global variables
  integer :: nx, ny, nz
  integer :: iJet1, iJet2, kJet1, kJet2
  real(WP) :: Lx, Ly, Lz, dx, dy, dz
  real(WP) :: jetDiameter, jetPosition
  logical :: includeSandpaper, conformToJet

end module jet_crossflow

subroutine jet_crossflow_grid

  ! Internal modules
  use jet_crossflow

  ! External modules
  use parallel

  implicit none

  ! Jet in crossflow schematic and the corresponding indeces 
  !
  !
  !                      _~~~~~~~~~~~~~~~~~~~~~~~~~~_
  ! j = 1  > ___________|                            |_________________|___:___|__________
  !
  !          ^          ^                            ^                 ^       ^         ^
  !         i=1        i=iTrip1                    i=iTrip2          i=iJet1  i=iJet2  i=nx

  ! Local variables
  integer :: i, j, k, n, nGrit, iTrip1, iTrip2
  real(WP) :: ytilde, r, tripLocation, tripWidth, tripHeight, gritHeight, gritWidth,         &
       sig, x, y, z, x0, z0, z12, amp, gauss, alpha, totalHeight, peakHeight, valleyHeight,  &
       shift, rand, norm(3)
  real(WP), dimension(:), allocatable :: deltaY
  logical :: stretchGrid, includeGrit

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

  ! Conform the mesh to the jet perimeter
  call parser_read('conform grid to jet', conformToJet, .false.)
  if (nz .eq. 1 .and. conformToJet) then
     write (*,'(A)') ''
     write (*,'(A)') "Warning, 'conform to jet' only enabled for 3D geometries"
     conformToJet = .false.
  end if

  ! Get jet parameters
  call parser_read('jet diameter', jetDiameter)
  call parser_read('jet position', jetPosition)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           coordinates(grid_index(i,j,k), 1) = Lx *  real(i - 1, WP) / real(nx - 1, WP)

           ! Create Y
           if (stretchGrid) then
              ytilde = real(ny-j, WP) / real(ny-1, WP)
              coordinates(grid_index(i,j,k), 2) = Ly * (1.0_WP - tanh(r * ytilde) / tanh(r))
           else
              coordinates(grid_index(i,j,k), 2) = Ly * real(j-1, WP) / real(ny-1, WP)
           end if

           ! Create Z
           if (nz .gt. 1) then
              coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /             &
                   real(nz - 1, WP) - 0.5_WP * Lz
           end if

        end do
     end do
  end do

  ! Smoothly adjust the grid near the jet
  if (conformToJet) then
     sig = jetDiameter / 4.0_WP
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              x = coordinates(grid_index(i,j,k), 1)
              y = coordinates(grid_index(i,j,k), 2)
              z = coordinates(grid_index(i,j,k), 3)
              norm(1) = x - jetPosition
              norm(2) = 0.0_WP
              norm(3) = z
              r = sqrt(sum(norm**2))
              norm = norm / (r + epsilon(1.0_WP))
              shift = r * exp(-0.5_WP * r**2 / sig**2)
              coordinates(grid_index(i,j,k), 1:3) = coordinates(grid_index(i,j,k), 1:3) +    &
                   shift * norm
           end do
        end do
     end do
  end if

  ! Generate the sandpaper
  call parser_read('include sandpaper', includeSandpaper, .false.)
  if (includeSandpaper) then
     call parser_read('sandpaper location', tripLocation)
     call parser_read('sandpaper width', tripWidth)
     call parser_read('sandpaper height', tripHeight)
     call parser_read('include sandpaper grit', includeGrit, .false.)
     if (includeGrit) then
        call parser_read('sandpaper grit height', gritHeight)
        call parser_read('sandpaper grit width', gritWidth)
        call parser_read('number of sandpaper particles', nGrit)
     end if

     ! Find the sandpaper extents
     iTrip1 = 1; iTrip2 = nx
     j = iStart(2); k = iStart(3)
     do i = iStart(1), iEnd(1)
        if (coordinates(grid_index(i,j,k), 1) .lt. tripLocation) iTrip1 = i + 1
     end do
     do i = iEnd(1), iStart(1), -1
        if (coordinates(grid_index(i,j,k),1) .gt. tripLocation + tripWidth) iTrip2 = i - 1
     end do
     call parallel_max(iTrip1)
     call parallel_min(iTrip2)
     if (iRank .eq. iRoot) then
        write (*,*)
        write (*,'(A, i0.0, A, i0.0, A)') 'Sandpaper extents: [', iTrip1, ', ', iTrip2, ']'
     end if

     ! Get unperturbed height variation
     allocate(deltaY(ny), source = 0.0_WP)
     i = iStart(1); k=iStart(3)
     do j = max(2, iStart(2)), iEnd(2)
        deltaY(j) = coordinates(grid_index(i,j,k), 2) - coordinates(grid_index(i,j-1,k), 2)
     end do

     ! Deform the mesh to the sandpaper height
     do k = iStart(3), iEnd(3)
        do i = iStart(1), iEnd(1)
           j = 1

           ! Create a smooth step
           sig = 20.0_WP
           x = coordinates(grid_index(i,j,k), 1)
           coordinates(grid_index(i,j,k), 2) = 0.5_WP * tripHeight *                         &
                (tanh(sig * (x - tripLocation)) - tanh(sig*(x - tripLocation - tripWidth)))

           ! Shift grid points above (smoothly)
           do j = 2, ny

              ! Get the current height
              y = coordinates(grid_index(i,j,k), 2)

              ! Adjust the current height
              ytilde = coordinates(grid_index(i,j-1,k),2) + deltaY(j)
              alpha = tanh(0.1_WP * real(j - 2, WP) / real(ny - 2, WP))
              coordinates(grid_index(i,j,k), 2) = ytilde * (1.0_wp - alpha) + y * alpha
           end do
        end do
     end do

     ! Embed the particles
     If (includeGrit) then

        ! Initialize the random number generator
        call random_setup

        ! Standard deviation
        sig = gritWidth

        ! Loop through number of particles.
        write (*,'(A)') ''
        write (*,'(A, i0.0, A)') 'Adding ', nGrit, ' particles...'
        do n = 1, nGrit

           ! Compute amplitude
           amp = gritHeight
           call random_number(rand)
           if (rand .lt. 0.5_WP) amp = -amp

           ! Get a random location
           call random_number(rand)
           x0 = tripLocation + 5.0_WP * sig + (tripWidth - 10.0_WP * sig) * rand
           call random_number(rand)
           z0 = (Lz - dz) * rand
           if (nz .eq. 1) z0 = 0.0_WP

           ! Modify the grid
           do k = iStart(3), iEnd(3)
              do i = max(iTrip1, iStart(1)), min(iTrip2, iEnd(1))
                 j = 1

                 ! Get the coordinates
                 x = coordinates(grid_index(i,j,k), 1)
                 y = coordinates(grid_index(i,j,k), 2)
                 z = 0.0_WP
                 if (nz .gt. 1) z = coordinates(grid_index(i,j,k), 3)

                 ! Represent sandpaper particles as Gaussian
                 gauss = amp * exp(-((x - x0)**2 / (2.0_WP * sig**2) +                       &
                      (z - z0)**2 / (2.0_WP * sig**2)))

                 ! Account for periodicity in z
                 z12 = Lz - abs(z - z0)
                 If (nz .gt. 1) gauss = gauss + amp *                                        &
                      exp(-((x - x0)**2 / (2.0_WP * sig**2) + z12**2 / (2.0_WP * sig**2)))

                 ! Update the vertical coordinate
                 coordinates(grid_index(i,j,k), 2) = y + gauss

                 ! Shift grid points above (smoothly)
                 do j = 2, ny

                    ! Get the current height
                    y = coordinates(grid_index(i,j,k), 2)

                    ! Adjust the current height
                    ytilde = coordinates(grid_index(i,j-1,k),2) + deltaY(j)
                    alpha = tanh(4.0_WP * real(j - 2, WP) / real(ny - 2, WP))
                    coordinates(grid_index(i,j,k), 2) = ytilde * (1.0_WP - alpha) + y * alpha
                 end do

              end do
           end do

           ! Output the progress
           if (mod(real(n, WP), real(nGrit, WP) / 10.0_WP) .le. 1.0E-9_WP .and.              &
                iRank .eq. iRoot) write(*,'(f5.1, A)') real(n, WP) / real(nGrit, WP) *       &
                100.0_WP, '% complete'

        end do
     end if

     ! Output surface roughness dimentions
     totalHeight = 0.0_WP
     peakHeight = 0.0_WP
     valleyHeight = huge(1.0_WP)
     j = 1
     do k = iStart(3), iEnd(3)
        do i = max(iTrip1, iStart(1)), min(iTrip2, iEnd(1))
           y = coordinates(grid_index(i,j,k), 2)
           totalHeight = max(totalHeight, y)
           peakHeight = max(peakHeight, y - tripHeight)
           valleyHeight = min(valleyHeight, y - tripHeight)
        end do
     end do
     call parallel_max(totalHeight)
     call parallel_max(peakHeight)
     call parallel_min(valleyHeight)
     if (iRank .eq. iRoot) then
        write (*,'(A)') ''
        print *, 'Max surface height:', real(totalHeight, 4)
        print *, 'Peak height:', real(peakHeight, 4)
        print *, 'Valley height:', real(valleyHeight, 4)
     end if

  end if !... if (includeSandpaper) then

  ! Find the jet extents
  iJet1 = nx; iJet2 = 1
  kJet1 = nz; kJet2 = 1
  j = 1
  do k = iStart(3), iEnd(3)
     do i = iStart(1), iEnd(1)
        x = coordinates(grid_index(i,j,k), 1)
        z = 0.0_WP
        if (nz .gt. 1) z = coordinates(grid_index(i,j,k), 3)
        r = sqrt((x - jetPosition)**2 + z**2)
        if (r .le. 0.5_WP * jetDiameter) then
           iJet1 = min(iJet1, i)
           iJet2 = max(iJet2, i)
           kJet1 = min(kJet1, k)
           kJet2 = max(kJet2, k)
        end if
     end do
  end do
  call parallel_min(iJet1)
  call parallel_min(kJet1)
  call parallel_max(iJet2)
  call parallel_max(kJet2)

  ! Stop if running in parallel
  if (conformToJet .and. nProcs .gt. 1)                                                      &
       call die('Conform grid to jet only implemented in serial!')

  if (iRank .eq. iRoot) then
     print *
     print *, 'Jet extents: [', iJet1, ',', iJet2, '] x [',kJet1, ',', kJet2,']'
  end if

  ! Setup the grid metrics
  call grid_metrics_setup

contains

  ! Mapping functional for grid stretching
  ! --------------------------------------
  subroutine mapping_function(s, b, c, sigma, x0, g)

    ! External modules
    use math, only : pi

    implicit none

    real(WP), intent(in) :: s(:), b, c, sigma, x0
    real(WP), intent(out) :: g(size(s))

    g = ((s - 0.5_WP + x0) * (1.0_WP + 2.0_WP * b) - b * sigma *                             &
         (exp(- ((s - 0.5_WP + x0 + c) / sigma) ** 2) / sqrt(pi) +                           &
         ((s - 0.5_WP + c) / sigma) * erf((s - 0.5_WP + x0 + c) / sigma) -                   &
         exp(- ((s - 0.5_WP + x0 - c) / sigma) ** 2) / sqrt(pi) -                            &
         ((s - 0.5_WP +x0 - c) / sigma) * erf((s - 0.5_WP + x0 - c) / sigma))) /             &
         (0.5_WP + b - b * sigma * (exp(- ((0.5_WP + c) / sigma) ** 2) /                     &
         sqrt(pi) + ((0.5_WP + c) / sigma) * erf((0.5_WP + c) / sigma) -                     &
         exp(- ((0.5_WP - c) / sigma) ** 2) / sqrt(pi) - ((0.5_WP - c) / sigma) *            &
         erf((0.5_WP - c) / sigma)))

    return
  end subroutine mapping_function

end subroutine jet_crossflow_grid


subroutine jet_crossflow_data

  ! Internal modules
  use jet_crossflow

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use random
  use solver_options
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, kk, H2, N2, O2
  real(WP) :: x0, U0, P0, x, y, z, xx, r, Re, delta, eta, sig, a, yDecay, jetVelocity,       &
       Yf0, Yo0, density, velocity(3), pressure, temperature, fuel, oxidizer, inert
  real(WP), dimension(:), allocatable :: Wi
  character(len = str_medium) :: velocityProfile, jetShape
  logical :: insideJet

  ! Solution to Blasius boundary layer
  real(WP) :: blasius0, blasius1
  real(WP) :: f2l, f2r, f0l, f0r
  real(WP), dimension(0:9) :: by0 = (/                                                       &
       0.0_8, 0.165571818583440_8, 0.650024518764203_8, 1.39680822972500_8,                  &
       2.30574664618049_8, 3.28327391871370_8, 4.27962110517696_8,                           &
       5.27923901129384_8, 6.27921363832835_8, 7.27921257797747_8 /)
  real(WP), dimension(0:9) :: by1 = (/                                                       &
       0.0_8, 0.329780063306651_8, 0.629765721178679_8, 0.84604458266019_8,                  &
       0.95551827831671_8, 0.99154183259084_8, 0.99897290050990_8,                           &
       0.9999216098795_8, 0.99999627301467_8, 0.99999989265063_8 /)
  real(WP), dimension(0:9) :: by2 = (/                                                       &
       0.332057384255589_8, 0.323007152241930_8, 0.266751564401387_8, 0.161360240845588_8,   &
       0.06423404047594_8, 0.01590689966410_8, 0.00240199722109_8,                           &
       0.00022016340923_8, 0.00001224984692_8, 0.00000041090325_8 /)

  ! Get the number of species
  call parser_read('number of species', nSpecies, 0)

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2 + nSpecies
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Setup solver_options to get species data
  call solver_options_setup

  ! Get species indices
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
           case default
              call die("shear_layer_data: unknown species: '" //                             &
                   trim(speciesName(k)) // "!")

           end select
        end do
     else
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

  ! Velocity properties
  call parser_read('free stream velocity', U0)
  call parser_read('blasius virtual origin', x0)
  call parser_read('Reynolds number', Re)
  call parser_read('jet velocity', jetVelocity)
  call parser_read('jet velocity profile', velocityProfile)
  call parser_read('jet shape', jetShape, 'square')

  ! Mixture properties
  call parser_read('initial fuel mass fraction', Yf0)
  call parser_read('initial oxidizer mass fraction', Yo0)

  ! Set constant pressure
  P0 = 1.0_WP / ratioOfSpecificHeats

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Local coordinates
           x = coordinates(grid_index(i,j,k), 1)
           y = coordinates(grid_index(i,j,k), 2)
           z = coordinates(grid_index(i,j,k), 3)

           ! Initialize the mixture
           velocity = 0.0_WP
           pressure = 1.0_WP / ratioOfSpecificHeats
           temperature = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
           fuel = 0.0_WP
           oxidizer = Yo0

           ! Offset by the virtual origin of the Blasius profile
           xx = x + x0

           ! Return the first derivative of the Blasius function
           delta = sqrt(xx / Re / U0)
           eta = y / delta
           if (eta .le. 0.0_WP) then
              blasius0 = 0.0_WP
              blasius1 = 0.0_WP
           else if (eta .ge. 9.0_WP) then
              blasius0 = by0(9) + (eta-9.0_WP)
              blasius1 = 1.0_WP
           else
              kk = int(eta)

              f0l = by0(kk)
              f0r = by0(kk + 1)
              f2l = by2(kk)
              f2r = by2(kk + 1)
              blasius0 = &
                   1.0_WP/6.0_WP*f2l*(real(kk+1, WP) - eta)**3 +                             &
                   1.0_WP/6.0_WP*f2r*(eta-real(kk, WP))**3 +                                 &
                   (f0l-1.0_WP/6.0_WP*f2l)*(real(kk+1, WP) - eta) +                          &
                   (f0r-1.0_WP/6.0_WP*f2r)*(eta-real(kk, WP))

              f0l = by1(kk)
              f0r = by1(kk+1)
              f2l = -0.5_WP * by0(kk) * by2(kk)
              f2r = -0.5_WP * by0(kk+1) * by2(kk+1)
              blasius1 =                                                                     &
                   1.0_WP/6.0_WP * f2l * (real(kk+1, WP)-eta)**3 +                           &
                   1.0_WP/6.0_WP * f2r * (eta-real(kk, WP))**3 +                             &
                   (f0l-1.0_WP/6.0_WP*f2l) * (real(kk+1, WP) - eta) +                        &
                   (f0r-1.0_WP/6.0_WP*f2r) * (eta-real(kk, WP))

           end if

           ! Update the velocity
           velocity(1) = U0 * blasius1
           velocity(2) = 0.5_WP * sqrt(U0 / xx / Re) * (eta * blasius1 - blasius0)

           ! Make sure wall velocities are zero
           if (j .eq. 1) velocity = 0.0_WP

           ! Jet conditions
           r = sqrt((x - jetPosition)**2 + z**2)
           insideJet = .false.
           select case(trim(jetShape))
           case('round')
              if (r .le. 1.0_WP * jetDiameter) insideJet = .true.
           case('square')
              if (i .ge. iJet1 .and. i .le. iJet2 .and. k .ge. kJet1 .and. k .le. kJet2)     &
                   insideJet = .true.
           case ('none')
              ! Nothing to do
           case default
              call die("Unknown jet shape '" // trim(jetShape) // "'")
           end select

           !if (insideJet) then
           yDecay = max(0.0_WP,                                                              &
                0.5_WP * (1.0_WP + tanh(10.0_WP * (0.5_WP - y / jetDiameter))))

           ! Fuel stream
           sig = 6.0_WP
           fuel = Yf0 * yDecay
           if (trim(jetShape) .eq. 'round') then
              fuel = fuel * 0.5_WP * (                                                       &
                   tanh(sig * (r + 0.5_WP * jetDiameter) / jetDiameter) -                    &
                   tanh(sig * (r - 0.5_WP * jetDiameter) / jetDiameter))
           else if (trim(jetShape) .eq. 'square') then
              fuel = fuel * 0.5_WP * (                                                       &
                   tanh(sig * (x - jetPosition + 0.5_WP * jetDiameter) / jetDiameter) -      &
                   tanh(sig * (x - jetPosition - 0.5_WP * jetDiameter) / jetDiameter))  *    &
                   0.5_WP * (                                                                &
                   tanh(sig * (z + 0.5_WP * jetDiameter) / jetDiameter) -                    &
                   tanh(sig * (z - 0.5_WP * jetDiameter) / jetDiameter))
           end if
           oxidizer = (1.0_WP - fuel) * Yo0

           ! Jet velocity
           sig = 6.0_WP
           velocity(2) = jetVelocity * yDecay
           select case(trim(velocityProfile))
           case('tanh')

              if (trim(jetShape) .eq. 'round') then
                 velocity(2) = velocity(2) * 0.5_WP * (                                      &
                      tanh(sig * (r + 0.5_WP * jetDiameter) / jetDiameter) -                 &
                      tanh(sig * (r - 0.5_WP * jetDiameter) / jetDiameter))
              else if (trim(jetShape) .eq. 'square') then
                 velocity(2) = velocity(2) * 0.5_WP * (                                      &
                      tanh(sig * (x - jetPosition + 0.5_WP * jetDiameter) / jetDiameter) -   &
                      tanh(sig * (x - jetPosition - 0.5_WP * jetDiameter) / jetDiameter))  * &
                      0.5_WP * (                                                             &
                      tanh(sig * (z + 0.5_WP * jetDiameter) / jetDiameter) -                 &
                      tanh(sig * (z - 0.5_WP * jetDiameter) / jetDiameter))
              end if

           case ('poiseuille')

              if (trim(jetShape) .eq. 'round') then
                 velocity(2) = max(0.0_WP,                                                   &
                      velocity(2) * 2.0_WP * (1.0_WP - (r / (0.5_WP * jetDiameter))**2))
              else if (trim(jetShape) .eq. 'SQUARE') then
                 velocity(2) = 2.0_WP * velocity(2) * max(0.0_WP, (1.0_WP -                  &
                      ((x - jetPosition) / (0.5_WP * jetDiameter))**2)) *                    &
                      max(0.0_WP, (1.0_WP - (z / (0.5_WP * jetDiameter))**2))
              end if

           case ('self similar')

              eta = (x - jetPosition)
              a = (sqrt(2.0_WP) - 1.0_WP) / (0.5_WP * jetDiameter)**2
              velocity(2) = velocity(2) * (1.0_WP + a * eta**2) ** (-2.0_WP)
              velocity(1) = 0.5_WP * jetVelocity * (eta - a * eta**3) /                      &
                   (1.0_WP + a * eta**2)**2 * yDecay

           case default

           end select

           !end if !... if (insideJet)

           ! Correct species mass fractions.
           inert = 1.0_WP - fuel - oxidizer
           if (inert .lt. 0.0_WP) then
              if (iRank .eq. iRoot)                                                          &
                   print *, 'Something is wrong with the species mass fraction!'
              oxidizer = oxidizer + inert
           end if

           ! Get density from the equation of state
           select case (equationOfState)
           case(IDEAL_GAS)
              density = ratioOfSpecificHeats * pressure /                                    &
                   (temperature * (ratioOfSpecificHeats - 1.0_WP))
           case (IDEAL_GAS_MIXTURE)
              density = ratioOfSpecificHeats * pressure /                                    &
                   ( temperature * (ratioOfSpecificHeats - 1.0_WP) *                         &
                   (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )
           end select

           ! Set the conserved variables
           conservedVariables(grid_index(i,j,k), 1) = density
           conservedVariables(grid_index(i,j,k), 2:nDimensions+1) =                          &
                density * velocity(1:nDimensions)
           conservedVariables(grid_index(i,j,k),nDimensions+2) = pressure /                  &
                (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP * density * sum(velocity ** 2)
           if (nSpecies .gt. 0) conservedVariables(grid_index(i,j,k), nDimensions+2+H2) =    &
                density * fuel
           if (nSpecies .gt. 1) conservedVariables(grid_index(i,j,k), nDimensions+2+O2) =    &
                density * oxidizer

        end do
     end do
  end do

  return
end subroutine jet_crossflow_data
