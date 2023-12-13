module boundary_layer

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

end module boundary_layer

subroutine boundary_layer_grid

  ! Internal modules
  use boundary_layer

  ! External modules
  use parallel

  implicit none

  ! Boundary layer schematic and the corresponding indeces 
  !
  !
  !                      _~~~~~~~~~~~~~~~~~~~~~~~~~~_
  ! j = 1  > ___________|                            |____________________________________
  !
  !          ^          ^                            ^                                   ^
  !         i=1        i=iTrip1                    i=iTrip2                            i=nx

  ! Local variables
  integer :: i, j, k, n, nGrit, iTrip1, iTrip2
  real(WP) :: ytilde, r, tripLocation, tripWidth, tripHeight, gritHeight, gritWidth,         &
       sig, x, y, z, x0, z0, z12, amp, gauss, alpha, rand,                                   &
       totalHeight, peakHeight, valleyHeight
  real(WP), allocatable :: delta(:)
  logical :: stretchGrid, includeSandpaper, includeGrit

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

  ! Generate the sandpaper
  call parser_read('include sandpaper', includeSandpaper, .false.)
  if (includeSandpaper) then
     if (nProcsDir(2) .gt. 1)                                                                &
          call die('Only 1 proc can be assigned in y when using sandpaper')
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
     allocate(delta(ny))
     i = iStart(1); k=iStart(3)
     do j = 2, ny
        delta(j) = coordinates(grid_index(i,j,k), 2) - coordinates(grid_index(i,j-1,k), 2)
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
              ytilde = coordinates(grid_index(i,j-1,k),2) + delta(j)
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

                 ! Get the coordinates.
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

                 ! Update the vertical coordinate.
                 coordinates(grid_index(i,j,k), 2) = y + gauss

                 ! Shift grid points above (smoothly)
                 do j = 2, ny

                    ! Get the current height
                    y = coordinates(grid_index(i,j,k), 2)

                    ! Adjust the current height.
                    ytilde = coordinates(grid_index(i,j-1,k),2) + delta(j)
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

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine boundary_layer_grid

subroutine boundary_layer_data

  ! Internal modules
  use boundary_layer

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, kk
  real(WP) :: gamma, x0, u0, p0, x, y, xx, Re, delta, eta

  ! Solution to Blasius boundary layer
  real(wp) :: blasius0, blasius1
  real(wp) :: f2l, f2r, f0l, f0r
  real(wp), dimension(0:9) :: by0 = (/                                                       &
       0.0_8, 0.165571818583440_8, 0.650024518764203_8, 1.39680822972500_8,                  &
       2.30574664618049_8, 3.28327391871370_8, 4.27962110517696_8,                           &
       5.27923901129384_8, 6.27921363832835_8, 7.27921257797747_8 /)
  real(wp), dimension(0:9) :: by1 = (/                                                       &
       0.0_8, 0.329780063306651_8, 0.629765721178679_8, 0.84604458266019_8,                  &
       0.95551827831671_8, 0.99154183259084_8, 0.99897290050990_8,                           &
       0.9999216098795_8, 0.99999627301467_8, 0.99999989265063_8 /)
  real(wp), dimension(0:9) :: by2 = (/                                                       &
       0.332057384255589_8, 0.323007152241930_8, 0.266751564401387_8, 0.161360240845588_8,   &
       0.06423404047594_8, 0.01590689966410_8, 0.00240199722109_8,                           &
       0.00022016340923_8, 0.00001224984692_8, 0.00000041090325_8 /)

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read the input
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('free stream velocity', u0)
  call parser_read('blasius virtual origin', x0)
  call parser_read('Reynolds number', Re)

  ! Set constant pressure
  p0 = 1.0_WP / gamma

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Local coordinates
           x  = coordinates(grid_index(i,j,k), 1)
           y  = coordinates(grid_index(i,j,k), 2)

           ! Offset by the virtual origin of the Blasius profile
           xx = x + x0

           ! Return the first derivative of the Blasius function
           delta = sqrt(xx / Re / u0)
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

           ! Set the density
           conservedVariables(grid_index(i,j,k), 1) = 1.0_WP

           ! Set the velocity
           conservedVariables(grid_index(i,j,k), 2) = u0 * blasius1
           conservedVariables(grid_index(i,j,k), 3) = 0.5_WP * sqrt(u0 / xx / Re) *          &
                (eta * blasius1 - blasius0)
           conservedVariables(grid_index(i,j,k), nDimensions+2) = 0.0_WP

           ! Make sure wall velocities are zero
           if (j .eq. 1) conservedVariables(grid_index(i,j,k), 2:nDimensions+1) = 0.0_WP

           ! Multiply by density
           conservedVariables(grid_index(i,j,k), 2:nDimensions+1) =                          &
                conservedVariables(grid_index(i,j,k), 2:nDimensions+1) *                     &
                conservedVariables(grid_index(i,j,k), 1)

           ! Set the energy
           conservedVariables(grid_index(i,j,k), nDimensions+2) =                            &
                p0 / (gamma - 1.0_wp) + 0.5_WP *                                             &
                sum(conservedVariables(grid_index(i,j,k), 2:nDimensions+1)**2) /             &
                conservedVariables(grid_index(i,j,k), 1)

        end do
     end do
  end do

  return
end subroutine boundary_layer_data
