module mixing_layer

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

end module mixing_layer

subroutine mixing_layer_grid

  ! Internal modules
  use mixing_layer

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

end subroutine mixing_layer_grid

subroutine mixing_layer_data

  ! Internal modules
  use mixing_layer

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, k, H2, O2, N2
  real(WP) :: density, temperature, pressure, velocity, upperVelocity, lowerVelocity,        &
       lowerTemperature, upperTemperature, T0, P0, x, y, x0, fuel, oxidizer, inert,          &
       mixtureFraction, Z0, YF0, Yo0, growthRate, minDensity, maxDensity
  real(WP), dimension(:), allocatable :: Wi

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
              call die("mixing_layer_data: unknown species: '" //                             &
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

  ! Zero-out the conserved variables
  conservedVariables = 0.0_WP

  ! Reference temperature and pressure
  T0 = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  P0 = 1.0_WP / ratioOfSpecificHeats

  ! Read in the shear layer parameters
  call parser_read('mixing layer upper velocity', upperVelocity)
  call parser_read('mixing layer lower velocity', lowerVelocity)
  call parser_read('mixing layer upper temperature', upperTemperature, T0)
  call parser_read('mixing layer lower temperature', lowerTemperature, T0)
  call parser_read('mixing layer pressure', pressure, P0)
  call parser_read('mixing layer growth rate', growthRate, 0.1_WP)
  call parser_read('mixing layer growth origin', x0, Lx)

  ! Mixture properties
  if (nSpecies .gt. 0) then
     call parser_read('initial mixture fraction', Z0, 1.0_WP)
     call parser_read('initial fuel mass fraction', Yf0)
     call parser_read('initial oxidizer mass fraction', Yo0)
  end if

  ! Mean profiles
  minDensity = huge(1.0_WP); maxDensity = -huge(1.0_WP)
  do i = 1, nGridPoints

     ! Coordinate
     x = coordinates(i, 1)
     y = coordinates(i, 2)

     ! Velocity profile
     velocity = lowerVelocity + 0.5_WP * (upperVelocity - lowerVelocity) *                   &
          (1.0_WP + tanh(0.5_WP * y / (1.0_WP + growthRate * max(0.0_WP, x - x0))))

     ! Temperature profile
     temperature = lowerTemperature + 0.5_WP * (upperTemperature - lowerTemperature) *       &
          (1.0_WP + tanh(2.0_WP * y / (1.0_WP + growthRate * max(0.0_WP, x - x0))))

     if (nSpecies .gt. 0) then
        ! Mixture fraction
        mixtureFraction = 0.5_WP * Z0 * ( 1.0_WP - erf(y / (1.0_WP + growthRate *            &
             max(0.0_WP, x - x0))))

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
        density = density / (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2))
     else
        density = 1.0_WP
     end if

     ! Keep track of min/max density
     minDensity = min(minDensity, density)
     maxDensity = max(maxDensity, density)

     ! Assign the density
     conservedVariables(i,1) = density

     ! Assign the momentum
     conservedVariables(i,2) = density * velocity

     ! Assign the energy
     conservedVariables(i, nDimensions+2) = pressure / (ratioOfSpecificHeats - 1.0_WP) +     &
          0.5_WP * conservedVariables(i,2)**2 / conservedVariables(i,1)
     if (nDimensions .gt. 1) conservedVariables(i, nDimensions+2) =                          &
          conservedVariables(i, nDimensions+2) + 0.5_WP * conservedVariables(i,3)**2 /       &
          conservedVariables(i,1)
     if (nDimensions .eq. 3) conservedVariables(i, nDimensions+2) =                          &
          conservedVariables(i, nDimensions+2) + 0.5_WP * conservedVariables(i,4)**2 /       &
          conservedVariables(i,1)

     ! Assign the species
     if (nSpecies .gt. 0) conservedVariables(i,nDimensions+2+H2) =                           &
          fuel * conservedVariables(i,1)
     if (nSpecies .gt. 1) conservedVariables(i,nDimensions+2+O2) =                           &
          oxidizer * conservedVariables(i,1)
  end do

  ! Output min/max density
  call parallel_min(minDensity)
  call parallel_max(maxDensity)
  if (iRank .eq. iRoot) then
     print *
     print *, 'min/max density:', minDensity, maxDensity
     print *
  end if

end subroutine mixing_layer_data
