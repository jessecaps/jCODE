module premixed

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

end module premixed

subroutine premixed_grid

  ! Internal modules
  use premixed

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

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           if (nx .gt. 1) then
              if (periodicityType(1) .eq. PLANE) then
                 coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /          &
                      real(nx - 1, WP) - 0.5_WP * Lx
              else
                 coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /                  &
                      real(nx - 1, WP) - 0.5_WP * Lx
              end if
           end if

           ! Create y
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

  return
end subroutine premixed_grid

subroutine premixed_data

  ! Internal modules
  use premixed

  ! External modules
  use parser
  use parallel
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: k, H2, O2, N2
  real(WP) :: density, T0, Tref, P0, Z0, Yf0, Yo0, fuel, oxidizer, inert
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

  ! Zero-out the conserved variables.
  conservedVariables = 0.0_WP

  ! Reference temperature and pressure
  Tref = 1.0_WP / (ratioOfSpecificHeats - 1.0_WP)
  call parser_read('initial temperature', T0, Tref)
  P0 = 1.0_WP / ratioOfSpecificHeats

  ! Read in the mixture parameters
  if (nSpecies .gt. 0) then
     call parser_read('initial mixture fraction', Z0, 1.0_WP)
     call parser_read('initial fuel mass fraction', YF0)
     call parser_read('initial oxidizer mass fraction', YO0)
  end if

  ! Components
  fuel = YF0 * Z0
  oxidizer = YO0 * (1.0_WP - Z0)

  ! Correct species mass fractions
  inert = 1.0_WP - fuel - oxidizer
  if (inert .lt. 0.0_WP) oxidizer = oxidizer + inert

  ! Get density from the equation of state
  select case (equationOfState)
  case(IDEAL_GAS)
     density = ratioOfSpecificHeats * P0 / (T0 * (ratioOfSpecificHeats - 1.0_WP))

  case (IDEAL_GAS_MIXTURE)
     density = ratioOfSpecificHeats * P0 / ( T0 * (ratioOfSpecificHeats - 1.0_WP) *          &
          (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )

  end select
  if (iRank .eq. iRoot) then
     print *
     print *, 'Mixture density = ', density
     print *
  end if

  ! Assign the state variables
  conservedVariables = 0.0_WP
  conservedVariables(:, 1) = density
  conservedVariables(:, 2:nDimensions+1) = 0.0_WP
  conservedVariables(:, nDimensions+2) = P0 / (ratioOfSpecificHeats - 1.0_WP)
  if (nSpecies .gt. 0) conservedVariables(:, nDimensions+2+H2) =                             &
       fuel * conservedVariables(:, 1)
  if (nSpecies .gt. 1) conservedVariables(:, nDimensions+2+O2) =                             &
       oxidizer * conservedVariables(:, 1)

end subroutine premixed_data
