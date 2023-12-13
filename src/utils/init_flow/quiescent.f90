module quiescent

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

end module quiescent

subroutine quiescent_grid

  ! Internal modules
  use quiescent

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

  return
end subroutine quiescent_grid

subroutine quiescent_data

  ! Internal modules
  use quiescent

  ! External modules
  use state, only : conservedVariables

  implicit none

  ! Local variables
  real(WP) :: gamma, P0, rho0

  ! Get the number of species
  call parser_read('number of species', nSpecies, 0)

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2 + nSpecies
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Ambient pressure
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  P0 = 1.0_WP / gamma

  ! Constant density
  rho0 = 1.0_WP

  ! Assign the state variables
  conservedVariables = 0.0_WP 
  conservedVariables(:, 1) = rho0 ! ... density
  conservedVariables(:, 2:nDimensions+1) = 0.0_WP ! ... momentum
  conservedVariables(:, nDimensions+2) = P0 / (gamma - 1.0_WP) ! ... energy

  return
end subroutine quiescent_data
