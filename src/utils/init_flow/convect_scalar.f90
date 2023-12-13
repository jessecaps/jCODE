module convect_scalar

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

end module convect_scalar

subroutine convect_scalar_grid

  ! Internal modules
  use convect_scalar

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
end subroutine convect_scalar_grid

subroutine convect_scalar_data

  ! Internal modules
  use convect_scalar

  ! External modules
  use parser
  use parallel
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, dim, direction, gridIndex
  real(WP) :: gamma, density, T0,P0, tau, Uc(3), xyz(3)
  character(len=str_short) :: orientation

  ! Set the number of conserved variables and allocate the array
  nSpecies = 1
  nUnknowns = nDimensions + 2 + nSpecies
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Zero-out the conserved variables.
  conservedVariables = 0.0_WP

  ! Reference quantities
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  density = 1.0_WP
  T0 = 1.0_WP / (gamma - 1.0_WP)
  P0 = 1.0_WP / gamma

  ! Set the convective velocity
  call parser_getsize('convective velocity', dim)
  if (dim .ne. 3) stop 'Convective velocity should be of size 3'
  call parser_read('convective velocity', Uc)

  ! Assign the state variables
  conservedVariables = 0.0_WP
  do i = 1, nGridPoints
     conservedVariables(i, 1) = density
     conservedVariables(i, 2:nDimensions+1) = Uc(1:nDimensions)
     conservedVariables(i, nDimensions+2) = P0 / (gamma - 1.0_WP) + 0.5_WP * density *       &
          sum(Uc ** 2)
  end do
  
  ! Create the scalar
  call parser_read('scalar size', tau)
  call parser_read('scalar orientation', orientation)
  select case (trim(orientation))
  case ('x')
     direction = 1
  case ('y')
     direction = 2
  case ('z')
     direction = 3
  end select
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           gridIndex = grid_index(i,j,k)
           xyz(1:nDimensions) = coordinates(gridIndex, 1:nDimensions)
           xyz(direction) = 0.0_WP
           conservedVariables(gridIndex, nDimensions+3) = density * exp(-sum(xyz**2)/tau**2)
        end do
     end do
  end do

  return
end subroutine convect_scalar_data
