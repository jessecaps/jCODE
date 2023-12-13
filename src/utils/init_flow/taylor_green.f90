module taylor_green

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
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, v0

end module taylor_green

subroutine taylor_green_grid

  ! Internal modules
  use taylor_green

  ! External modules
  use parallel
  use math, only : pi

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
  dx = Lx / real(nx, WP)
  dy = Ly / real(ny, WP)
  dz = Lz / real(nz, WP)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           if (nx .gt. 1) then
                 coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /          &
                      real(nx - 1, WP) - pi
           end if

           ! Create X
           if (ny .gt. 1) then
                 coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) /          &
                      real(ny - 1, WP) - pi
           end if


           ! Create Z
           if (nz .gt. 1) then
                 coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /          &
                      real(nz - 1, WP) - pi
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine taylor_green_grid

subroutine taylor_green_data

  ! Internal modules
  use taylor_green

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: gamma, pressure, velocity(3), T0, P0, rho0, x, y, z

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Zero-out the conserved variables
  conservedVariables = 0.0_WP

  ! Reference quantities
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  T0 = 1.0_WP / (gamma - 1.0_WP)
  P0 = 1.0_WP / gamma
  rho0 = 1.0_WP
  velocity = 0.0_WP

  ! Read in the shear layer parameters
  call parser_read('initial velocity', v0)
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Ly', Ly, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Assign the initial conditions
  do i = 1, nGridPoints

     ! Local cordinates
     x = coordinates(i, 1)
     y = coordinates(i, 2)
     if (nz.gt.1) then
        z = coordinates(i, 3)
     else
        z = 0.5_WP * pi
     end if

     ! Assign the velocity
     velocity(1) =  v0 * sin(x) * cos(y) * sin(z)
     velocity(2) = -v0 * cos(x) * sin(y) * sin(z)

     ! Assign the pressure
     pressure = P0 + rho0 * v0**2 / 16.0_WP * (cos(2.0_WP * x) + cos(2.0_WP * y)) *          &
          (cos(2.0_WP * z) + 2.0_WP)

     ! Set the density
     conservedVariables(i,1) = pressure / T0 / (gamma - 1.0_WP) * gamma

     ! Set the momentum
     conservedVariables(i,2:nDimensions+1) = conservedVariables(i,1) * velocity(1:nDimensions)

     ! Set the energy
     conservedVariables(i, nDimensions+2) = pressure / (gamma - 1.0_WP) +                    &
          0.5_WP * conservedVariables(i,1) * sum(velocity(1:nDimensions))**2

  end do

end subroutine taylor_green_data
