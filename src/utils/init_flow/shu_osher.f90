module shu_osher

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

end module shu_osher


subroutine shu_osher_grid

  ! Internal modules
  use shu_osher

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
  call parser_read('Lx', Lx)
  call parser_read('Ly', Ly)
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

end subroutine shu_osher_grid


subroutine shu_osher_data

  ! Internal modules
  use shu_osher

  ! External modules
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, dir
  real(WP) :: gamma, x, x0, p1, rho1, u1, density, pressure, velocity

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read in shock parameters
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('shock position', x0, -4.0_WP)
  call parser_read('shock direction', dir, 1)
  call parser_read('post-shock pressure', p1, 31.0_WP/3.0_WP)
  call parser_read('post-shock density', rho1, 27.0_WP/7.0_WP)
  call parser_read('post-shock velocity', u1, 4.0_WP*sqrt(35.0_WP)/9.0_WP)

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Get local coordinate
           x = coordinates(grid_index(i,j,k), dir)

           ! Set the density & pressure
           if (x.le.x0) then
              density  = rho1
              pressure = p1
              velocity = u1
           else
              density  = 1.0_WP + 0.2_WP * sin(5.0_WP * x)
              pressure = 1.0_WP
              velocity = 0.0_WP
           end if

           ! Set the conserved variables
           conservedVariables(grid_index(i,j,k), 1) = density
           conservedVariables(grid_index(i,j,k), 2:nDimensions+1) = 0.0_WP
           conservedVariables(grid_index(i,j,k), nDimensions+2) = pressure /                 &
                (gamma - 1.0_WP) + 0.5_WP * density * velocity ** 2

           conservedVariables(grid_index(i,j,k), 1+dir) = density * velocity
        end do
     end do
  end do

  return
end subroutine shu_osher_data
