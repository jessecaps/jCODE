module pipe

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
  real(WP) :: Lx, Ly, dx, dy

end module pipe

subroutine pipe_grid

  ! Internal modules
  use pipe

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: rtilde, s
  logical :: stretchGrid

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Only implemented in 2D (for now)
  if (nz.gt.1) call die('pipe flow must be 2D!')

  ! Read in the grid size
  call parser_read('Lx', Lx)
  call parser_read('Ly', Ly)

  ! Compute the grid spacing
  dx = Lx / real(nx-1, WP)
  dy = Ly / real(ny, WP)

  ! Should we stretch the mesh?
  call parser_read('grid stretching', stretchGrid, .false.)
  if (stretchGrid) s = 2.0_WP

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create R
           if (stretchGrid) then
              rtilde = real(i-1, WP) / real(nx-1, WP)
              coordinates(grid_index(i,j,k), 1) = Lx * tanh(s * rtilde) / tanh(s)
           else
              coordinates(grid_index(i,j,k), 1) = Lx * real(i-1, WP) / real(nx-1, WP)
           end if

           ! Create L
           coordinates(grid_index(i,j,k), 2) = Ly *  real(j - 1, WP) / real(ny - 1, WP)

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine pipe_grid

subroutine pipe_data

  ! Internal modules
  use pipe

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, u0, p0, r, dpdx(3), Re
  logical :: laminar

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read the input
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('pipe velocity', u0)
  call parser_read('Reynolds number', Re, 0.0_WP)
  call parser_read('poiseuille profile', laminar, .false.)
  dpdx = 0.0_WP
  if (laminar) call parser_read('pressure gradient', dpdx)

  ! Set constant pressure
  p0 = 1.0_WP / gamma

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Set the density
           conservedVariables(grid_index(i,j,k), 1) = 1.0_WP

           ! Set the velocity
           conservedVariables(grid_index(i,j,k), 2:nDimensions+1) = 0.0_WP
           if (laminar) then
              ! Initialize a Piseuille profile
              r = abs(coordinates(grid_index(i,j,k), 1))
              conservedVariables(grid_index(i,j,k), 2:nDimensions+1) = 0.25_WP * Re *        &
                   dpdx(1:nDimensions) * (Lx**2 - r**2)
           else
              ! Uniform velocity
              conservedVariables(grid_index(i,j,k), 3) = u0
           end if

           ! Make sure wall velocities are zero
           if (i .eq. nx) conservedVariables(grid_index(i,j,k), 2:nDimensions+1) = 0.0_WP

           ! Set the energy
           conservedVariables(grid_index(i,j,k), nDimensions+2) =                            &
                p0 / (gamma - 1.0_wp) + 0.5_WP *                                             &
                sum(conservedVariables(grid_index(i,j,k), 2:nDimensions+1)**2) /             &
                conservedVariables(grid_index(i,j,k), 1)

        end do
     end do
  end do

  return
end subroutine pipe_data
