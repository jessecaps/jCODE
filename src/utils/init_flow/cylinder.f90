module cylinder

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
  real(WP) :: Lx, Ly, Lz, dz

end module cylinder


! ============================================= !
! Generate the grid: x=> radius, y=>theta, z=>z !
! ============================================= !
subroutine cylinder_grid

  ! Internal modules
  use cylinder

  ! External modules
  use math, only : pi

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: theta, r, rpole, s, ytilde
  logical :: stretchR

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)
  Ly = 2.0_WP * pi

  ! Compute the grid spacing
  dz = Lz / real(nz, WP)

  ! Get the cylinder radius
  call parser_read('cylinder radius', rpole)

  ! Grid Stretching
  call parser_read('stretch radius', stretchR, .false.)
  if (stretchR) s = 2.0_WP

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        theta = real(j-1, WP) / real(ny-1, WP) * Ly
        do i = iStart(1), iEnd(1)

           ! Radius
           if (stretchR) then
              ytilde = real(nx-i, WP) / real(nx-1, WP)
              r = rpole + (0.5_WP * Lx - rpole) * (1.0_WP - tanh(s * ytilde) / tanh(s))
           else
              r = real(i-1, WP) / real(nx-1, WP) * (0.5_WP * Lx - rpole) + rpole
           end if

           ! Lumpy cylinder
           !r = r + 0.02_WP * sin(real(j-1,WP) / real(ny-1,WP) * 2.0_WP * pi * 8.0_WP)

           if (j.eq.1 .or. j.eq.ny) then

              ! Create X
              coordinates(grid_index(i,j,k), 1) = r

              ! Create Y
              coordinates(grid_index(i,j,k), 2) = 0.0_WP

           else

              ! Create X
              coordinates(grid_index(i,j,k), 1) = r * cos(theta)

              ! Create Y
              coordinates(grid_index(i,j,k), 2) = r * sin(theta)

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
end subroutine cylinder_grid


! =============================== !
! Generate the initial conditions !
! =============================== !
subroutine cylinder_data

  ! Internal modules
  use cylinder

  ! External modules
  use math
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: gamma, P0, rho0, u0, v, sigma, freq, amp, r
  logical ::  perturbFlow

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Read from input
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('crossflow mach number', u0, 0.0_WP)
  call parser_read('cylinder perturbations', perturbFlow, .false.)
  if (perturbFlow) then
     call parser_read('cylinder perturbation amplitude', amp, 0.1_WP*u0)
     call parser_read('cylinder perturbation frequency', freq, 10.0_WP)
     call parser_read('cylinder perturbation sigma', sigma, 1.0_WP)
  end if

  ! Ambient pressure
  P0 = 1.0_WP / gamma

  ! Constant density
  rho0 = 1.0_WP

  ! Assign the state variables
  conservedVariables(:, 1) = rho0 ! ... density
  conservedVariables(:, 2) = rho0 * u0 ! ... streamwise momentum
  conservedVariables(:, 3:nDimensions+1) = 0.0_WP ! ... momentum

  ! Add perturbations to vertical velocity
  if (perturbFlow) then
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Avoid boundary points
              if (i.eq.1 .or. i.eq.globalGridSize(1)) cycle

              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))
              r = sqrt(sum(coordinates(gridIndex,1:2)**2))
              v = amp * sin(r * freq / Lx * 4.0_WP * pi)*exp(-r**2 / (2.0_WP * sigma**2))
              conservedVariables(gridIndex, 3) = rho0 * v
           end do
        end do
     end do
  end if

  conservedVariables(:, nDimensions+2) = p0 / (gamma - 1.0_wp) + 0.5_WP *                    &
       sum(conservedVariables(:, 2:nDimensions+1)**2, 2) / conservedVariables(:, 1) ! ... energy

  return
end subroutine cylinder_data
