module sphere

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
  real(WP) :: L

end module sphere


! ============================================= !
! Generate the grid: x=> radius, y=>theta, z=>z !
! ============================================= !
subroutine sphere_grid

  ! Internal modules
  use sphere

  ! External modules
  use math, only : pi

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: theta, phi, phi0, r, rpole, ytilde, s
  logical :: stretchGrid

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  if (nz.eq.1) call die('sphere_grid: spherical geometry must be 3D!')

  ! Get the sphere radius
  call parser_read('domain length', L)
  call parser_read('sphere radius', rpole)

  ! Grid Stretching
  call parser_read('use stretching', stretchGrid, .false.)
  if (stretchGrid) s = 2.0_WP

  ! Avoid singularity of polar angle
  phi0 = 0.5_WP*pi/real(nz,WP)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     ! Phi (avoid pole singularity!)
     phi = real(k-1, WP) / real(nz-1, WP) * (pi-2.0_WP*phi0) + phi0
     do j = iStart(2), iEnd(2)
        ! Theta
        theta = real(j-1, WP) / real(ny-1, WP) * 2.0_WP * pi
        do i = iStart(1), iEnd(1)

           ! Radius
           if (stretchGrid) then
              ytilde = real(nx-i, WP) / real(nx-1, WP)
              r = rpole + (0.5_WP*L - rpole) * (1.0_WP - tanh(s * ytilde) / tanh(s))
           else
              r = real(i-1, WP) / real(nx-1, WP) * (0.5_WP * L - rpole) + rpole
           end if

           if (j .eq. 1 .or. j .eq. ny) then

              ! Create X
              coordinates(grid_index(i,j,k), 1) = r * sin(phi)

              ! Create Z
              coordinates(grid_index(i,j,k), 3) = 0.0_WP

           else

              ! Create X
              coordinates(grid_index(i,j,k), 1) = r * cos(theta) * sin(phi)

              ! Create Z
              coordinates(grid_index(i,j,k), 3) = r * sin(theta) * sin(phi)

           end if

           ! Create Y
           coordinates(grid_index(i,j,k), 2) = r * cos(phi)

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

  return
end subroutine sphere_grid


! =============================== !
! Generate the initial conditions !
! =============================== !
subroutine sphere_data

  ! Internal modules
  use sphere

  ! External modules
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, P0, rho0, u0, pmax, pos(3), r0, p, buf
  logical :: usePulse

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Ambient conditions
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  P0 = 1.0_WP / gamma
  rho0 = 1.0_WP

  ! Read from input
  call parser_read('crossflow mach number', u0, 0.0_WP)
  call parser_read('use pressure pulse', usePulse, .false.)
  if (usePulse) then
     call parser_read('peak pressure', pmax)
     call parser_read('pressure x', pos(1))
     call parser_read('pressure y', pos(2))
     call parser_read('pressure z', pos(3))
     call parser_read('pressure r', r0)
  end if

  ! Assign the state variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Specify conditions
           p = p0
           if (usePulse) then
              buf = sum((coordinates(grid_index(i,j,k),1:nDimensions) -                      &
                   pos(1:nDimensions))**2)
              buf = exp(-0.5_WP * 9.0_WP * buf / r0**2)
              p = p + pmax * buf
           end if
           
           conservedVariables(grid_index(i,j,k), 1) = rho0 ! ... density
           conservedVariables(grid_index(i,j,k), 2) = rho0 * u0 ! ... streamwise momentum
           conservedVariables(grid_index(i,j,k), 3:nDimensions+1) = 0.0_WP ! ... momentum
           conservedVariables(grid_index(i,j,k), nDimensions+2) = p / (gamma - 1.0_WP) +     &
                0.5_WP * sum(conservedVariables(grid_index(i,j,k), 2:nDimensions+1)**2) /    &
                conservedVariables(grid_index(i,j,k), 1) ! ... energy
        end do
     end do
  end do

  return
end subroutine sphere_data
