module channel

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

end module channel


! =============== !
! Create the grid !
! =============== !
subroutine channel_grid

  ! Internal modules
  use channel

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: ytilde, stretchParameter

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

  ! Should we stretch the mesh?
  call parser_read('stretching', stretchParameter, 2.0_WP)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) / real(nx - 1, WP)

           ! Create Y
           ytilde = 2.0_WP * real(j-1,WP) / real(ny-1,WP) - 1.0_WP
           coordinates(grid_index(i,j,k), 2) = 0.5_WP * Ly *                                 &
                tanh(stretchParameter * ytilde) / tanh(stretchParameter)

           ! Create Z
           if (nz .gt. 1) then
              coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /             &
                   real(nz - 1, WP) - 0.5_WP * Lz
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine channel_grid


! ========================= !
! Create the variable array !
! ========================= !
subroutine channel_data

  ! Internal modules
  use channel

  ! External modules
  use math
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, u0, p0, x, y, z, y1, y2, amp, u, v, w, coeff, rnd
  logical :: laminar

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Initialize the random number generator
  call random_setup

  ! Read the input
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('mean velocity', u0)
  call parser_read('fluctuation amplitude', amp, 0.0_WP)
  call parser_read('laminar profile', laminar, .true.)

  ! Set constant pressure
  p0 = 1.0_WP / gamma

  ! Get min/max y coordinates
  i = iStart(1); j = iStart(2)+1; k = iStart(3)
  y = coordinates(grid_index(i,j,k), 2)
  call parallel_min(y, y1)
  j = iEnd(2)-1
  y = coordinates(grid_index(i,j,k), 2)
  call parallel_max(y, y2)

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Local coordinates
           x  = coordinates(grid_index(i,j,k), 1)
           y  = coordinates(grid_index(i,j,k), 2)
           if (nz.gt.1) z = coordinates(grid_index(i,j,k), 3)

           ! Laminar profile
           coeff = 6.0_WP
           u = u0 * coeff * (y-y1) * (y2-y) / (y2-y1)**2
           w = 0.0_WP

           ! For faster transition
           if (.not.laminar .and. nz.gt.1) then
              ! Fluctuations in X for W
              w = w + amp * abs(u0) * cos(8.0_WP * twoPi * x / Lx)
              ! Fluctuations in Z for U
              u = u + amp * abs(u0) * cos(8.0_WP * twoPi * z / Lz)
           end if

           ! Random values
           call random_number(rnd)
           u = u + (rnd - 0.5_WP) * amp * u
           call random_number(rnd)
           v = v + (rnd - 0.5_WP) * amp * v
           call random_number(rnd)
           w = w + (rnd - 0.5_WP) * amp * w

           ! Zero-out the walls
           if (j.eq.1 .or. j.eq.ny) then
              u = 0.0_WP
              v = 0.0_WP
              w = 0.0_WP
           end if

           ! Set the density
           conservedVariables(grid_index(i,j,k), 1) = 1.0_WP

           ! Set the velocity
           conservedVariables(grid_index(i,j,k), 2) = u
           conservedVariables(grid_index(i,j,k), 3) = v
           if (nz.gt.1) conservedVariables(grid_index(i,j,k), 4) = w

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
end subroutine channel_data
