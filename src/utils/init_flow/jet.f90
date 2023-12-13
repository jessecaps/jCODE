module jet

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
  real(WP) :: Lx, Ly, Lz

end module jet

subroutine jet_grid

  ! Internal modules
  use jet

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k, nxs
  real(WP) :: sigma, b, c, Rj, dxu, xStretchRate
  real(WP), allocatable, dimension(:) :: s, g
  logical :: stretchGridX, stretchGridR

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Ly', Ly, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Should we stretch the mesh?
  call parser_read('axial stretching', stretchGridX, .false.)
  call parser_read('radial stretching', stretchGridR, .false.)

  ! Read in the jet radius
  call parser_read('jet radius', Rj, 0.5_WP)
  
  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           if (.not. stretchGridX) then
              coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) / real(nx - 1, WP)
           end if

           ! Create Y
           if (ny .gt. 1 .and. .not. stretchGridR) then
              coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /                     &
                   real(ny - 1, WP) - 0.5_WP * Ly
           end if

           ! Create Z
           if (nz .gt. 1 .and. .not. stretchGridR) then
              coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /                     &
                   real(nz - 1, WP) - 0.5_WP * Lz
           end if

        end do
     end do
  end do

  ! Stretch in x
  !--------------------------------------------------
  if (stretchGridX) then

     call parser_read('x stretching grid', nxs)
     call parser_read('x stretching base', dxu)
     call parser_read('x stretching rate', xStretchRate)

     allocate(g(nx))
     do i = 1, nx-nxs
        ! Uniform portion
        g(i) = real(i - 1, WP) * dxu
     end do
     do i = 1+nx-nxs, nx
        ! Stretched portion
        g(i) = g(i-1) + dxu * xStretchRate ** (i - (1+nx-nxs))
     end do
     ! Hard code the end boundary grid
     g(nx) = Lx

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1),iEnd(1)
              coordinates(grid_index(i,j,k), 1) = g(i)
           end do
        end do
     end do
     deallocate(g)
  end if
  
  if (stretchGridR) then
     ! Grid stretching parameters
     sigma = 0.21_WP
     b = 24.0_WP
     c = 0.6_WP

     ! Stretch in y
     !--------------------------------------------------
     if (ny.gt.1) then

        ! Create uniform spacing
        allocate(s(ny))
        do j = 1, ny
           s(j) = real(j - 1, WP) / real(ny - 1, WP)
        end do

        ! Compute mapping g(s)
        allocate(g(ny))
        call mapping_function(s, b, c, sigma, g)

        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              do i = iStart(1), iEnd(1)
                 ! Create y
                 coordinates(grid_index(i,j,k), 2) = 0.5_WP * Ly * (1.0_WP + g(j))        &
                      - 0.5_WP * Ly
              end do
           end do
        end do
        deallocate(s)
        deallocate(g)
     end if

     ! Stretch in z
     !--------------------------------------------------
     if (nz.gt.1) then
        ! Create uniform spacing
        allocate(s(nz))
        do k = 1, nz
           s(k) = real(k - 1, WP) / real(nz - 1, WP)
        end do

        ! Compute mapping g(s)
        allocate(g(nz))
        call mapping_function(s, b, c, sigma, g)

        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              do i = iStart(1), iEnd(1)
                 ! Create z
                 coordinates(grid_index(i,j,k), 3) = 0.5_WP * Lz * (1.0_WP + g(k)) - 0.5_WP * Lz
              end do
           end do
        end do
        deallocate(s)
        deallocate(g)
     end if
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

end subroutine jet_grid

subroutine jet_data

  ! Internal modules
  use jet

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: gamma, density, pressure, velocity(3), mixtureFraction, uj, uc, P0, rho0,      &
       x, y, z, r, Rj, RjI, theta0, slope, growthRate, x0

  ! Initialize the random number generator
  ! call random_init
  
  ! Get the number of species
  call parser_read('number of species', nSpecies, 0)
  if (nSpecies .gt. 1)call die("jet_data: nSpecies <= 1 for now!")
  
  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2 + nSpecies
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Setup solver_options to get species data
  call solver_options_setup

  if (equationOfState .ne. IDEAL_GAS) call die("jet_data: only use ideal gas law for now!")

  ! Zero-out the conserved variables
  conservedVariables = 0.0_WP

  ! Reference density and pressure
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  rho0 = 1.0_WP
  P0 = 1.0_WP / ratioOfSpecificHeats

  ! Read jet parameters
  call parser_read('jet radius', Rj, 0.5_WP)
  call parser_read('jet velocity', Uj)
  call parser_read('co-flow velocity', Uc)
  call parser_read('jet momentum thickness', theta0, 0.05_WP*Rj)
  call parser_read('jet growth rate', growthRate, 1.0_WP)
  call parser_read('jet growth origin', x0, Lx)

  slope = Rj / (4.0_WP * theta0)
  RjI = 1.0_WP / Rj

  do i = 1, nGridPoints

     ! Coordinate
     x = 0.0_WP; y = 0.0_WP; z = 0.0_WP
     if (nx.gt.1) x = coordinates(i, 1)
     if (ny.gt.1) y = coordinates(i, 2)
     if (nz.gt.1) z = coordinates(i, 3)

     ! Determine conditions
     r = sqrt(y**2+z**2)
     density = rho0
     velocity = 0.0_WP
     pressure = P0
     
     ! Initialize mixture fraction and velocity profile
     mixtureFraction = 0.5_WP * ( 1.0_WP + tanh(slope * (Rj/r - RjI*r)                       &
             / (1.0_WP + growthRate * max(0.0_WP, x - x0))) )
     velocity(1) = (uj - uc) * mixturefraction + uc

     ! Assign the density
     conservedVariables(i,1) = density

     ! Assign the momentum
     conservedVariables(i,2:nDimensions+1) = density * velocity(1:nDimensions)
     !conservedVariables(i,3:nDimensions+1) = 0.0_WP

     ! Assign the energy
     conservedVariables(i, nDimensions+2) = pressure / (gamma - 1.0_WP) + 0.5_WP *           &
                sum(conservedVariables(i, 2:nDimensions+1)**2) / conservedVariables(i, 1)

     if (nSpecies .gt. 0) then
        ! Assign mixture fraction profile
        conservedVariables(i, nDimensions+3) = density * mixtureFraction
     end if

  end do

end subroutine jet_data
