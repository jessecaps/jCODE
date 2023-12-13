module ibm_init

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
  real(WP) :: Lx, Ly, dx, dy

end module ibm_init

subroutine ibm_init_grid

  ! Internal modules
  use ibm_init

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: domainHeight, x0, y0, xtilde, ytilde, alpha, minAlpha, maxAlpha, stretchingFactor
  logical :: stretchX, stretchY

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Only implemented in 2D and serial for now
  if (nz.gt.1) call die('ibm cylinder only implemented for 2D!')
  if (nProcs.gt.1) call die('ibm cylinder only implemented in serial!')

  ! Read in the grid size
  call parser_read('Lx', Lx)
  call parser_read('Ly', Ly)

  ! Compute the grid spacing
  dx = Lx / real(nx-1, WP)
  if (periodicityType(2) .eq. PLANE) then
     dy = Ly / real(ny, WP)
     domainHeight = Ly - dy
  else
     dy = Ly / real(ny-1, WP)
     domainHeight = Ly
  end if
  
  ! Should we stretch the mesh?
  call parser_is_defined('stretch x position', stretchX)
  if (stretchX) call parser_read('stretch x position', x0); x0 = x0 - 0.5_WP * Lx
  call parser_is_defined('stretch y position', stretchY)
  if (stretchY) call parser_read('stretch y position', y0)
  if (stretchX .or. stretchY) then
     call parser_read('stretching factor', stretchingFactor, 2.0_WP)
  end if

  ! Create X
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           if (stretchX) then
              xtilde = 2.0_WP * real(i-1,WP) / real(nx-1,WP) - 1.0_WP
              alpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (xtilde - x0 / Lx)))
              minAlpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (-1.0_WP - x0 / Lx)))
              maxAlpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (1.0_WP - x0 / Lx)))
              alpha = (alpha - minAlpha) / (maxAlpha - minAlpha)
              coordinates(grid_index(i,j,k), 1) = Lx * alpha
           else
              coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) / real(nx - 1, WP)
           end if
        end do
     end do
  end do

  ! Create Y
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           if (stretchY) then
              ytilde = 2.0_WP * real(j-1,WP) / real(ny-1,WP) - 1.0_WP
              alpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (ytilde - y0 / Ly)))
              minAlpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (-1.0_WP - y0 / Ly)))
              maxAlpha = 0.5_WP * (1.0_WP + sinh(stretchingFactor * (1.0_WP - y0 / Ly)))
              alpha = (alpha - minAlpha) / (maxAlpha - minAlpha)
              coordinates(grid_index(i,j,k), 2) = domainHeight * alpha - 0.5_WP * Ly
           else
              coordinates(grid_index(i,j,k), 2) = domainHeight * real(j - 1, WP)             &
                   / real(ny - 1, WP) -  0.5_WP * Ly
           end if
        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine ibm_init_grid


subroutine ibm_init_data

  ! Internal modules
  use ibm_init

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use random
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, x, x0, p1, p2, rho1, rho2, u1, u2, density, pressure, velocity
  logical :: useShock

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Uniform flow or shock?
  call parser_read('include shock', useShock, .false.)

  ! Specific heat ratio
  call parser_read('ratio of specific heats', gamma, 1.4_WP)

  if (useShock) then
     ! Read in shock parameters
     call parser_read('shock position', x0)
     call parser_read('post-shock pressure', p1)
     call parser_read('pre-shock pressure', p2)
     call parser_read('post-shock density', rho1)
     call parser_read('pre-shock density', rho2)
     call parser_read('post-shock velocity', u1, 0.0_WP)
     call parser_read('pre-shock velocity', u2, 0.0_WP)
  else
     call parser_read('crossflow mach number', velocity)
     density = 1.0_WP
     pressure = 1.0_WP / gamma
  end if

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           if (useShock) then
              ! Get local coordinate
              x = coordinates(grid_index(i,j,k), 1)

              ! Set the density & pressure
              if (x.le.x0) then
                 density = rho1
                 pressure = p1
                 velocity = u1
              else
                 density = rho2
                 pressure = p2
                 velocity = u2
              end if
           end if

           ! Set the conserved variables
           conservedVariables(grid_index(i,j,k), 1) = density
           conservedVariables(grid_index(i,j,k), 2) = density * velocity
           conservedVariables(grid_index(i,j,k), 3:nDimensions+1) = 0.0_WP
           conservedVariables(grid_index(i,j,k), nDimensions+2) = pressure /                 &
                (gamma - 1.0_WP) + 0.5_WP * density * velocity ** 2
        end do
     end do
  end do

  return
end subroutine ibm_init_data


subroutine ibm_init_levelset

  ! Internal modules
  use ibm_init

  ! External modules
  use string
  use math
  use grid_levelset

  implicit none

  ! Local variables
  integer :: i, j
  real(WP) :: R, L, width, height, x0, y0, distance
  real(WP), dimension(nDimensions) :: pos0, pos
  character(len = str_medium) :: filename, case

  ! Return if not writing a levelset file
  call parser_read('init levelset file', filename, '')
  if (len_trim(filename).eq.0) return

  ! Allocate distance array
  allocate(levelset(nGridPoints,1))
  levelset=0.0_WP

  ! Get IBM case
  call parser_read('ibm case', case)

  select case (trim(case))

  case ('cylinder')

     ! Read in placement parameters
     call parser_read('cylinder radius', R)
     call parser_read('cylinder x position', x0)
     call parser_read('cylinder y position', y0, 0.0_WP)
     pos = 0.0_WP; pos(1) = x0; pos(2) = y0

     do i = 1, nGridpoints
        ! Initialize distance       
        distance = huge(1.0_WP)

        ! Get distance to particle center
        pos0 = coordinates(i,:) - pos(1:nDimensions)

        ! Apply periodicity
        do j = 1, nDimensions
           if (periodicityType(j).eq.PLANE) then
              if (abs(pos0(j)).gt.abs(pos0(j)-periodicLength(j))) then
                 pos0(j) = pos0(j) - periodicLength(j)
              end if
              if (abs(pos0(j)).gt.abs(pos0(j)+periodicLength(j))) then
                 pos0(j) = pos0(j) + periodicLength(j)
              end if
           end if
        end do

        ! Get signed distance to particle interface and compare to current
        distance = min(distance, sqrt(sum(pos0**2)) - R)

        ! Update the levelset
        levelset(i,1) = distance

     end do

  case ('plate')

     ! Read in placement parameters
     call parser_read('plate length', L, 1.0_WP)
     call parser_read('plate width', width)
     call parser_read('plate x position', x0)
     call parser_read('plate y position', y0, 0.0_WP)

     call die('Levelset for plate not yet implemented')

  case ('wedge')

     call die('Levelset for wedge not yet implemented')

     ! Read in placement parameters
     call parser_read('wedge width', width)
     call parser_read('wedge height', height)
     call parser_read('wedge x position', x0)
     call parser_read('wedge y position', y0, 0.0_WP)

     do i = 1, nGridpoints
        ! Initialize distance       
        distance = huge(1.0_WP)

        ! Get distance to right side
        pos0 = coordinates(i,:) - pos(1:nDimensions)

        ! Apply periodicity
        do j = 1, nDimensions
           if (periodicityType(j).eq.PLANE) then
              if (abs(pos0(j)).gt.abs(pos0(j)-periodicLength(j))) then
                 pos0(j) = pos0(j) - periodicLength(j)
              end if
              if (abs(pos0(j)).gt.abs(pos0(j)+periodicLength(j))) then
                 pos0(j) = pos0(j) + periodicLength(j)
              end if
           end if
        end do

        ! Get signed distance to particle interface and compare to current
        distance = min(distance, sqrt(sum(pos0**2)) - R)

        ! Get distance to top

        ! Get distance to bottom

        ! Update the levelset
        levelset(i,1) = distance

     end do

    
  case default

     call die('unknown case: "' // trim(case) // '"')

  end select

  return
end subroutine ibm_init_levelset
