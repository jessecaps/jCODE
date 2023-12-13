module jet_impingement

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
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, maxDxDyDz
  real(WP) :: L, ID1, ID2,OD1, OD2, steepness, L0, D0, wallThickness, innerDiameter,         &
       outerDiameter, gap, inletDiameter, sig, outerInletDiameter, outerExitDiameter, L_br

  ! Nozzle typess
  integer :: nozzleType
  integer, parameter ::                                                                      &
       LINEAR_NOZZLE = 1,                                                                    &
       TANH_NOZZLE   = 2,                                                                    &
       TANH_NOZZLE2  = 3

end module jet_impingement

subroutine jet_impingement_grid

  ! Internal modules
  use jet_impingement

  ! External modules
  use parallel

  ! Local variables
  integer :: i, j, k
  real(WP) :: xmin, xmax
  real(WP) :: ymin, ymax, ymin_uni, ymax_uni, dy_min, dy_max
  real(WP) :: zmin, zmax, zmin_uni, zmax_uni, dz_min, dz_max
  real(WP) :: cornerMin(3), cornerMax(3), cornerMinU(3), cornerMaxU(3), hMin(3), hMax(3)
  real(WP), dimension(:), allocatable :: y1D, z1D
  logical :: stretch


  call parser_read('stretch grid', stretch, .false.)

  if (stretch) then

     ! Should not be used for periodic cases
     do i = 1, nDimensions
        if (periodicityType(i) .ne. NONE)                                                    &
             call die('Stretched grid only for non-periodic domains')
     end do

     call parser_read('corner min', cornerMin)
     call parser_read('corner max', cornerMax)
     call parser_read('corner uniform min', cornerMinU)
     call parser_read('corner uniform max', cornerMaxU)
     call parser_read('h min', hMax)
     call parser_read('h max', hMin)

     nx = globalGridSize(1)
     ! Create y
     ny = 1
     if (globalGridSize(2) .gt. 1) then
        ymin = min(cornerMin(2), cornerMax(2))
        ymax = max(cornerMin(2), cornerMax(2))
        ymin_uni = min(cornerMinU(2), cornerMaxU(2))
        ymax_uni = max(cornerMinU(2), cornerMaxU(2))
        dy_min = min(hMin(2), hMax(2))
        dy_max = max(hMin(2), hMax(2))
        call grid_stretch(ymin, ymax, ymin_uni, ymax_uni, dy_min, dy_max, y1D)
        ny = size(y1d)
        Ly = ymax - ymin
     end if

     ! Create z
     nz = 1
     if (globalGridSize(3) .gt. 1) then
        zmin = min(cornerMin(3), cornerMax(3))
        zmax = max(cornerMin(3), cornerMax(3))
        zmin_uni = min(cornerMinU(3), cornerMaxU(3))
        zmax_uni = max(cornerMinU(3), cornerMaxU(3))
        dz_min = min(hMin(3), hMax(3))
        dz_max = max(hMin(3), hMax(3))
        call grid_stretch(zmin, zmax, zmin_uni, zmax_uni, dz_min, dz_max, z1D)
        nz = size(z1d)
        Lz = zmax - zmin
     end if

     ! Update grid metrics
     globalGridSize(1) = nx
     globalGridSize(2) = ny
     globalGridSize(3) = nz
     call grid_cleanup
     call operator_cleanup
     call geometry_setup
     call grid_setup
     call operator_setup
     
     ! Create x
     
     call parser_read('Lx', Lx, 0.0_WP)
     ! Compute the grid spacing
     if (periodicityType(1) .eq. PLANE) then
        dx = Lx / real(nx, WP)
     else
        dx = Lx / real(nx-1, WP)
     end if
     
     xmin = huge(1.0_WP)
     xmax = -huge(1.0_WP)

     ! Generate the grid
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              if (nx .gt. 1) then
                 if (periodicityType(1) .eq. PLANE) then
                    coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1,WP) /        &
                         real(nx - 1, WP)
                 else
                    coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /               &
                         real(nx - 1, WP)
                 end if
                 xmin = min(xmin, coordinates(grid_index(i,j,k), 1))
                 xmax = max(xmax, coordinates(grid_index(i,j,k), 1))
              end if
              if (ny.gt.1) coordinates(grid_index(i,j,k), 2) = y1D(j)
              if (nz.gt.1) coordinates(grid_index(i,j,k), 3) = z1D(k)
           end do
        end do
     end do
     call parallel_max(xmax)
     call parallel_min(xmin)
     
     if (nDimensions .gt. 1) dx = min(dx, dy_min)
     if (nDimensions .gt. 2) dx = min(dx, dy_min, dz_min)

  else

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

     maxDxDyDz = dx
     if (nDimensions .gt. 1) maxDxDyDz = max(maxDxDyDz, dy)
     if (nDimensions .gt. 2) maxDxDyDz = max(maxDxDyDz, dz)

     ! Generate the grid
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)

              ! Create X
              if (nx .gt. 1) then
                 if (periodicityType(1) .eq. PLANE) then
                    coordinates(grid_index(i,j,k), 1) = (Lx - dx) *  real(i - 1, WP) /       &
                         real(nx - 1, WP)
                 else
                    coordinates(grid_index(i,j,k), 1) = Lx * real(i - 1, WP) /               &
                         real(nx - 1, WP)
                 end if
              end if

              ! Create Y
              if (ny .gt. 1) then
                 if (periodicityType(2) .eq. PLANE) then
                    coordinates(grid_index(i,j,k), 2) = (Ly - dy) *  real(j - 1, WP) /       &
                         real(ny - 1, WP) - 0.5_WP * Ly
                 else
                    coordinates(grid_index(i,j,k), 2) = Ly * real(j - 1, WP) /               &
                         real(ny - 1, WP) - 0.5_WP * Ly
                 end if
              end if

              ! Create Z
              if (nz .gt. 1) then
                 if (periodicityType(3) .eq. PLANE) then
                    coordinates(grid_index(i,j,k), 3) = (Lz - dz) *  real(k - 1, WP) /       &
                         real(nz - 1, WP) - 0.5_WP * Lz
                 else
                    coordinates(grid_index(i,j,k), 3) = Lz * real(k - 1, WP) /               &
                         real(nz - 1, WP) - 0.5_WP * Lz
                 end if
              end if

           end do
        end do
     end do

  end if

  ! Setup the grid metrics
  call grid_metrics_setup

contains

  ! --------------------------------------------------------- !
  ! Generate a stretched grid (in 1D)                         !
  !                                                           !
  !   dx_max                 dx_min                   dx_max  !
  !  |-----|---|--|-|-|-|-|-|-|-|-|-|-|-|-|-|-|--|---|-----|  !
  ! xmin        xmin_uni                   xmax_uni      xmax !
  ! --------------------------------------------------------- !

  subroutine grid_stretch(xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max, x)

    implicit none

    ! Arguments
    real(WP), intent(in) :: xmin, xmax, xmin_uni, xmax_uni, dx_min, dx_max
    real(WP), dimension(:), allocatable, intent(out) :: x

    ! Local variables
    integer :: nu, ns1, ns2, nx
    real(WP) :: ratio1, ratio2, ratio, sc

    ! sanity check
    if (dx_min .gt. dx_max) call die('dx_min is larger than dx_max')
    if (xmin_uni .gt. xmax_uni) call die('xmin_uni is larger than xmax_uni')

    sc = min(dx_min, abs(xmin_uni-xmin), abs(xmax-xmax_uni), abs(xmax_uni-xmin_uni))
    if (abs(sc-dx_min) .gt. epsilon(1.0_WP)) call die('dx_min too large')

    sc = min(dx_max, abs(xmin_uni-xmin), abs(xmax-xmax_uni), abs(xmax_uni-xmin_uni))
    if (abs(sc-dx_max) .gt. epsilon(1.0_WP)) call die('dx_max too large')

    ! Determine number of grid points and stretching rate
    nu = ceiling((xmax_uni-xmin_uni) / dx_min) ! Number of uniform grid
    if (abs(dx_min-dx_max) .lt. epsilon(1.0_WP)) then
       ratio1 = 1.0_WP; ratio2 = 1.0_WP
       ns1 = ceiling((xmin_uni-xmin) / dx_min) + 1
       ns2 = ceiling((xmax-xmax_uni) / dx_min) + 1
    else
       ! Stretching ratio on the left
       ratio1 = (xmin_uni-xmin)/(xmin_uni-xmin-dx_max+dx_min)
       ! Stretching ratio on the right
       ratio2 = (xmax-xmax_uni)/(xmax-xmax_uni-dx_max+dx_min)

       ! Preserve the same ratio around the uniform zone
       ratio = min(ratio1,ratio2)

       !nb of points from xmin to xmin_uni (xmin_uni included)
       ns1 = ceiling(1  + log(1+(ratio-1)*(xmin_uni-xmin)/(dx_min*ratio))/log(ratio))
       !nb of points from xmax to xmax_uni included (xmax_uni included)
       ns2 = ceiling(1  + log(1+(ratio-1)*(xmax-xmax_uni)/(dx_min*ratio))/log(ratio))

       ! Update (instead of using a newton algorithm to get the exact ratio)
       ratio = (dx_max/dx_min)**(1/(real(max(ns1,ns2)-1,WP)))
       ns1 = ceiling(1  + log(1+(ratio-1)*(xmin_uni-xmin)/(dx_min*ratio))/log(ratio))
       ns2 = ceiling(1  + log(1+(ratio-1)*(xmax-xmax_uni)/(dx_min*ratio))/log(ratio))

    end if

    ! Total number of grid points
    nx = nu + ns1 + ns2 - 1

    ! Generate the grid
    allocate(x(nx))

    ! left
    x(ns1) = xmin_uni
    do i = 1,ns1-1
       x(ns1-i) = x(ns1-i+1) - dx_min * ratio**real(i,WP)
    end do

    ! uniform
    do i=ns1+1,ns1+nu-1
       x(i) = x(i-1) + dx_min
    end do
    x(ns1+nu) = xmax_uni

    !right
    do i = ns1+nu+1,nx
       x(i) = x(i-1) + dx_min * ratio**real(i-(ns1+nu),WP)
    end do

    !adjust the boundaries
    x(1) = xmin
    x(nx) = xmax

    !prevent the first and last interval to be too small or larger than dx_max
    if (x(2)-xmin < x(3)-x(2)) then
       x(2) = (xmin + x(3)) / 2.0_WP
    end if

    if (xmax-x(nx-1) < x(nx-1)-x(nx-2)) then
       x(nx-1) = (xmax + x(nx-2)) / 2.0_WP
    end if

    return
  end subroutine grid_stretch

end subroutine jet_impingement_grid


subroutine jet_impingement_data

  ! Internal modules
  use jet_impingement

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: gamma, density, pressure, velocity, x, r, x0, U0, p0, T0, rho0, D, buf
  character(len = str_medium) :: inputString

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  ! Select the nozzle type and read in geometric parameters
  call parser_read('nozzle type', inputString, 'linear')
  select case (trim(inputString))
  case ('linear', 'Linear', 'LINEAR')
     nozzleType = LINEAR_NOZZLE

     call parser_read('inner diameter 1', ID1)
     call parser_read('inner diameter 2', ID2)
     call parser_read('length', L)
     call parser_read('wall thickness', wallThickness)
     call parser_read('straight length', L0)

  case ('tanh', 'TANH')
     nozzleType = TANH_NOZZLE

     call parser_read('inner diameter 1', ID1)
     call parser_read('inner diameter 2', ID2)
     call parser_read('outer diameter 1', OD1)
     call parser_read('outer diameter 2', OD2)
     call parser_read('length', L)
     call parser_read('steepness', steepness, 10.0_WP)
     call parser_read('convergence location', L0, 0.5_WP * L)
     
  case ('tanh2', 'TANH2')
     nozzleType = TANH_NOZZLE2
     
     ! Read in geometric parameters
     call parser_read('nozzle inner exit diameter', innerDiameter)
     call parser_read('nozzle outer inlet diameter', outerInletDiameter)
     call parser_read('nozzle outer exit diameter', outerExitDiameter)
     call parser_read('nozzle inner inlet diameter', inletDiameter)
     ! call parser_read('nozzle gap', gap)
     call parser_read('straight section length', L_br)
     call parser_read('nozzle tanh constant', sig)
     call parser_read('nozzle length', L)
     call parser_read('nozzle break length', L0, 0.9_WP*L)
     call parser_read('nozzle break diameter', D0, 1.1_WP*innerDiameter)

  case default
     call die("nozzle_data: invalid nozzle type '" // trim(inputString) // "'!")

  end select

  ! Read in fluid parameters
  call parser_read('ratio of specific heats', gamma, 1.4_WP)
  call parser_read('nozzle temperature', T0, 1.0_WP / (gamma-1.0_WP))
  call parser_read('nozzle pressure', p0)
  call parser_read('nozzle velocity', U0)

  ! Apply ideal gas law
  rho0 = gamma / (gamma - 1.0_WP) * p0 / T0

  ! Set the conserved variables
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Get position
           x = coordinates(grid_index(i,j,k), 1)
           r = sqrt(sum(coordinates(grid_index(i,j,k), 2:nDimensions)**2))

           ! Set state variables based on nozzle geometry
           select case (nozzleType)
           case (LINEAR_NOZZLE)
              x0 = L0
              if ((r.le.0.5_WP*ID1 .and. x.le.x0)) then
                 buf = 1.0_WP
              else
                 buf = 0.0_WP
              end if
              velocity = (U0 - 0.0_WP) * buf + 0.0_WP
              density = (rho0 - 1.0_WP) * buf + 1.0_WP
              pressure = (p0 - 1.0_WP / gamma) * buf + 1.0_WP/gamma


           case (TANH_NOZZLE)
              x0 = 1.0_WP*L
              D = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (x - L0)))
              if (r.le.0.5_WP * D .and. x.le.L) then
                 buf = 1.0_WP
              else
                 buf = 0.0_WP
              end if
              velocity = (U0 - 0.0_WP) * buf + 0.0_WP
              density = (rho0 - 1.0_WP) * buf + 1.0_WP
              pressure = (p0 - 1.0_WP / gamma) * buf + 1.0_WP/gamma

           case (TANH_NOZZLE2)
              x0 = 1.0_WP*L_br
              gap = 0.5_WP*(outerExitDiameter-innerDiameter)
              outerdiameter = outerInletDiameter-2.0_WP*gap
              ! Get position
              x = coordinates(grid_index(i,j,k), 1)
              r = sqrt(sum(coordinates(grid_index(i,j,k), 2:nDimensions)**2))

              ! Get local nozzle diameter and set state variables
              D = D0 + (D0-inletDiameter) * tanh(sig*20.0_WP*(x-L_br-L0))
              if (r-0.025_WP*D.le.0.5_WP*D) then
                 velocity = 0.5_WP*(0.0_WP+U0) + 0.5_WP*(0.0_WP-U0) * tanh(5.0_WP*(x-L_br-x0))
                 density = 0.5_WP*(1.0_WP+rho0) + 0.5_WP*(1.0_WP-rho0) * tanh(5.0_WP*(x-L_br-x0))
                 pressure = 0.5_WP*(1.0_WP/gamma+p0) + 0.5_WP*(1.0_WP/gamma-p0) * tanh(5.0_WP*(x-L_br-x0))
                 if (x .le. L) then
                    velocity = U0
                    density = rho0
                    pressure = p0
                 end if
              else
                 velocity = 0.0_WP
                 density  = 1.0_WP
                 pressure = 1.0_WP / gamma
              end if
           end select


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
end subroutine jet_impingement_data


subroutine jet_impingement_levelset

  ! Internal modules
  use jet_impingement

  ! External modules
  use math
  use grid_levelset
  use parallel

  implicit none

  ! General
  integer :: gridIndex, i, j, k, ip
  real(WP), parameter :: eps = 1.0e-9_WP
  real(WP) :: r, x, y, z, rp, xp, rin, rout, ID1e, ID2e, xComp, rComp, straightHeightLevel
  real(WP) :: sig2, rtgt
  integer, dimension(nGridPoints) :: normIndex
  real(WP), dimension(nGridPoints) :: thetaGrid, normIndicator
  real(WP), dimension(2) :: tempVec

  ! Levelset
  integer :: count
  real(WP) :: mydist

  ! Progress monitoring
  integer :: togo, prog
  integer :: iratio_new, iratio_old

  ! Levelset file
  real(WP) :: tmp

  ! Allocate distance array
  allocate(levelset(nGridPoints,1)); levelset=huge(1.0_WP)

  select case (nozzleType)

  case (LINEAR_NOZZLE)

     ! Define the nozzle levelset using straight lines
     ! -----------------------------------------------

     select case (nDimensions)   

     case(2)
        xComp = sin(atan((ID1*0.5_WP - ID2*0.5_WP)/L))
        rComp = sqrt(1.0_WP**2 - xComp**2)

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = 0.0_WP
           if (y .le. -1e-15_WP) thetaGrid(i) = 1.5_WP*pi
           if (y .ge. 1e-15) thetaGrid(i) = 0.5_WP*pi

           if (abs(y) .lt. 1e-15_WP) then
              thetaGrid(i) = 0.0_WP
           end if

           levelset(i,1) = y**2 + z**2
        end do

        levelset(:,1) = sqrt(levelset(:,1))

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)

           if (x .le. L0) then
              rin = 0.5_WP * ID1
           else
              rin = (ID2*0.5_WP - ID1*0.5_WP)/ L * x + ID1 * 0.5_WP                          &
                   - (ID2*0.5_WP - ID1*0.5_WP) / L * L0
           end if
           rout = rin + wallThickness

           normIndicator(i) = sign(1.0_WP,levelset(i,1) - (rin + rout)*0.5_WP)

           if(abs(levelset(i,1) - (rin + rout)*0.5_WP) .lt. 1e-15_WP) then 
              normIndicator(i) = 1.0_WP
           end if

           tempVec(1) = abs(levelset(i,1) - (rin + rout)*0.5_WP) - (rout - rin)*0.5_WP
           tempVec(2) = x - (L + L0) 
           straightHeightLevel = x - L0

           normIndex(i) = maxloc(tempVec(:), DIM = 1)
           levelset(i,1) = tempVec(normIndex(i))
           if(abs(tempVec(1) - tempVec(2)) .lt. 1e-15_WP) then
              levelset(i,1) = tempVec(2)
              normIndex(i) = 2    
           end if

        end do

     case(3)
        xComp = sin(atan((ID1*0.5_WP - ID2*0.5_WP)/L))
        rComp = sqrt(1.0_WP**2 - xComp**2)   

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = coordinates(i,3)

           thetaGrid(i) = pi - 0.5_WP*pi*(1.0_WP + sign(1.0_WP,z))                           &
                * (1.0_WP - sign(1.0_WP,y**2)) - 0.25_WP * pi * (2.0_WP + sign(1.0_WP, z))   &
                * (sign(1.0_WP,y)) - sign(1.0_WP, y*z)                                       &
                * atan((abs(z) - abs(y))/(abs(z) + abs(y)))
           if (abs(y) .lt. 1e-14_WP .and. abs(z) .lt. 1e-14_WP) then
              thetaGrid(i) = 0.0_WP
           end if

           levelset(i,1) = y**2 + z**2
        end do

        levelset(:,1) = sqrt(levelset(:,1))

        do i = 1, nGridPoints
           x = coordinates(i,1)
           y = coordinates(i,2)
           z = coordinates(i,3)

           if (x .le. L0) then
              rin = ID1*0.5_WP
              rout = rin + wallThickness
           else
              rin = (ID2*0.5_WP - ID1*0.5_WP) / L * x + ID1*0.5_WP                           &
                   - (ID2*0.5_WP - ID1*0.5_WP)/L * L0
              rout = rin + wallThickness
           end if

           normIndicator(i) = sign(1.0_WP,levelset(i,1) - (rin + rout)*0.5_WP)

           if(abs(levelset(i,1) - (rin + rout)*0.5_WP) .lt. 4e-14_WP) then 
              normIndicator(i) = 1.0_WP
           end if

           tempVec(1) = abs(levelset(i,1) - (rin + rout)*0.5_WP) - (rout - rin)*0.5_WP
           tempVec(2) = x - (L + L0)
           straightHeightLevel = x - L0

           normIndex(i) = maxloc(tempVec(:), DIM = 1)
           levelset(i,1) = tempVec(normIndex(i))
           if(abs(tempVec(1) - tempVec(2)) .lt. 4e-14_WP) then
              levelset(i,1) = tempVec(2)
              normIndex(i) = 2    
           end if

        end do
     end select

  case (TANH_NOZZLE)

     ! Define the nozzle levelset as a tanh function
     ! ---------------------------------------------

     ! Prepare counter
     prog = 0
     togo = nGridPoints
     iratio_old = 0

     ! Store effective inner diameters
     ID1e = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (0.0_WP - L0)))
     ID2e = ID2 + (ID1 - ID2) * 0.5_WP * (1.0_WP - tanh(steepness / L * (L - L0)))

     ! Compute projection
     if (iRank .eq. iRoot) print *,'Computing levelset analytically...'
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              ! Prepare projections
              count = 0
              mydist = huge(1.0_WP)

              ! Get position
              x = coordinates(gridIndex, 1)
              r = sqrt(sum(coordinates(gridIndex, 2:nDimensions)**2))

              ! Determine radius to inner and outer nozzle geometry
              rin = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                                   &
                   (1.0_WP - tanh(steepness / L * (x - L0))))
              rout = 0.5_WP * OD1 + 0.5_WP * (OD2 - OD1) * x / L

              ! Calculate distance
              if (r.lt.rin .and. x.le.L) then
                 ! Inner nozzle
                 do ip = 1,1000
                    xp = L * real(ip - 1,WP) / real(1000 - 1, WP)
                    rp = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                              &
                         (1.0_WP - tanh(steepness / L * (xp - L0))))
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp .lt. mydist) mydist = tmp
                 end do
              elseif (r .gt. rout .and. x.le.L) then
                 ! Cone
                 mydist = (r-rout) * cos(atan(0.5_WP*(OD2 - ID2e) / L))
              elseif (r.ge.rin .and. r.le.rout .and. x.le.L) then
                 ! In the solid
                 do ip = 1,1000
                    xp = L * real(ip - 1,WP) / real(1000 - 1, WP)
                    rp = 0.5_WP * (ID2 + (ID1 - ID2) * 0.5_WP *                              &
                         (1.0_WP - tanh(steepness / L * (xp - L0))))
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp.lt.mydist) mydist = tmp
                 end do
                 mydist = min(mydist, abs((rout-r) * cos(atan(0.5_WP * (OD2 - ID2e) / L))))
                 mydist = min(mydist, abs(x-L))
                 mydist = -mydist
              elseif (r.gt.0.5_WP * OD2 .and. x.gt.L) then
                 mydist = sqrt((x - L)**2 + (r - 0.5_WP * OD2)**2)
              elseif (r.le.0.5_WP * OD2 .and. r.ge.0.5_WP * ID2e .and. x.gt.L) then
                 mydist = x - L
              else
                 mydist = sqrt((x - L)**2 + (r - 0.5_WP * ID2e)**2)
              end if

              ! Store the minimum distance
              levelset(gridIndex,1) = mydist

              ! Add point to counter
              if (irank.eq.iroot) then
                 prog=prog+1
                 iratio_new=int(real(prog,WP)/real(togo,WP)*100.0_WP)
                 if (iratio_new.ge.iratio_old+20) then
                    iratio_old=iratio_new
                    write(*,'(i3,x,a1)') iratio_new,'%'
                 end if
              end if

           end do
        end do
     end do
     
  case(TANH_NOZZLE2)
     !L = L + L_br 
     ! Prepare counter
     prog=0
     togo=nGridPoints
     iratio_old=0
     if (irank .eq. iroot) then
        print*, 'Length of nozzle', L
        L = L + L_br
     end if
     ! Compute projection
     if (irank.eq.iroot) print *,'Computing levelset analytically...'
     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              gridIndex = i - gridOffset(1) + localGridSize(1) *                             &
                   (j - 1 - gridOffset(2) + localGridSize(2) *                               &
                   (k - 1 - gridOffset(3)))

              ! Prepare projections
              count = 0
              mydist = huge(1.0_WP)

              ! Get position
              x = coordinates(gridIndex, 1)
              r = sqrt(sum(coordinates(gridIndex, 2:nDimensions)**2))
              ! Calculate levelset analytically
              sig2 = sig*(D0-inletDiameter)/(innerDiameter-D0)
              rin = 0.5_WP*max(D0 + (D0-inletDiameter) * tanh(sig*20.0_WP*(x-L_br-L0)),D0 +  &
                   (innerDiameter-D0) * tanh(sig2*25.0_WP*(x-L0)))
              if (x .le. L_br) then
                 rin =  0.5_WP*inletDiameter
              end if
              rout = gap + 0.5_WP * outerDiameter + 0.5_WP *                                 &
                   (innerDiameter - outerDiameter) * x / L
              rtgt = gap + 0.5_WP * innerDiameter +  L * (x-L) /                             &
                   (0.5_WP * (outerDiameter - innerDiameter))

              ! Calculate normal and distance
              if (r.lt.rin .and. x.le.L) then
                 do ip = 1,1000
                    xp = L*real(ip-1,WP)/real(1000-1,WP)
                    rp = 0.5_WP*max(D0 + (D0-inletDiameter) * tanh(sig*20.0_WP*(xp-L_br-L0)),&
                         D0 + (innerDiameter-D0) * tanh(sig2*25.0_WP*(xp-L_br-L0)))
                    if (xp .le. L_br) then
                       !print*, 'rp: ', rp
                       rp =  0.5_WP*inletDiameter
                    end if
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp.lt.mydist) mydist = tmp
                 end do
              elseif (r.gt.max(rout,rtgt)) then
                 ! Cone
                 mydist = (r-rout)*cos(atan(0.5_WP*(outerDiameter - innerDiameter) / L))
              elseif (r.ge.rin .and. r.le.rout .and. x.le.L) then
                 ! Iterior
                 do ip = 1,1000
                    xp = L*real(ip-1,WP)/real(1000-1,WP)
                    rp = 0.5_WP*max(D0 + (D0-inletDiameter) * tanh(sig*20.0_WP*(xp-L_br-L0)),&
                         D0 + (innerDiameter-D0) * tanh(sig2*25.0_WP*(xp-L_br-L0)))
                    if (xp .le. L_br) then
                       rp =  0.5_WP*inletDiameter
                    end if
                    tmp = sqrt((x-xp)**2+(r-rp)**2)
                    if (tmp.lt.mydist) mydist = tmp
                 end do
                 mydist = min(mydist, abs((rout-r) * cos(atan(0.5_WP *                       &
                      (outerDiameter - innerDiameter) / L))))
                 mydist = min(mydist, abs(x-L))
                 mydist = -mydist 
              elseif (r.le.rtgt .and. r.gt.(gap + 0.5_WP * innerDiameter) .and. x.gt.L) then
                 mydist = sqrt((x-L)**2+(r-(gap + 0.5_WP * innerDiameter))**2)
              elseif (r.le.(gap + 0.5_WP * innerDiameter) .and. r.ge.0.5_WP *                &
                   innerDiameter .and. x.gt.L) then
                 mydist = x-L
              else
                 mydist = sqrt((x-L)**2+(r-0.5_WP * innerDiameter)**2)
              end if

              ! Postprocess distance
              levelset(gridIndex,1) = mydist

              ! Add point to counter
              if (irank.eq.iroot) then
                 prog=prog+1
                 iratio_new=int(real(prog,WP)/real(togo,WP)*100.0_WP)
                 if (iratio_new.ge.iratio_old+20) then
                    iratio_old=iratio_new
                    write(*,'(i3,x,a1)') iratio_new,'%'
                 end if
              end if

           end do
        end do
     end do
  

  end select

  return
end subroutine jet_impingement_levelset


subroutine jet_impingement_particles

  ! Internal modules
  use jet_impingement

  ! External modules
  use random
  use math
  use parallel
  use particle

  implicit none

  ! Local variables
  integer :: i, particleOffset, npartx, nparty, npartz, ix, iy, iz
  real(WP) :: diameter, domainVolume, volumeFraction, particleVolume, volumeFactor,          &
       sumParticleVolume, Lp, Lpx, Lpz, rand, hbed, Tp
  character(len = str_medium) :: filename, particleDistribution

  ! Return if not writing a particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename) .eq. 0) return 

  ! Initialize the random number generator
  call random_setup

  ! Read in the particle parameters
  call parser_read('particle diameter', diameter)
  call parser_read('particle volume fraction', volumeFraction)
  call parser_read('bed height', hbed)
  call parser_read('particle temperature', Tp, 2.5_WP)

  ! Compute domain volume
  domainVolume = sum(gridNorm(:,1))
  call parallel_sum(domainVolume)

  ! Get particle volume
  if (nDimensions .eq. 2) then
     volumeFactor = 1.0_WP / 4.0_WP
  else if (nDimensions .eq. 3) then
     volumeFactor = 1.0_WP / 6.0_WP
  end if
  particleVolume = pi * volumeFactor * diameter ** nDimensions

  ! Get the particle distribution type
  call parser_read('particle distribution', particleDistribution, 'random')

  ! Determine number of particles
  select case (trim(particleDistribution))
  case ('random')

     nParticlesGlobal = int(volumeFraction * domainVolume / particleVolume)
     
  case ('uniform')

     ! Mean interparticle distance
     npartx = 1; nparty = 1; npartz = 1
     Lp = (particleVolume / volumeFraction)**(1.0_WP / real(nDimensions, WP))
     nparty = int(Ly / Lp)
     Lp = Ly / real(nparty, WP)
     npartx = int(hbed / Lp)
     Lpx = hbed / real(npartx, WP)
     if (nz .gt. 1) npartz = int(Lz / Lp)
     Lpz = Lz / real(npartz, WP)
     nParticlesGlobal = npartx * nparty * npartz

  case default

     call die("Unknown particle distribution '" // trim(particleDistribution) // "'")

  end select

  ! Distribute particles to processors
  call pigeon_hole(nParticlesGlobal, nProcs, iRank, particleOffset, nParticles)

  ! Allocate the particle vector
  if (allocated(particles)) deallocate(particles)
  allocate(particles(nParticles))

  ! Initialize the particles
  sumParticleVolume = 0.0_WP
  do i = 1, nParticles

     ! Unique id
     particles(i)%id = int(particleOffset + i, kind = 8)

     ! Particle diameter
     particles(i)%diameter = diameter

     ! Particle temperature
     particles(i)%temperature = Tp

     ! Other parameters
     particles(i)%velocity = 0.0_WP
     particles(i)%angularVelocity = 0.0_WP
     particles(i)%collision = 0.0_WP
     particles(i)%torque = 0.0_WP

     ! Distribute the particles
     particles(i)%position = 0.0_WP
     select case (trim(particleDistribution))
     case ('random')
        if (nx .gt. 1) then
           call random_number(rand)
           particles(i)%position(1) = Lx - hbed * rand
        end if
        if (ny .gt. 1) then
           call random_number(rand)
           particles(i)%position(2) = 0.5_WP * Ly * (2.0_WP * rand - 1.0_WP)
        end if
        if (nz .gt. 1) then
           call random_number(rand)
           particles(i)%position(3) = 0.5_WP * Lz * (2.0_WP * rand - 1.0_WP)
        end if

     case ('uniform')
        iy = (particleOffset + i - 1) / (npartx * npartz)
        ix = (particleOffset + i - 1 - npartx * npartz * iy) / npartz
        iz = particleOffset + i - 1 - npartx * npartz * iy - npartz * ix
        call random_number(rand)        
        particles(i)%position(2) = (real(iy, WP) + 0.5_WP) * Lp - 0.5_WP * Ly                &
                 + 0.8_WP * (Lp - diameter) * (rand - 0.5_WP)
        if (particles(i)%position(2) .gt. 0.5_WP*Ly - 0.5_WP*diameter .or.                   &
             particles(i)%position(2) .lt.  -0.5_WP*Ly + 0.5_WP*diameter) then
           particles(i)%position(2) =  (real(iy, WP) + 0.5_WP) * Lp - 0.5_WP * Ly
        end if

        particles(i)%position(1) = Lx - hbed - 0.5_WP * diameter +                           &
             (real(ix, WP) + 0.5_WP) * Lpx
        if (nz .gt. 1) then
           particles(i)%position(3) = (real(iz, WP) + 0.5_WP) * Lpz - 0.5_WP * Lz
           call random_number(rand)        
           particles(i)%position(3) = (real(iz, WP) + 0.5_WP) * Lpz - 0.5_WP * Lz            &
                + 0.8_WP*(Lpz-diameter)*(rand-0.5_WP)
           if (particles(i)%position(3) .gt. 0.5_WP*Lz - 0.5_WP*diameter .or.                &
                particles(i)%position(3) .lt.  -0.5_WP*Lz + 0.5_WP*diameter) then
              particles(i)%position(3) =  (iz + 0.5_WP) * Lpz - 0.5_WP * Lz
           end if
        end if
     end select

     ! The particle exists
     particles(i)%stop = 0

     ! Sum particle volume
     sumParticleVolume = sumParticleVolume + pi * volumeFactor *                             &
          particles(i)%diameter ** nDimensions
  end do

  ! Compute effective volume fraction
  call parallel_sum(sumParticleVolume)
  volumeFraction = sumParticleVolume / domainVolume

  ! Output stuff to the screen.
  if (iRank .eq. iRoot) then
     print *
     print *, 'Number of particles: ', nParticlesGlobal
     print *, 'Volume fraction: ', volumeFraction
  end if

  return
end subroutine jet_impingement_particles

 
