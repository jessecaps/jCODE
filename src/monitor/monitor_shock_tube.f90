module monitor_shock_tube

  ! External modules
  use monitor

  implicit none

  ! Global variables
  integer :: i1, i2
  real(WP) :: threshold

end module monitor_shock_tube


! ========================================= !
! Setup the routine to monitor a Shock Tube !
! ========================================= !
subroutine monitor_shock_tube_setup

  ! Internal modules
  use monitor_shock_tube

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use grid_functions
  use particle

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP) :: x1, x2
  real(WP), dimension(:), allocatable :: x

  if (trim(simulationName).ne.'shock tube' .or. .not. useParticles) return

  ! Read from input
  call parser_read('particle threshold', threshold, 0.1_WP)
  call parser_read('pressure position 1',x1, 0.0_WP)
  call parser_read('pressure position 2',x2, 0.0_WP)

  ! Get index of pressure probes
  allocate(x(globalGridSize(1))); x = 0.0_WP
  j = iStart(2); k = iStart(3)
  do i = iStart(1), iEnd(1)  
     x(i) = coordinates(grid_index(i,j,k), 1)
  end do
  call parallel_sum_dir(x, 1)
  i1 = globalGridSize(1); i2 = globalGridSize(1)
  do i = 1, globalGridSize(1)-1
     if (x1.ge.x(i) .and. x1.lt.x(i+1)) i1 = i
     if (x2.ge.x(i) .and. x2.lt.x(i+1)) i2 = i
  end do

  ! Set the monitor names
  call monitor_create('shock_tube', 8)
  call monitor_set_header(1, 'p1', 'r')
  call monitor_set_header(2, 'p2', 'r')
  call monitor_set_header(3, 'xp_left', 'r')
  call monitor_set_header(4, 'xp_right', 'r')
  call monitor_set_header(5, 'xp_1%', 'r')
  call monitor_set_header(6, 'xp_99%', 'r')
  call monitor_set_header(7, 'min_xp', 'r')
  call monitor_set_header(8, 'max_xp', 'r')
  
  return
end subroutine monitor_shock_tube_setup

! ================================== !
! Compute Shock Tube statistics      !
! ================================== !
subroutine monitor_shock_tube_timestep
  ! Internal modules
  use monitor_shock_tube

  ! External modules
  use parallel
  use solver_options
  use grid_functions
  use state
  use particle

  implicit none

  ! Local variables
  integer, parameter :: nbin = 100
  integer :: i, j, k, nx, iref, gridIndex
  real(WP) :: minValue, X1Dold, Xmix, minPart, maxPart, buf, xleft(2), xright(2), pdf(nbin), &
       meanP1, meanP2, vol1, vol2
  real(WP), dimension(:), allocatable :: x, X1D, normSum

  if (trim(simulationName).ne.'shock tube' .or. .not. useParticles) return

  ! Set up vectors
  nx = globalGridSize(1)
  pdf = 0.0_WP
  allocate(x(nx)); x = 0.0_WP
  allocate(X1D(nx))
  allocate(normSum(nx)); normSum = 0.0_WP

  ! Get the horizontal coordinate
  j = iStart(2); k = iStart(3)
  do i = iStart(1), iEnd(1)  
     x(i) = coordinates(grid_index(i,j,k), 1)
  end do
  call parallel_sum_dir(x, 1)

  xleft = -huge(1.0_WP)
  xright = -huge(1.0_WP)

  X1D = 0.0_WP
  select case (nDimensions)
  case (2)
     do i = iStart(1), iEnd(1)
        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              gridIndex = grid_Index(i,j,k)
              X1D(i) = X1D(i) + (1.0_WP - volumeFraction(gridIndex,1)) *                     &
                   gridSpacing(gridIndex,2)
              normSum(i) = normSum(i) + gridSpacing(gridIndex,2)
           end do
        end do
     end do
  case (3)
     do i = iStart(1), iEnd(1)
        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              gridIndex = grid_Index(i,j,k)
              X1D(i) = X1D(i) + (1.0_WP - volumeFraction(gridIndex,1)) *                     &
                   gridSpacing(gridIndex,2) * gridSpacing(gridIndex,3)
              normSum(i) = normSum(i) + gridSpacing(gridIndex,2) *                           &
                   gridSpacing(gridIndex,3)
           end do
        end do
     end do
  end select
  call parallel_sum(X1D)
  call parallel_sum(normSum)
  X1D = X1D / normSum

  ! Calculating minimum from left to right
  minValue = huge(1.0_WP)
  X1Dold = -huge(1.0_WP)
  do i = 1, nx     
     Xmix = threshold - X1D(i)
     if ((Xmix .lt. minValue) .and. (X1D(i) .gt. X1Dold)) then
        minValue = Xmix
        iref = i
        if (Xmix .lt. 0) then
           exit
        end if
     end if
     X1Dold = X1D(i)
  end do

  if (X1D(iref) .gt. threshold) then
     xleft(1) = (x(iref) - x(iref-1)) / (X1D(iref) - X1D(iref-1)) *                          &
          (threshold - X1D(iref)) + x(iref)
  else if (X1D(iref) .lt. threshold) then
     xleft(1) = (x(iref) - x(iref+1)) / (X1D(iref) - X1D(iref+1)) *                          &
          (threshold - X1D(iref)) + x(iref)
  else
     xleft(1) = x(iref)
  end if

  ! Calculating minimum from right to left
  minValue = huge(1.0_WP)
  X1Dold = -huge(1.0_WP)
  do i = nx, 1, -1     
     Xmix = threshold - X1D(i)
     if ((Xmix .lt. minValue) .and. (X1D(i) .gt. X1Dold)) then
        minValue = Xmix
        iref = i
        if (Xmix .lt. 0) then 
           exit
        end if
     end if
     X1Dold = X1D(i)
  end do

  if (X1D(iref) .gt. threshold) then
     xright(1) = (x(iref) - x(iref+1)) / (X1D(iref) - X1D(iref+1)) *                         &
          (threshold - X1D(iref)) + x(iref)
  else if (X1D(iref) .lt. threshold) then
     xright(1) = (x(iref) - x(iref-1)) / (X1D(iref) - X1D(iref-1)) *                         &
          (threshold - X1D(iref)) + x(iref)
  else
     xright(1) = x(iref)
  end if

  ! Determine minimum and maximum particle x-position
  minPart =  huge(1.0_WP)
  maxPart = -huge(1.0_WP)
  do i = 1, nParticles
     minPart = min(minPart, particles(i)%position(1))
     maxPart = max(maxPart, particles(i)%position(1))
  end do
  call parallel_min(minPart)
  call parallel_max(maxPart)

  ! Generate PDF of particle position
  do i = 1, nParticles
     j = min(floor((particles(i)%position(1)-minPart) / (maxPart-minPart) * nbin) + 1, nbin)
     pdf(j) = pdf(j) + 1.0_WP
  end do
  call parallel_sum(pdf)

  if (nParticlesGlobal.gt.0) pdf = pdf / real(nParticlesGlobal, WP)

  ! Find position of 1st percentile
  i=1; buf=pdf(1)
  do while (buf.lt.0.01_WP .and. nParticlesGlobal.ne.0)
     i=i+1
     buf=buf+pdf(i)
  end do
  xleft(2) = real(i,WP) / real(nbin,WP) * (maxPart - minPart) + minPart

  ! Find position of 99th percentile
  i=1; buf=pdf(1)
  do while (buf.lt.0.99_WP .and. nParticlesGlobal.ne.0)
     i=i+1
     buf=buf+pdf(i)
  end do
  xright(2) = real(i,WP) / real(nbin,WP) * (maxPart - minPart) + minPart

  ! Get mean pressures
  meanP1 = 0.0_WP
  meanP2 = 0.0_WP
  vol1 = 0.0_WP
  vol2 = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        ! Average pressure at first probe
        do i = max(iStart(1), i1), min(iEnd(1), i1)
           gridIndex = grid_Index(i,j,k)
           meanP1 = meanP1 + pressure(gridIndex, 1) * gridNorm(gridIndex, 1)
           vol1 = vol1 + gridNorm(gridIndex, 1)
        end do
        ! Average pressure at second probe
        do i = max(iStart(1), i2), min(iEnd(1), i2)
           gridIndex = grid_Index(i,j,k)
           meanP2 = meanP2 + pressure(gridIndex, 1) * gridNorm(gridIndex, 1)
           vol2 = vol2 + gridNorm(gridIndex, 1)
        end do
     end do
  end do
  call parallel_sum(vol1)
  call parallel_sum(meanP1); meanP1 = meanP1 / vol1
  call parallel_sum(vol2)
  call parallel_sum(meanP2); meanP2 = meanP2 / vol2

  ! Set the shock tube parameters
  call monitor_select('shock_tube')
  call monitor_set_single_value(1, meanP1)
  call monitor_set_single_value(2, meanP2)
  call monitor_set_single_value(3, xleft(1))
  call monitor_set_single_value(4, xright(1))
  call monitor_set_single_value(5, xleft(2))
  call monitor_set_single_value(6, xright(2))
  call monitor_set_single_value(7, minPart)
  call monitor_set_single_value(8, maxPart)

  ! Cleanup
  deallocate(x, X1D, normSum)

  return
end subroutine monitor_shock_tube_timestep
