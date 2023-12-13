module monitor_impulse

  ! External modules
  use monitor

  implicit none

  ! Global variables
  integer :: jLoc
  real(WP) :: U0i

end module monitor_impulse


! ========================================================= !
! Setup the routine to monitor an impulsively started plate !
! ========================================================= !
subroutine monitor_impulse_setup

  ! Internal modules
  use monitor_impulse

  ! External modules
  use parallel
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use grid

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex, ierror
  real(WP) :: yLoc, minY(2)

  if (trim(simulationName) .ne. 'impulsive plate' .or. isDomainCurvilinear) return

  ! Read in the free-stream velocity
  call parser_read('free stream velocity', U0i)
  U0i = 1.0_WP / U0i

  ! Determine y-value to compute max pressure
  call parser_read('monitor impulse y value', yLoc, 1.0_WP)

  ! Find vertical coordinate associated with yStar
  minY = huge(1.0_WP)
  i = iStart(1); k = iStart(3)
  do j = iStart(2), iEnd(2)
     gridIndex = i - gridOffset(1) + localGridSize(1) *                                      &
          (j - 1 - gridOffset(2) + localGridSize(2) *                                        &
          (k - 1 - gridOffset(3)))
     if (abs(coordinates(gridIndex, 2) - yLoc) .lt. minY(1)) then
        minY(1) = abs(coordinates(gridIndex, 2) - yLoc)
        minY(2) = real(j, WP)
     end if
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE, minY, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC,               &
       comm, ierror)
  jLoc = int(minY(2))

  ! Set the monitor names
  call monitor_create('impulse', 3)
  call monitor_set_header(1, 'delta', 'r')
  call monitor_set_header(2, 'meanP_y', 'r')
  call monitor_set_header(3, 'maxP_y', 'r')

  return
end subroutine monitor_impulse_setup


! ========================== !
! Compute impulse statistics !
! ========================== !
subroutine monitor_impulse_timestep

  ! Internal modules
  use monitor_impulse

  ! External modules
  use parallel
  use solver_options
  use geometry
  use grid
  use state

  implicit none

  ! Local variables
  integer :: i, j, k, gridIndex
  real(WP) :: maxP, meanP, delta, vol

  if (trim(simulationName) .ne. 'impulsive plate' .or. isDomainCurvilinear) return

  ! Compute the displacement thickness
  delta = 0.0_WP
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           delta = delta + (1.0_WP - velocity(gridIndex, 1) * U0i) * gridSpacing(gridIndex, 2)
        end do
     end do
  end do
  call parallel_sum(delta)
  delta = delta / real(globalGridSize(1) * globalGridSize(3), WP)

  ! Compute the max/mean pressure at some distance above the plate
  vol = 0.0_WP
  meanP = 0.0_WP
  maxP=-huge(1.0_WP)
  j = jLoc
  if (j .ge. iStart(2) .and. j .le. iEnd(2)) then
     do k = iStart(3), iEnd(3)
        do i = iStart(1), iEnd(1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           vol = vol + gridNorm(gridIndex, 1)
           meanP = meanP + pressure(gridIndex,1) * gridNorm(gridIndex, 1)
           maxP = max(maxP, pressure(gridIndex, 1))
        end do
     end do
  end if
  call parallel_max(maxP)
  call parallel_sum(meanP)
  call parallel_sum(vol)
  meanP = meanP / vol

  ! Set the impulse parameters
  call monitor_select('impulse')
  call monitor_set_single_value(1, delta)
  call monitor_set_single_value(2, meanP)
  call monitor_set_single_value(3, maxP)

  return
end subroutine monitor_impulse_timestep
