module monitor_drag

  ! External modules
  use monitor
  use functional

  implicit none

  ! Global variables
  logical :: monitorDrag

end module monitor_drag


! ========================================================== !
! Setup the routine to monitor Drag                          !
! ========================================================== !
subroutine monitor_drag_setup

  ! Internal modules
  use monitor_drag

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use grid_functions
  use grid_metrics

  implicit none

  ! Check if we should monitor drag based on periodicity type
  monitorDrag = .false.
  if (nDimensions .gt. 1 .and. isDomainCurvilinear) then
     if (.not.isPeriodic(1) .and. periodicityType(2).eq.OVERLAP) monitorDrag = .true.
  end if

  ! Return if we shouldn't monitor drag
  if (.not. monitorDrag) return

  ! Set the monitor names
  call monitor_create('drag', 4)
  call monitor_set_header(1, 'Area', 'r')
  call monitor_set_header(2, 'Fx', 'r')
  call monitor_set_header(3, 'Fy', 'r')
  call monitor_set_header(4, 'Fz', 'r')

  return
end subroutine monitor_drag_setup

! ================================== !
! Compute Drag statistics            !
! ================================== !
subroutine monitor_drag_timestep
  ! Internal modules
  use monitor_drag

  ! External modules
  use parallel
  use geometry
  use solver_options
  use grid
  use grid_functions
  use state
  use first_derivative

  implicit none

  ! Local variables
  integer :: i, j, k, ii, gridIndex
  real(WP) :: normBoundaryFactor, surfaceArea, F(3)

  ! Return if we shouldn't monitor drag
  if (.not. monitorDrag) return

  normBoundaryFactor = 1.0_WP! / firstDerivative(1)%normBoundary(1)
  surfaceArea = 0.0_WP
  F = 0.0_WP

  if (iRankDir(1) .eq. 0 .and. useViscosity) then
     i = 1
     do k = iStart(3), iEnd(3)
        do j = iStart(2), min(iEnd(2), globalGridSize(2) - 1)
           gridIndex = i - gridOffset(1) + localGridSize(1) *                                &
                (j - 1 - gridOffset(2) + localGridSize(2) *                                  &
                (k - 1 - gridOffset(3)))
           surfaceArea = surfaceArea + sqrt(sum(metrics(gridIndex,1:nDimensions) *           &
                metrics(gridIndex,1:nDimensions)))
           do ii = 1, nDimensions
              ! Viscous stress
              F(ii) = F(ii) + sum(metrics(gridIndex,1:nDimensions) * normBoundaryFactor *    &
                   stressTensor(gridIndex, 1+nDimensions*(ii-1) : ii*nDimensions))
              ! Pressure
              F(ii) = F(ii) - metrics(gridIndex,ii) * normBoundaryFactor *                   &
                   pressure(gridIndex,1)
           end do
        end do
     end do
  end if

  call parallel_sum(surfaceArea)
  do i = 1, nDimensions
     call parallel_sum(F(i))
  end do

  ! Set the drag parameters
  call monitor_select('drag')
  call monitor_set_single_value(1, surfaceArea)
  call monitor_set_single_value(2, F(1))
  call monitor_set_single_value(3, F(2))
  call monitor_set_single_value(4, F(3))

  return
end subroutine monitor_drag_timestep
