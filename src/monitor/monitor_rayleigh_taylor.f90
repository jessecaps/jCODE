module monitor_rayleigh_taylor

  ! External modules
  use monitor

  implicit none

end module monitor_rayleigh_taylor


! ========================================================== !
! Setup the routine to monitor a Rayleigh-Taylor instability !
! ========================================================== !
subroutine monitor_rayleigh_taylor_setup

  ! Internal modules
  use monitor_rayleigh_taylor

  ! External modules
  use parser
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  if (trim(simulationName) .ne. 'rayleigh taylor' .or. nSpecies .eq. 0) return

  ! Set the monitor names
  call monitor_create('rayleigh_taylor', 2)
  call monitor_set_header(1, 'height', 'r')
  call monitor_set_header(2, 'amplitude', 'r')

  return
end subroutine monitor_rayleigh_taylor_setup


! ================================== !
! Compute Rayleigh-Taylor statistics !
! ================================== !
subroutine monitor_rayleigh_taylor_timestep

  ! Internal modules
  use monitor_rayleigh_taylor

  ! External modules
  use parallel
  use solver_options
  use grid_functions
  use state

  implicit none

  ! Local variables
  integer :: i, j, k, ny, jmin, ierror
  real(WP) :: minValues(2), Ymix, ymax, ymin, yInterp, midPlane, eta, height
  real(WP), dimension(:), allocatable :: y, Y1D

  if (trim(simulationName) .ne. 'rayleigh taylor' .or. nSpecies .eq. 0) return

  ! Set up vectors
  ny = globalGridSize(2)
  allocate(y(ny)); y = 0.0_WP
  allocate(Y1D(ny))

  ! Get the vertical coordinate
  i = iStart(1); k = iStart(3)
  do j = iStart(2), iEnd(2)  
     y(j) = coordinates(grid_index(i,j,k), 2)
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE, y, ny, MPI_REAL_WP, MPI_SUM, commDir(2), ierror)

  ymin = huge(1.0_WP)
  ymax = -huge(1.0_WP)
  do k = iStart(3), iEnd(3)
     do i = iStart(1), iEnd(1)

        minValues(1) = huge(1.0_WP)

        ! Find interface at this location (i,k)
        Y1D = 0.0_WP
        do j = iStart(2), iEnd(2)
           Ymix = abs(massFraction(grid_index(i,j,k), 1) - 0.5_WP)
           if (Ymix .lt. minValues(1)) then
              minValues(1) = Ymix
              minValues(2) = real(j, WP)
           end if

           ! Store mass fraction
           Y1D(j) = massFraction(grid_index(i,j,k), 1)

        end do
        call MPI_ALLREDUCE(MPI_IN_PLACE, minValues, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC,    &
             commDir(2), ierror)
        call MPI_ALLREDUCE(MPI_IN_PLACE, Y1D, ny, MPI_REAL_WP, MPI_SUM, commDir(2), ierror)

        ! Interpolate the vertical coordinate
        jmin = int(minValues(2))
        jmin = max(jmin, 2); jmin = min(jmin, ny-1)
        if (Y1D(jmin) .gt. 0.5_WP) then
           yInterp = (y(jmin) - y(jmin-1)) / (Y1D(jmin) - Y1D(jmin-1)) *                     &
                (0.5_WP - Y1D(jmin)) + y(jmin)
        else if (Y1D(jmin) .lt. 0.5_WP) then
           yInterp = (y(jmin) - y(jmin+1)) / (Y1D(jmin) - Y1D(jmin+1)) *                     &
                (0.5_WP - Y1D(jmin)) + y(jmin)
        else
           yInterp = y(jmin);
        end if
        ymin = min(ymin, yInterp)
        ymax = max(ymax, yInterp)
     end do
  end do

  ! Get global min/max
  call parallel_min(ymin)
  call parallel_max(ymax)

  ! Get the mid plane and mixin height
  height = ymax - ymin
  midPlane = ymin + 0.5_WP * height

  ! Get the amplitude at x=0
  eta = 0.0_WP
  i = ceiling(0.5_WP * real(globalGridSize(1), WP))
  k = ceiling(0.5_WP * real(globalGridSize(3), WP) )
  if (i.ge.iStart(1) .and. i.le.iEnd(1) .and. k.ge.iStart(3) .and. k.le.iEnd(3)) then

     minValues(1) = huge(1.0_WP)

     ! Find interface at this location (i,k)
     Y1D = 0.0_WP
     do j = iStart(2), iEnd(2)
        Ymix = abs(massFraction(grid_index(i,j,k), 1) - 0.5_WP)
        if (Ymix .lt. minValues(1)) then
           minValues(1) = Ymix
           minValues(2) = real(j, WP)
        end if

        ! Store mass fraction
        Y1D(j) = massFraction(grid_index(i,j,k), 1)

     end do
     call MPI_ALLREDUCE(MPI_IN_PLACE, minValues, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC,    &
          commDir(2), ierror)
     call MPI_ALLREDUCE(MPI_IN_PLACE, Y1D, ny, MPI_REAL_WP, MPI_SUM, commDir(2), ierror)

     ! Interpolate the vertical coordinate
     jmin = int(minValues(2))
     jmin = max(jmin, 2); jmin = min(jmin, ny-2)
     if (Y1D(jmin) .gt. 0.5_WP) then
        yInterp = (y(jmin) - y(jmin-1)) / (Y1D(jmin) - Y1D(jmin-1)) *                     &
             (0.5_WP - Y1D(jmin)) + y(jmin)
     else if (Y1D(jmin) .lt. 0.5_WP) then
        yInterp = (y(jmin) - y(jmin+1)) / (Y1D(jmin) - Y1D(jmin+1)) *                     &
             (0.5_WP - Y1D(jmin)) + y(jmin)
     else
        yInterp = y(jmin);
     end if
     eta = abs(yInterp - midPlane)
  end if
  call parallel_max(eta)

  ! Set the Rayleigh-Taylor parameters
  call monitor_select('rayleigh_taylor')
  call monitor_set_single_value(1, height)
  call monitor_set_single_value(2, eta)

  ! Cleanup
  deallocate(y, Y1D)

  return
end subroutine monitor_rayleigh_taylor_timestep
