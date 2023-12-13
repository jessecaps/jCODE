module ignition_source

  ! External modules
  use precision

  implicit none

  integer :: nModes
  real(WP) :: location(3), radius(3), energy, timeStart, timeDuration
  real(WP), dimension(:), allocatable :: amplitudes, frequencies, phases
  logical :: useIgnition

end module ignition_source


! ========================== !
! Setup the ignition sources !
! ========================== !
subroutine ignition_source_setup

  ! Internal modules
  use ignition_source

  ! External modules
  use parser
  use string
  use math, only : twoPi
  use geometry

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: rand, wmax, wmin, geometricGrowthFactor

  call parser_read('use ignition source', useIgnition, .false.)

  if (.not. useIgnition) return

  ! Read ignition parameters
  call parser_read('ignition source x', location(1), 0.0_WP)
  call parser_read('ignition source y', location(2), 0.0_WP)
  call parser_read('ignition source z', location(3), 0.0_WP)
  call parser_read('ignition source radius x', radius(1))
  call parser_read('ignition source radius y', radius(2), radius(1))
  call parser_read('ignition source radius z', radius(3), radius(1))
  call parser_read('ignition source energy',  energy, 0.0_WP)
  call parser_read('ignition source time start',  timeStart, 0.0_WP)
  call parser_read('ignition source time duration',  timeDuration, 0.0_WP)
  call parser_read('ignition source time modes', nModes, 0)
  if (nModes .gt. 0) then
     if (timeDuration .gt. 0.0_WP)                                                           &
          call die('ignition_source_setup: cannot specify nModes > 0 and timeDuration > 0!')
     call parser_read('ignition source max frequency', wmax)
     call parser_read('ignition source min frequency', wmin)
     geometricGrowthFactor = 1.0_WP / (10.0_WP**0.25_WP)
     allocate(amplitudes(nModes))
     allocate(frequencies(nModes))
     allocate(phases(nModes))
     do i = 1, nModes
        amplitudes(i) = 1.0_WP
        !frequencies(i) = wmax * geometricGrowthFactor ** real(i - 1, WP)
        frequencies(i) = wmin + real(i-1,WP) * (wmax-wmin) / real(nModes-1,WP)
        call random_number(rand)
        phases(i) = twoPi * rand
     end do
  end if

  return
end subroutine ignition_source_setup


! ============================ !
! Cleanup the ignition sources !
! ============================ !
subroutine ignition_source_cleanup

  ! Internal modules
  use ignition_source

  implicit none

  if (allocated(amplitudes)) deallocate(amplitudes)
  if (allocated(frequencies)) deallocate(frequencies)
  if (allocated(phases)) deallocate(phases)

  return
end subroutine ignition_source_cleanup


! =============================================== !
! Add the ignition sources during the forward run !
! =============================================== !
subroutine ignition_source_forward(time, source)

  ! Internal modules
  use ignition_source

  ! External modules
  use math
  use solver_options
  use geometry
  use grid

  implicit none

  ! Arguments
  real(WP), intent(in) :: time
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j
  real(WP) :: power, timePortion, gaussianFactor(3), distance(3)

  if (.not. useIgnition) return

  ! Get power
  power = energy
  if (minval(radius(1:nDimensions)) .gt. 0.0_WP) power = power /                             &
       sqrt(2.0_WP * pi) ** nDimensions / product(radius(1:nDimensions))

  ! Get time portion
  if (timeDuration .gt. 0.0_WP) then
     timePortion = exp( -0.5_WP * (time - timeStart)**2 / timeDuration**2)
     power = power / timeDuration / sqrt(2.0_WP * pi)
  else
     timePortion = 1.0_WP
  end if

  ! Laser pulse
  if (nModes .gt. 0) then
     timePortion = 0.0_WP
     do i = 1, nModes
        timePortion = timePortion + 0.5_WP * amplitudes(i) / real(nModes, WP) *              &
             (1.0_WP + cos(frequencies(i) * time + phases(i)))
     end do
  end if

  ! Prepare spatial extent
  gaussianFactor = 0.0_WP
  if (minval(radius(1:nDimensions)) .gt. 0.0_WP)                                             &
       gaussianFactor(1:nDimensions) = 0.5_WP / radius(1:nDimensions)**2
  distance = 0.0_WP

  do i = 1, nGridPoints
     ! Get minimal distance based on periodicity
     distance(1:nDimensions) =  coordinates(i,:) - location(1:nDimensions)
     do j = 1, nDimensions
        if (isPeriodic(j)) then
           if (abs(coordinates(i,j) - periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) - periodicLength(j)         &
                - location(j)
           if (abs(coordinates(i,j) + periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) + periodicLength(j)         &
                - location(j)
        end if
     end do

     ! Add the ignition source
     source(i, nDimensions+2) = source(i, nDimensions+2) +                                   &
          power * timePortion * exp(- sum(gaussianFactor(1:nDimensions) *                    &
          distance(1:nDimensions)**2) )
  end do

  return
end subroutine ignition_source_forward
