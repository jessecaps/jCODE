module monitor_state

  ! External modules
  use monitor

  implicit none

end module monitor_state


! ====================================== !
! Setup the routine to monitor the state !
! ====================================== !
subroutine monitor_state_setup(mode, controlIteration)

  ! Internal modules
  use monitor_state

  ! External modules
  use simulation_flags
  use solver_options
  use geometry

  implicit none

  ! Arguments
  integer, intent(in) :: mode, controlIteration

  ! Local variables
  integer :: i
  character(len = 1) :: vel(3)
  character(len = 2) :: speciesIndex

  select case (mode)

  case (FORWARD)

     if (controlIteration .eq. 0) then

        ! -------------------------------
        ! Monitor the baseline prediction

        call monitor_create('density', 3)
        call monitor_set_header(1, 'min', 'r')
        call monitor_set_header(2, 'max', 'r')
        call monitor_set_header(3, 'mean', 'r')

        call monitor_create('velocity', nDimensions * 3)
        vel =  (/ 'U', 'V', 'W' /)
        do i = 1, nDimensions
           call monitor_set_header((i-1) * 3 + 1, 'min'//vel(i), 'r')
           call monitor_set_header((i-1) * 3 + 2, 'max'//vel(i), 'r')
           call monitor_set_header((i-1) * 3 + 3, 'mean'//vel(i), 'r')
        end do

        call monitor_create('temperature', 3)
        call monitor_set_header(1, 'min', 'r')
        call monitor_set_header(2, 'max', 'r')
        call monitor_set_header(3, 'mean', 'r')

        call monitor_create('pressure', 3)
        call monitor_set_header(1, 'min', 'r')
        call monitor_set_header(2, 'max', 'r')
        call monitor_set_header(3, 'mean', 'r')

        if (nSpecies .gt. 0) then
           call monitor_create('mass_fraction', nSpecies * 3)
           if (allocated(speciesName)) then
              do i = 1, nSpecies
                 call monitor_set_header((i-1) * 3 + 1, 'min_'//trim(speciesName(i)), 'r')
                 call monitor_set_header((i-1) * 3 + 2, 'max_'//trim(speciesName(i)), 'r')
                 call monitor_set_header((i-1) * 3 + 3, 'mean_'//trim(speciesName(i)), 'r')
              end do
           else
              do i = 1, nSpecies
                 write(speciesIndex, '(I2.2)') i 
                 call monitor_set_header((i-1) * 3 + 1, 'minY'//speciesIndex, 'r')
                 call monitor_set_header((i-1) * 3 + 2, 'maxY'//speciesIndex, 'r')
                 call monitor_set_header((i-1) * 3 + 3, 'meanY'//speciesIndex, 'r')
              end do
           end if
        end if

        if (useParticles) then
           call monitor_create('volume_fraction', 3)
           call monitor_set_header(1, 'min', 'r')
           call monitor_set_header(2, 'max', 'r')
           call monitor_set_header(3, 'mean', 'r')
        end if

     else

        ! -------------------------------------
        ! Monitor the forward control iteration

     end if

  case (ADJOINT)

     if (controlIteration .eq. 0) then

        ! -------------------------------------
        ! Monitor the baseline adjoint solution

     else

        ! -------------------------------------
        ! Monitor the adjoint control iteration

     end if

  end select

  return
end subroutine monitor_state_setup


! ============================================================ !
! Compute various statistics of the state during a forward run !
! ============================================================ !
subroutine monitor_state_forward(controlIteration)

  ! Internal modules
  use monitor_state

  ! External modules
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use state

  implicit none

  ! Arguments
  integer, intent(in) :: controlIteration

  ! Local variables
  integer :: i, j
  real(WP) :: minRho, maxRho, meanRho, minU(3), maxU(3), meanU(3),                           &
       minT, maxT, meanT, minP, maxP, meanP, minVF, maxVF, meanVF, Tref, volumeInverse
  real(WP), allocatable, dimension(:) :: minY, maxY, meanY

  if (controlIteration .eq. 0) then

     ! -------------------------------
     ! Monitor the baseline prediction

     ! Prepare the stats
     volumeInverse = 0.0_WP
     minRho = huge(1.0_WP); maxRho = -huge(1.0_WP); meanRho = 0.0_WP
     minU = huge(1.0_WP); maxU = -huge(1.0_WP); meanU = 0.0_WP
     minT = huge(1.0_WP); maxT = -huge(1.0_WP); meanT = 0.0_WP
     minP = huge(1.0_WP); maxP = -huge(1.0_WP); meanP = 0.0_WP
     if (nSpecies .gt. 0) then
        allocate(minY(nSpecies))
        allocate(maxY(nSpecies))
        allocate(meanY(nSpecies))
        minY = huge(1.0_WP); maxY = -huge(1.0_WP); meanY = 0.0_WP
     end if

     ! Compute min, max, and Favre-mean
     do i = 1, nGridPoints
        volumeInverse = volumeInverse + gridNorm(i, 1)
        minRho = min(minRho, conservedVariables(i, 1))
        maxRho = max(maxRho, conservedVariables(i, 1))
        meanRho = meanRho + conservedVariables(i, 1) * gridNorm(i, 1)
        do j = 1, nDimensions
           minU(j) = min(minU(j), velocity(i, j))
           maxU(j) = max(maxU(j), velocity(i, j))
           meanU(j) = meanU(j) + conservedVariables(i, j+1) * gridNorm(i, 1)
        end do
        minT = min(minT, temperature(i, 1))
        maxT = max(maxT, temperature(i, 1))
        meanT = meanT + temperature(i, 1) * gridNorm(i, 1)
        minP = min(minP, pressure(i, 1))
        maxP = max(maxP, pressure(i, 1))
        meanP = meanP + pressure(i, 1) * gridNorm(i, 1)
        do j = 1, nSpecies
           minY(j) = min(minY(j), massFraction(i, j))
           maxY(j) = max(maxY(j), massFraction(i, j))
           meanY(j) = meanY(j) + conservedVariables(i, nDimensions+2+j) * gridNorm(i, 1)
        end do
     end do
     call parallel_sum(volumeInverse)
     volumeInverse = 1.0_WP / volumeInverse
     call parallel_min(minRho)
     call parallel_max(maxRho)
     call parallel_sum(meanRho)
     meanRho = meanRho * volumeInverse
     do i = 1, nDimensions
        call parallel_min(minU(i))
        call parallel_max(maxU(i))
        call parallel_sum(meanU(i))
        meanU(i) = meanU(i) * volumeInverse / meanRho
     end do
     call parallel_min(minT)
     call parallel_max(maxT)
     call parallel_sum(meanT)
     meanT = meanT * volumeInverse
     call parallel_min(minP)
     call parallel_max(maxP)
     call parallel_sum(meanP)
     meanP = meanP * volumeInverse
     do i = 1, nSpecies
        call parallel_min(minY(i))
        call parallel_max(maxY(i))
        call parallel_sum(meanY(i))
        meanY(i) = meanY(i) * volumeInverse / meanRho
     end do

     ! Dimensionalize the temperature
     Tref = (ratioOfSpecificHeats - 1.0_WP) * 293.15_WP
     minT = minT * Tref
     maxT = maxT * Tref
     meanT = meanT * Tref

     ! Set the density values
     call monitor_select('density')
     call monitor_set_single_value(1, minRho)
     call monitor_set_single_value(2, maxRho)
     call monitor_set_single_value(3, meanRho)

     ! Set the velocity values
     call monitor_select('velocity')
     do i = 1, nDimensions
        call monitor_set_single_value((i-1) * 3 + 1, minU(i))
        call monitor_set_single_value((i-1) * 3 + 2, maxU(i))
        call monitor_set_single_value((i-1) * 3 + 3, meanU(i))
     end do

     ! Set the temperature values
     call monitor_select('temperature')
     call monitor_set_single_value(1, minT)
     call monitor_set_single_value(2, maxT)
     call monitor_set_single_value(3, meanT)

     ! Set the pressure values
     call monitor_select('pressure')
     call monitor_set_single_value(1, minP)
     call monitor_set_single_value(2, maxP)
     call monitor_set_single_value(3, meanP)

     ! Set the mass fraction values
     if (nSpecies .gt. 0) then
        call monitor_select('mass_fraction')
        do i = 1, nSpecies
           call monitor_set_single_value((i-1) * 3 + 1, minY(i))
           call monitor_set_single_value((i-1) * 3 + 2, maxY(i))
           call monitor_set_single_value((i-1) * 3 + 3, meanY(i))
        end do
     end if

     ! Set the volume fraction values
     if (useParticles) then
        minVF = huge(1.0_WP); maxVF = -huge(1.0_WP); meanVF = 0.0_WP
        do i = 1, nGridPoints
           minVF = min(minVF, volumeFraction(i, 1))
           maxVF = max(maxVF, volumeFraction(i, 1))
           meanVF = meanVF + volumeFraction(i, 1) * gridNorm(i, 1)
        end do
        call parallel_min(minVF)
        call parallel_max(maxVF)
        call parallel_sum(meanVF)
        meanVF = meanVF * volumeInverse

        call monitor_select('volume_fraction')
        call monitor_set_single_value(1, minVF)
        call monitor_set_single_value(2, maxVF)
        call monitor_set_single_value(3, meanVF)
     end if

  else

     ! --------------------------------------
     ! Monitor the forward control iteration


  end if

  return
end subroutine monitor_state_forward
