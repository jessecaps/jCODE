module ignition_controller

  ! External modules
  use controller
  use precision

  implicit none

end module ignition_controller


! =========================== !
! Setup the ignition actuator !
! =========================== !
subroutine ignition_controller_setup

  ! Internal modules
  use ignition_controller

  ! External modules
  use string
  use parser
  use simulation_flags
  use ignition_source

  implicit none

  ! Local variables
  integer :: i, j
  character(len = str_medium) :: name

  if (predictionOnly) return

  ! Ignition actuator requires ignition source to be used
  if (.not. useIgnition)                                                                     &
       call die('ignition_controller_setup: ignition actuator requires an ignition source!')

  ! Initialize the baseline ignition parameters
  allocate(baselineValue(nControlParameters))
  baselineValue = 0.0_WP
  j = 1
  do i = 1, nControlParameters

     select case (trim(sensitivityParameter(i)))

     case ('energy')
        baselineValue(i) = energy

     case ('x0')
        baselineValue(i) = location(1)

     case ('y0')
        baselineValue(i) = location(2)

     case ('z0')
        baselineValue(i) = location(3)

     case ('rx')
        baselineValue(i) = radius(1)

     case ('ry')
        baselineValue(i) = radius(2)

     case ('rz')
        baselineValue(i) = radius(3)

     case ('t0')
        baselineValue(i) = timeStart

     case ('duration')
        baselineValue(i) = timeDuration

     case default

        if (j.le.nModes) then
           ! Set parameter name
           write(name, '(A,I3.3)') "a_", j
           sensitivityParameter(i) = trim(name)

           ! Set the baseline amplitude
           baselineValue(i) = amplitudes(j)
           j = j + 1
        else if (j.gt.nModes .and. j.le.2*nModes) then
           ! Set parameter name
           write(name, '(A,I3.3)') "theta_", j-nModes
           sensitivityParameter(i) = trim(name)

           ! Set the baseline amplitude
           baselineValue(i) = phases(j-nModes)
           j = j + 1
        else
           call die("ignition_controller_setup: unknown sensitivity parameter '" //          &
                trim(sensitivityParameter(i)) // "'!")
        end if

     end select

  end do

  return
end subroutine ignition_controller_setup


! ============================= !
! Cleanup the ignition actuator !
! ============================= !
subroutine ignition_controller_cleanup

  ! Internal modules
  use ignition_controller

  implicit none

  return
end subroutine ignition_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine ignition_controller_update_parameters(actuationAmount)

  ! Internal modules
  use ignition_controller

  ! External modules
  use geometry
  use solver_options
  use time_info
  use ignition_source

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Local variables
  integer :: i, j
  real(WP), parameter :: smallNumber = 1.0E-9_WP

  j = 1
  do i = 1, nControlParameters
     select case (trim(sensitivityParameter(i)))

     case ('energy')
        energy = max( smallNumber, baselineValue(i) +                                        &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case ('x0')
        location(1) = baselineValue(i) +                                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i)

     case ('y0')
        location(2) = baselineValue(i) +                                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i)

     case ('z0')
        location(3) = baselineValue(i) +                                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i)

     case ('rx')
        radius(1) = max( smallNumber, baselineValue(i) +                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case ('ry')
        radius(2) = max( smallNumber, baselineValue(i) +                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case ('rz')
        radius(3) = max( smallNumber, baselineValue(i) +                                     &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case ('t0')
        timeStart = max( 0.0_WP, baselineValue(i) +                                          &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case ('duration')
        timeDuration = max( smallNumber, baselineValue(i) +                                  &
             real(gradientDirection, WP) * actuationAmount * controlGradient(i) )

     case default

        if (j.le.nModes) then
           amplitudes(j) = baselineValue(i) +                                                &
                real(gradientDirection, WP) * actuationAmount * controlGradient(i)
           j = j + 1
        else if (j.gt.nModes .and. j.le.2*nModes) then
           phases(j-nModes) = baselineValue(i) +                                             &
                real(gradientDirection, WP) * actuationAmount * controlGradient(i)
           j = j + 1
        end if

     end select
  end do

  return
end subroutine ignition_controller_update_parameters


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine ignition_controller_compute_sensitivity

  ! Internal modules
  use ignition_controller

  ! External modules
  use geometry
  use solver_options
  use grid, only : coordinates
  use grid_functions, only : inner_product
  use state, only : adjointVariables
  use time_info
  use ignition_source

  implicit none

  ! Local variables
  integer :: i, j
  real(WP) ::  currentTime, dimFactor, timePortion, temp
  real(WP), dimension(:,:), allocatable :: F, ignitionSource, distance

  ! Get relevant quantities
  currentTime = adjointCoefficientTime
  dimFactor = 1.0_WP
  if (nDimensions .eq. 3) dimFactor = 2.0_WP

  ! Allocate temporary arrays
  allocate(F(nGridPoints, 2))
  allocate(ignitionSource(nGridPoints, nUnknowns))
  allocate(distance(nGridPoints, 3))

  ! Get the minimal distance to the ignition source (account for periodicity)
  distance = 0.0_WP
  do i = 1, nGridPoints
     distance(i,1:nDimensions) = coordinates(i,:) - location(1:nDimensions)
     do j = 1, nDimensions
        if (isPeriodic(j)) then
           if (abs(coordinates(i,j) - periodicLength(j) - location(j)) .lt.                  &
                abs(distance(i,j))) distance(i,j) = coordinates(i,j) - periodicLength(j)     &
                - location(j)
           if (abs(coordinates(i,j) + periodicLength(j) - location(j)) .lt.                  &
                abs(distance(i,j))) distance(i,j) = coordinates(i,j) + periodicLength(j)     &
                - location(j)
        end if
     end do
  end do

  ! Recompute time portion if using laser pulse
  if (nModes .gt. 0) then
     timePortion = 0.0_WP
     do i = 1, nModes
        timePortion = timePortion + 0.5_WP * amplitudes(i) / real(nModes, WP) *              &
             (1.0_WP + cos(frequencies(i) * currentTime + phases(i)))
     end do
  end if

  ! Compute the source terms
  ignitionSource = 0.0_WP
  call ignition_source_forward(currentTime, ignitionSource)

  ! Compute sensitivities of ignition parameters
  instantaneousSensitivity = 0.0_WP
  F(:,1) = adjointVariables(:,nDimensions+2)

  j = 1
  do i = 1, nControlParameters
     select case (trim(sensitivityParameter(i)))

     case ('energy')
        F(:,2) = ignitionSource(:,nDimensions+2) / energy

     case ('x0')
        F(:,2) = ignitionSource(:,nDimensions+2) * distance(:,1) / radius(1)**2

     case ('y0')
        F(:,2) = ignitionSource(:,nDimensions+2) * distance(:,2) / radius(2)**2

     case ('z0')
        F(:,2) = ignitionSource(:,nDimensions+2) * distance(:,3) / radius(3)**2

     case ('rx')
        F(:,2) = ignitionSource(:,nDimensions+2) / radius(1)**3 *                            &
             (distance(:,1)**2 - radius(1)**2)

     case ('ry')
        F(:,2) = ignitionSource(:,nDimensions+2) / radius(2)**3 *                            &
             (distance(:,2)**2 - radius(2)**2)

     case ('rz')
        F(:,2) = ignitionSource(:,nDimensions+2) / radius(3)**3 *                            &
             (distance(:,3)**2 - radius(3)**2)

     case ('t0')
        F(:,2) = ignitionSource(:,nDimensions+2) * (currentTime - timeStart) / timeDuration**2

     case ('duration')
        F(:,2) = ignitionSource(:,nDimensions+2) / timeDuration**3 *                         &
             ((currentTime - timeStart) ** 2 - timeDuration**2)

     case default

        if (j.le.nModes) then

           ! Amplitude sensitivity
           temp = 0.5_WP / real(nModes, WP) *                                                &
                (1.0_WP + cos(frequencies(j) * currentTime + phases(j)))

           F(:,2) = ignitionSource(:,nDimensions+2) / timePortion * temp

           j = j + 1
        else if (j.gt.nModes .and. j.le.2*nModes) then

           ! Phase sensitivity
           temp = - 0.5_WP * amplitudes(j-nModes) / real(nModes, WP) *                      &
                sin(frequencies(j-nModes) * currentTime + phases(j-nModes))

           F(:,2) = ignitionSource(:,nDimensions+2) / timePortion * temp

           j = j + 1
        end if

     end select

     instantaneousSensitivity(i) = instantaneousSensitivity(i) + inner_product(F(:,1), F(:,2))

  end do

  ! Clean up
  deallocate(F)
  deallocate(ignitionSource)
  deallocate(distance)

return
end subroutine ignition_controller_compute_sensitivity
