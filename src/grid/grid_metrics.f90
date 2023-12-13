module grid_metrics

  ! External modules
  use grid

  implicit none

end module grid_metrics


! ====================== !
! Setup the grid metrics !
! ====================== !
subroutine grid_metrics_setup

  ! Internal modules
  use grid_metrics

  ! External modules
  use string
  use parallel
  use simulation_flags
  use geometry
  use first_derivative
  use grid_functions

  implicit none

  ! Local variables
  integer :: i, j, k
  real(WP), allocatable :: jacobianMatrixInverse(:,:), F(:,:)
  logical :: hasNegativeJacobian, hasNegativeJacobian_
  character(len = str_medium) :: message
  real(WP) :: jacobianOutsideRange

  allocate(jacobianMatrixInverse(nGridPoints, nDimensions ** 2))
  do i = 1, nDimensions
     call compute_coordinate_derivative(i,                                                   &
          jacobianMatrixInverse(:,(i-1)*nDimensions+1:i*nDimensions))
  end do

  select case (nDimensions)

  case (1)
     jacobian(:,1) = jacobianMatrixInverse(:,1)
     metrics(:,1) = 1.0_WP
     gridSpacing(:,1) = abs(jacobianMatrixInverse(:,1))
     if (allocated(arcLengths)) arcLengths(:,1) = abs(metrics(:,1))

  case (2)

     if (isDomainCurvilinear) then
        jacobian(:,1) = (jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,4) -           &
             jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,3))
        metrics(:,1) = jacobianMatrixInverse(:,4)
        metrics(:,2) = - jacobianMatrixInverse(:,3)
        metrics(:,3) = - jacobianMatrixInverse(:,2)
        metrics(:,4) = jacobianMatrixInverse(:,1)
        gridSpacing(:,1) = abs(jacobianMatrixInverse(:,1) + jacobianMatrixInverse(:,3))
        gridSpacing(:,2) = abs(jacobianMatrixInverse(:,2) + jacobianMatrixInverse(:,4))
        if (allocated(arcLengths)) then
           arcLengths(:,1) = sqrt(metrics(:,1) ** 2 + metrics(:,2) ** 2)
           arcLengths(:,2) = sqrt(metrics(:,3) ** 2 + metrics(:,4) ** 2)
        end if
     else
        jacobian(:,1) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,4)
        metrics(:,1) = jacobianMatrixInverse(:,4)
        metrics(:,2) = 0.0_WP
        metrics(:,3) = 0.0_WP
        metrics(:,4) = jacobianMatrixInverse(:,1)
        gridSpacing(:,1) = abs(jacobianMatrixInverse(:,1))
        gridSpacing(:,2) = abs(jacobianMatrixInverse(:,4))
        if (allocated(arcLengths)) then
           arcLengths(:,1) = abs(metrics(:,1))
           arcLengths(:,2) = abs(metrics(:,4))
        end if
     end if

  case (3)

     if (isDomainCurvilinear) then
        jacobian(:,1) =                                                                      &
             jacobianMatrixInverse(:,1) *                                                    &
             (jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9) -                      &
             jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,6)) +                      &
             jacobianMatrixInverse(:,4) *                                                    &
             (jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,3) -                      &
             jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,9)) +                      &
             jacobianMatrixInverse(:,7) *                                                    &
             (jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,6) -                      &
             jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,3))
        gridSpacing(:,1) = abs(jacobianMatrixInverse(:,1) + jacobianMatrixInverse(:,4) +     &
             jacobianMatrixInverse(:,7))
        gridSpacing(:,2) = abs(jacobianMatrixInverse(:,2) + jacobianMatrixInverse(:,5) +     &
             jacobianMatrixInverse(:,8))
        gridSpacing(:,3) = abs(jacobianMatrixInverse(:,3) + jacobianMatrixInverse(:,6) +     &
             jacobianMatrixInverse(:,9))
     else
        jacobian(:,1) = jacobianMatrixInverse(:,1) *                                         &
             jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9)
        gridSpacing(:,1) = abs(jacobianMatrixInverse(:,1))
        gridSpacing(:,2) = abs(jacobianMatrixInverse(:,5))
        gridSpacing(:,3) = abs(jacobianMatrixInverse(:,9))
     end if

     if (any(periodicityType .eq. PLANE)) then

        if (isDomainCurvilinear) then
           metrics(:,1) = jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9) -          &
                jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,6)
           metrics(:,2) = jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,6) -          &
                jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,9)
           metrics(:,3) = jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,8) -          &
                jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,5)
           metrics(:,4) = jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,3) -          &
                jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,9)
           metrics(:,5) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,9) -          &
                jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,3)
           metrics(:,6) = jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,2) -          &
                jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,8)
           metrics(:,7) = jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,6) -          &
                jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,3)
           metrics(:,8) = jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,3) -          &
                jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,6)
           metrics(:,9) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,5) -          &
                jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,2)
        else
           metrics(:,1) = jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9)
           metrics(:,2) = 0.0_WP
           metrics(:,3) = 0.0_WP
           metrics(:,4) = 0.0_WP
           metrics(:,5) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,9)
           metrics(:,6) = 0.0_WP
           metrics(:,7) = 0.0_WP
           metrics(:,8) = 0.0_WP
           metrics(:,9) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,5)
        end if

     else

        allocate(F(nGridPoints, 1))

        F(:,1) = jacobianMatrixInverse(:,5) * coordinates(:,3)
        call first_derivative_apply(3, F)
        metrics(:,1) = F(:,1)
        if (isDomainCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,8) * coordinates(:,3)
           call first_derivative_apply(2, F)
           metrics(:,1) = metrics(:,1) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,2) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,6) * coordinates(:,1)
           call first_derivative_apply(3, F)
           metrics(:,2) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,9) * coordinates(:,1)
           call first_derivative_apply(2, F)
           metrics(:,2) = metrics(:,2) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,3) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,4) * coordinates(:,2)
           call first_derivative_apply(3, F)
           metrics(:,3) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,7) * coordinates(:,2)
           call first_derivative_apply(2, F)
           metrics(:,3) = metrics(:,3) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,4) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,8) * coordinates(:,3)
           call first_derivative_apply(1, F)
           metrics(:,4) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,2) * coordinates(:,3)
           call first_derivative_apply(3, F)
           metrics(:,4) = metrics(:,4) - F(:,1)
        end if

        F(:,1) = jacobianMatrixInverse(:,9) * coordinates(:,1)
        call first_derivative_apply(1, F)
        metrics(:,5) = F(:,1)
        if (isDomainCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,3) * coordinates(:,1)
           call first_derivative_apply(3, F)
           metrics(:,5) = metrics(:,5) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,6) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,7) * coordinates(:,2)
           call first_derivative_apply(1, F)
           metrics(:,6) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,1) * coordinates(:,2)
           call first_derivative_apply(3, F)
           metrics(:,6) = metrics(:,6) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,7) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,2) * coordinates(:,3)
           call first_derivative_apply(2, F)
           metrics(:,7) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,5) * coordinates(:,3)
           call first_derivative_apply(1, F)
           metrics(:,7) = metrics(:,7) - F(:,1)
        end if

        if (.not. isDomainCurvilinear) then
           metrics(:,8) = 0.0_WP
        else
           F(:,1) = jacobianMatrixInverse(:,3) * coordinates(:,1)
           call first_derivative_apply(2, F)
           metrics(:,8) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,6) * coordinates(:,1)
           call first_derivative_apply(1, F)
           metrics(:,8) = metrics(:,8) - F(:,1)
        end if

        F(:,1) = jacobianMatrixInverse(:,1) * coordinates(:,2)
        call first_derivative_apply(2, F)
        metrics(:,9) = F(:,1)
        if (isDomainCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,4) * coordinates(:,2)
           call first_derivative_apply(1, F)
           metrics(:,9) = metrics(:,9) - F(:,1)
        end if

        deallocate(F)

     end if

     if (allocated(arcLengths)) then
        if (isDomainCurvilinear) then
           arcLengths(:,1) = sqrt(metrics(:,1) ** 2 + metrics(:,2) ** 2 + metrics(:,3) ** 2)
           arcLengths(:,2) = sqrt(metrics(:,4) ** 2 + metrics(:,5) ** 2 + metrics(:,6) ** 2)
           arcLengths(:,3) = sqrt(metrics(:,7) ** 2 + metrics(:,8) ** 2 + metrics(:,9) ** 2)
        else
           arcLengths(:,1) = abs(metrics(:,1))
           arcLengths(:,2) = abs(metrics(:,5))
           arcLengths(:,3) = abs(metrics(:,9))
        end if
     end if

  end select

  deallocate(jacobianMatrixInverse)

  hasNegativeJacobian_ = .not. isVariableWithinRange(jacobian(:,1), jacobianOutsideRange,    &
       i, j, k, minValue = 0.0_WP)

  call parallel_lor(hasNegativeJacobian_, hasNegativeJacobian)

  if (hasNegativeJacobian .and. iRank .eq. iRoot) then
     write(message, '(A,3(I0.0,A),(ES11.4,A))') "Jacobian at (", i, ", ", j, ", ", k, "): ", &
          jacobianOutsideRange, " is non-positive!"
     call die(trim(message))
  end if

  ! Update the norm matrix, using the norm from the first derivative operator
  gridNorm = 1.0_WP
  do i = 1, nDimensions
     call first_derivative_apply_norm(i, gridNorm)
  end do
  gridNorm = gridNorm * jacobian

  jacobian = 1.0_WP / jacobian

  ! Determine minimum grid spacing
  minGridSpacing = minval(gridSpacing)
  call parallel_min(minGridSpacing)

  ! Determine min/max domain extent
  domainExtent = 0.0_WP
  do i = 1, nDimensions
     domainExtent(i,1) = minval(coordinates(:,i))
     call parallel_min(domainExtent(i,1))
     domainExtent(i,2) = maxval(coordinates(:,i))
     call parallel_max(domainExtent(i,2))
  end do

  ! Compute the grid volume
  call compute_grid_volume(globalGridVolume)

  return
end subroutine grid_metrics_setup
