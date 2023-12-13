module grid

  ! External modules
  use precision
  use geometry

  implicit none

  integer, allocatable :: iblank(:)
  real(WP), dimension(:,:), allocatable :: coordinates, jacobian, metrics, gridNorm,         &
       gridSpacing, arcLengths

contains

  ! Allocate grid data
  ! ------------------
  subroutine allocate_grid

    ! External modules
    use simulation_flags
    use dissipation, only : compositeDissipation

    implicit none

    allocate(iblank(nGridPoints)); iblank = 1
    allocate(coordinates(nGridPoints, nDimensions))
    allocate(jacobian(nGridPoints, 1))
    allocate(metrics(nGridPoints, nDimensions ** 2))
    allocate(gridNorm(nGridPoints, 1)); gridNorm = 1.0_WP
    allocate(gridSpacing(nGridPoints, nDimensions))

    if (.not. compositeDissipation) allocate(arcLengths(nGridPoints, nDimensions))

    return
  end subroutine allocate_grid


  ! Make unit cube (computational coordinates)
  ! ------------------------------------------
  subroutine make_unit_cube

    implicit none

    ! Local variables
    integer :: i, j, k, l
    real(wp) :: h
    real(wp), allocatable :: unitInterval(:)

    do l = 1, nDimensions

       if (allocated(unitInterval)) deallocate(unitInterval)
       allocate(unitInterval(globalGridSize(l)))

       if (periodicityType(l).eq.PLANE .or. periodicityType(l).eq.POLAR) then
          h = 1.0_WP / real(globalGridSize(l), WP)
       else
          h = 1.0_WP / real(globalGridSize(l) - 1, WP)
       end if

       do i = 1, globalGridSize(l)
          unitInterval(i) = real(i - 1, WP) * h
       end do

       select case (l)

       case (1)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   coordinates(i + localGridSize(1) * (j - 1 +                               &
                        localGridSize(2) * (k - 1)), l) = unitInterval(i + gridOffset(l))
                end do
             end do
          end do

       case (2)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   coordinates(i + localGridSize(1) * (j - 1 +                               &
                        localGridSize(2) * (k - 1)), l) = unitInterval(j + gridOffset(l))
                end do
             end do
          end do

       case (3)
          do k = 1, localGridSize(3)
             do j = 1, localGridSize(2)
                do i = 1, localGridSize(1)
                   coordinates(i + localGridSize(1) * (j - 1 +                               &
                        localGridSize(2) * (k - 1)), l) = unitInterval(k + gridOffset(l))
                end do
             end do
          end do

       end select

       if (allocated(unitInterval)) deallocate(unitInterval)

    end do

    return
  end subroutine make_unit_cube

end module grid


! ============== !
! Setup the grid !
! ============== !
subroutine grid_setup

  ! Internal modules
  use grid

  implicit none

  call allocate_grid
  call make_unit_cube

  return
end subroutine grid_setup


! ================ !
! Cleanup the grid !
! ================ !
subroutine grid_cleanup

  ! Internal modules
  use grid

  implicit none

  call grid_levelset_cleanup

  if (allocated(iblank)) deallocate(iblank)
  if (allocated(coordinates)) deallocate(coordinates)
  if (allocated(jacobian)) deallocate(jacobian)
  if (allocated(metrics)) deallocate(metrics)
  if (allocated(gridNorm)) deallocate(gridNorm)
  if (allocated(gridSpacing)) deallocate(gridSpacing)
  if (allocated(arcLengths)) deallocate(arcLengths)

  return
end subroutine grid_cleanup
