module grid_levelset

  ! External modules
  use grid
  use precision
  use geometry

  implicit none

  real(WP), dimension(:,:), allocatable :: levelset, levelsetNormal, indicatorFunction,      &
       primitiveGridNorm
  real(WP) :: globalVolumeFraction
  logical useLevelset

contains

  ! Allocate levelset arrays
  ! ------------------------
  subroutine allocate_levelset

    ! Allocate levelset arrays
    allocate(levelset(nGridPoints, 1))
    allocate(levelsetNormal(nGridPoints, nDimensions))
    allocate(indicatorFunction(nGridPoints, 1))

    ! Store the unmodified grid norm
    allocate(primitiveGridNorm(nGridPoints, 1))
    primitiveGridNorm = gridNorm

    return
  end subroutine allocate_levelset

end module grid_levelset


! ====================== !
! Setup levelset routine !
! ====================== !
subroutine grid_levelset_setup

  ! Internal modules
  use grid_levelset

  ! External modules
  use parser
  use string
  use solver_options

  implicit none

  ! Local variables
  logical :: readLevelset
  character(len=str_medium) :: levelsetFile

  ! Read levelset
  call parser_is_defined('levelset file to read', readLevelset)
  if (readLevelset) then
     call parser_read('levelset file to read', levelsetFile)
     call allocate_levelset
     call simulation_read(IO_LEVELSET, trim(levelsetFile))
     call grid_levelset_dependent_variables
  end if

  return
end subroutine grid_levelset_setup


! ======================== ! 
! Cleanup levelset routine !
! ======================== !
subroutine grid_levelset_cleanup

  ! Internal modules
  use grid_levelset

  implicit none

  if (allocated(levelset)) deallocate(levelset)
  if (allocated(levelsetNormal)) deallocate(levelsetNormal)
  if (allocated(indicatorFunction)) deallocate(indicatorFunction)
  if (allocated(primitiveGridNorm)) deallocate(primitiveGridNorm)

  return
end subroutine grid_levelset_cleanup


! ==================================== !
! Compute levelset dependent variables !
! ==================================== !
subroutine grid_levelset_dependent_variables

  ! Internal modules
  use grid_levelset

  ! External modules
  use parallel
  use simulation_flags
  use grid
  use grid_functions
  use filter

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: buf

  if (.not. allocated(levelset)) return

  ! Compute normal vector to surface via gradient of levelset
  call gradient(levelset(:,1), levelsetNormal)

  ! Filter the levelset normal to avoid discontinuous functions
  !if (useIBM) then
  !   do i = 1, nDimensions
  !      call gaussian_filter_apply(i, levelsetNormal(:,i:i))
  !   end do
  !end if
  
  do i = 1, nGridPoints

     ! Make the levelset normal a unit norm
     buf = sqrt(sum(levelsetNormal(i,:)**2))
     if (buf .gt. 0.0_WP) levelsetNormal(i,:) = levelsetNormal(i,:) / buf

     ! Compute indicator function
     if (levelset(i, 1) .le. 0.0_WP) then
        indicatorFunction(i, 1) = 0.0_WP
     else
        indicatorFunction(i, 1) = 1.0_WP
     end if

     ! Update the grid norm
     gridNorm(i, 1) = primitiveGridNorm(i, 1) * indicatorFunction(i, 1)
  end do

  ! Update the grid volume and volume fraction
  call compute_grid_volume(globalGridVolume)
  globalVolumeFraction = sum(primitiveGridNorm(:,1))
  call parallel_sum(globalVolumeFraction)
  globalVolumeFraction = globalGridVolume / globalVolumeFraction

  return
end subroutine grid_levelset_dependent_variables
