module stat

  ! External modules
  use precision
  use string

  implicit none
  
  ! Type of outputs
  integer :: nStat
  integer, allocatable :: statType(:)
  integer, parameter :: COMPRESS_0D    = 0,                                                  &
                        COMPRESS_1D    = 1,                                                  &
                        COMPRESS_2D    = 2,                                                  &
                        COMPRESS_PATCH = 3,                                                  &
                        COMPRESS_IBM   = 4,                                                  &
                        COMPRESS_IBM_1D= 5,                                                  &
                        COMPRESS_PART  = 6

  ! Frequency of output
  real(WP) :: statFrequency
  
end module stat


! ======================================== !
! Initialize statistics dumping procedures !
! ======================================== !
subroutine stat_setup

  ! Internal modules
  use stat

  ! External modules
  use parser
  use parallel

  implicit none

  ! Local variables
  integer :: i
  logical :: isDefined
  character(len = str_medium), allocatable :: inputType(:)
  
  ! Read the input file
  call parser_is_defined('stat type', isDefined)
  nStat = 0

  ! Return if nothing to do
  if (.not. isDefined) return

  ! Create the directory
  if (iRank .eq. iRoot) call CREATE_FOLDER("stat")

  ! Determine number of stats
  call parser_getsize('stat type', nStat)
  allocate(statType(nStat))
  allocate(inputType(nStat))
  call parser_read('stat type', inputType)
  
  ! State frequency
  call parser_read('stat frequency', statFrequency, -1.0_WP)

  ! Setup the output type
  do i = 1, nStat
     select case (trim(inputType(i)))

     case ('0D')
        statType(i) = COMPRESS_0D

     case ('1D')
        statType(i) = COMPRESS_1D
        call stat_1D_setup

     case ('2D')
        statType(i) = COMPRESS_2D

     case ('IBM', 'ibm')
        statType(i) = COMPRESS_IBM
        call stat_ibm_setup

     case ('IBM-1D', 'ibm-1d')
        statType(i) = COMPRESS_IBM_1D
        call stat_ibm_1d_setup

     case ('PATCH', 'patch', 'Patch')
        statType(i) = COMPRESS_PATCH
        call stat_patch_setup

     case ('particle')
        statType(i) = COMPRESS_PART
        call stat_part_setup

     case default
        call die("stat_setup: unknown output type '" // trim(inputType(i)) // "'")

     end select
  end do

  ! Cleanup
  deallocate(inputType)
  
  return
end subroutine stat_setup


! ===================================== !
! Cleanup statistics dumping procedures !
! ===================================== !
subroutine stat_cleanup

  ! Internal modules
  use stat

  implicit none

  ! Local variables
  integer :: i

  do i = 1, nStat
     select case (statType(i))

     case (COMPRESS_0D)

     case (COMPRESS_1D)
        call stat_1d_cleanup

     case (COMPRESS_2D)

     case (COMPRESS_IBM)
        call stat_ibm_cleanup

     case (COMPRESS_IBM_1D)
        call stat_ibm_1d_cleanup

     case (COMPRESS_PATCH)
        call stat_patch_cleanup

     case (COMPRESS_PART)
        call stat_part_cleanup

     end select
  end do

  if (allocated(statType)) deallocate(statType)

  return
end subroutine stat_cleanup


! ========================================= !
! Dump stats at each time step if necessary !
! ========================================= !
subroutine dump_statistics

  ! Internal modules
  use stat

  ! External modules
  use time_info

  implicit none

  ! Local variables
  integer :: i

  ! Return if not used
  if (nStat .eq. 0) return

  ! Compute statistics
  do i = 1, nStat
     select case (statType(i))

     case (COMPRESS_0D)

     case (COMPRESS_1D)
        call stat_1d_compute

     case (COMPRESS_2D)

     case (COMPRESS_IBM)
        call stat_ibm_compute

     case (COMPRESS_IBM_1D)
        call stat_ibm_1d_compute

     case (COMPRESS_PATCH)
        call stat_patch_compute

     case (COMPRESS_PART)
        call stat_part_compute

     end select
  end do

  ! Dump stats if needed
  do i = 1, nStat
     if (statFrequency .gt. 0.0_WP) then
        if (mod(time, statFrequency) .lt. 0.5_WP * timeStepSize .or.                         &
             mod(time, statFrequency) .ge. statFrequency - 0.5_WP * timeStepSize) then

           ! Select statistic format
           select case (statType(i))

           case (COMPRESS_0D)

           case (COMPRESS_1D)
              call stat_1d_write

           case (COMPRESS_2D)

           case (COMPRESS_IBM)
              call stat_ibm_write

           case (COMPRESS_IBM_1D)
              call stat_ibm_1d_write

           case (COMPRESS_PATCH)
              call stat_patch_write

           case (COMPRESS_PART)
              call stat_part_write

           end select

        end if
     end if
  end do

  return
end subroutine dump_statistics
