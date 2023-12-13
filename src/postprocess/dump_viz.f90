module dump_viz

  ! External modules
  use precision
  use string
  use grid_patch

  implicit none
  
  ! Type of outputs
  integer :: outputType
  integer, parameter ::                                                                      &
       NO_OUTPUT      = 0,                                                                   &
       ENSIGHT_OUTPUT = 1,                                                                   &
       VISIT_OUTPUT   = 2

  ! Visualization patch to dump a portion of the solution
  integer :: nVizPatches
  type(t_Patch), pointer :: vizPatch(:)

  ! Frequency of output
  real(WP) :: outputFrequency

contains

  ! Setup the visualization patch
  ! -----------------------------
  subroutine setup_viz_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch

    ! If no viz patch is specified define the patch over the entire domain
    if (nVizPatches .eq. 0) then
       patch%name = 'jcode'
       patch%patchType = VISUALIZATION
       patch%normalDirection = 0
       patch%iMin = 1; patch%iMax = -1
       patch%jMin = 1; patch%jMax = -1
       patch%kMin = 1; patch%kMax = -1
       call patch_setup(patch)
       nVizPatches = 1
    end if

    ! Verify the patch type
    call verify_viz_patch(patch)

    return
  end subroutine setup_viz_patch
  

  ! Verify that the visualization patch is correct
  ! ----------------------------------------------
  subroutine verify_viz_patch(patch)

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    if (patch%normalDirection .ne. 0) call die ('verify_viz_patch: normal direction /= 0!')

    return
  end subroutine verify_viz_patch
  
end module dump_viz


! ================================== !
! Initialize data dumping procedures !
! ================================== !
subroutine dump_setup(mode)

  ! Internal modules
  use dump_viz

  ! External modules
  use parser

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Local variables
  integer :: i, j
  character(len = str_medium) :: type
  
  ! Read the input file
  call parser_read('output type', type, '')
  call parser_read('output frequency', outputFrequency, -1.0_WP)

  ! Setup the output type
  select case (trim(type))

  case ('')
     ! No output type
     return

  case ('ensight')
     outputType = ENSIGHT_OUTPUT

  case ('visit')
     outputType = VISIT_OUTPUT

  case default
     call die("dump_setup: unknown output type '" // trim(type) // "'")

  end select

  ! Find the number of viz patches
  nVizPatches = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. VISUALIZATION) then
        nVizPatches = nVizPatches + 1
        if (j .eq. 0) j = i
     end if
  end do

  ! If no viz patch is specified make one
  if (nVizPatches .gt. 0) then
     vizPatch => patches(j:j+nVizPatches-1)
  else
     allocate(vizPatch(1))
  end if
  
  ! Setup the viz patches
  do i = 1, max(nVizPatches, 1)
     call setup_viz_patch(vizPatch(i))
  end do

  ! Setup the specific visualization routine
  select case (outputType)

  case (ENSIGHT_OUTPUT)
     call dump_ensight_setup(mode)

  case (VISIT_OUTPUT)
     ! Not yet implemented

  end select
  
  return
end subroutine dump_setup


! =============================== !
! Cleanup data dumping procedures !
! =============================== !
subroutine dump_cleanup

  ! Internal modules
  use dump_viz

  implicit none

  call dump_ensight_cleanup

  if (nVizPatches .gt. 0) nullify(vizPatch)

  nVizPatches = 0

  return
end subroutine dump_cleanup


! ======================================== !
! Dump data at each time step if necessary !
! ======================================== !
subroutine dump_result(mode)

  ! Internal modules
  use dump_viz

  ! External modules
  use time_info

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Check if time to dump
  if (outputFrequency .gt. 0.0_WP) then
     if (mod(time, outputFrequency) .lt. 0.5_WP * timeStepSize .or.                          &
          mod(time, outputFrequency) .ge. outputFrequency - 0.5_WP * timeStepSize) then

        ! Select dump format
        select case (outputType)

        case (ENSIGHT_OUTPUT)
           call dump_ensight_data(mode)
           call monitor_log('EnSight file written')

        case (VISIT_OUTPUT)
           call monitor_log('VisIt file written')

        end select

     end if
  end if

  return
end subroutine dump_result
