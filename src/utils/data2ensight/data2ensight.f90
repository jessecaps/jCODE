program data2ensight

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use geometry
  use grid
  use grid_patch
  use state
  use equation_of_state
  use time_info
  use dump_viz
  use dump_ensight

  implicit none

  ! Local variables
  integer :: i, j, mode
  character(len = str_medium) :: input, filename, partname, gridFile
  logical :: usePart

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  ! Read information from input file

  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '==========================================='
     write (*,*) '| jCODE - Data to ENSIGHT GOLD converter  |'
     write (*,*) '==========================================='
     write (*,*)
  end if

  ! Read file data
  call parser_read('data file to convert', filename)
  call parser_is_defined('particle file to convert', usePart)
  if (usePart) call parser_read('particle file to convert', partname)
  call parser_read('mode (+1 for forward or -1 for adjoint)', mode, 1)

  ! Check for errors
  if (mode .ne. FORWARD .and. mode .ne. ADJOINT) call die('Unknown mode')

  ! Set up the grid and stencil operators
  call simulation_flags_setup
  call get_dimensions(gridFile)
  call geometry_setup
  call get_nUnknowns
  call solver_options_setup
  call operator_setup
  call grid_setup
  call simulation_read(IO_GRID, trim(gridFile))
  call grid_patch_setup

  ! Setup the metrics and prepare the state
  call grid_metrics_setup
  call grid_levelset_setup
  call state_setup
  if (usePart) call particle_setup

  if (.not. predictionOnly) then
     select case (mode)
     case (FORWARD)
        ! Setup the cost functional for computing the target mollfier
        call functional_setup
     case (ADJOINT)
        ! Setup the controller for computing the controller mollfier
        call controller_setup
     end select
  end if

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

  ! Setup EnSight and write the geometry
  call dump_ensight_setup(mode)

  select case (mode)

  case (FORWARD)
    
     call simulation_read(IO_FORWARD_STATE, trim(filename))

     ! Read in the particle file (if you want)
     if (dumpParticles .and. usePart) then
        call simulation_read(IO_PARTICLE, trim(partname))
     end if
     
     ! Update the state and write the EnSight files
     call update_state

  case (ADJOINT)

     call simulation_read(IO_ADJOINT_STATE, trim(filename))

  end select
     
  ! Write the EnSight files
  call dump_ensight_data(mode)

  ! Finalize the parallel environment
  call parallel_finalize

end program data2ensight



! -------------------------------------------
subroutine die(errorText)

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: errorText

  call parallel_kill(errorText)
  
  return
end subroutine die
