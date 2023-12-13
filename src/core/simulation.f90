module simulation

  ! External modules
  use precision
  use string
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use state
  use time_info
  use particle

  implicit none

end module simulation


! ==================== !
! Setup the simulation !
! ==================== !
subroutine simulation_setup

  ! Internal modules
  use simulation

  ! External modules
  use geometry
  use grid
  use grid_io

  implicit none

  ! Local variables
  character(len = str_medium) :: gridFile

  ! Clean slate
  call simulation_cleanup

  ! Initialize simulation flags
  call simulation_flags_setup

  ! Get the number of dimensions from the grid file header
  call get_dimensions(gridFile)

  ! Initialize the geometric parameters
  call geometry_setup

  ! Determine the total number of conserved variables from the solution file
  call get_nUnknowns

  ! Initialize solver options
  call solver_options_setup

  ! Setup the stecil operators
  call operator_setup

  ! Setup the grid
  call grid_setup

  ! Read the grid file
  call simulation_read(IO_GRID, gridFile)

    ! Setup patches
  call grid_patch_setup

  ! Compute the grid metrics
  call grid_metrics_setup

  ! Setup the levelset
  call grid_levelset_setup

  ! Setup the state
  call state_setup

  ! Setup the particles
  call particle_setup

  ! Setup the solver
  call solver_setup

  return
end subroutine simulation_setup


! ====================== !
! Cleanup the simulation !
! ====================== !
subroutine simulation_cleanup

  ! Internal modules
  use simulation

  implicit none

  call solver_cleanup
  call boundary_cleanup
  call grid_patch_cleanup
  call grid_cleanup
  call operator_cleanup
  call state_cleanup
  call particle_cleanup
  call solver_options_cleanup

  return
end subroutine simulation_cleanup


! ==================== !
! Read simulation data !
! ==================== !
subroutine simulation_read(quantityOfInterest, filename)

  ! Internal modules
  use simulation

  implicit none

  ! Arguments
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename

  ! Start the i/o timer
  call timing_start('i/o')

  select case (quantityOfInterest)

  case (IO_GRID)

     ! Read the grid
     if (useSerialIO) then
        call grid_read_serial(trim(filename))
     else
        call grid_read_parallel(trim(filename))
     end if

  case (IO_FORWARD_STATE)

     ! Read the solution
     if (useSerialIO) then
        call state_read_serial(conservedVariables, trim(filename))
     else
        call state_read_parallel(conservedVariables, trim(filename))
     end if

  case (IO_TARGET_STATE)

     ! Read target data
     if (useSerialIO) then
        call state_read_serial(targetState, trim(filename))
     else
        call state_read_parallel(targetState, trim(filename))
     end if

  case (IO_ADJOINT_STATE)

     ! Read adjoint variables
     if (useSerialIO) then
        call state_read_serial(adjointVariables, trim(filename))
     else
        call state_read_parallel(adjointVariables, trim(filename))
     end if

  case (IO_PARTICLE)

     ! Read particle data
     if (useSerialIO) then
        call particle_read_serial(trim(filename))
     else
        call particle_read_parallel(trim(filename))
     end if

  case (IO_IBM)

     ! Read immersed boundary data
     if (useSerialIO) then
        call ibm_read_serial(trim(filename))
     else
        call ibm_read_parallel(trim(filename))
     end if

  case (IO_LEVELSET)

     ! Read immersed boundary data
     if (useSerialIO) then
        call levelset_read_serial(trim(filename))
     else
        call levelset_read_parallel(trim(filename))
     end if

  end select

  ! Log the file read
  call monitor_log('File "'//trim(filename)//'" read')

  ! Stop the i/o timer
  call timing_stop('i/o')

  return
end subroutine simulation_read


! ===================== !
! Write simulation data !
! ===================== !
subroutine simulation_write(quantityOfInterest, filename)

  ! Internal modules
  use simulation

  implicit none

  ! Arguments
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename

  ! Start the i/o timer
  call timing_start('i/o')

  select case (quantityOfInterest)

  case (IO_GRID)

     ! Write the grid
     if (useSerialIO) then
        call grid_write_serial(trim(filename))
     else
        call grid_write_parallel(trim(filename))
     end if

  case (IO_FORWARD_STATE)

     ! Write the solution
     if (useSerialIO) then
        call state_write_serial(conservedVariables, trim(filename))
     else
        call state_write_parallel(conservedVariables, trim(filename))
     end if

  case (IO_TARGET_STATE)

     ! Write the target state
     if (useSerialIO) then
        call state_write_serial(targetState, trim(filename))
     else
        call state_write_parallel(targetState, trim(filename))
     end if

  case (IO_ADJOINT_STATE)

     ! Write the adjoint state
     if (useSerialIO) then
        call state_write_serial(adjointVariables, trim(filename))
     else
        call state_write_parallel(adjointVariables, trim(filename))
     end if

  case (IO_PARTICLE)

     ! Write particle data
     if (useSerialIO) then
        call particle_write_serial(trim(filename))
     else
        call particle_write_parallel(trim(filename))
     end if

  case (IO_IBM)

     ! Write immersed boundary data
     if (useSerialIO) then
        call ibm_write_serial(trim(filename))
     else
        call ibm_write_parallel(trim(filename))
     end if

  case (IO_LEVELSET)

     ! Write immersed boundary data
     if (useSerialIO) then
        call levelset_write_serial(trim(filename))
     else
        call levelset_write_parallel(trim(filename))
     end if

  end select

  ! Log the file write
  call monitor_log('File "'//trim(filename)//'" written')

  ! Stop the i/o timer
  call timing_stop('i/o')

  return
end subroutine simulation_write
