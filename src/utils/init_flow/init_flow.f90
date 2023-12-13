program init_flow

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use simulation_flags
  use geometry
  use solver_options
  use grid
  use time_info
  use particle

  implicit none

  ! Local variables
  integer :: i
  character(len = str_medium) :: input, filename, val
  character(len = 1) :: xyz(3)

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)

  ! Disable manual processor decomposition
  disableManualDecomp = .true.

  ! Get the grid dimensions from input
  xyz(1) = 'x'; xyz(2) = 'y'; xyz(3) = 'z'
  nDimensions = 0
  do i = 1, 3
     call parser_read('n'//xyz(i), globalGridSize(i), 1)
     if (globalGridSize(i) .gt. 1) nDimensions = nDimensions + 1
  end do

  ! Get periodicity information from input
  periodicityType = NONE
  periodicLength = 0.0_WP
  do i = 1, nDimensions
     call parser_read('periodicity type in ' // xyz(i), val, '')
     periodicityType(i) = NONE
     if (trim(val) .eq. 'plane') then
        periodicityType(i) = PLANE
     else if (trim(val) .eq. 'overlap') then
        periodicityType(i) = OVERLAP
     else if (trim(val) .eq. 'polar') then
        periodicityType(i) = POLAR
     end if
     if (periodicityType(i) == PLANE)                                                        &
          call parser_read('L' // xyz(i), periodicLength(i))
  end do
  isPeriodic = (periodicityType .ne. NONE)

  ! Is the domain curvilinear?
  call parser_read('curvilinear domain', isDomainCurvilinear, .false.)

  ! Initialize some routines to access later
  call simulation_flags_setup
  call geometry_setup
  call operator_setup
  call grid_setup

  ! Initialize the simulation
  call get_simulation

  ! Report grid information
  call monitor_grid_diagnostics

  ! Initialize the time
  time = 0.0_WP

  ! Write the grid
  call parser_read('init grid file', filename)
  call simulation_write(IO_GRID, filename)

  ! Write initial conditions
  call parser_read('init solution file', filename)
  call simulation_write(IO_FORWARD_STATE, filename)

  ! Write the particle file
  call parser_read('init particle file', filename, '')
  if (len_trim(filename).ne.0) call simulation_write(IO_PARTICLE, filename)

  ! Write the IBM file
  call parser_read('init ibm file', filename, '')
  if (len_trim(filename).ne.0) call simulation_write(IO_IBM, filename)

  ! Write the levelset file
  call parser_read('init levelset file', filename, '')
  if (len_trim(filename).ne.0) call simulation_write(IO_LEVELSET, filename)

  ! Write the target solution if necessary
  call parser_read('init target state file', filename, '')
  if (len_trim(filename).ne.0) call simulation_write(IO_TARGET_STATE, filename)

  ! Finalize the parallel environment
  call parallel_finalize

contains

  subroutine get_simulation

    implicit none

    ! Local variables
    character(len = str_medium) :: simulation

    ! Detect the simulation type
    call parser_read('simulation name', simulation)

    select case (trim(simulation))

    case ('boundary layer')
       call boundary_layer_grid
       call boundary_layer_data

    case ('channel')
       call channel_grid
       call channel_data

    case ('cylinder')
       call cylinder_grid
       call cylinder_data

    case ('el ibm')
       call el_ibm_grid
       call el_ibm_data
       call el_ibm_objects
       call el_ibm_particles

    case ('erosion')
       call erosion_grid
       call erosion_data
       call erosion_particles

    case ('flame')
       call flame_grid
       call flame_data

    case ('fluidized bed')
       call fluidized_bed_grid
       call fluidized_bed_data
       call fluidized_bed_particles

    case ('HIT')
       call hit_grid
       call hit_data
       call hit_ibm

    case ('ibm')
       call ibm_init_grid
       call ibm_init_data
       call ibm_init_levelset

    case ('ibm particle')
       call ibm_part_grid
       call ibm_part_data
       call ibm_part_objects

    case ('impulsive plate')
       call impulsive_plate_grid
       call impulsive_plate_data
       call impulsive_plate_particles

    case ('jet')
       call jet_grid
       call jet_data

    case ('jet impingement')
       call jet_impingement_grid
       call jet_impingement_data
       call jet_impingement_particles
       call jet_impingement_levelset

    case ('jet in crossflow')
       call jet_crossflow_grid
       call jet_crossflow_data

    case ('laser induced breakdown')
       call lib_grid
       call lib_data

    case ('mixing layer')
       call mixing_layer_grid
       call mixing_layer_data

    case ('nozzle')
       call nozzle_grid
       call nozzle_data
       call nozzle_levelset
       call nozzle_particles
       call nozzle_ibm_particles

    case ('particle box')
       call particle_box_grid
       call particle_box_data
       call particle_box_particles

    case ('pipe')
       call pipe_grid
       call pipe_data

    case ('premixed')
       call premixed_grid
       call premixed_data

    case ('quiescent')
       call quiescent_grid
       call quiescent_data

    case ('rayleigh taylor')
       call rayleigh_taylor_grid
       call rayleigh_taylor_data

    case ('scalar convection')
       call convect_scalar_grid
       call convect_scalar_data

    case ('shear layer')
       call shear_layer_grid
       call shear_layer_data
       call shear_layer_particles

    case ('shock tube')
       call shock_tube_grid
       call shock_tube_data
       call shock_tube_particles

    case ('shu osher')
       call shu_osher_grid
       call shu_osher_data

    case ('sphere')
       call sphere_grid
       call sphere_data

    case ('taylor green vortex')
       call taylor_green_grid
       call taylor_green_data

    case default
       call die("Unknown simulation '" // trim(simulation) // "'!")

    end select

    return
  end subroutine get_simulation

end program init_flow



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
