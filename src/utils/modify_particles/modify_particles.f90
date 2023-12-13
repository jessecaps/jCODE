program modify_particles

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use geometry
  use grid
  use state
  use equation_of_state
  use time_info
  use particle
  use particle_exchange, only : interpolate_fluid_to_particle

  implicit none

  ! Program choices
  integer, parameter :: RESET_VELOCITY    = 1
  integer, parameter :: INTERPOLATE_FLUID = 2
  integer, parameter :: GLUE_PARTICLES    = 3
  integer, parameter :: GLUE_BOTTOM       = 4

  ! Local variables
  integer :: choice, p, i
  real(WP) :: fluidVelocity(3)
  character(len = str_medium) :: input, gridFile, solutionFile, partFile1, partFile2, message

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  ! Read information from standard input
  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '=================================='
     write (*,*) '| jCODE - Modify a particle file |'
     write (*,*) '=================================='
     write (*,*)
     write (*,*) '1. Zero out velocity'
     write (*,*) '2. Interpolate fluid velocity to particle'
     write (*,*) '3. Glue particles in place'
     write (*,*)
     write (*,"(a10)", advance = 'no') ' Choice : '
     read "(i1)", choice
     write (*,*)
     write (*,"(a29)", advance = 'no') ' Particle filename to read : '
     read "(A)", partFile1
     write (*,"(a30)", advance = 'no') ' Particle filename to write : '
     read "(A)", partFile2
  end if
  call parallel_bc(choice)
  call parallel_bc(partFile1)
  call parallel_bc(partFile2)
  !call parser_read('analyzer choice',choice)
  !call parser_read('particle filename to read',partFile1)
  !call parser_read('particle filename to write',partFile2)

  ! Set up the grid and stencil operators
  call simulation_flags_setup
  call get_dimensions(gridFile)
  call geometry_setup
  call get_nUnknowns
  call solver_options_setup
  call operator_setup
  call grid_setup
  call simulation_read(IO_GRID, trim(gridFile))

  select case (choice)

  case (RESET_VELOCITY)

     ! Zero out the particle velocity
     ! ------------------------------

     ! Read the particle file
     call simulation_read(IO_PARTICLE, trim(partFile1))

     ! Modify the particles
     do p = 1, nParticles
        particles(p)%velocity = 0.0_WP
     end do

     ! Write the new particle file
     call simulation_write(IO_PARTICLE, trim(partFile2))

  case (INTERPOLATE_FLUID)

     ! Interpolate the fluid velocity to particle position
     ! ---------------------------------------------------

     ! Read the grid file
     call simulation_read(IO_GRID, trim(gridFile))

     ! Compute the grid metrics
     call grid_metrics_setup

     ! Get the solution filename and reset the solver options
     if (iRank .eq. iRoot) then
        write (*,"(a25)", advance = 'no') ' Solution file to read : '
        read "(A)", solutionFile
     end if
     call parallel_bc(solutionFile)
     !call parser_read('solution filename to read', solutionFile)
     call get_nUnknowns
     call solver_options_setup

     ! Allocate the state arrays
     call allocate_state

     ! Read the solution file
     call simulation_read(IO_FORWARD_STATE, trim(solutionFile))

     ! Get the fluid velocity
     call compute_dependent_variables(conservedVariables, velocity = velocity)

     ! Setup the particle solver
     call particle_solver_setup

     ! Setup the Cartesian arrays for interphase exchange
     call particle_exchange_setup

     ! Prepare fluid arrays for interpolation
     call prepare_cartesian_arrays

     ! Read the particle file
     call simulation_read(IO_PARTICLE, trim(partFile1))

     ! Modify the particles
     do p = 1, nParticles
        fluidVelocity = 0.0_WP
        do i = 1, nDimensions
           call interpolate_fluid_to_particle(particles(p)%gridIndex, particles(p)%position, &
                velocity = fluidVelocity)
        end do
        particles(p)%velocity = fluidVelocity
     end do

     ! Write the new particle file
     call simulation_write(IO_PARTICLE, trim(partFile2))

  case (GLUE_PARTICLES)

     ! Set particle id to 0 to prevent them from moving during run time
     ! ----------------------------------------------------------------

     ! Read the particle file
     call simulation_read(IO_PARTICLE, trim(partFile1))

     ! Modify the particles
     do p = 1, nParticles
        particles(p)%id = int(0, kind = 8)
     end do

     ! Write the new particle file
     call simulation_write(IO_PARTICLE, trim(partFile2))

  case (GLUE_BOTTOM)

     ! Glue particles located near y=0
     ! -------------------------------

     ! Read the particle file
     call simulation_read(IO_PARTICLE, trim(partFile1))

     ! Modify the particles
     do p = 1, nParticles
        if (particles(p)%position(2).lt.particles(p)%diameter)                               &
             particles(p)%id = int(0, kind = 8)
     end do

     ! Write the new particle file
     call simulation_write(IO_PARTICLE, trim(partFile2))

  case default

     write(message, '(A,I1,A)') "Unknown choice: ", choice, "!"
     call die(trim(message))

  end select

  ! Finalize the parallel environment
  call parallel_finalize

end program modify_particles


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
