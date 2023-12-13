program ibm2part

  ! External modules
  use precision
  use string
  use math
  use parser
  use parallel
  use geometry
  use grid
  use state
  use equation_of_state
  use time_info
  use ibm
  use particle

  implicit none

  ! Local variables
  integer :: i, stationary
  character(len = str_medium) :: input, ibmFile, partFile, gridFile

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  ! Read information from standard input
  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '========================================='
     write (*,*) '| jCODE - convert IBM file to part file |'
     write (*,*) '========================================='
     write (*,*)
     write (*,"(a20)", advance = 'no') ' IBM file to read : '
     read "(A)", ibmFile
     write (*,"(a26)", advance = 'no') ' Particle file to write : '
     read "(A)", partFile
     write (*,"(a40)", advance = 'no') ' Stationary particles (0: no, 1: yes) : '
     read "(i1)", stationary
  end if
  call parallel_bc(ibmFile)
  call parallel_bc(partFile)

  ! Set up the grid and stencil operators
  call simulation_flags_setup
  useParticles = .true. ! Force the code to use particles
  call get_dimensions(gridFile)
  call geometry_setup
  call get_nUnknowns
  call solver_options_setup
  call operator_setup
  call grid_setup
  call simulation_read(IO_GRID, trim(gridFile))

  ! Read the IBM file
  call simulation_read(IO_IBM, trim(ibmFile))

  ! Create the particles
  nParticlesGlobal = nObjects
  nParticles = nParticlesGlobal
  if (iRank .ne. iRoot) nParticles = 0
  allocate(particles(nParticles))
  do i = 1, nParticles
     if (stationary.eq.1) then
        particles(i)%id = 0
     else
        particles(i)%id = i
     end if
     particles(i)%diameter = (2.0_WP * real(nDimensions, WP) * object(i)%volume / pi)        &
          ** (1.0_WP / real(nDimensions, WP))
     particles(i)%temperature = 2.5_WP
     particles(i)%position = object(i)%position
     particles(i)%velocity = object(i)%velocity
     particles(i)%angularVelocity = object(i)%angularVelocity
     particles(i)%collision = object(i)%cForce
     particles(i)%torque = object(i)%cTorque
     particles(i)%stop = 0
  end do
  
  ! Write the new particle file
  call simulation_write(IO_PARTICLE, trim(partFile))

  ! Finalize the parallel environment
  call parallel_finalize

end program ibm2part


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
