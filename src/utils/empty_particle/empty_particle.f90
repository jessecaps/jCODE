program empty_particle

  ! External modules
  use precision
  use string
  use fileio

  implicit none
  
  ! Particle data
  integer :: npart, iunit, ierror
  character(len=str_medium) :: filename
  real(WP) :: time

  ! Define the particle type (must match t_Particle in particle.f90)
  type :: t_Particle
     integer(KIND=8) :: id
     real(WP) :: diameter
     real(WP) :: temperature
     real(WP), dimension(3) :: position
     real(WP), dimension(3) :: velocity
     real(WP), dimension(3) :: angularVelocity
     real(WP), dimension(3) :: collision
     real(WP), dimension(3) :: torque
     integer :: gridIndex(3)
     integer :: stop
  end type t_Particle
  type(t_Particle), dimension(:), allocatable :: particles
  
  ! Type size
  integer, parameter :: particleSize = 160
  
  ! Read file name from standard input
  print*,'================================='
  print*,'| jCODE - create empty part file |'
  print*,'================================='
  print*
  print "(a28,$)", " part file to write : "
  read "(a)", filename

  ! Set zero particles
  npart = 0

  ! Allocate part
  allocate(particles(1:npart))
  
  ! Set time info
  time = 0.0_WP

  ! WRITE EMPTY PARTICLE FILE
  ! ** Open the part file to write **
  call BINARY_FILE_OPEN(iunit, trim(filename),"w",ierror)
  call BINARY_FILE_WRITE(iunit, npart, 1, kind(npart), ierror)
  call BINARY_FILE_WRITE(iunit, particleSize,1, kind(particleSize), ierror)
  call BINARY_FILE_WRITE(iunit, time, 1, kind(time), ierror)
  call BINARY_FILE_WRITE(iunit, particles, npart, particleSize, ierror)
  call BINARY_FILE_CLOSE(iunit, ierror)
  
end program empty_particle
