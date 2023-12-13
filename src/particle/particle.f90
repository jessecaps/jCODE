module particle

  ! External modules
  use precision

  implicit none

  ! Define the particle type
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
  type(t_Particle), dimension(:), allocatable :: particles, ghostParticles

  ! Particle parameters
  integer :: nParticles=0, nParticlesGlobal=0, nGhostParticles=0
  integer :: nParticlesIn, nParticlesPrevious
  integer :: MPI_PARTICLE, SIZE_MPI_PARTICLE
  real(WP) :: particleDensity, particleSpecificHeat
  logical :: useSaffmanLift, useAddedMass, useParticleHeat, usePhaseChange, useFriction

  ! Monitor particle info
  type :: t_MonitorParticle
     real(WP) :: drag
     real(WP) :: stress
     real(WP) :: gravity
     real(WP) :: collision
     real(WP) :: lift
     real(WP) :: addedMass
  end type t_MonitorParticle
  type(t_MonitorParticle) :: monitorParticle(3)
  integer :: nParticleCollisions

contains


  ! Send `stopped` particles to the end and resize
  ! ----------------------------------------------
  subroutine recycle_particles(particleVector)

    implicit none

    ! Arguments
    type(t_Particle), allocatable, intent(inout) :: particleVector(:)

    ! Local variables
    integer :: i, newSize

    ! Compact real particles at the beginning of the array
    newSize = 0
    if (allocated(particleVector)) then
       do i = 1, size(particleVector)
          if (particleVector(i)%stop .eq. 0) then
             newSize = newSize + 1
             if (i .ne. newSize) then
                particleVector(newSize) = particleVector(i)
                particleVector(i)%stop = 1
             end if
          end if
       end do
    end if

    ! Resize array to newSize
    call resize_particles(particleVector, newSize)

    return
  end subroutine recycle_particles


  ! Resize the particle vector based on new size `n`
  ! -----------------------------------------------
  subroutine resize_particles(particleVector, n)

    implicit none

    ! Arguments
    type(t_Particle), allocatable, intent(inout) :: particleVector(:)
    integer, intent(in) :: n

    ! Local variables
    integer :: i, nParticlesOld
    type(t_Particle), allocatable :: tempParticles(:)

    ! Resize part array to size n
    if (.not. allocated(particleVector)) then
       ! Particle is of size 0
       if (n .eq. 0) then
          ! Nothing to do, that's what we want
       else
          ! Allocate directly of size n
          allocate(particleVector(n))
          particleVector(1:n)%stop = 1
       end if
    else if (n .eq. 0) then
       ! Empty the particle vector
       deallocate(particleVector); allocate(particleVector(0))
    else
       ! Update non zero size to another non zero size
       nParticlesOld = size(particleVector)
       if (n .gt. nParticlesOld) then
          ! Increase from `nParticlesOld` to `n`
          allocate(tempParticles(n))
          do i = 1, nParticlesOld
             tempParticles(i) = particleVector(i)
          end do
          deallocate(particleVector)
          allocate(particleVector(size(tempParticles)))
          particleVector = tempParticles
          particleVector(nParticlesOld + 1:n)%stop = 1
       else if (n .lt. nParticlesOld) then
          ! Decrease from `nParticlesOld` to `n`
          allocate(tempParticles(n))
          do i = 1, n
             tempParticles(i) = particleVector(i)
          end do
          deallocate(particleVector)
          allocate(particleVector(size(tempParticles)))
          particleVector = tempParticles
       end if
    end if

    return
  end subroutine resize_particles

end module particle


! ========================= !
! Setup the particle module !
! ========================= !
subroutine particle_setup

  ! Internal modules
  use particle

  ! External modules
  use simulation_flags

  implicit none

  if (.not. useParticles) return

  ! Setup the particle time integrator
  call particle_integrator_setup

  ! Setup the particle solver
  call particle_solver_setup

  ! Setup the Cartesian arrays for interphase exchange
  call particle_exchange_setup

  ! Setup the particle injector routine
  call particle_injector_setup

  return
end subroutine particle_setup


! =========================== !
! Cleanup the particle module !
! =========================== !
subroutine particle_cleanup

  ! Internal modules
  use particle

  implicit none

  call particle_injector_cleanup
  call particle_exchange_cleanup
  call particle_solver_cleanup
  call particle_integrator_cleanup

  if (allocated(particles)) deallocate(particles)
  if (allocated(ghostParticles)) deallocate(ghostParticles)
  nParticles = 0
  nParticlesGlobal = 0

  return
end subroutine particle_cleanup


! ============================================ !
! Setup MPI particle for communication and i/o !
! ============================================ !
subroutine prepare_mpi_particle

  ! Internal modules
  use particle

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer, dimension(23) :: types, lengths, displacement ! ... size based on t_Particle
  integer :: ierror

  ! Create the MPI structure to send particles
  types( 1: 1) = MPI_INTEGER8
  types( 2:18) = MPI_REAL_WP
  types(19:22) = MPI_INTEGER4
  types(23)    = MPI_UB
  lengths(:)   = 1

  ! Hard-code particle type displacement here
  displacement( 1) = 0                     ! ... ID
  displacement( 2) = displacement( 1) + 8  ! ... Diameter
  displacement( 3) = displacement( 2) + WP ! ... Temperature
  displacement( 4) = displacement( 3) + WP ! ... Position(1)
  displacement( 5) = displacement( 4) + WP ! ... Position(2)
  displacement( 6) = displacement( 5) + WP ! ... Position(3)
  displacement( 7) = displacement( 6) + WP ! ... Velocity(1)
  displacement( 8) = displacement( 7) + WP ! ... Velocity(2)
  displacement( 9) = displacement( 8) + WP ! ... Velocity(3)
  displacement(10) = displacement( 9) + WP ! ... AngularVelocity(1)
  displacement(11) = displacement(10) + WP ! ... AngularVelocity(2)
  displacement(12) = displacement(11) + WP ! ... AngularVelocity(3)
  displacement(13) = displacement(12) + WP ! ... Collision(1)
  displacement(14) = displacement(13) + WP ! ... Collision(2)
  displacement(15) = displacement(14) + WP ! ... Collision(3)
  displacement(16) = displacement(15) + WP ! ... Torque(1)
  displacement(17) = displacement(16) + WP ! ... Torque(2)
  displacement(18) = displacement(17) + WP ! ... Torque(3)
  displacement(19) = displacement(18) + WP ! ... GridIndex(1)
  displacement(20) = displacement(19) + 4  ! ... GridIndex(2)
  displacement(21) = displacement(20) + 4  ! ... GridIndex(3)
  displacement(22) = displacement(21) + 4  ! ... Stop
  displacement(23) = displacement(22) + 4  ! ... UpperBound

  ! Finalize by creating and commiting the new type
  call MPI_Type_struct(23, lengths, displacement, types, MPI_PARTICLE, ierror)
  call MPI_Type_commit(MPI_PARTICLE, ierror)

  ! If problem, say it
  if (ierror .ne. 0) call die('Problem with MPI_PARTICLE!')

  ! Get the size of this type
  call MPI_type_size(MPI_PARTICLE, SIZE_MPI_PARTICLE, ierror)

  return
end subroutine prepare_mpi_particle


! ======================== !
! Update the particle size !
! ======================== !
subroutine update_particle_size

  ! Internal modules
  use particle

  ! External modules
  use parallel

  implicit none

  ! Update number of particles
  nParticles = size(particles)
  call parallel_sum(nParticles, nParticlesGlobal)

  return
end subroutine update_particle_size
