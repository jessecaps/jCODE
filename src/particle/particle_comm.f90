module particle_comm

  ! External modules
  use particle
  use particle_exchange

  implicit none

contains

  ! Communicate particles in a given direction
  ! ------------------------------------------
  subroutine communicate_direction(dir)

    ! External modules
    use parallel

    implicit none

    ! Arguments
    integer, intent(in) :: dir

    ! Local variables
    integer :: i, i1, i2, rankLeft, rankRight, nSendLeft, nSendRight, nRecv, ierror
    type(t_Particle), dimension(:), allocatable :: particlesToRecv, particlesToSendLeft,     &
         particlesToSendRight

    ! Initialize MPI variables
    call MPI_Cart_shift(commDir(dir), 0, 1, rankLeft, rankRight, ierror)

    ! Get processor extents
    i1 = iStart(dir)
    i2 = iEnd(dir)

    ! Determine the rank the particle should belong to
    nSendLeft = 0; nSendRight = 0
    do i = 1, nParticles
       if (particles(i)%gridIndex(dir) .lt. i1) then
          nSendLeft = nSendLeft + 1
       else if (particles(i)%gridIndex(dir) .gt. i2) then
          nSendRight = nSendRight + 1
       end if
    end do

    ! Allocate the send buffers
    allocate(particlesToSendLeft(nSendLeft))
    allocate(particlesToSendRight(nSendRight))

    ! Copy particles to buffer
    nSendLeft = 0; nSendRight = 0
    do i = 1, nParticles
       if (particles(i)%gridIndex(dir) .lt. i1) then
          nSendLeft = nSendLeft + 1
          particlesToSendLeft(nSendLeft) = particles(i)
          particles(i)%stop = 1 ! ... remove the particle
       else if (particles(i)%gridIndex(dir) .gt. i2) then
          nSendRight = nSendRight + 1
          particlesToSendRight(nSendRight) = particles(i)
          particles(i)%stop = 1 ! ... remove the particle
       end if
    end do

     ! <<<<<<<< Communicate LEFT <<<<<<<<<<<<

     ! Send to left, receive from right
     nRecv = 0
     call MPI_SendRecv(nSendLeft, 1, MPI_INTEGER, rankLeft, 0, nRecv, 1, MPI_INTEGER,        &
          rankRight, 0, commDir(dir), MPI_STATUS_IGNORE, ierror)
     allocate(particlesToRecv(nRecv))
     call MPI_SendRecv(particlesToSendLeft, nSendLeft, MPI_PARTICLE, rankLeft, 1,            &
          particlesToRecv, nRecv, MPI_PARTICLE, rankRight,1, commDir(dir),  &
          MPI_STATUS_IGNORE, ierror)
     deallocate(particlesToSendLeft)

     ! Add to the particle vector
     if (nRecv .gt. 0) then
        call resize_particles(particles, nParticles + nRecv)
        particles(nParticles + 1 : nParticles + nRecv) = particlesToRecv(1:nRecv)
        nParticles = nParticles + nRecv
     end if
     deallocate(particlesToRecv)

     ! >>>>>>>> Communicate RIGHT >>>>>>>>>

     ! Send to right, receive from left
     nRecv = 0
     call MPI_SendRecv(nSendRight, 1, MPI_INTEGER, rankRight, 2, nRecv, 1, MPI_INTEGER,      &
          rankLeft, 2, commDir(dir), MPI_STATUS_IGNORE, ierror)
     allocate(particlesToRecv(nRecv))
     call MPI_SendRecv(particlesToSendRight, nSendRight, MPI_PARTICLE, rankRight, 3,         &
          particlesToRecv, nRecv, MPI_PARTICLE, rankLeft, 3, commDir(dir),                   &
          MPI_STATUS_IGNORE, ierror)
     deallocate(particlesToSendRight)

     ! Add to the particle vector
     if (nRecv .gt. 0) then
        call resize_particles(particles, nParticles + nRecv)
        particles(nParticles + 1:nParticles + nRecv) = particlesToRecv(1:nRecv)
        nParticles = nParticles + nRecv
     end if
     deallocate(particlesToRecv)

     ! Recycle
     call recycle_particles(particles)

     ! Update number of particles
     call update_particle_size

    return
  end subroutine communicate_direction


  ! Communicate ghost particles in a given direction
  ! ------------------------------------------------
  subroutine communicate_border(dir)

    ! External modules
    use parallel
    use math, only : bisection

    implicit none

    ! Arguments
    integer, intent(in) :: dir

    ! Local variables
    integer :: n, nSend, nRecv, rankLeft, rankRight, tag(3,2), ierror
    integer :: iStatus(MPI_STATUS_SIZE)
    type(t_Particle), dimension(:), allocatable :: particlesToSend, particlesToRecv

    ! Define unitque tags for send/ecv
    tag(1,1) = 0
    tag(1,2) = 2
    tag(2,1) = 4
    tag(2,2) = 6
    tag(3,1) = 8
    tag(3,2) = 10

    ! <<<<<<<< Communicate LEFT <<<<<<<<<<<<

    ! Count number of particles to send left
    nSend = 0
    do n = 1, nParticles
       if (particles(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1))           &
            nSend = nSend + 1
    end do
    do n = 1, nGhostParticles
       if (ghostParticles(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1))      &
            nSend = nSend + 1
    end do

    ! Allocate send buffer
    allocate(particlesToSend(nSend))

    ! Copy particles in buffer
    nSend = 0
    do n = 1, nParticles
       if (particles(n)%gridIndex(dir) .lt. gridOffset(dir) + 1 + nOverlap(dir,1)) then
          nSend = nSend + 1
          particlesToSend(nSend) = particles(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. 0) then
             particlesToSend(nSend)%position(dir) = particlesToSend(nSend)%position(dir) +   &
                  periodicLength(dir)
             particlesToSend(nSend)%gridIndex(dir) = particlesToSend(nSend)%gridIndex(dir) + &
                  globalGridSize(dir)
          end if
       end if
    end do
    do n = 1, nGhostParticles
       if (ghostParticles(n)%gridIndex(dir) .lt. gridOffset(dir) +1 + nOverlap(dir,1)) then
          nSend = nSend + 1
          particlesToSend(nSend) = ghostParticles(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. 0) then
             particlesToSend(nSend)%position(dir) = particlesToSend(nSend)%position(dir) +   &
                  periodicLength(dir)
             particlesToSend(nSend)%gridIndex(dir) = particlesToSend(nSend)%gridIndex(dir) + &
                  globalGridSize(dir)
          end if
       end if
    end do

    ! Communicate sizes
    nRecv = 0
    call MPI_Cart_shift(commDir(dir), 0, -1, rankRight, rankLeft, ierror)
    call MPI_SendRecv(nSend, 1, MPI_INTEGER, rankLeft, tag(dir,1), nRecv, 1, MPI_INTEGER,    &
         rankRight, tag(dir,1), commDir(dir), MPI_STATUS_IGNORE, ierror)

    ! Allocate recv buffer
    allocate(particlesToRecv(nRecv))

    ! Communicate
    call MPI_SendRecv(particlesToSend, nSend, MPI_PARTICLE, rankLeft, 0, particlesToRecv,    &
         nRecv, MPI_PARTICLE, rankRight, 0, commDir(dir), MPI_STATUS_IGNORE, ierror)

    ! Add to ghost particles
    if (nRecv .gt. 0) then
       call resize_particles(ghostParticles, nGhostParticles+nRecv)
       ghostParticles(nGhostParticles+1:nGhostParticles + nRecv) = particlesToRecv(1:nRecv)
       nGhostParticles = nGhostParticles + nRecv
    end if

    ! Clean up
    deallocate(particlesToSend)
    deallocate(particlesToRecv)

    ! >>>>>>>> Communicate RIGHT >>>>>>>>>

    ! Count number of particles to send right
    nSend = 0
    do n = 1, nParticles
       if (particles(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -           &
            nOverlap(dir,2)) nSend = nSend + 1
    end do
    do n = 1, nGhostParticles - nRecv
       if (ghostParticles(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -      &
            nOverlap(dir,2)) nSend = nSend + 1
    end do

    ! Allocate send buffer
    allocate(particlesToSend(nSend))

    ! Copy particles in buffer
    nSend = 0
    do n = 1, nParticles
       if (particles(n)%gridIndex(dir) .gt. localGridSize(dir) + gridOffset(dir) -           &
            nOverlap(dir,2)) then
          nSend = nSend + 1
          particlesToSend(nSend) = particles(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir)-1) then
             particlesToSend(nSend)%position(dir) = particlesToSend(nSend)%position(dir) -   &
                  periodicLength(dir)
             particlesToSend(nSend)%gridIndex(dir) = particlesToSend(nSend)%gridIndex(dir) - &
                  globalGridSize(dir)
          end if
       end if
    end do
    do n = 1, nGhostParticles - nRecv
       if (ghostParticles(n)%gridIndex(dir) .gt. localGridSize(dir)  + gridOffset(dir) -     &
            nOverlap(dir,2)) then
          nSend = nSend + 1
          particlesToSend(nSend) = ghostParticles(n)
          if (isPeriodic(dir) .and. procCoords(dir) .eq. nProcsDir(dir)-1) then
             particlesToSend(nSend)%position(dir) = particlesToSend(nSend)%position(dir) -   &
                  periodicLength(dir)
             particlesToSend(nSend)%gridIndex(dir) = particlesToSend(nSend)%gridIndex(dir) - &
                  globalGridSize(dir)
          end if
       end if
    end do

    ! Communicate sizes
    nRecv = 0
    call MPI_Cart_shift(commDir(dir), 0, +1, rankLeft, rankRight, ierror)
    call MPI_SendRecv(nSend, 1, MPI_INTEGER, rankRight, tag(dir,2), nRecv, 1, MPI_INTEGER,   &
         rankLeft, tag(dir,2), commDir(dir), iStatus, ierror)

    ! Allocate recv buffer
    allocate(particlesToRecv(nRecv))

    ! Communicate
    call MPI_SendRecv(particlesToSend, nSend, MPI_PARTICLE, rankRight, 0, particlesToRecv,   &
         nRecv, MPI_PARTICLE, rankLeft, 0, commDir(dir), iStatus, ierror)

    ! Add to ghost particles
    if (nRecv .gt. 0) then
       call resize_particles(ghostParticles, nGhostParticles+nRecv)
       ghostParticles(nGhostParticles+1:nGhostParticles + nRecv) = particlesToRecv(1:nRecv)
       nGhostParticles = nGhostParticles + nRecv
    end if

    ! Clean up
    deallocate(particlesToSend)
    deallocate(particlesToRecv)

    return
  end subroutine communicate_border


  ! Check if a particle left the domain
  ! -----------------------------------
  subroutine check_particle_bounds(part)

    implicit none

    ! Arguments
    type(t_Particle), intent(inout) :: part

    if (part%stop .eq. 1) return

    select case (nDimensions)

    case(1)

       if (.not. isPeriodic(1)) then
          if (part%position(1) .lt. cartesianCoordinates(1,1) .or. part%position(1) .gt.     &
               cartesianCoordinates(globalGridSize(1),1)) part%stop = 1
       end if

    case(2)

       if (.not. isPeriodic(1)) then
          if (part%position(1) .lt. cartesianCoordinates(1,1) .or. part%position(1) .gt.     &
               cartesianCoordinates(globalGridSize(1),1)) part%stop = 1
       end if

       if (.not. isPeriodic(2)) then
          if (part%position(2) .lt. cartesianCoordinates(1,2) .or. part%position(2) .gt.     &
               cartesianCoordinates(globalGridSize(2),2)) part%stop = 1
       end if

    case(3)

       if (.not. isPeriodic(1)) then
          if (part%position(1) .lt. cartesianCoordinates(1,1) .or. part%position(1) .gt.     &
               cartesianCoordinates(globalGridSize(1),1)) part%stop = 1
       end if

       if (.not. isPeriodic(2)) then
          if (part%position(2) .lt. cartesianCoordinates(1,2) .or. part%position(2) .gt.     &
               cartesianCoordinates(globalGridSize(2),2)) part%stop = 1
       end if

       if (.not. isPeriodic(3)) then
          if (part%position(3) .lt. cartesianCoordinates(1,3) .or. part%position(3) .gt.     &
               cartesianCoordinates(globalGridSize(3),3)) part%stop = 1
       end if

    end select

    return
  end subroutine check_particle_bounds


  ! Localize a particle to the grid
  ! -------------------------------
  subroutine localize_particle(part)

    implicit none

    ! Arguments
    type(t_Particle), intent(inout) :: part

    ! Local variables
    integer :: i, ip

    if (part%stop .eq. 1) return

    do i = 1, nDimensions
       ip = part%gridIndex(i)
       do while (part%position(i) - cartesianCoordinates(ip, i) .lt. 0.0_WP)
          ip = ip - 1
       end do
       do while (cartesianCoordinates(ip + 1, i) - part%position(i) .le. 0.0_WP)
          ip = ip + 1
       end do
       part%gridIndex(i) = ip
    end do

    return
  end subroutine localize_particle


  ! Correct a particle position for periodicity
  ! -------------------------------------------
  subroutine correct_particle_position(part)

    implicit none

    ! Arguments
    type(t_Particle), intent(inout) :: part

    if (part%stop .eq. 1) return

    select case (nDimensions)

    case(1)

       if (isPeriodic(1)) part%position(1) = cartesianCoordinates(1,1) +                     &
            modulo(part%position(1) - cartesianCoordinates(1,1), periodicLength(1))

    case(2)

       if (isPeriodic(1)) part%position(1) = cartesianCoordinates(1,1) +                     &
            modulo(part%position(1) - cartesianCoordinates(1,1), periodicLength(1))

       if (isPeriodic(2)) part%position(2) = cartesianCoordinates(1,2) +                     &
            modulo(part%position(2) - cartesianCoordinates(1,2), periodicLength(2))

    case(3)

       if (isPeriodic(1)) part%position(1) = cartesianCoordinates(1,1) +                     &
            modulo(part%position(1) - cartesianCoordinates(1,1), periodicLength(1))

       if (isPeriodic(2)) part%position(2) = cartesianCoordinates(1,2) +                     &
            modulo(part%position(2) - cartesianCoordinates(1,2), periodicLength(2))

       if (isPeriodic(3)) part%position(3) = cartesianCoordinates(1,3) +                     &
            modulo(part%position(3) - cartesianCoordinates(1,3), periodicLength(3))

    end select

    call localize_particle(part)

    return
  end subroutine correct_particle_position

end module particle_comm


! =============================================== !
! Communicate particles to neighboring processors !
! =============================================== !
subroutine communicate_particles

  ! Internal modules
  use particle_comm

  ! External modules
  use parallel

  implicit none

  ! Local variables
  integer :: i

  ! Recycle
  call recycle_particles(particles)

  ! Update number of particles
  call update_particle_size

  select case (nDimensions)

  case (1)

     call communicate_direction(1)

  case (2)

     call communicate_direction(1)
     call communicate_direction(2)

  case (3)

     call communicate_direction(1)
     call communicate_direction(2)
     call communicate_direction(3)

  end select

  ! Update number of particles
  nParticles = size(particles)
  call parallel_sum(nParticles, nParticlesGlobal)

  ! Correct particle position
  do i = 1, nParticles
     call correct_particle_position(particles(i))
  end do

  return
end subroutine communicate_particles


! ========================================== !
! Communicate ghost particles for collisions !
! ========================================== !
subroutine communicate_ghost_particles

  ! Internal modules
  use particle_comm

  ! External modules
  use parallel

  implicit none

  ! Clean ghost particles structure
  call resize_particles(ghostParticles, 0)
  nGhostParticles = 0

  select case (nDimensions)

  case (1)

     call communicate_border(1)

  case (2)

     call communicate_border(1)
     call communicate_border(2)

  case (3)

     call communicate_border(1)
     call communicate_border(2)
     call communicate_border(3)

  end select

  ! Update number of particles
  if (allocated(ghostParticles)) nGhostParticles = size(ghostParticles)

  return
end subroutine communicate_ghost_particles
