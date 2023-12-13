module particle_io

  ! External modules
  use particle
  use parallel

  implicit none

contains

  ! Communicate particles globally - only need to do this during initial read
  ! -------------------------------------------------------------------------
  subroutine communicate_globally

    ! External modules
    use precision
    use geometry
    use grid
    use math, only : bisection

    implicit none

    ! Local variables
    type(t_Particle), allocatable :: particleMatrix(:,:)
    integer :: i, j, p, gridIndex, request(2), rankSend, rankRecv, nParticlesSend,           &
         nParticlesRecv, nParticlesOld, localExtents(3,2), ierror
    integer, dimension(MPI_STATUS_SIZE, 2) :: statuses
    integer, dimension(3) :: ijk
    integer, allocatable :: getProcRank(:,:,:), particlesToSend(:), particlesToRecv(:),      &
         counter(:), particleRank(:), procExtents(:,:,:)
    real(WP), allocatable :: cartesianCoordinates(:,:)

    ! Prepare Cartesian vectors for particle localization
    allocate(cartesianCoordinates(maxval(globalGridSize), nDimensions))
    cartesianCoordinates = -huge(0.0_WP)
    do j = 1, nDimensions
       do i = iStart(j), iEnd(j)
          ijk = iStart
          ijk(j) = i
          gridIndex = ijk(1) - gridOffset(1) + localGridSize(1) *                            &
               (ijk(2) - 1 - gridOffset(2) + localGridSize(2) *                              &
               (ijk(3) - 1 - gridOffset(3)))
          cartesianCoordinates(i, j) = coordinates(gridIndex, j)
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, cartesianCoordinates, maxval(globalGridSize) *          &
         nDimensions, MPI_REAL_WP, MPI_MAX, comm, ierror)

    ! Locate the closest grid point
    do p = 1, size(particles)
       ijk = iStart
       do j = 1, nDimensions
          call bisection(particles(p)%position(j), ijk(j),                                   &
               cartesianCoordinates(1:globalGridSize(j), j), 1, globalGridSize(j))
       end do
       particles(p)%gridIndex = ijk
    end do

    ! Get the local Cartesian extents
    localExtents(:, 1) = iStart
    localExtents(:, 2) = iEnd

    ! Get the Cartesian processor extents
    allocate(procExtents(maxval(nProcsDir), 3, 2))
    procExtents = 0
    do i = 1, nDimensions
       call MPI_Allgather(localExtents(i,1), 1, MPI_INTEGER,                                 &
            procExtents(1:nProcsDir(i), i, 1), 1, MPI_INTEGER, commDir(i), ierror)
       call MPI_Allgather(localExtents(i,2), 1, MPI_INTEGER,                                 &
            procExtents(1:nProcsDir(i), i, 2), 1, MPI_INTEGER, commDir(i), ierror)
    end do

    ! Get processor rank associated with the Cartesian coordinates
    allocate(getProcRank(nProcsDir(1), nProcsDir(2), nProcsDir(3))); getProcRank = -1
    ijk = 1
    ijk(1:nDimensions) = procCoords(1:nDimensions) + 1
    getProcRank(ijk(1), ijk(2), ijk(3)) = iRank
    call MPI_Allreduce(MPI_IN_PLACE, getProcRank, product(nProcsDir), MPI_INTEGER, MPI_MAX,  &
         comm, ierror)

    ! Recycle
    call recycle_particles(particles)

    ! Prepare information about who sends what to whom
    allocate(particlesToSend(nProcs))
    allocate(particleRank(size(particles)))
    particlesToSend = 0
    do p = 1, size(particles)
       ijk = 1
       do j = 1, nDimensions
          innerLoop: do i = 1, nProcsDir(j)
             if (particles(p)%gridIndex(j) .ge. procExtents(i,j,1) .and.                     &
                  particles(p)%gridIndex(j) .le. procExtents(i,j,2)) then
                ijk(j) = i
                exit innerLoop
             end if
          end do innerLoop
       end do
       particleRank(p) = getProcRank(ijk(1), ijk(2), ijk(3))
       particlesToSend(particleRank(p)+1) = particlesToSend(particleRank(p)+1) + 1

    end do
    particlesToSend(iRank+1) = 0 ! ... remove the diagonal

    ! Prepare information about who receives what from whom
    allocate(particlesToRecv(nProcs))
    do rankRecv = 0, nProcs - 1
       call MPI_Gather(particlesToSend(rankRecv+1), 1, MPI_INTEGER, particlesToRecv, 1,      &
            MPI_INTEGER, rankRecv, comm, ierror)
    end do

    ! Prepare the buffers
    nParticlesSend = maxval(particlesToSend)
    nParticlesRecv = sum(particlesToRecv)

    ! Allocate buffers to send particles
    allocate(particleMatrix(nParticlesSend, nProcs))

    ! Prepare the counter
    allocate(counter(nProcs))
    counter = 0
    do p = 1, size(particles)
       if (particleRank(p) .ne. iRank) then
          ! Prepare for sending.
          counter(particleRank(p)+1) = counter(particleRank(p)+1) + 1
          particleMatrix(counter(particleRank(p)+1), particleRank(p)+1) = particles(p)
          ! Need to remove the particle later.
          particles(p)%stop = 1
       end if
    end do

    ! Everybody resizes
    nParticlesOld = size(particles)
    call resize_particles(particles, nParticlesOld + nParticlesRecv)

    ! Loop through the CPUs, pack particles in `particleMatrix`, send, unpack
    request = MPI_REQUEST_NULL
    do rankSend = 0, nProcs - 1
       if (iRank .eq. rankSend) then
          ! Non-blocking send.
          do rankRecv = 0, nProcs - 1
             if (particlesToSend(rankRecv+1) .gt. 0) then
                call MPI_Isend(particleMatrix(:,rankRecv+1), particlesToSend(rankRecv+1),    &
                     MPI_PARTICLE, rankRecv, 0, comm, request(1), ierror)
             end if
          end do
       else
          ! Non-blocking receive.
          if (particlesToRecv(rankSend+1) .gt. 0) then
             call MPI_Irecv(particles(nParticlesOld + 1 : nParticlesOld +                    &
                  particlesToRecv(rankSend+1)),  particlesToRecv(rankSend+1),                &
                  MPI_PARTICLE, rankSend, 0, comm, request(2), ierror)
             nParticlesOld = nParticlesOld + particlesToRecv(rankSend+1)
          end if
       end if
       call MPI_Waitall(2, request, statuses, ierror)
    end do

    ! Recycle
    call recycle_particles(particles)

    ! Update the particle size
    nParticles = size(particles)
    call parallel_sum(nParticles, nParticlesGlobal)

    ! Cleanup
    deallocate(cartesianCoordinates)
    deallocate(counter)
    deallocate(particleRank)
    deallocate(particleMatrix)
    deallocate(particlesToSend)
    deallocate(particlesToRecv)
    deallocate(getProcRank)
    deallocate(procExtents)

    return
  end subroutine communicate_globally

end module particle_io


! ============================== !
! Read particle data in parallel !
! ============================== !
subroutine particle_read_parallel(filename_)

  ! Internal modules
  use particle_io

  ! External modules
  use math
  use geometry
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, particleSize, ifile, ierror
  integer (kind=MPI_OFFSET_KIND) :: offset
  real(WP) :: buf
  character(len=str_medium) :: filename, message

  ! Prefix the filename to prevent potential file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Open the file
  call MPI_File_open(comm, trim(filename), MPI_MODE_RDONLY, mpiInfo, ifile, ierror)

  ! Check if there was an error opening the file
  if (ierror.ne.0) call die('Error reading ' // trim(filename_))

  ! Number of particles
  call MPI_File_read_all(ifile, nParticlesGlobal, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)

  ! Particle size
  call MPI_File_read_all(ifile, particleSize, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)

  if (particleSize .ne. SIZE_MPI_PARTICLE) then
     write(message, '(3A,I0.0,A,I0.0)') "Problem reading '", trim(filename),                 &
          "': Particle size = ", particleSize, ", expected size = ", SIZE_MPI_PARTICLE
     call die(trim(message))
  end if

  ! Time
  call MPI_File_read_all(ifile, buf, 1, MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)

  ! Distribute the particles evenly among the processors
  call pigeon_hole(nParticlesGlobal, nProcs, iRank, i, nParticles)
  if (allocated(particles)) deallocate(particles)
  allocate(particles(nParticles))

  ! Create offset
  offset = 4*2 + int(WP, MPI_OFFSET_KIND) + int(i, MPI_OFFSET_KIND) *                        &
       int(SIZE_MPI_PARTICLE, MPI_OFFSET_KIND)

  if (nParticles .gt. 0) call MPI_File_read_at(ifile, offset, particles, nParticles,         &
       MPI_PARTICLE, MPI_STATUS_IGNORE, ierror)

  ! Close the file
  call MPI_File_close(ifile, ierror)

  ! Communicate particles among the processors
  call communicate_globally

  return
end subroutine particle_read_parallel


! ============================ !
! Read particle data in serial !
! ============================ !
subroutine particle_read_serial(filename_)

  ! Internal modules
  use particle_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ibuf, ifile, ierror
  integer, dimension(3) :: dims
  character(len=str_medium) :: filename, message
  logical :: fileIsThere

  ! Get the name of the header file to read in
  filename = trim(filename_) // '/part.header'

  ! If file is not there revert to parallel read
  inquire(file = filename, exist = fileIsThere)
  if (.not. fileIsThere) then
     call particle_read_parallel(trim(filename_))
     return
  end if

  ! Root reads the header file
  if (irank.eq.iroot) then
     ! Open the file
     call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
     ! Read the particle header
     call BINARY_FILE_READ(ifile, dims, 3, kind(dims), ierror)
     ! Check number of processors is consistent
     if (dims(1).ne.nProcs) then
        print*, 'Expected ', dims(1), ' processor(s)'
        print*, 'Currently ', nProcs, ' processor(s)'
        call die('Number of processors incompatible with serial particle file')
     end if

     ! Set number of particles
     nParticlesGlobal = dims(2)

     ! Check particle size
     if (dims(3) .ne. SIZE_MPI_PARTICLE) then
        write(message, '(3A,I0.0,A,I0.0)') "Problem reading '", trim(filename),               &
             "': Particle size = ", dims(3), ", expected size = ", SIZE_MPI_PARTICLE
        call die(trim(message))
     end if
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Broadcast information
  call parallel_bc(nParticlesGlobal)

  ! Build actual filename to read for each processor
  filename = trim(filename_) // '/part.'
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open each file and the number of particles for this processor
  call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)
  call BINARY_FILE_READ(ifile, nParticles, 1, kind(nParticles), ierror)

  ! Allocate particle vector and read
  if (allocated(particles)) deallocate(particles)
  allocate(particles(nParticles))
  if (nParticles.gt.0) call binary_file_read(ifile, particles(1:nParticles), nParticles,     &
       SIZE_MPI_PARTICLE, ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  ! Safety check
  call parallel_sum(nParticles, ibuf)
  if (ibuf.ne.nParticlesGlobal)                                                              &
       call die('Total number of particles read from ' // trim(filename_) // ' is incorrect.')

  ! Communicate particles among the processors
  call communicate_globally

  return
end subroutine particle_read_serial


! =============================== !
! Write particle data in parallel !
! =============================== !
subroutine particle_write_parallel(filename_)

  ! Internal modules
  use particle_io

  ! External modules
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: i, ifile, ierror
  integer, allocatable :: nParticlesProc(:)
  integer (kind=MPI_OFFSET_KIND) :: offset
  character(len=str_medium) :: filename
  logical :: fileIsThere

  ! Prefix the filename to prevent potentials file locking
  filename = trim(mpiiofs) // trim(filename_)

  ! Get number of particles per processor in `comm`
  allocate(nParticlesProc(nProcs))
  call MPI_Allgather(nParticles, 1, MPI_INTEGER, nParticlesProc, 1, MPI_INTEGER,             &
       comm, ierror)

  ! Open the file to write
  inquire(file = filename, exist = fileIsThere)
  if (fileIsThere .and. irank.eq.iroot) call MPI_FILE_DELETE(filename, mpiInfo, ierror)
  call MPI_FILE_OPEN(comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo,         &
       ifile, ierror)

  ! Root process writes the header
  if (iRank .eq. iRoot) then
     call MPI_FILE_WRITE(ifile, nParticlesGlobal, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
     call MPI_FILE_WRITE(ifile, SIZE_MPI_PARTICLE, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
     call MPI_FILE_WRITE(ifile, time, 1, MPI_REAL_WP, MPI_STATUS_IGNORE, ierror)
  end if

  ! Get the local particle offset
  offset = 4*2 + int(WP, MPI_OFFSET_KIND)
  do i = 1, iRank
     offset = offset + int(nParticlesProc(i), MPI_OFFSET_KIND) *                             &
          int(SIZE_MPI_PARTICLE, MPI_OFFSET_KIND)
  end do

  ! Write local particle data
  if (nParticles .gt. 0) call MPI_file_write_at(ifile, offset, particles, nParticles,        &
       MPI_PARTICLE, MPI_STATUS_IGNORE, ierror)

  ! Close the file
  call MPI_File_close(ifile, ierror)

  ! Cleanup
  deallocate(nParticlesProc)

  return
end subroutine particle_write_parallel


! ============================= !
! Write particle data in serial !
! ============================= !
subroutine particle_write_serial(filename_)

  ! Internal modules
  use particle_io

  ! External modules
  use fileio
  use geometry
  use solver_options
  use time_info

  implicit none

  ! Arguments
  character(len = *), intent(in) :: filename_

  ! Local variables
  integer :: ifile, ierror
  integer, dimension(3) :: dims
  character(len=str_medium) :: filename, dirname

  ! Set the name of the directory
  dirname = trim(adjustl(filename_))

  ! Create directory
  if (irank.eq.iroot) call CREATE_FOLDER(trim(dirname))
  call MPI_BARRIER(comm, ierror)

  ! Root dumps header file
  if (irank.eq.iroot) then
     ! Create header name
     filename = trim(adjustl(dirname)) // "/" // "part.header"
     ! Open file to write
     call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)
     ! Write common integer info
     dims(1) = nProcs
     dims(2) = nParticlesGlobal
     dims(3) = SIZE_MPI_PARTICLE
     call BINARY_FILE_WRITE(ifile, dims, 3, kind(dims), ierror)
     ! Add time info
     call BINARY_FILE_WRITE(ifile, time, 1, kind(time), ierror)
     ! Close the file
     call BINARY_FILE_CLOSE(ifile, ierror)
  end if

  ! Create filename
  filename = trim(adjustl(dirname)) // "/part."
  write(filename(len_trim(filename)+1:len_trim(filename)+6),'(i6.6)') irank

  ! Open the file to write
  call BINARY_FILE_OPEN(ifile, trim(adjustl(filename)), "w", ierror)

  ! Write header
  call BINARY_FILE_WRITE(ifile, nParticles, 1, kind(nParticles), ierror)

  ! Write the particle data
  if (nParticles.gt.0) call BINARY_FILE_WRITE(ifile, particles(1:nParticles), nParticles,    &
       SIZE_MPI_PARTICLE, ierror)

  ! Close the file
  call BINARY_FILE_CLOSE(ifile, ierror)

  return
end subroutine particle_write_serial
