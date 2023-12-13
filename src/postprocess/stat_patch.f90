module stat_patch

  ! External modules
  use stat
  use string
  use grid_patch

  implicit none

    integer, parameter ::                                                                    &
       STAT_PATCH_RHO   = 0,                                                                 &
       STAT_PATCH_U     = 1,                                                                 &
       STAT_PATCH_V     = 2,                                                                 &
       STAT_PATCH_W     = 3,                                                                 &
       STAT_PATCH_P     = 4,                                                                 &
       STAT_PATCH_T     = 5,                                                                 &
       STAT_PATCH_ALPHA = 10,                                                                &
       STAT_PATCH_Y     = 100

  ! Global variables
  integer :: nStats, stat_nvar
  integer, allocatable :: stat_index(:)
  character(len=str_medium), allocatable :: stat_name(:)

  type, private :: t_Stat
     real(WP) :: Delta_t
     real(WP), allocatable :: statVar(:,:)
  end type t_Stat

  type(t_Stat), allocatable :: statData(:)
  type(t_Patch), pointer :: statPatch(:)

contains

  subroutine setup_stat_patch(patch, data)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Stat), intent(inout) :: data

    ! Local variables
    integer :: i, ifile, ierror
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer(kind=MPI_Offset_kind) :: disp
    integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
    integer(kind=MPI_Offset_kind) :: int_MOK, WP_MOK
    integer(kind=MPI_Offset_kind) :: NVARS_MOK
    character(len=str_medium) :: filename, filename_
    logical :: fileIsThere
    real(WP), allocatable :: buffer(:,:)

    ! Verify the patch type
    call verify_stat_patch(patch)

    data%Delta_t = 0.0_WP

    if (patch%comm .eq. MPI_COMM_NULL) return

    allocate(data%statVar(patch%nPatchPoints, 2*stat_nvar))
    data%statVar = 0.0_WP

    ! Create the directory to save the stat patch data
    if (iRank .eq. patch%masterRank) then
       call CREATE_FOLDER('stat/' // trim(patch%name))
    end if
    ! call MPI_BARRIER(patch%comm, ierror)

    ! Write the stat patch coodinates
    filename_ = 'stat/' // trim(patch%name) // '/' // 'grid' !... set the file name

    ! Prefix the filename to prevent potential file locking
    filename = trim(mpiiofs) // trim(filename_)

    ! Open the file to write
    inquire(file = filename, exist = fileIsThere)
    if (fileIsThere) return
    call MPI_FILE_OPEN(patch%comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
         ifile, ierror)
    ! Check if there was an error opening the file
    if (ierror .ne. 0) call die('Error writing ' // trim(filename_))

    ! Write dimensions in the header
    if (iRank .eq. patch%masterRank) then
       call MPI_FILE_WRITE(ifile, patch%globalSize, 3, MPI_INTEGER, status, ierror)
    end if

    ! Resize some integers so MPI can read even the biggest files
    nx_MOK    = int(patch%globalSize(1), MPI_Offset_kind)
    ny_MOK    = int(patch%globalSize(2), MPI_Offset_kind)
    nz_MOK    = int(patch%globalSize(3), MPI_Offset_kind)
    int_MOK   = int(4,                   MPI_Offset_kind)
    WP_MOK    = int(WP,                  MPI_Offset_kind)
    NVARS_MOK = int(nDimensions,         MPI_Offset_kind)

    ! Compute the header displacement
    disp = 3*int_MOK + 0*WP_MOK

    ! Write the coordinates for each dimension
    allocate(buffer(patch%nPatchPoints, nDimensions))
    call patch_collect(patch, coordinates, buffer)
    do i = 1, nDimensions
       call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, patch%mpiRealSubarrayType, "native", &
            mpiInfo, ierror)
       call MPI_FILE_WRITE_ALL(ifile, buffer(:,i), patch%nPatchPoints, MPI_REAL_WP, status,  &
            ierror)

       ! Update the displacement
       disp = disp + nx_MOK*ny_MOK*nz_MOK*WP_MOK
    end do
    deallocate(buffer)

    ! Close the file
    call MPI_FILE_CLOSE(ifile, ierror)

    return
  end subroutine setup_stat_patch


  subroutine cleanup_stat_patch(patch, data)

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Stat), intent(inout) :: data

    if (allocated(data%statVar)) deallocate(data%statVar)

    return
  end subroutine cleanup_stat_patch


  subroutine verify_stat_patch(patch)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch

    ! Local variables
    integer :: i, n, extent(6)
    character(len = str_long) :: message

    if (patch%patchType .ne. STATISTICS_PATCH)                                               &
         call die('verify_stat_patch: Patch type mismatch!')

    if (patch%normalDirection .gt. nDimensions .or. patch%normalDirection .lt. 0)            &
         call die("verify_stat_patch: Normal direction is invalid for '" //                  &
         trim(patch%name) // "'!")

    extent(1) = patch%iMin; extent(2) = patch%iMax
    extent(3) = patch%jMin; extent(4) = patch%jMax
    extent(5) = patch%kMin; extent(6) = patch%kMax

    n = nDimensions
    do i = 1, nDimensions
       if (extent((i-1)*2+1) .lt. 0 .or. extent((i-1)*2+2) .gt. globalGridSize(i) .or.       &
            extent((i-1)*2+1) .gt. extent((i-1)*2+2))                                        &
            call die("verify_stat_patch: Invalid extent on '" // trim(patch%name) // "'!")
       if (extent((i-1)*2+1) .eq. extent((i-1)*2+2)) n = n - 1
    end do

    i = abs(patch%normalDirection)
    if (i .eq. 0) then
       if (n .ne. nDimensions) then
          write(message, '(2(A,I0.0),3A)') "verify_stat_patch: Expected a ", nDimensions,    &
               "D patch, but extent represents a ", n, "D patch on '", trim(patch%name), "'!"
          call die(trim(message))
       end if
    else
       if ((patch%normalDirection .gt. 0 .and. extent((i-1)*2+1) .ne. 1) .or.                &
            (patch%normalDirection .lt. 0 .and. extent((i-1)*2+2) .ne. globalGridSize(i)))   &
            call die("verify_stat_patch: '" // trim(patch%name) //                           &
            "' not aligned with a computational boundary!")
    end if

    return
  end subroutine verify_stat_patch


  subroutine compute_stat_patch(patch, data, timeStepSize)

    ! External madules
    use state

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Stat), intent(inout) :: data
    real(WP), intent(in) :: timeStepSize

    ! Local variables
    integer :: i, j, k, l, gridIndex, patchIndex
    real(WP) :: localVar

     if (patch%comm .eq. MPI_COMM_NULL) return

    data%Delta_t = data%Delta_t + timeStepSize

    do k = patch%iStart(3), patch%iEnd(3)
       do j = patch%iStart(2), patch%iEnd(2)
          do i = patch%iStart(1), patch%iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             patchIndex = i - patch%offset(1) + patch%localSize(1) *                         &
                  (j - 1 - patch%offset(2) + patch%localSize(2) *                            &
                  (k - 1 - patch%offset(3)))

             do l = 1, stat_nvar
                select case (stat_index(l))
                case (STAT_PATCH_RHO)
                   localVar = conservedVariables(gridIndex, 1)
                case (STAT_PATCH_U, STAT_PATCH_V, STAT_PATCH_W)
                   localVar = velocity(gridIndex, stat_index(l))
                case (STAT_PATCH_P)
                   localVar = pressure(gridIndex, 1)
                case (STAT_PATCH_T)
                   localVar = temperature(gridIndex, 1)
                case (STAT_PATCH_ALPHA)
                   localVar = volumeFraction(gridIndex, 1)
                case (STAT_PATCH_Y+1:)
                   localVar = massFraction(gridIndex, stat_index(l) - STAT_PATCH_Y)
                end select

                data%statVar(patchIndex,l) = data%statVar(patchIndex,l)                      &
                     + timeStepSize * localVar
                data%statVar(patchIndex,l+stat_nvar) = data%statVar(patchIndex,l+stat_nvar)  &
                     + timeStepSize * localVar ** 2
             end do
          end do
       end do
    end do

    return
  end subroutine compute_stat_patch


  subroutine read_stat_patch(patch, data, time)

    implicit NONE

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Stat), intent(inout) :: data
    real(WP), intent(in) :: time

    ! Local variables
    integer :: ifile, ierror, nvars, var
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(4) :: dims
    integer(kind=MPI_Offset_kind) :: disp
    integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
    integer(kind=MPI_Offset_kind) :: int_MOK, WP_MOK
    integer(kind=MPI_Offset_kind) :: NVARS_MOK
    character(len=str_medium) :: filename, filename_
    logical :: fileIsThere
    real(WP) :: time_
    real(WP), allocatable :: buffer(:)

    if (patch%comm .eq. MPI_COMM_NULL) return

    ! Set the file name
    write(filename_, '(5A,ES10.4E2)') 'stat/', trim(patch%name), '/', trim(patch%name),      &
         '-T-', time

    ! Check if file exists
    inquire(file = filename_, exist = fileIsThere)
    if (.not. fileIsThere) then
       data%Delta_t = 0.0_WP
       return
    end if

    ! Prefix the filename to prevent potential file locking
    filename = trim(mpiiofs) // trim(filename_) 

    ! Open the file
    call MPI_FILE_OPEN(patch%comm, filename, MPI_MODE_RDONLY, mpiInfo, ifile, ierror)

    ! Check if there was an error opening the file
    if (ierror .ne. 0) call die('Error reading ' // trim(filename_))

    ! Read dimensions from header
    call MPI_FILE_READ_ALL(ifile, dims, 4, MPI_INTEGER, status, ierror)
    if ((dims(1).ne.patch%globalSize(1)) .or. (dims(2).ne.patch%globalSize(2)) .or.          &
         (dims(3).ne.patch%globalSize(3))) then
       print*, 'grid = ', patch%globalSize(1), patch%globalSize(2), patch%globalSize(3)
       print*, trim(filename_), ' = ',dims(1), dims(2), dims(3)
       call die('The size of the' // trim(filename_) //                                      &
            'file does not correspond to the patch grid file')
    end if

    ! Check that number of variables is consistent
    nvars = dims(4)
    if (nvars .ne. 2*stat_nvar) then
       print *, 'Expected ', nvars, ' stat patch variables (mean and squared mean)'
       print *, 'Currently ', 2*stat_nvar, ' stat patch variables (mean and squared mean)'
       call die('The size of the' // trim(filename_) //                                      &
            'file is not consistent with the simulation')
    end if

    ! Read the current Delta_t and time
    call MPI_FILE_READ_ALL(ifile, data%Delta_t, 1, MPI_REAL_WP, status, ierror)
    call MPI_FILE_READ_ALL(ifile, time_, 1, MPI_REAL_WP, status, ierror)

    ! Resize some integers so MPI can read even the biggest files
    nx_MOK    = int(patch%globalSize(1), MPI_Offset_kind)
    ny_MOK    = int(patch%globalSize(2), MPI_Offset_kind)
    nz_MOK    = int(patch%globalSize(3), MPI_Offset_kind)
    int_MOK   = int(4,                   MPI_Offset_kind)
    WP_MOK    = int(WP,                  MPI_Offset_kind)
    NVARS_MOK = int(2*stat_nvar,         MPI_Offset_kind)

    ! Compute the header displacement
    disp = 4*int_MOK + 2*WP_MOK

    ! Read the data for each variable
    allocate(buffer(patch%nPatchPoints))
    do var = 1, 2*stat_nvar
       ! Read the mean value
       call MPI_FILE_SET_VIEW(ifile, disp, MPI_REAL_WP, patch%mpiRealSubarrayType, "native", &
            mpiInfo, ierror)
       call MPI_FILE_READ_ALL(ifile, buffer(:), patch%nPatchPoints, MPI_REAL_WP, status,     &
            ierror)

       ! Update the displacement
       disp = disp + nx_MOK*ny_MOK*nz_MOK*WP_MOK

       ! Compute data%statVar
       data%statVar(:,var) = buffer(:) * data%Delta_t
    end do
    deallocate(buffer)

    ! Close the file
    call MPI_FILE_CLOSE(ifile, ierror)

    return
  end subroutine read_stat_patch


  subroutine write_stat_patch(patch, data, time)

    implicit NONE

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_Stat), intent(inout) :: data
    real(WP), intent(in) :: time

    ! Local variables
    integer :: iunit, ierror, var, nBytes
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(4) :: dims
    integer(kind=MPI_Offset_kind) :: disp
    integer(kind=MPI_Offset_kind) :: nx_MOK, ny_MOK, nz_MOK
    integer(kind=MPI_Offset_kind) :: int_MOK, WP_MOK
    integer(kind=MPI_Offset_kind) :: NVARS_MOK
    character(len=str_medium) :: filename, filename_
    logical :: fileIsThere
    real(WP), allocatable :: buffer(:)

    ! Nothing to write
    if (patch%comm .eq. MPI_COMM_NULL .or. data%Delta_t .le. epsilon(1.0_WP)) return

    ! Set the file name
    write(filename_, '(2A)') 'stat/', trim(patch%name)

    ! Prefix the filename to prevent potential file locking
    filename = trim(mpiiofs) // trim(filename_) 

    ! Open the file to write
    inquire(file = filename, exist = fileIsThere)
    if (fileIsThere .and. iRank .eq. patch%masterRank)                                       &
         call MPI_FILE_DELETE(filename, mpiInfo, ierror)
    call MPI_FILE_OPEN(patch%comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
         iunit, ierror)

    ! Check if there was an error opening the file
    if (ierror .ne. 0) call die('Error writing ' // trim(filename_))

    ! Write header
    if (iRank .eq. patch%masterRank) then
       ! Write dimensions
       dims(1:3) = patch%globalSize
       dims(4) = 2 * stat_nvar
       call MPI_FILE_WRITE(iunit, dims, 4, MPI_INTEGER, status, ierror)
       ! Write the current Delta_t
       call MPI_FILE_WRITE(iunit, data%Delta_t, 1, MPI_REAL_WP, status, ierror)
       ! Write the current time
       call MPI_FILE_WRITE(iunit, time, 1, MPI_REAL_WP, status, ierror)
    end if

    ! Resize some integers so MPI can read even the biggest files
    nx_MOK    = int(patch%globalSize(1), MPI_Offset_kind)
    ny_MOK    = int(patch%globalSize(2), MPI_Offset_kind)
    nz_MOK    = int(patch%globalSize(3), MPI_Offset_kind)
    int_MOK   = int(4,                   MPI_Offset_kind)
    WP_MOK    = int(WP,                  MPI_Offset_kind)
    NVARS_MOK = int(2*stat_nvar,         MPI_Offset_kind)

    ! Compute the header displacement
    disp = 4*int_MOK + 2*WP_MOK

    ! Write the data for each variable
    allocate(buffer(patch%nPatchPoints))
    do var = 1, 2*stat_nvar
       ! Compute the mean value
       buffer(:) = data%statVar(:,var) / data%Delta_t

       ! Write the mean value
       call MPI_FILE_SET_VIEW(iunit, disp, MPI_REAL_WP, patch%mpiRealSubarrayType, "native", &
            mpiInfo, ierror)
       call MPI_FILE_WRITE_ALL(iunit, buffer(:), patch%nPatchPoints, MPI_REAL_WP, status,    &
            ierror)

       ! Update the displacement
       disp = disp + nx_MOK*ny_MOK*nz_MOK*WP_MOK
    end do
    deallocate(buffer)

    ! Close the file
    call MPI_FILE_CLOSE(iunit, ierror)

    ! Write ACII file if under 100 mb
!!$    nBytes = 4*2*stat_nvar*patch%nPatchPoints + 4*nDimensions*patch%nPatchPoints
!!$    if (nBytes .lt. 100.0e6_WP) then
!!$
!!$       ! Root processor writes
!!$       if (iRank .eq. patch%masterRank) then
!!$          iunit = iopen()
!!$          write(filename_, '(3A)') 'stat/', trim(patch%name), '.txt'
!!$          open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
!!$          write(iunit,'(10000a20)') (trim(adjustl(stat_name(var))),var=1,stat_nvar)
!!$          do i = 1, nStatPoints
!!$             write(iunit,'(10000ES20.12)') grid1D(i), buf1D(i,:)
!!$          end do
!!$          close(iclose(iunit))
!!$       end if
!!$
!!$    end if

    return
  end subroutine write_stat_patch

end module stat_patch


! =========================== !
! Initialize patch statistics !
! =========================== !
subroutine stat_patch_setup

  ! Internal modules
  use stat_patch

  ! External modules
  use parser
  use simulation_flags, only : twoWayCoupling
  use solver_options, only: nSpecies, speciesName

  implicit none

  ! Local variables
  integer :: i, j, k, ibuf(3)
  character(len=1) :: U3D(3)

  ! Find the number of stat patches
  nStats = 0; j = 0
  do i = 1, nPatches
     if (patches(i)%patchType .eq. STATISTICS_PATCH) then
        nStats = nStats + 1
        if (j .eq. 0) j = i
     end if
  end do
  if (nStats .eq. 0) return

  ! Allocate the stat data
  allocate(statData(nStats))

  ! Connect the stat patch
  statPatch => patches(j:j+nStats-1)

  ! Compute statistics on the primitive variables (rho, Ui, P, T, Yk, alpha)
  stat_nvar = 3 + nDimensions + nSpecies
  if (twoWayCoupling) stat_nvar = stat_nvar + 1
  allocate(stat_name(2*stat_nvar)) !... for their mean and variance
  allocate(stat_index(stat_nvar))
  stat_name(1) = 'rho'; stat_index(1) = STAT_PATCH_RHO
  U3D = (/ 'U', 'V', 'W' /); ibuf = (/ STAT_PATCH_U, STAT_PATCH_V, STAT_PATCH_W /)
  do i = 1, nDimensions
     stat_name(i+1) = U3D(i)
     stat_index(i+1) = ibuf(i)
  end do
  stat_name(nDimensions+2:nDimensions+3) = (/ 'P', 'T' /)
  stat_index(nDimensions+2:nDimensions+3) = (/ STAT_PATCH_P, STAT_PATCH_T /)
  do k = 1, nSpecies
     stat_name(nDimensions+3+k) = trim(speciesName(k))
     stat_index(nDimensions+3+k) = STAT_PATCH_Y + k
  end do
  if (twoWayCoupling) then
     stat_name(nDimensions+nSpecies+4) = 'alpha'
     stat_index(nDimensions+nSpecies+4) = STAT_PATCH_ALPHA
  end if

  ! Update stat name for mean and variance
  do i = 1, stat_nvar
     stat_name(i) = trim(stat_name(i)) // '_MEAN'
     stat_name(i+stat_nvar) = trim(stat_name(i)) // '_SQ_MEAN'
  end do

  ! Setup the stat conditions
  do i = 1, nStats
     call setup_stat_patch(statPatch(i), statData(i))
  end do

  ! Read initial stat if exists
  call stat_patch_read

  return
end subroutine stat_patch_setup


! ======================== !
! Cleanup patch statistics !
! ======================== !
subroutine stat_patch_cleanup

  ! Internal modules
  use stat_patch

  implicit none

  ! Local variables
  integer :: i

  if (nStats .gt. 0) then
     deallocate(stat_name)
     deallocate(stat_index)
     do i = 1, nStats
        call cleanup_stat_patch(statPatch(i), statData(i))
     end do
     deallocate(statData)
     nullify(statPatch)
  end if

  nStats = 0

  return
end subroutine stat_patch_cleanup


! ======================== !
! Compute patch statistics !
! ======================== !
subroutine stat_patch_compute

  ! Internal modules
  use stat_patch

  ! External modules
  use time_info

  implicit none

  ! Local variables
  integer  :: i

  do i = 1, nStats
     call compute_stat_patch(statPatch(i), statData(i), timeStepSize)
  end do

  return
end subroutine stat_patch_compute


! ===================== !
! Read patch statistics !
! ===================== !
subroutine stat_patch_read

  ! Internal modules
  use stat_patch

  ! External modules
  use time_info

  implicit none

  ! Local variables
  integer :: i

  do i = 1, nStats
     call read_stat_patch(statPatch(i), statData(i), time)
  end do

  return
end subroutine stat_patch_read


! ====================== !
! Write patch statistics !
! ====================== !
subroutine stat_patch_write

  ! Internal modules
  use stat_patch

  ! External madules
  use time_info

  implicit none

  ! Local variables
  integer  :: i

  do i = 1, nStats
     call write_stat_patch(statPatch(i), statData(i), time)
  end do

  ! Log
  call monitor_log("STATISTICS PATCH FILE WRITTEN")

  return
end subroutine stat_patch_write
