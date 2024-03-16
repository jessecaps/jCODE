module dump_ensight

  ! External modules
  use string
  use precision
  use dump_viz

  implicit none

  ! Time info
  integer :: nOutputTimes
  real(WP), allocatable :: outputTimes(:)
  
  ! Patch-specific EnSight data
  type, private :: t_EnsightData
     character(len = str_medium) :: directory
     integer :: subarrayType
     real(WP), allocatable :: buffer_WP(:,:)
     real(SP), allocatable :: buffer_SP(:,:)
  end type t_EnsightData
  type(t_EnsightData), pointer :: ensightData(:)

  ! Optional variables to dump
  logical :: dumpPrimitiveVariables, dumpQcrit, dumpVorticity, dumpCPU, dumpHeatRelease,     &
       dumpIgnition, dumpParticles, dumpAllParticles, dumpSensitivity, dumpTargetMollifier,  &
       dumpControlMollifier, dumpDilatation, dumpJacobian, dumpGridNorm, dumpGridSpacing,    &
       dumpViscosity, dumpArcLength, dumpMachNumber, dumpSchlieren, dumpGridIndex, dumpIBM,  &
       dumpParticleAcceleration, dumpParticleCollision, dumpPressureStrain, dumpShadowgraph, &
       dumpBtorque, dumpIBMparticles

  ! Additional info
  integer :: schlierenDirection, schlierenComponent
  logical, allocatable :: dumpParticleToPatch(:)

contains

  ! Setup the individual visualization patch
  ! ----------------------------------------
  subroutine setup_ensight_patch(patch, data)

    ! External modules
    use parser
    use parallel
    use geometry
    use grid
    use time_info

    implicit none

    ! Arguments
    type(t_Patch), intent(inout) :: patch
    type(t_EnsightData), intent(inout) :: data

    ! Local variables
    integer :: extent(6), offset(3), ierror

    ! Return if current process does not own visualization patch
    if (patch%nPatchPoints .le. 0) return

    ! Allocate buffer arrays
    if (.not. allocated(data%buffer_WP)) allocate(data%buffer_WP(patch%nPatchPoints, 3))
    if (.not. allocated(data%buffer_SP)) allocate(data%buffer_SP(patch%nPatchPoints, 3))

    ! Setup single precision MPI subararry type
    extent = (/ patch%iMin, patch%iMax, patch%jMin, patch%jMax, patch%kMin, patch%kMax /)
    offset = patch%offset - extent(1::2) + 1
    call MPI_Type_create_subarray(3, patch%globalSize, patch%localSize, offset,              &
         MPI_ORDER_FORTRAN, MPI_REAL, data%subarrayType, ierror)
    call MPI_Type_commit(data%subarrayType, ierror)

    ! Store the directory to write to
    if (nVizPatches .eq. 1) then
       data%directory = "ensight-3D"
    else
       data%directory = "ensight-3D/" // trim(adjustl(patch%name))
    end if

    return
  end subroutine setup_ensight_patch


  ! Cleanup the EnSight patch data
  ! ------------------------------
  subroutine cleanup_ensight_patch(data)

    ! Arguments
    type(t_EnsightData), intent(inout) :: data

    if (allocated(data%buffer_WP)) deallocate(data%buffer_WP)
    if (allocated(data%buffer_SP)) deallocate(data%buffer_SP)

    return
  end subroutine cleanup_ensight_patch


  ! Dump binary EnSight gold data - case file
  ! -----------------------------------------
  subroutine dump_ensight_case

    ! External modules
    use fileio
    use parallel
    use simulation_flags
    use time_info
    use dissipation
    use particle_exchange
    use combustion
    use ignition_source
    use ibm

    implicit none

    ! Local variables
    integer :: i, j, iunit, ierror
    real(WP), dimension(:), allocatable :: buffer
    character(len=80) :: str
    character(len = str_medium) :: name
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data

    ! Update the time info
    if (nOutputTimes .gt. 0) then
       allocate(buffer(nOutputTimes))
       buffer(1:nOutputTimes) = outputTimes(1:nOutputTimes)
       deallocate(outputTimes)
       nOutputTimes = nOutputTimes + 1
       allocate(outputTimes(nOutputTimes))
       outputTimes(1:nOutputTimes - 1) = buffer(1:nOutputTimes - 1)
       deallocate(buffer)
       outputTimes(nOutputTimes) = time
    else
       deallocate(outputTimes)
       nOutputTimes = 1
       allocate(outputTimes(nOutputTimes))
       outputTimes(nOutputTimes) = time
    end if

    ! Update the case file for each patch
    do j = 1, nVizPatches

       ! Pull off the data
       patch => vizPatch(j)
       data => ensightData(j)

       ! Master rank associated with the patch writes (in ASCII)
       if (iRank .ne. patch%masterRank) cycle

       ! Open the file
       iunit = iopen()
       open(iunit, file = trim(adjustl(data%directory))//"/"//trim(adjustl(patch%name))//    &
            ".case", form = "formatted", iostat = ierror, status="REPLACE")

       ! Write the case
       str = 'FORMAT'
       write(iunit, '(a80)') str
       str = 'type: ensight gold'
       write(iunit, '(a80)') str
       str = 'GEOMETRY'
       write(iunit, '(a80)') str
       str= 'model: geometry'
       write(iunit, '(a80)') str
       str = 'measured: 1 PARTICLES/PARTICLES.******'
       if (dumpParticleToPatch(j) .or. dumpIBMparticles.and.j.eq.1) write(iunit,'(a80)') str
       str = 'VARIABLE'
       write(iunit, '(a80)') str

       ! CPU
       if (dumpCPU) then
          str='scalar per node: 1 CPU CPU/CPU.******'
          write(iunit,'(a80)') str
       end if

       ! Primitive variables
       if (dumpPrimitiveVariables) then
          ! Density
          str='scalar per node: 1 DENSITY DENSITY/DENSITY.******'
          write(iunit,'(a80)') str

          ! Velocity
          str='vector per node: 1 VELOCITY VELOCITY/VELOCITY.******'
          write(iunit,'(a80)') str

          ! Temperature
          str='scalar per node: 1 TEMPERATURE TEMPERATURE/TEMPERATURE.******'
          write(iunit,'(a80)') str

          ! Pressure
          str='scalar per node: 1 PRESSURE PRESSURE/PRESSURE.******'
          write(iunit,'(a80)') str

          ! Mass fraction
          if (nSpecies .gt. 0) then
             do i = 1, nSpecies
                write(name, "(A,I2.2)") "MASS_FRACTION_", i
                str = 'scalar per node: 1 ' // trim(adjustl(name)) // ' ' //                 &
                     trim(adjustl(name)) // '/' // trim(adjustl(name))//'.******'
                write(iunit,'(a80)') str
             end do
          end if
       end if

       ! Viscosity
       if (dumpViscosity) then
          str='scalar per node: 1 VISC VISC/VISC.******'
          write(iunit,'(a80)') str

          ! Artificial (and turbulent) viscosities
          if (useLES) then
             str='scalar per node: 1 VISC_T VISC_T/VISC_T.******'
             write(iunit,'(a80)') str
          end if
          if (useShockCapturing) then
             str='scalar per node: 1 VISC_SHOCK VISC_SHOCK/VISC_SHOCK.******'
             write(iunit,'(a80)') str
             str='scalar per node: 1 BULK_VISC_SHOCK BULK_VISC_SHOCK/BULK_VISC_SHOCK.******'
             write(iunit,'(a80)') str
             str='scalar per node: 1 DIFF_SHOCK DIFF_SHOCK/DIFF_SHOCK.******'
             write(iunit,'(a80)') str       
          end if
       end if

       ! Dissipation sensor
       if (hybridDissipation) then
          str='scalar per node: 1 SENSOR SENSOR/SENSOR.******'
          write(iunit,'(a80)') str
       end if

       ! Q-criterion
       if (dumpQcrit) then
          str='scalar per node: 1 QCRIT QCRIT/QCRIT.******'
          write(iunit,'(a80)') str
       end if

       ! Vorticity
       if (dumpVorticity) then
          if (nDimensions.eq.3) then
             str='vector per node: 1 VORT VORT/VORT.******'
          else
             str='scalar per node: 1 VORT VORT/VORT.******'
          end if
          write(iunit,'(a80)') str
       end if

       ! Dilatation
       if (dumpDilatation) then
          str='scalar per node: 1 DIVU DIVU/DIVU.******'
          write(iunit,'(a80)') str
       end if

       ! Baroclinic torque
       if (dumpBtorque) then
          if (nDimensions.eq.3) then
             str='vector per node: 1 B_TORQUE B_TORQUE/B_TORQUE.******'
          else
             str='scalar per node: 1 B_TORQUE B_TORQUE/B_TORQUE.******'
          end if
          write(iunit,'(a80)') str
       end if
       
       ! Mach number
       if (dumpMachNumber) then
          str='scalar per node: 1 MACH MACH/MACH.******'
          write(iunit,'(a80)') str
       end if

       ! Numerical schlieren
       if (dumpSchlieren) then
          str='scalar per node: 1 SCHLIEREN SCHLIEREN/SCHLIEREN.******'
          write(iunit,'(a80)') str
       end if

       ! ShadowGraph
       if (dumpShadowgraph) then
          str='scalar per node: 1 SHADOWGRAPH SHADOWGRAPH/SHADOWGRAPH.******'
          write(iunit,'(a80)') str
       end if

       ! Pressure-strain correlation <p' div(u')>
       if (dumpPressureStrain) then
          str='scalar per node: 1 PSTRAIN PSTRAIN/PSTRAIN.******'
          write(iunit,'(a80)') str
       end if

       ! Combustion
       if (dumpHeatRelease .and. nReactions .gt. 0) then
          str='scalar per node: 1 HEAT_RELEASE HEAT_RELEASE/HEAT_RELEASE.******'
          write(iunit,'(a80)') str
       end if

       ! Ignition
       if (dumpIgnition) then
          str='scalar per node: 1 IGNITION_SOURCE IGNITION_SOURCE/IGNITION_SOURCE.******'
          write(iunit,'(a80)') str
       end if

       ! Immersed boundary
       if (dumpIBM) then
          ! Levelset
          str='scalar per node: 1 LEVELSET LEVELSET/LEVELSET.******'
          write(iunit,'(a80)') str

          ! Object index
          str='scalar per node: 1 OBJ_ID OBJ_ID/OBJ_ID.******'
          write(iunit,'(a80)') str
          
          if (ibm_move) then
             str='vector per node: 1 IBM_VEL IBM_VEL/IBM_VEL.******'
             write(iunit,'(a80)') str
          end if
       end if

       ! Jacobian
       if (dumpJacobian) then
          str='scalar per node: 1 JAC JAC/JAC.******'
          write(iunit,'(a80)') str
       end if

       ! Grid index
       if (dumpGridIndex) then
          str='vector per node: 1 IJK IJK/IJK.******'
          write(iunit,'(a80)') str
       end if

       ! Grid spacing
       if (dumpGridSpacing) then
          str='vector per node: 1 DX DX/DX.******'
          write(iunit,'(a80)') str
       end if

       ! Grid norm
       if (dumpGridNorm) then
          str='scalar per node: 1 NORM NORM/NORM.******'
          write(iunit,'(a80)') str
       end if

       ! Arc length
       if (dumpArcLength .and. allocated(arcLengths)) then
          str='vector per node: 1 ARC ARC/ARC.******'
          write(iunit,'(a80)') str
       end if

       ! Target mollifier
       if (dumpTargetMollifier .and. .not. predictionOnly) then
          str='scalar per node: 1 TARGET_MOLLIFIER TARGET_MOLLIFIER/TARGET_MOLLIFIER.******'
          write(iunit,'(a80)') str
       end if

       ! Particles
       if (useParticles) then
          str='scalar per node: 1 VOLUME_FRACTION VOLUME_FRACTION/VOLUME_FRACTION.******'
          write(iunit,'(a80)') str
          if (useGranularTemperature) then
             str='scalar per node: 1 GRANULAR_TEMP GRANULAR_TEMP/GRANULAR_TEMP.******'
             write(iunit,'(a80)') str
          end if
       end if

       if (dumpParticleToPatch(j)) then
          str='scalar per measured node: 1 DIAMETER PARTICLES/DIAMETER.******'
          write(iunit,'(a80)') str
          str='scalar per measured node: 1 Tp PARTICLES/Tp.******'
          write(iunit,'(a80)') str
          str='vector per measured node: 1 Vp PARTICLES/Vp.******'
          write(iunit,'(a80)') str
          str='scalar per measured node: 1 ID PARTICLES/ID.******'
          write(iunit,'(a80)') str
          if (useFriction) then
             str='vector per measured node: 1 OMEGAp PARTICLES/OMEGAp.******'
             write(iunit,'(a80)') str
          end if
          if (dumpParticleAcceleration) then
             str='vector per measured node: 1 Ap PARTICLES/Ap.******'
             write(iunit,'(a80)') str
          end if
          if (dumpParticleCollision) then
             str='vector per measured node: 1 COLp PARTICLES/COLp.******'
             write(iunit,'(a80)') str
          end if
       end if

       if (dumpIBMparticles .and. j.eq.1) then
          str='scalar per measured node: 1 DIAMETER PARTICLES/DIAMETER.******'
          write(iunit,'(a80)') str
          str='vector per measured node: 1 Vp PARTICLES/Vp.******'
          write(iunit,'(a80)') str
          str='vector per measured node: 1 OMEGAp PARTICLES/OMEGAp.******'
          write(iunit,'(a80)') str
       end if

       ! Time section
       str='TIME'
       write(iunit,'(a80)') str
       str='time set: 1'
       write(iunit,'(a80)') str
       str='number of steps:'
       write(iunit,'(a16,x,i12)') str, nOutputTimes
       str='filename start number: 1'
       write(iunit,'(a80)') str
       str='filename increment: 1'
       write(iunit,'(a80)') str
       str='time values:'
       write(iunit,'(a12,x,10000000(3(ES14.7,x),/))') str, outputTimes

       ! Close the file
       close(iclose(iunit))

    end do

    return
  end subroutine dump_ensight_case


  ! Dump binary EnSight gold data - adjoint case file
  ! -------------------------------------------------
  subroutine dump_ensight_case_adjoint

    ! External modules
    use fileio
    use parallel
    use simulation_flags
    use solver_options
    use time_info

    implicit none

    ! Local variables
    integer :: i, j, iunit, ierror
    real(WP), dimension(:), allocatable :: buffer
    character(len=80) :: str
    character(len = str_medium) :: name
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data

    ! Update the time info
    if (nOutputTimes .gt. 0) then
       allocate(buffer(nOutputTimes))
       buffer(1:nOutputTimes) = outputTimes(1:nOutputTimes)
       deallocate(outputTimes)
       nOutputTimes = nOutputTimes + 1
       allocate(outputTimes(nOutputTimes))
       outputTimes(1) = time
       outputTimes(2:nOutputTimes) = buffer(1:nOutputTimes - 1)
       deallocate(buffer)
    else
       deallocate(outputTimes)
       nOutputTimes = 1
       allocate(outputTimes(nOutputTimes))
       outputTimes(nOutputTimes) = time
    end if

    ! Update the case file for each patch
    do j = 1, nVizPatches

       ! Pull off the data
       patch => vizPatch(j)
       data => ensightData(j)

       ! Master rank associated with the patch writes (in ASCII)
       if (iRank .ne. patch%masterRank) cycle

       ! Open the file
       iunit = iopen()
       open(iunit, file = trim(adjustl(data%directory))//"/"//trim(adjustl(patch%name))//    &
            "_adjoint.case", form = "formatted", iostat = ierror, status="REPLACE")

       ! Write the case
       str='FORMAT'
       write(iunit, '(a80)') str
       str='type: ensight gold'
       write(iunit, '(a80)') str
       str='GEOMETRY'
       write(iunit, '(a80)') str
       str= 'model: geometry'
       write(iunit, '(a80)') str
       str='VARIABLE'
       write(iunit, '(a80)') str

       ! Adjoint variables
       do i = 1, nUnknowns
          write(name, "(A,I2.2)") 'ADJOINT_', i
          write(str, "(7A)") 'scalar per node: 1 ', trim(name), ' ', trim(name), '/',        &
               trim(name), '.******'
          write(iunit,'(a80)') str
       end do

       ! Sensitivity gradient
       if (dumpSensitivity) then
          str='scalar per node: 1 SENSITIVITY SENSITIVITY/SENSITIVITY.******'
          write(iunit,'(a80)') str
       end if

       ! Control mollifier
       if (dumpControlMollifier) then
          str='scalar per node: 1 CONTROL_MOLLIFIER CONTROL_MOLLIFIER/CONTROL_MOLLIFIER.******'
          write(iunit,'(a80)') str
       end if

       ! Time section
       str='TIME'
       write(iunit,'(a80)') str
       str='time set: 1'
       write(iunit,'(a80)') str
       str='number of steps:'
       write(iunit,'(a16, i12)') str, nOutputTimes
       str='filename start number: '
       write(iunit,'(a23, i5)') str, nOutputTimes
       str='filename increment: -1'
       write(iunit,'(a80)') str
       str='time values:'
       write(iunit,'(a12,x,10000000(3(ES12.5,x),/))') str, outputTimes

       ! Close the file
       close(iclose(iunit))

    end do

    return
  end subroutine dump_ensight_case_adjoint


  ! Dump binary EnSight gold data - geometry
  ! ----------------------------------------
  subroutine dump_ensight_geometry(patch, data)

    ! External modules
    use parallel
    use grid

    implicit none

    ! Arguments
    type(t_Patch), intent(in) :: patch
    type(t_EnsightData), intent(inout) :: data  

    ! Local variables
    integer :: i, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = 80) :: filename, buffer
    logical :: fileExists

    ! Return if current process does not own visualization patch
    if (patch%nPatchPoints .le. 0) return

    ! Open the file
    filename = trim(adjustl(data%directory)) // "/geometry"
    inquire(file = filename, exist = fileExists)
    if (fileExists .and. iRank .eq. patch%masterRank)                                        &
         call MPI_FILE_DELETE(filename, mpiInfo, ierror)
    call MPI_FILE_OPEN(patch%comm, trim(mpiiofs)//trim(adjustl(filename)),                      &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, iunit, ierror)

    ! Write the geometry header
    if (iRank .eq. patch%masterRank) then
       buffer = 'C Binary'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'Ensight Gold Geometry File'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'Curvilinear Geometry from jCODE'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'node id off'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'element id off'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'part'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       i = 1
       call MPI_FILE_WRITE(iunit, i, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
       buffer = 'jCODE Grid'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'block'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_WRITE(iunit, patch%globalSize(1), 1, MPI_INTEGER, MPI_STATUS_IGNORE,    &
            ierror)
       call MPI_FILE_WRITE(iunit, patch%globalSize(2), 1, MPI_INTEGER, MPI_STATUS_IGNORE,    &
            ierror)
       call MPI_FILE_WRITE(iunit, patch%globalSize(3), 1, MPI_INTEGER, MPI_STATUS_IGNORE,    &
            ierror)
    end if

    ! Store coordinates in single-precision array
    data%buffer_SP = 0.0_SP
    call patch_collect(patch, coordinates, data%buffer_WP(:,1:nDimensions))
    do i = 1, nDimensions
       data%buffer_SP(:, i) = real(data%buffer_WP(:, i), SP)
    end do

    ! Write the grid coordinates
    offset = 80 * 8 + 4 * 4
    do i = 1, 3
       call MPI_File_set_view(iunit, offset, MPI_REAL, data%subarrayType, "native",          &
            mpiInfo, ierror)
       call MPI_File_write_all(iunit, data%buffer_SP(:, i), patch%nPatchPoints, MPI_REAL,    &
            MPI_STATUS_IGNORE, ierror)
       offset = offset + 4 * int(product(patch%globalSize), MPI_OFFSET_KIND)
    end do

    ! Close the file
    call MPI_File_close(iunit, ierror)

    return
  end subroutine dump_ensight_geometry


  ! Dump binary EnSight gold data - scalar
  ! --------------------------------------
  subroutine dump_ensight_scalar(name, scalar)

    ! External modules
    use parallel
    use fileio
    use grid_patch

    implicit none

    ! Arguments
    character(len = str_medium), intent(in) :: name
    real(WP), dimension(:), intent(in) :: scalar

    ! Local variables
    integer :: i, j, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = str_long) :: filename
    character(len = 80) :: buffer
    logical :: fileExists
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data

    ! Start the timer
    call timing_start('i/o')

    ! Extract data from each patch and write
    do j = 1, nVizPatches

       ! Pull off the data
       patch => vizPatch(j)
       data => ensightData(j)

       if (patch%nPatchPoints .le. 0) cycle

       call patch_collect(patch, scalar, data%buffer_WP(:,1))
       data%buffer_SP(:,1) = real(data%buffer_WP(:,1), SP)

       ! Generate the file
       if (iRank .eq. patch%masterRank) call CREATE_FOLDER(trim(adjustl(data%directory)) //  &
            "/" // trim(adjustl(name)))
       filename = trim(adjustl(data%directory)) // "/" // trim(adjustl(name)) // "/" //      &
            trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes

       ! Open the file
       inquire(file = trim(filename), exist = fileExists)
       if (fileExists .and. iRank .eq. patch%masterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(patch%comm, trim(mpiiofs)//trim(adjustl(filename)),                &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, iunit, ierror)

       ! Write header (master rank only)
       if (iRank .eq. patch%masterRank) then
          buffer = trim(adjustl(name))
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          buffer = 'part'
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          i = 1
          call MPI_FILE_WRITE(iunit, i, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
          buffer = 'block'
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       end if

       ! Write the file
       offset = 3 * 80 + 4
       call MPI_FILE_SET_VIEW(iunit, offset, MPI_REAL, data%subarrayType, "native",          &
            mpiInfo, ierror)
       call MPI_FILE_WRITE_ALL(iunit, data%buffer_SP(:,1), patch%nPatchPoints, MPI_REAL,     &
            MPI_STATUS_IGNORE, ierror)

       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

    end do
    
    ! Stop the timer
    call timing_stop('i/o')

    return
  end subroutine dump_ensight_scalar


  ! Dump binary EnSight gold data - vector
  ! --------------------------------------
  subroutine dump_ensight_vector(name, vector)

    ! External modules
    use parallel
    use fileio
    use grid_patch

    implicit none

    ! Arguments
    character(len = str_medium), intent(in) :: name
    real(WP), dimension(:,:), intent(in) :: vector

    ! Local variables
    integer :: i, j, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = str_long) :: filename
    character(len = 80) :: buffer
    logical :: fileExists
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data

    ! Start the timer
    call timing_start('i/o')

    ! Extract data from each patch and write
    do j = 1, nVizPatches

       ! Pull off the data
       patch => vizPatch(j)
       data => ensightData(j)

       if (patch%nPatchPoints .le. 0) cycle

       call patch_collect(patch, vector(:,1:nDimensions), data%buffer_WP(:,1:nDimensions))

       data%buffer_SP = 0.0_SP
       do i = 1, nDimensions
          data%buffer_SP(:, i) = real(data%buffer_WP(:, i), SP)
       end do

       ! Generate the file
       if (iRank .eq. patch%masterRank) call CREATE_FOLDER(trim(adjustl(data%directory)) //  &
            "/" // trim(adjustl(name)))
       filename = trim(adjustl(data%directory)) // "/" // trim(adjustl(name)) // "/" //      &
            trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes

       ! Open the file
       inquire(file = trim(filename), exist = fileExists)
       if (fileExists .and. iRank .eq. patch%masterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(patch%comm, trim(mpiiofs)//trim(adjustl(filename)),                &
            IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, iunit, ierror)

       ! Write header (root process only)
       if (iRank .eq. patch%masterRank) then
          buffer = trim(adjustl(name))
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          buffer = 'part'
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          i = 1
          call MPI_FILE_WRITE(iunit, i, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
          buffer = 'block'
          call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       end if

       ! Write the file
       offset = 3 * 80 + 4
       do i = 1, 3
          call MPI_File_set_view(iunit, offset, MPI_REAL, data%subarrayType, "native",       &
               mpiInfo, ierror)
          call MPI_File_write_all(iunit, data%buffer_SP(:, i), patch%nPatchPoints,           &
               MPI_REAL, MPI_STATUS_IGNORE, ierror)
          offset = offset + 4 * int(product(patch%globalSize), MPI_OFFSET_KIND)
       end do

       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

    end do

    ! Stop the timer
    call timing_stop('i/o')

    return
  end subroutine dump_ensight_vector


  ! Dump binary EnSight gold data - particles
  ! -----------------------------------------
  subroutine dump_ensight_particle

    ! External modules
    use parallel
    use fileio
    use simulation_flags
    use particle
    use math

    implicit none

    ! Local variables
    integer :: i, j, k, iunit, ierror, nVizParticles, nVizParticlesGlobal, vizMasterRank,    &
         vizRank, nVizProcs, vizComm, nVizPatches_
    integer, allocatable :: nParticlesProc(:)
    integer(kind=MPI_OFFSET_KIND) :: headerOffset, offset, particleOffset
    integer, allocatable :: partIndex(:), partInteger(:)
    real(SP), allocatable :: partScalar(:), partVector(:)
    real(WP) :: dxdt(3), dudt(3), dwdt(3), dTdt, dmdt
    character(len = str_long) :: particleDirectory, name, filename
    character(len = str_medium) :: vizDirectory
    character(len = 80) :: cbuffer
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data
    logical :: fileIsThere

    ! Start the timer
    call timing_start('i/o')

    if (nParticles.gt.0) allocate(partIndex(nParticles))

    ! Only write to one patch if all particles are dumped
    nVizPatches_ = nVizPatches
    if (dumpAllParticles) nVizPatches_ = 1

    ! Extract data from each patch and write
    do k = 1, nVizPatches

       if (.not.dumpParticleToPatch(k)) cycle

       ! Pull of the data
       patch => vizPatch(k)
       data => ensightData(k)

       if (.not.dumpAllParticles .and. patch%nPatchPoints.le.0) cycle

       ! Set number of visualization particles to zero for now
       nVizParticles = 0

       if (dumpAllParticles .or. product(patch%globalSize).eq.product(globalGridSize)) then
          if (k.eq.1) then
             ! Visualize all particles within the domain
             nVizParticles = nParticles
             nVizParticlesGlobal = nParticlesGlobal
             do i = 1, nParticles
                partIndex(i) = i
             end do
             vizComm = comm
             nVizProcs = nProcs
             vizRank = irank
             vizMasterRank = huge(1)
             if (patch%nPatchPoints .gt. 0) then
                vizDirectory = trim(adjustl(data%directory))
                vizMasterRank = patch%masterRank
             end if
             call parallel_min(vizMasterRank)
             call MPI_BCAST(vizDirectory, len(vizDirectory), MPI_CHARACTER, vizMasterRank, comm,&
                  ierror)
          end if
       else
          ! Visualize a subset of particles
          do i = 1, nParticles
             if ( particles(i)%gridIndex(1) .ge. patch%iMin .and.                            &
                  particles(i)%gridIndex(1) .le. patch%iMax .and.                            &
                  particles(i)%gridIndex(2) .ge. patch%jMin .and.                            &
                  particles(i)%gridIndex(2) .le. patch%jMax .and.                            &
                  particles(i)%gridIndex(3) .ge. patch%kMin .and.                            &
                  particles(i)%gridIndex(3) .le. patch%kMax ) then
                nVizParticles = nVizParticles + 1
                partIndex(nVizParticles) = i
             end if
          end do
          vizComm = patch%comm
          vizDirectory = trim(adjustl(data%directory))
          call MPI_ALLREDUCE(nVizParticles, nVizParticlesGlobal, 1, MPI_INTEGER, MPI_SUM,    &
               vizComm, ierror)
          call MPI_Comm_Size(vizComm, nVizProcs, ierror)
          call MPI_Comm_rank(vizComm, vizRank, ierror)
          call MPI_ALLREDUCE(vizRank, vizMasterRank, 1, MPI_INTEGER, MPI_MIN, vizComm, ierror)
       end if

       ! Allocate buffer arrays for writing
       if (nVizParticles.gt.0) then
          allocate(partInteger(nVizParticles))
          allocate(partScalar(nVizParticles))
          allocate(partVector(nVizParticles*3))
       end if

       ! Get number of viz particles per processor in `vizComm`
       allocate(nParticlesProc(nVizProcs)); nParticlesProc = 0
       call MPI_Allgather(nVizParticles, 1, MPI_INTEGER, nParticlesProc, 1, MPI_INTEGER,     &
            vizComm, ierror)

       ! Generate the file - Particles
       particleDirectory = trim(adjustl(vizDirectory)) // '/PARTICLES'
       name = 'PARTICLES'
       if (vizRank .eq. vizMasterRank) call CREATE_FOLDER(trim(adjustl(particleDirectory)))
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       filename = trim(mpiiofs) // trim(filename)
       ! Open the file
       inquire(file=filename, exist = fileIsThere)
       if (fileIsThere .and. vizRank.eq.vizMasterRank)                                       &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
            iunit, ierror)
       ! Write header
       if (vizRank .eq. vizMasterRank) then
          cbuffer = 'C Binary'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          cbuffer = 'Particle coordinates from jCODE'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          cbuffer = 'particle coordinates'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          call MPI_FILE_WRITE(iunit, max(nVizParticlesGlobal, 1), 1, MPI_INTEGER,            &
               MPI_STATUS_IGNORE, ierror)
          ! Zero particle case
          if (nVizParticlesGlobal .eq. 0) then
             call MPI_FILE_WRITE(iunit, 1, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
          end if
       end if
       if (nVizParticles .gt. 0) then
          headerOffset = 80 * 3 + 4
          particleOffset = 0
          do i = 1, vizRank
             particleOffset = particleOffset + int(nParticlesProc(i), MPI_OFFSET_KIND)
          end do
          offset = headerOffset + particleOffset * 4
          do i = 1, nVizParticles
             partInteger(i) = i + int(particleOffset)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partInteger, nVizParticles, MPI_INTEGER,     &
               MPI_STATUS_IGNORE, ierror)
          ! Write the particle positions
          offset = headerOffset + nVizParticlesGlobal*4 + particleOffset*3*4
          do i = 1, nVizParticles
             j = partIndex(i)
             partVector((i-1)*3+1) = real(particles(j)%position(1), SP)
             partVector((i-1)*3+2) = real(particles(j)%position(2), SP)
             partVector((i-1)*3+3) = real(particles(j)%position(3), SP)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partVector, nVizParticles*3, MPI_REAL_SP,    &
               MPI_STATUS_IGNORE, ierror)
       end if
       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle ID
       name = 'ID'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       filename = trim(mpiiofs) // trim(filename)
       ! Open the file
       inquire(file=filename, exist = fileIsThere)
       if (fileIsThere .and. vizRank .eq. vizMasterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
            iunit, ierror)
       ! Write header
       if (vizRank .eq. vizMasterRank) then
          cbuffer = 'Particle ID'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          ! Zero particle case
          if (nVizParticlesGlobal .eq. 0) then
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
          end if
       end if
       ! Compute the offset and write IDs
       if (nVizParticles .gt. 0) then
          offset = 80 + particleOffset*4
          do i = 1, nVizParticles
             j = partIndex(i)
             partScalar(i) = real(particles(j)%id, SP)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partScalar, nVizParticles, MPI_REAL_SP,      &
               MPI_STATUS_IGNORE, ierror)
       end if
       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle diameter
       name = 'DIAMETER'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       filename = trim(mpiiofs) // trim(filename)
       ! Open the file
       inquire(file=filename, exist = fileIsThere)
       if (fileIsThere .and. vizRank .eq. vizMasterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
            iunit, ierror)
       ! Write header
       if (vizRank .eq. vizMasterRank) then
          cbuffer = 'Particle diameter'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          ! Zero particle case
          if (nVizParticlesGlobal .eq. 0) then
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
          end if
       end if
       ! Compute the offset and write diameters
       if (nVizParticles .gt. 0) then
          offset = 80 + particleOffset*4
          do i = 1, nVizParticles
             j = partIndex(i)
             partScalar(i) = real(particles(j)%diameter, SP)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partScalar, nVizParticles, MPI_REAL_SP,      &
               MPI_STATUS_IGNORE, ierror)
       end if
       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle temperature
       name = 'Tp'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       filename = trim(mpiiofs) // trim(filename)
       ! Open the file
       inquire(file=filename, exist = fileIsThere)
       if (fileIsThere .and. vizRank .eq. vizMasterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
            iunit, ierror)
       ! Write header
       if (vizRank .eq. vizMasterRank) then
          cbuffer = 'Particle temperature'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          ! Zero particle case
          if (nVizParticlesGlobal .eq. 0) then
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
          end if
       end if
       ! Compute the offset and write temperatures
       if (nVizParticles .gt. 0) then
          offset = 80 + particleOffset*4
          do i = 1, nVizParticles
             j = partIndex(i)
             partScalar(i) = real(particles(j)%temperature, SP)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partScalar, nVizParticles, MPI_REAL_SP,      &
               MPI_STATUS_IGNORE, ierror)
       end if
       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle velocity
       name = 'Vp'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       filename = trim(mpiiofs) // trim(filename)
       ! Open the file
       inquire(file=filename, exist = fileIsThere)
       if (fileIsThere .and. vizRank .eq. vizMasterRank)                                     &
            call MPI_FILE_DELETE(filename, mpiInfo, ierror)
       call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiInfo, &
            iunit, ierror)
       ! Write header
       if (vizRank .eq. vizMasterRank) then
          cbuffer = 'Particle velocity'
          call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
          ! Zero particle case
          if (nVizParticlesGlobal .eq. 0) then
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
          end if
       end if
       ! Compute the offset and write velocities
       if (nVizParticles .gt. 0) then
          offset = 80 + particleOffset*4*3
          do i = 1, nVizParticles
             j = partIndex(i)
             partVector((i-1)*3+1) = real(particles(j)%velocity(1), SP)
             partVector((i-1)*3+2) = real(particles(j)%velocity(2), SP)
             partVector((i-1)*3+3) = real(particles(j)%velocity(3), SP)
          end do
          call MPI_FILE_WRITE_AT(iunit, offset, partVector, nVizParticles*3, MPI_REAL_SP,    &
               MPI_STATUS_IGNORE, ierror)
       end if
       ! Close the file
       call MPI_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle angular velocity
       if (useFriction) then
          name = 'OMEGAp'
          filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
          write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')         &
               nOutputTimes
          filename = trim(mpiiofs) // trim(filename)
          ! Open the file
          inquire(file=filename, exist = fileIsThere)
          if (fileIsThere .and. vizRank .eq. vizMasterRank)                                  &
               call MPI_FILE_DELETE(filename, mpiInfo, ierror)
          call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE),       &
               mpiInfo, iunit, ierror)
          ! Write header
          if (vizRank .eq. vizMasterRank) then
             cbuffer = 'Particle angular velocity'
             call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
             ! Zero particle case
             if (nVizParticlesGlobal .eq. 0) then
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             end if
          end if
          ! Compute the offset and write IDs
          if (nVizParticles .gt. 0) then
             offset = 80 + particleOffset*4*3
             do i = 1, nVizParticles
                j = partIndex(i)
                partVector((i-1)*3+1) = real(particles(j)%angularVelocity(1), SP)
                partVector((i-1)*3+2) = real(particles(j)%angularVelocity(2), SP)
                partVector((i-1)*3+3) = real(particles(j)%angularVelocity(3), SP)
             end do
             call MPI_FILE_WRITE_AT(iunit, offset, partVector, nVizParticles*3, MPI_REAL_SP, &
                  MPI_STATUS_IGNORE, ierror)
          end if
          ! Close the file
          call MPI_FILE_CLOSE(iunit, ierror)
       end if

       ! Generate the file - Particle acceleration
       if (dumpParticleAcceleration) then
          name = 'Ap'
          filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
          write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')         &
               nOutputTimes
          filename = trim(mpiiofs) // trim(filename)
          ! Open the file
          inquire(file=filename, exist = fileIsThere)
          if (fileIsThere .and. vizRank .eq. vizMasterRank)                                  &
               call MPI_FILE_DELETE(filename, mpiInfo, ierror)
          call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE),       &
               mpiInfo, iunit, ierror)
          ! Write header
          if (vizRank .eq. vizMasterRank) then
             cbuffer = 'Particle acceleration'
             call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
             ! Zero particle case
             if (nVizParticlesGlobal .eq. 0) then
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             end if
          end if
          ! Compute the offset and write acceleration
          if (nVizParticles .gt. 0) then
             offset = 80 + particleOffset*4*3
             do i = 1, nVizParticles
                j = partIndex(i)
                call particle_rhs(particles(j), dxdt, dudt, dwdt, dTdt, dmdt)
                partVector((i-1)*3+1) = real(dudt(1), SP)
                partVector((i-1)*3+2) = real(dudt(2), SP)
                partVector((i-1)*3+3) = real(dudt(3), SP)
             end do
             call MPI_FILE_WRITE_AT(iunit, offset, partVector, nVizParticles*3, MPI_REAL_SP, &
                  MPI_STATUS_IGNORE, ierror)
          end if
          ! Close the file
          call MPI_FILE_CLOSE(iunit, ierror)
       end if

       ! Generate the file - Particle collision force
       if (dumpParticleCollision) then
          name = 'COLp'
          filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
          write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')         &
               nOutputTimes
          filename = trim(mpiiofs) // trim(filename)
          ! Open the file
          inquire(file=filename, exist = fileIsThere)
          if (fileIsThere .and. vizRank .eq. vizMasterRank)                                  &
               call MPI_FILE_DELETE(filename, mpiInfo, ierror)
          call MPI_FILE_OPEN(vizComm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE),       &
               mpiInfo, iunit, ierror)
          ! Write header
          if (vizRank .eq. vizMasterRank) then
             cbuffer = 'Particle collision'
             call MPI_FILE_WRITE(iunit, cbuffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
             ! Zero particle case
             if (nVizParticlesGlobal .eq. 0) then
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
                call MPI_FILE_WRITE(iunit, 0.0_SP, 1, MPI_REAL_SP, MPI_STATUS_IGNORE, ierror)
             end if
          end if
          ! Compute the offset and write collisions
          if (nVizParticles .gt. 0) then
             offset = 80 + particleOffset*4*3
             do i = 1, nVizParticles
                j = partIndex(i)
                partVector((i-1)*3+1) = real(particles(j)%collision(1), SP)
                partVector((i-1)*3+2) = real(particles(j)%collision(2), SP)
                partVector((i-1)*3+3) = real(particles(j)%collision(3), SP)
             end do
             call MPI_FILE_WRITE_AT(iunit, offset, partVector, nVizParticles*3, MPI_REAL_SP, &
                  MPI_STATUS_IGNORE, ierror)
          end if
          ! Close the file
          call MPI_FILE_CLOSE(iunit, ierror)
       end if

       ! Cleanup patch data
       if (allocated(partInteger)) deallocate(partInteger)
       if (allocated(partScalar)) deallocate(partScalar)
       if (allocated(partVector)) deallocate(partVector)
       if (allocated(nParticlesProc)) deallocate(nParticlesProc)

    end do

    ! Cleanup
    if (allocated(partIndex)) deallocate(partIndex)

    ! Stop the timer
    call timing_stop('i/o')

    return
  end subroutine dump_ensight_particle


  ! Dump binary EnSight gold data - IBM particles
  ! This will write all IBM particles to one patch
  ! ----------------------------------------------
  subroutine dump_ensight_ibm_particle

    ! External modules
    use parallel
    use fileio
    use simulation_flags
    use ibm
    use math

    implicit none

    ! Local variables
    integer :: i, k, ibuffer, iunit, ierror, vizMasterRank
    real(SP) :: rbuffer
    character(len = str_long) :: particleDirectory, name, filename
    character(len = str_medium) :: vizDirectory
    character(len = 80) :: cbuffer
    type(t_Patch), pointer :: patch
    type(t_EnsightData), pointer :: data

    ! Start the timer
    call timing_start('i/o')

    ! Write to file associated with patch 1
    k = 1
    patch => vizPatch(k)
    data => ensightData(k)
    vizMasterRank = patch%masterRank
    vizDirectory = trim(adjustl(data%directory))

    ! Only master rank writes
    if (irank .eq. vizMasterRank) then

       ! Generate the file - Particles
       particleDirectory = trim(adjustl(vizDirectory)) // '/PARTICLES'
       name = 'PARTICLES'
       call CREATE_FOLDER(trim(adjustl(particleDirectory)))
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       ! Open the file
       call BINARY_FILE_OPEN(iunit, trim(filename), "w", ierror)
       ! Write header
       cbuffer = 'C Binary'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       cbuffer = 'Particle coordinates from jCODE'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       cbuffer = 'particle coordinates'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       ibuffer = max(nObjects, 1)
       call BINARY_FILE_WRITE(iunit, ibuffer, 1, kind(ibuffer), ierror)
       do i = 1, nObjects
          ibuffer = i
          call BINARY_FILE_WRITE(iunit, ibuffer, 1, kind(ibuffer), ierror)
       end do
       ! Particle coordinates
       do i = 1, nObjects
          rbuffer = real(object(i)%position(1), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%position(2), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%position(3), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end do
       ! Zero particle case
       if (nObjects .eq. 0) then
          ibuffer = 1
          call BINARY_FILE_WRITE(iunit, ibuffer, 1, kind(ibuffer), ierror)
          rbuffer = 0.0_SP
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end if
       ! Close the file
       call BINARY_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle diameter
       name = 'DIAMETER'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       ! Open the file
       call BINARY_FILE_OPEN(iunit, trim(filename), "w", ierror)
       ! Write header
       cbuffer = 'Particle diameter'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       ! Particle diameter
       do i = 1, nObjects
          rbuffer = real((2.0_WP * real(nDimensions, WP) * object(i)%volume / pi)            &
               ** (1.0_WP / real(nDimensions, WP)), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end do
       ! Single particle case (not sure why this is needed)
       if (nObjects .eq. 1) call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       ! Zero particle case
       if (nObjects .eq. 0) then
          rbuffer = 0.0_SP
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end if
       ! Close the file
       call BINARY_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle velocity
       name = 'Vp'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       call BINARY_FILE_OPEN(iunit, trim(filename), "w", ierror)
       ! Write header
       cbuffer = 'Particle velocity'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       ! Particle velocity
       do i = 1, nObjects
          rbuffer = real(object(i)%velocity(1), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%velocity(2), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%velocity(3), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end do
       ! Zero particle case
       if (nObjects .eq. 0) then
          rbuffer = 0.0_SP
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer) ,ierror)
       end if
       ! Close the file
       call BINARY_FILE_CLOSE(iunit, ierror)

       ! Generate the file - Particle angular velocity
       name = 'OMEGAp'
       filename = trim(adjustl(particleDirectory)) // "/" // trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            nOutputTimes
       call BINARY_FILE_OPEN(iunit, trim(filename), "w", ierror)
       ! Write header
       cbuffer = 'Particle velocity'
       call BINARY_FILE_WRITE(iunit, cbuffer, 80, kind(cbuffer), ierror)
       ! Particle velocity
       do i = 1, nObjects
          rbuffer = real(object(i)%angularVelocity(1), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%angularVelocity(2), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          rbuffer = real(object(i)%angularVelocity(3), SP)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
       end do
       ! Zero particle case
       if (nObjects .eq. 0) then
          rbuffer = 0.0_SP
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer), ierror)
          call BINARY_FILE_WRITE(iunit, rbuffer, 1, kind(rbuffer) ,ierror)
       end if
       ! Close the file
       call BINARY_FILE_CLOSE(iunit, ierror)

    end if

    ! Stop the timer
    call timing_stop('i/o')

    return
  end subroutine dump_ensight_ibm_particle

end module dump_ensight


! ============================================== !
! Dump binary EnSight gold data - Initialization !
! ============================================== !
subroutine dump_ensight_setup(mode)

  ! Internal modules
  use dump_ensight

  ! External modules
  use parser
  use fileio
  use parallel
  use simulation_flags
  use solver_options
  use geometry
  use time_info
  use controller
  use dissipation
  use ignition_source
  use ibm

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Local variables
  integer :: i, ierror
  logical :: fileIsThere
  character(len = 80) :: filename

  ! Setup the individual visualization patches
  allocate(ensightData(nVizPatches))
  do i = 1, nVizPatches
     call setup_ensight_patch(vizPatch(i), ensightData(i))
  end do

  select case (mode)

  case (FORWARD)

     ! Decide what to dump
     call parser_read('dump primitive variables', dumpPrimitiveVariables, .true.)
     call parser_read('dump processor decomp', dumpCPU, .true.)
     call parser_read('dump viscosity', dumpViscosity, .true.)
     call parser_read('dump Q criterion', dumpQcrit, .true.)
     call parser_read('dump vorticity', dumpVorticity, .false.)
     call parser_read('dump dilatation', dumpDilatation, .false.)
     call parser_read('dump baroclinic torque', dumpBtorque, .false.)
     call parser_read('dump Mach number', dumpMachNumber, .true.)
     call parser_read('dump pressure strain', dumpPressureStrain, .false.)
     call parser_read('dump schlieren', dumpSchlieren, .false.)
     call parser_read('dump shadowgraph', dumpShadowgraph, .false.)
     call parser_read('dump heat release', dumpHeatRelease, .true.)
     call parser_read('dump ignition source', dumpIgnition, .true.)
     call parser_read('dump particles', dumpParticles, .true.)
     call parser_read('dump jacobian', dumpJacobian, .false.)
     call parser_read('dump grid index', dumpGridIndex, .false.)
     call parser_read('dump grid spacing', dumpGridSpacing, .false.)
     call parser_read('dump grid norm', dumpGridNorm, .false.)
     call parser_read('dump arc length', dumpArcLength, .false.)
     call parser_read('dump target mollifier', dumpTargetMollifier, .true.)
     call parser_read('dump ibm', dumpIBM, useIBM)
     call parser_read('dump ibm particles', dumpIBMparticles, .true.)
     call parser_read('dump particle acceleration', dumpParticleAcceleration, .false.)
     call parser_read('dump particle collision', dumpParticleCollision, .false.)

     ! Turn off if not used
     if (.not.useParticles) dumpParticles = .false.
     if (.not.useViscosity) then
        dumpQcrit = .false.
        if (.not.useShockCapturing .and. .not.useLES) dumpViscosity = .false.
     end if
     if (.not.useIgnition) dumpIgnition = .false.
     if (.not. useIBM) then
        dumpIBM = .false.
        dumpIBMparticles = .false.
     else
        if (ibmType .ne. IBM_PARTICLE) dumpIBMparticles = .false.
     end if

     ! Avoid dumping Lagrangian particles and IBM particles
     if (dumpParticles .and. dumpIBMparticles)                                               &
          call die('dump_ensight_setup: cannot dump Lagrangian particles and IBM particles')

     ! Determine if all particles should be dumped and which patch should write
     allocate(dumpParticleToPatch(nVizPatches))
     dumpParticleToPatch = .false.
     if (dumpParticles) then
        call parser_read('dump all particles', dumpAllParticles, .false.)
        if (dumpAllParticles) then
           dumpParticleToPatch(1) = .true.
        else
           dumpParticleToPatch = .true.
        end if
     end if

     ! Determine additional information for numerical schlieren
     if (dumpSchlieren) then
        call parser_read('schlieren line of sight', schlierenDirection, 0)
        if (schlierenDirection.lt.0 .or. schlierenDirection.gt.3)                            &
             call die('dump_ensight_setup: schlieren direction must be between 0 & 3!')
        call parser_read('schlieren component', schlierenComponent, 0)
        if (schlierenComponent.lt.0 .or. schlierenComponent.gt.3)                            &
             call die('dump_ensight_setup: schlieren component must be between 0 & 3!')
     end if

     ! Create the directory
     call CREATE_FOLDER("ensight-3D")

     ! Create sub-directories if multiple patches are used
     if (nVizPatches .gt. 1) then
        do i = 1, nVizPatches
           if (iRank .eq. vizPatch(i)%masterRank) then
              call CREATE_FOLDER(trim(adjustl(ensightData(i)%directory)))
           end if
        end do
     end if

     ! Write the geometry
     do i = 1, nVizPatches
        call dump_ensight_geometry(vizPatch(i), ensightData(i))
     end do
     
     ! Read time info from the first patch
     if (iRank .eq. vizPatch(1)%masterRank) then
        filename = trim(adjustl(ensightData(1)%directory)) // "/" //                         &
             trim(adjustl(vizPatch(1)%name)) // ".case"
     end if
     call MPI_BCAST(filename, len(filename), MPI_CHARACTER, vizPatch(1)%masterRank, comm,    &
          ierror)
     inquire(file=trim(adjustl(filename)), exist = fileIsThere)
     if (fileIsThere) then
        ! Read the file
        call parser_parsefile(trim(adjustl(filename)))
        ! Get the time
        call parser_getsize('time values', nOutputTimes)
        allocate(outputTimes(nOutputTimes))
        call parser_read('time values', outputTimes)
        ! Remove future time
        future: do i = 1, size(outputTimes)
           if (outputTimes(i) .ge. time*0.99999_WP) then
              nOutputTimes = i - 1
              exit future
           end if
        end do future
     else
        ! Set the time
        nOutputTimes = 0
        allocate(outPutTimes(1))
     end if
     call MPI_BARRIER(comm, ierror)

  case (ADJOINT)

     ! Decide what to dump
     call parser_read('dump adjoint sensitivity', dumpSensitivity, .true.)
     if (.not. spaceTimeGradient) dumpSensitivity = .false.
     call parser_read('dump control mollifier', dumpControlMollifier, .true.)

     ! Assume ensight directory was created during the forward run

     ! Reset output times
     nOutputTimes = 0
     if (allocated(outputTimes)) deallocate(outputTimes)
     allocate(outputTimes(1))

  end select

  return
end subroutine dump_ensight_setup


! ====================== !
! Cleanup EnSight module !
! ====================== !
subroutine dump_ensight_cleanup

  ! Internal modules
  use dump_ensight

  implicit none

  ! Local variables
  integer :: i

  do i = 1, nVizPatches
     call cleanup_ensight_patch(ensightData(i))
  end do

  if (allocated(outputTimes)) deallocate(outputTimes)
  if (nVizPatches .gt. 0) nullify(vizPatch)
  if (associated(ensightData)) nullify(ensightData)

  return
end subroutine dump_ensight_cleanup


! ============================= !
! Dump binary EnSight gold data !
! ============================= !
subroutine dump_ensight_data(mode)

  ! Internal modules
  use dump_ensight

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use dissipation
  use grid_levelset
  use grid_patch
  use grid_functions
  use operator
  use state
  use state_functions
  use particle_exchange
  use ibm
  use time_info
  use combustion
  use ignition_source
  use functional
  use controller

  implicit none

  ! Arguments
  integer, intent(in) :: mode

  ! Local variables
  character(len = str_medium) :: name
  integer :: i, j, k, gridIndex
  real(WP), dimension(nGridPoints, nUnknowns) :: temp
  real(WP), allocatable, dimension(:,:) :: temp2

  select case(mode)

  case (FORWARD)

     ! Update the case file
     call dump_ensight_case

     ! CPU
     if (dumpCPU) then
        name = 'CPU'
        temp(:,1) = real(iRank + 1, SP)
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Primitive variables
     if (dumpPrimitiveVariables) then
        ! Density
        name = 'DENSITY'
        if (twoWayCoupling) then
           call dump_ensight_scalar(name, conservedVariables(:,1) / volumeFraction(:,1))
        else
           call dump_ensight_scalar(name, conservedVariables(:,1))
        end if

        ! Velocity
        name = 'VELOCITY'
        call dump_ensight_vector(name, velocity)

        ! Temperature
        name = 'TEMPERATURE'
        call dump_ensight_scalar(name, temperature(:,1))

        ! Pressure
        name = 'PRESSURE'
        call dump_ensight_scalar(name, pressure(:,1))

        ! Mass fractions
        do i = 1, nSpecies
           write(name, "(A,I2.2)") "MASS_FRACTION_", i
           call dump_ensight_scalar(name, massFraction(:,i))
        end do

     end if

     ! Viscosity
     if (dumpViscosity) then

        ! Store total viscosity
        temp(:,1) = dynamicViscosity(:,1)

        ! Turbulent viscosity
        if (useLES) then
           name = 'VISC_T'
           call dump_ensight_scalar(name, turbulentViscosity(:,1))
           
           ! Remove turbulent viscosity contribution
           temp(:,1) = temp(:,1) - turbulentViscosity(:,1)
        end if

        ! Artificial viscosities
        if (useShockCapturing) then
           name = 'VISC_SHOCK'
           call dump_ensight_scalar(name, artificialShearViscosity(:,1))
           name = 'BULK_VISC_SHOCK'
           call dump_ensight_scalar(name, artificialBulkViscosity(:,1))
           name = 'DIFF_SHOCK'
           call dump_ensight_scalar(name, artificialThermalDiffusivity(:,1))

           ! Remove artificial viscosity contribution
           temp(:,1) = temp(:,1) - artificialShearViscosity(:,1)
        end if

        ! Molecular viscosity
        name = 'VISC'
        call dump_ensight_scalar(name, temp(:,1))
        
     end if

     ! Hybrid dissipation
     if (hybridDissipation) then
        name = 'SENSOR'
        call dump_ensight_scalar(name, dissipationSensor(:,1))
     end if

     ! Q-criterion
     if (dumpQcrit) then
        name = 'QCRIT'
        call compute_Q_criterion(velocityGradient, temp(:,1))
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Vorticity
     if (dumpVorticity) then
        name = 'VORT'
        call compute_vorticity(velocityGradient, temp(:,1:3))
        if (nDimensions.eq.3) then
           call dump_ensight_vector(name, temp(:,1:3))
        else
           call dump_ensight_scalar(name, temp(:,3))
        end if
     end if

     ! Dilatation
     if (dumpDilatation) then
        name = 'DIVU'
        call compute_dilatation(velocityGradient, temp(:,1))
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Baroclinic torque: grad(rho) x grad(p) / rho^2
     if (dumpBtorque) then
        allocate(temp2(nGridPoints, 2*nDimensions))
        name = 'B_TORQUE'
        call gradient(conservedVariables(:,1), temp(:,1:nDimensions))
        call gradient(pressure(:,1), temp2(:,1:nDimensions))
        if (nDimensions.eq.3) then
           temp2(:,4) = (temp(:,2) * temp2(:,3) - temp(:,3) * temp2(:,2)) /                  &
                conservedVariables(:,1)**2
           temp2(:,5) = (temp(:,3) * temp2(:,1) - temp(:,1) * temp2(:,3)) /                  &
                conservedVariables(:,1)**2
           temp2(:,6) = (temp(:,1) * temp2(:,2) - temp(:,2) * temp2(:,1)) /                  &
                conservedVariables(:,1)**2
           call dump_ensight_vector(name, temp2(:,4:6))
        else
           temp(:,3) = (temp(:,1) * temp2(:,2) - temp(:,2) * temp2(:,1)) /                   &
                conservedVariables(:,1)**2
           call dump_ensight_scalar(name, temp(:,3))
        end if
        deallocate(temp2)
     end if

     ! Mach number
     if (dumpMachNumber) then
        name = 'MACH'
        temp(:,1) = sqrt(sum(velocity(:,1:nDimensions)**2, 2)) /                             &
        sqrt(ratioOfSpecificHeats * specificVolume(:,1) * pressure(:,1))
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Pressure strain <p' div(u')>
     if (dumpPressureStrain) then
        name = 'PSTRAIN'        
        do i = 1, nGridPoints
           select case (nDimensions)
           case (1)
              temp(i,1) = velocityGradient(i,1)
           case (2)
              temp(i,1) = velocityGradient(i,1) + velocityGradient(i,4)
           case (3)
              temp(i,1) = velocityGradient(i,1) + velocityGradient(i,5) +                    &
                   velocityGradient(i,9)
           end select
           temp(i,1) = temp(i,1) * pressure(i,1)
        end do
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Numerical schlieren
     if (dumpSchlieren) then
        name = 'SCHLIEREN'        
        temp(:,1) = conservedVariables(:,1)
        if (twoWayCoupling) temp(:,1) = temp(:,1) / volumeFraction(:,1)
        call compute_schlieren(temp(:,1), schlierenDirection, schlierenComponent)
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Shadowgraph
     if (dumpShadowgraph) then
        name = 'SHADOWGRAPH'        
        temp(:,1) = conservedVariables(:,1)
        if (twoWayCoupling) temp(:,1) = temp(:,1) / volumeFraction(:,1)
        call laplacian(temp(:,1), temp(:,2))
        call dump_ensight_scalar(name, temp(:,2))
     end if

     ! Combustion
     if (dumpHeatRelease .and. nReactions .gt. 0) then
        name = 'HEAT_RELEASE'
        temp = 0.0_WP
        call combustion_source(FORWARD, temp)
        call dump_ensight_scalar(name, temp(:,nDimensions+2))
     end if

     ! Ignition
     if (dumpIgnition) then
        name = 'IGNITION_SOURCE'
        temp = 0.0_WP
        call ignition_source_forward(time, temp)
        call dump_ensight_scalar(name, temp(:,nDimensions+2))
     end if

     ! Immersed boundary
     if (dumpIBM) then
        ! Levelset
        name = 'LEVELSET'
        call dump_ensight_scalar(name, levelset(:,1))

        ! Object index
        name = 'OBJ_ID'
        call dump_ensight_scalar(name, real(objectIndex,WP))

        if (ibm_move) then
           name = 'IBM_VEL'
           call dump_ensight_vector(name, ibmVelocity)
        end if

        ! IBM particles
        if (dumpIBMparticles) call dump_ensight_ibm_particle
     end if

     if (useParticles) then
        ! Volume fraction
        name = 'VOLUME_FRACTION'
        call dump_ensight_scalar(name, volumeFraction(:,1))

        ! Granular temperature
        if (useGranularTemperature) then
           name = 'GRANULAR_TEMP'
           call dump_ensight_scalar(name, granularTemperature(:,1))
        end if
     end if

     ! Jacobian
     if (dumpJacobian) then
        name = 'JAC'
        call dump_ensight_scalar(name, jacobian(:,1))
     end if

     ! Grid spacing
     if (dumpGridIndex) then
        do k = iStart(3), iEnd(3)
           do j = iStart(2), iEnd(2)
              do i = iStart(1), iEnd(1)
                 gridIndex = i - gridOffset(1) + localGridSize(1) *                          &
                      (j - 1 - gridOffset(2) + localGridSize(2) *                            &
                      (k - 1 - gridOffset(3)))
                 temp(gridIndex, 1) = real(i, WP)
                 temp(gridIndex, 2) = real(j, WP)
                 temp(gridIndex, 3) = real(k, WP)
              end do
           end do
        end do
        name = 'IJK'
        call dump_ensight_vector(name, temp(:,1:nDimensions))
     end if

     ! Grid spacing
     if (dumpGridSpacing) then
        name = 'DX'
        call dump_ensight_vector(name, gridSpacing)
     end if

     ! Grid norm
     if (dumpGridNorm) then
        name = 'NORM'
        call dump_ensight_scalar(name, gridNorm(:,1))
     end if

     ! Arc length
     if (dumpArcLength .and. allocated(arcLengths)) then
        name = 'ARC'
        call dump_ensight_vector(name, arcLengths)
     end if

     ! Target mollifier
     if (dumpTargetMollifier .and. .not. predictionOnly) then
        name = 'TARGET_MOLLIFIER'
        call dump_ensight_scalar(name, targetMollifier(:,1))
     end if

     ! Particle data
     if (dumpParticles) call dump_ensight_particle

  case (ADJOINT)

     ! Update the case file
     call dump_ensight_case_adjoint

     ! Write the adjoint variables
     do i = 1, nUnknowns
        write(name, "(A,I2.2)") "ADJOINT_", i
        call dump_ensight_scalar(name, adjointVariables(:,i))
     end do

     ! Sensitivity
     if (dumpSensitivity) then
        name = 'SENSITIVITY'
        if (iGradientBuffer .eq. 0) then
           i = size(gradientBuffer, 3)
        else
           i = iGradientBuffer
        end if
        call patch_disperse(controllerPatch, gradientBuffer(:,1,i), temp(:,1))
        call dump_ensight_scalar(name, temp(:,1))
     end if

     ! Control mollifier
     if (dumpControlMollifier) then
        name = 'CONTROL_MOLLIFIER'
        call dump_ensight_scalar(name, controlMollifier(:,1))
     end if

  end select

  return
end subroutine dump_ensight_data
