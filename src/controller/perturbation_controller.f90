module perturbation_controller

  ! External modules
  use math, only : pi
  use controller
  use precision
  use solver_options

  implicit none

  ! Perturbation parameters
  integer  :: nModesx, nModesz, iModex, iModez, iMode
  real(WP) :: P0, T0(1), minAmplitude, maxAmplitude,initialThicknessInverse,                 &
       perturbationEnergy, maxControlMollifierInverse
  real(WP), allocatable, dimension(:) :: waveNumberx, waveNumberz,                           &
       amplitude1, amplitude2, amplitude3, amplitude4, phase1, phase2
  logical  :: isothermalInitialConditions, sineTerms, constantEnergy, velocityCorrection,    &
       waveNumberSensitivity, amplitudeSensitivity, phaseSensitivity

contains

  subroutine write_initial_perturbation(controlIteration)

    ! External modules
    use parallel
    use fileio

    implicit none

    ! Arguments
    integer, intent(in)  :: controlIteration

    ! Local variables
    integer :: i, j, m, fileUnit, iostat
    character(len = str_medium) :: filename

    ! Set the name
    filename = 'perturbation_init.txt'
    if (controlIteration .gt. 0)                                                             &
         write(filename, '(A,I3.3,A)') 'perturbation_init_', controlIteration, '.txt'

    ! Create the directory
    call CREATE_FOLDER("perturbation_init")
    filename = 'perturbation_init/' // trim(adjustl(filename))

    ! Root process writes
    if (iRank .eq. iRoot) then

       ! Open the file
       fileUnit = iOpen()
       open(unit = fileUnit, file = trim(filename), action = 'write', status = 'unknown',    &
            iostat = iostat)

       ! Write the file
       if (nDimensions .eq. 3) then

          if (sineTerms) then
             write(fileUnit, '(A4,2X,A4,6(2X,A24))') 'i', 'j', 'Wave number x',              &
                  'Wave number z', 'Amplitude 1', 'Amplitude 2', 'Amplitude 3', 'Amplitude 4'
             do i = 1, nModesx
                do j = 1, nModesz
                   m = (j-1)*nModesx + i
                   write(fileUnit, '(I4,2X,I4,6(2X,ES24.15E3))') i, j, waveNumberx(i),       &
                        waveNumberz(j), amplitude1(m), amplitude2(m), amplitude3(m),         &
                        amplitude4(m)
                end do
             end do
          else
             write(fileUnit, '(A4,2X,A4,5(2X,A24))') 'i', 'j', 'Wave number x',              &
                  'Wave number z', 'Amplitude', 'Phase x', 'Phase z'
             do i = 1, nModesx
                do j = 1, nModesz
                   m = (j-1)*nModesx + i
                   write(fileUnit, '(I4,2X,I4,5(2X,ES24.15E3))') i, j, waveNumberx(i),       &
                        waveNumberz(j), amplitude1(m), phase1(i), phase2(j)
                end do
             end do
          end if

       else

          if (sineTerms) then
             write(fileUnit, '(A4,3(2X,A24))') 'i', 'Wave number', 'Amplitude 1',            &
                  'Amplitude 2'
             do i = 1, nModesx
                write(fileUnit, '(I4,3(2X,ES24.15E3))') i, waveNumberx(i), amplitude1(i),    &
                     amplitude2(i)
             end do
          else
             write(fileUnit, '(A4,3(2X,A24))') 'i', 'Wave number', 'Amplitude', 'Phase'          
             do i = 1, nModesx
                write(fileUnit, '(I4,3(2X,ES24.15E3))') i, waveNumberx(i), amplitude1(i),    &
                     phase1(i)
             end do
          end if

       end if

       ! Close the file
       close(iclose(fileUnit))
    end if

    return
  end subroutine write_initial_perturbation

  subroutine read_initial_perturbation

    ! External modules
    use parallel
    use fileio
    use parser

    implicit none

    ! Local variables
    integer :: i, j, m, fileUnit, iostat, iBuffer, jBuffer
    character(len = str_medium) :: filename

    ! Read in the file name
    call parser_read('initial perturbation file', filename, 'perturbation_init_read.txt')    

    ! Root process reads
    if (iRank .eq. iRoot) then
       ! Open the file
       fileUnit = iOpen()
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'unknown',     &
            iostat = iostat)
       if (iostat .ne. 0) call die('read_initial_perturbation: Failed to open initial &
            &perturbation file!')

       ! Read the file header
       read(fileUnit, *, iostat = iostat)

       if (nDimensions .eq. 3) then

          if (sineTerms) then
             do i = 1, nModesx
                do j = 1, nModesz
                   m = (j-1)*nModesx + i
                   read(fileUnit, *, iostat = iostat) iBuffer, jBuffer, waveNumberx(i),      &
                        waveNumberz(j), amplitude1(m), amplitude2(m), amplitude3(m),         &
                        amplitude4(m)
                end do
             end do
          else
             do i = 1, nModesx
                do j = 1, nModesz
                   m = (j-1)*nModesx + i
                   read(fileUnit, *, iostat = iostat) iBuffer, jBuffer, waveNumberx(i),      &
                        waveNumberz(j), amplitude1(m), phase1(i), phase2(j)
                end do
             end do
          end if

       else

          if (sineTerms) then
             do i = 1, nModesx
                read(fileUnit, *, iostat = iostat) iBuffer, waveNumberx(i), amplitude1(i),   &
                     amplitude2(i)
             end do
          else          
             do i = 1, nModesx
                read(fileUnit, *, iostat = iostat) iBuffer, waveNumberx(i), amplitude1(i),   &
                     phase1(i)
             end do
          end if

       end if

       ! Close the file
       close(iclose(fileUnit))
    end if

    ! Broadcast read data with all processors
    call parallel_bc(waveNumberx)
    if (nDimensions .eq. 3) then
       call parallel_bc(waveNumberz)
    else
       waveNumberz = 0.0_WP
    end if
    call parallel_bc(amplitude1)
    if (sineTerms) then
       call parallel_bc(amplitude2)
       if (nDimensions .eq. 3) then
          call parallel_bc(amplitude3)
          call parallel_bc(amplitude4)
       else
          amplitude3 = 0.0_WP
          amplitude4 = 0.0_WP
       end if
    else
       call parallel_bc(phase1)
       if (nDimensions .eq. 3) then
          call parallel_bc(phase2)
       else
          phase2 = 0.0_WP
       end if
    end if

    return
  end subroutine read_initial_perturbation

  ! Compute the perturbation
  ! --------------------------
  subroutine calculate_eta_2d(eta_2d, cosX, sinX, cosZ, sinZ)

    ! External modules
    use parallel
    use grid_functions

    implicit none

    ! Arguments
    real(WP), dimension(:,:), intent(out) :: eta_2d, cosX, sinX, cosZ, sinZ

    ! Local variables
    integer :: i, j, k, m, n, p, nx, nz, ierror
    real(WP) :: x, z

    ! number of grid points in the horizontal directions
    nx = globalGridSize(1)
    nz = globalGridSize(3)

    eta_2d = 0.0_WP
    cosX = 1.0_WP; sinX = 0.0_WP
    cosZ = 1.0_WP; sinZ = 0.0_WP
    if (iRankDir(2).eq.iRoot) then
       j = iStart(2)
       do k = iStart(3), iEnd(3)
          do i = iStart(1), iEnd(1)

             x = coordinates(grid_index(i,j,k), 1)
             z = 0.0_WP
             if (nz .gt. 1) z = coordinates(grid_index(i,j,k), 3)

             if (sineTerms) then
                do n = 1, nModesz
                   do m = 1, nModesx
                      p = (n-1)*nModesx + m
                      if (k.eq.iStart(3) .and. n.eq.1) then
                         cosX(i,m) = cos( waveNumberx(m) * x )
                         sinX(i,m) = sin( waveNumberx(m) * x )
                      end if
                      if (nDimensions.eq.3 .and. i.eq.iStart(1) .and. m.eq.1) then
                         cosZ(k,n) = cos( waveNumberz(n) * z )
                         sinZ(k,n) = sin( waveNumberz(n) * z )
                      end if
                      eta_2d(i,k) = eta_2d(i,k)                                              &
                           + amplitude1(p) * cosX(i,m) * cosZ(k,n)                           &
                           + amplitude2(p) * sinX(i,m) * cosZ(k,n)                           &
                           + amplitude3(p) * cosX(i,m) * sinZ(k,n)                           &
                           + amplitude4(p) * sinX(i,m) * sinZ(k,n)
                   end do
                end do
             else
                do n = 1, nModesz
                   do m = 1, nModesx
                      p = (n-1)*nModesx + m
                      if (k.eq.iStart(3) .and. n.eq.1) then
                         cosX(i,m) = cos( waveNumberx(m) * x + phase1(m) )
                         sinX(i,m) = sin( waveNumberx(m) * x + phase1(m) )
                      end if
                      if (nDimensions.eq.3 .and. i.eq.iStart(1) .and. m.eq.1) then
                         cosZ(k,n) = cos( waveNumberz(n) * z + phase2(n) )
                         sinZ(k,n) = sin( waveNumberz(n) * z + phase2(n) )
                      end if
                      eta_2d(i,k) = eta_2d(i,k) + amplitude1(p) * cosX(i,m) * cosZ(k,n)
                   end do
                end do
             end if

          end do
       end do
    end if ! (iRankDir(2).eq.iRoot)

    ! Broadcast eta to other cores in the vertical direction
    call MPI_BCAST(eta_2d, nx*nz,      MPI_REAL_WP, iRoot, commDir(2), ierror)
    call MPI_BCAST(cosX,   nx*nModesx, MPI_REAL_WP, iRoot, commDir(2), ierror)
    call MPI_BCAST(sinX,   nx*nModesx, MPI_REAL_WP, iRoot, commDir(2), ierror)
    call MPI_BCAST(cosZ,   nz*nModesz, MPI_REAL_WP, iRoot, commDir(2), ierror)
    call MPI_BCAST(sinZ,   nz*nModesz, MPI_REAL_WP, iRoot, commDir(2), ierror)

  end subroutine calculate_eta_2d

end module perturbation_controller


! =============================== !
! Setup the perturbation actuator !
! =============================== !
subroutine perturbation_controller_setup

  ! Internal modules
  use perturbation_controller

  ! External modules
  use string
  use parallel
  use parser
  use equation_of_state
  use random

  implicit none

  ! Local variables
  integer  :: i, j, m, n
  real(WP) :: meanAmplitude, rmsAmplitude, rand, temp
  logical  :: readPerturbation, randomAmplitude, randomPhase
  character(len = str_medium) :: key

  ! Initialize the random number generator
  call random_setup

  ! Check in if reads an initial perturbation parameters from a file
  call parser_read('read initial perturbation', readPerturbation, .false.)

  ! Read in perturbation properties 
  call parser_read('number of modes of x', nModesx, 1)
  nModesz = 1
  if (nDimensions .eq. 3)  call parser_read('number of modes of z', nModesz, 1)

  ! Check in being constant energy
  call parser_read('constant perturbation energy', constantEnergy, .false.)

  ! Check if add constraint has been imposed for the optimizer library
  call parser_read('number of constraints', nConstraint, 0)
  !if (nConstraint .gt. 1) call die('perturbation_controller_setup: Having more than 1 &
  !     &constraint is not available for now!')
  ! Check equality of constraints
  allocate(equalityConstraint(nConstraint))
  do i = 1, nConstraint
     write(key, '(A,I1.1)') 'equality for constraint ', i 
     call parser_read(key, equalityConstraint(i), .true.)
  end do

  ! No constantEnergy is necessary if nConstraint>0 for optimization
  if (nConstraint.gt.0) constantEnergy = .false.

  ! Check in sensitivity types
  call parser_read('compute wave number sensitivity', waveNumberSensitivity, .false.)
  call parser_read('compute amplitude sensitivity', amplitudeSensitivity, .false.)
  call parser_read('compute phase sensitivity', phaseSensitivity, .false.)

  if (constantEnergy .and. .not.amplitudeSensitivity)                                        &
       call die('perturbation_controller_setup: constant perturbation energy only works when &
       &amplitude sensitivity is requested!')

  randomAmplitude = .false.
  randomPhase     = .false.
  if (nModesx .gt. 1 .or. nModesz .gt. 1) then
     call parser_read('random perturbation amplitude', randomAmplitude, .false.)
     call parser_read('random perturbation phase', randomPhase, .false.)
  end if

  sineTerms = .false.
  if (randomAmplitude .and. randomPhase) call parser_read('perturbation sine term',         &
       sineTerms, .false.)

  if (sineTerms .and. constantEnergy) call die('perturbation_controller_setup: Not yet &
       &implemented for including both constant perturbation energy and sine terms')

  ! Check to make sure this is consistent
  n = 0
  if (waveNumberSensitivity) n = n + nModesx + (nDimensions - 2) * nModesz
  if (amplitudeSensitivity) n = n + nModesx * nModesz
  if (phaseSensitivity) then
     if (sineTerms) then
        n = n + (2 * (nDimensions - 2) + 1) * nModesx * nModesz
     else
        n = n + nModesx + (nDimensions - 2) * nModesz
     end if
  end if
  if (constantEnergy) n = n - 1 ! ... one of the amplitudes depends on the others!
  if (n .eq. 0) call die('perturbation_controller_setup: must compute sensitivity of &
       &something!')
  if (nControlParameters .ne. n) call die("perturbation_controller_setup: number &
       &of control parametes must be " // char(n) // "!")

  ! Allocate the parameters
  allocate(waveNumberx(nModesx), waveNumberz(nModesz))
  allocate(amplitude1(nModesx*nModesz))
  if (sineTerms) then
     allocate(amplitude2(nModesx*nModesz), amplitude3(nModesx*nModesz),                      &
          amplitude4(nModesx*nModesz))
  else
     allocate(phase1(nModesx), phase2(nModesz))
  end if

  ! Set wave numbers of x direction
  do i = 1, nModesx
     write(key, '(A,I4.4)') 'perturbation wave number of x ', i
     call parser_read(trim(key), waveNumberx(i), real(i, WP))
  end do

  ! Set wave numbers of z direction
  if (nDimensions.eq.3) then
     do i = 1, nModesz
        write(key, '(A,I4.4)') 'perturbation wave number of z ', i
        call parser_read(trim(key), waveNumberz(i), real(i, WP))
     end do
  else
     waveNumberz = 0.0_WP
  end if

  ! Read in amplitude constraints
  call parser_read('maximum perturbation amplitude', maxAmplitude, huge(1.0_WP))
  if (maxAmplitude .lt. 0.0_WP) call die("perturbation_controller_setup: maximum &
       &perturbation amplitude must NOT be negative!")

  call parser_read('minimum perturbation amplitude', minAmplitude, -huge(1.0_WP))
  if (minAmplitude .gt. maxAmplitude) call die("perturbation_controller_setup: minimum &
       &perturbation amplitude must NOT be greater than maximum perturbation amplitude!")

  call parser_read('mean perturbation amplitude', meanAmplitude,                          &
       (maxAmplitude - minAmplitude) / 2.0_WP)
  if (meanAmplitude .lt. 0.0_WP) call die("perturbation_controller_setup: mean &
       &perturbation amplitude must NOT be negative!") 

  call parser_read('rms perturbation amplitude', rmsAmplitude, 0.0_WP)
  if (rmsAmplitude .lt. 0.0_WP) call die("perturbation_controller_setup: rms &
       &perturbation amplitude must NOT be negative!")

  if (.not. readPerturbation) then
     ! Set amplitude(s) and phase(s)
     if (iRank .eq. iRoot) then
        do j = 1, nModesz
           do i = 1, nModesx
              m = (j-1)*nModesx + i

              if (randomAmplitude) then
                 amplitude1(m) = max(minAmplitude, random_normal(meanAmplitude,rmsAmplitude))
              else
                 amplitude1(m) = maxAmplitude
              end if

              if (sineTerms) then
                 amplitude2(m) = max(minAmplitude, random_normal(meanAmplitude,rmsAmplitude))
                 if (nDimensions .eq. 3) then
                    amplitude3(m)=max(minAmplitude,random_normal(meanAmplitude,rmsAmplitude))
                    amplitude4(m)=max(minAmplitude,random_normal(meanAmplitude,rmsAmplitude))
                 else
                    amplitude3(m) = 0.0_WP
                    amplitude4(m) = 0.0_WP
                 end if
              elseif (randomPhase) then
                 call random_number(rand)
                 phase1(i) = rand * 2.0_WP * pi
                 if (nDimensions .eq. 3) then
                    call random_number(rand)
                    phase2(j) = rand * 2.0_WP * pi
                 else
                    phase2(j) = 0.0_WP
                 end if
              else
                 phase1(i) = 0.0_WP
                 phase2(j) = 0.0_WP
              end if

           end do
        end do
     end if

     ! Broadcast parameters with all processors
     call parallel_bc(amplitude1)
     if (sineTerms) then
        call parallel_bc(amplitude2)
        call parallel_bc(amplitude3)
        call parallel_bc(amplitude4)
     else
        call parallel_bc(phase1)
        call parallel_bc(phase2)
     end if

  else
     call read_initial_perturbation
  end if

  ! Store maximum magnitude of control mollifier
  maxControlMollifierInverse = maxVal(controlMollifier(:,1))
  call parallel_max(maxControlMollifierInverse)
  maxControlMollifierInverse = 1.0_WP / maxControlMollifierInverse

  ! Compute the initial perturbation energy
  temp = sum(amplitude1**2)
  if (sineTerms .and. phaseSensitivity) then ! if not phase sens, not be stored at baseline
     temp = temp + sum(amplitude2**2) + sum(amplitude3**2) + sum(amplitude4**2)
  end if

  ! Read in the initial perturbation energy
  call parser_read('initial perturbation energy', perturbationEnergy, temp)
  ! if (abs(temp - perturbationEnergy).gt.0.0_WP) call die('perturbation_controller_setup: &
  !      &the provided initial perturbation energy perturbation is wrong!')

  ! Read in the dependent mode wavenumbers for constant energy problem 
  iModex = 1
  iModez = 1
  if (constantEnergy) then
     call parser_read('Dependent mode index in x', iModex, 1)
     iModex = max(1,       iModex)
     iModex = min(nModesx, iModex)
     if (nDimensions.eq.3) call parser_read('Dependent mode index in z', iModez, 1)
     iModez = max(1,       iModez)
     iModez = min(nModesz, iModez)
  end if
  iMode = (iModez-1)*nModesx + iModex

  ! Check in if the amplitude of the dependend mode is not zero
  if (constantEnergy .and. abs(amplitude1(iMode)).le.0.0_WP)                                 &
       call die('perturbation_controller_setup: The amplitude of the dependent mode is zero, &
             &infinite sensitivities may be obtained')

  ! Initialize the baseline perturbation
  allocate(baselineValue(nControlParameters))
  baselineValue = 0.0_WP

  ! Set the baseline values and rename sensitivityParameter
  n = 0
  if (waveNumberSensitivity) then
     do i = 1, nModesx
        baselineValue(i) = waveNumberx(i)
        write(sensitivityParameter(i),'(A,I4.4)') "kx_", i
     end do
     n = n + nModesx

     if (nDimensions .eq. 3) then
        do i = 1, nModesz
           baselineValue(n+i) = waveNumberz(i)
           write(sensitivityParameter(n+i),'(A,I4.4)') "kz_", i
        end do
        n = n + nModesz
     end if
  end if

  if (amplitudeSensitivity) then
     do j = 1, nModesz
        do i = 1, nModesx
           m = (j-1)*nModesx + i
           if (constantEnergy .and. m.eq.iMode) then
              ! Skip amplitude1(iMode)
              n = n - 1
              cycle 
           end if

           baselineValue(n + m) = amplitude1(m)
           if (nDimensions .eq. 2) then
              write(sensitivityParameter(n + m), '(A,(I4.4))') "a_", i
           else
              write(sensitivityParameter(n + m), '(2(A,(I4.4)))') "a_", i, '_', j
           end if
        end do
     end do
     n = n + nModesx*nModesz
  end if

  if (phaseSensitivity) then

     if (sineTerms) then
        do j = 1, nModesz
           do i = 1, nModesx

              m = (j-1)*nModesx + i
              baselineValue(n + m) = amplitude2(m)
              if (nDimensions .eq. 2) then
                 write(sensitivityParameter(n + m), '(A,(I4.4))') "b_", i
              else
                 write(sensitivityParameter(n + m), '(2(A,(I4.4)))') "b_", i, '_', j
                 baselineValue(n + nModesx*nModesz + m) = amplitude3(m)
                 write(sensitivityParameter(n + nModesx*nModesz + m), '(2(A,(I4.4)))')       &
                      "c_", i, '_', j
                 baselineValue(n + 2*nModesx*nModesz + m) = amplitude4(m)
                 write(sensitivityParameter(n + 2*nModesx*nModesz + m), '(2(A,(I4.4)))')     &
                      "d_", i, '_', j
              end if

           end do
        end do
     else
        do i = 1, nModesx
           baselineValue(n+i) = phase1(i)
           if (nDimensions .eq. 2) then
              write(sensitivityParameter(n+i),'(A,I4.4)') "theta_", i
           else
              write(sensitivityParameter(n+i),'(A,I4.4)') "theta_x_", i
           end if
        end do
        n = n + nModesx

        if (nDimensions .eq. 3) then
           do i = 1, nModesz
              baselineValue(n+i) = phase2(i)
              write(sensitivityParameter(n+i),'(A,I4.4)') "theta_z_", i
           end do
        end if
     end if

  end if

  ! Check in isothermal initial conditions Rayleigh-Taylor initialazing
  call parser_read('isothermal initial conditions', isothermalInitialConditions, .true.)

  ! Read in reference pressure for Rayleigh-Taylor initialazing
  call parser_read('reference pressure', P0, 1.0_WP / ratioOfSpecificHeats)

  if (isothermalInitialConditions) then
     ! Read in reference temperature
     call parser_read('reference temperature', T0(1), 1.0_WP/(ratioOfSpecificHeats - 1.0_WP))
  else
     ! Calculate the reference temperature, reference density equals to the molecular weight
     T0 = ratioOfSpecificHeats * P0 / (ratioOfSpecificHeats - 1.0_WP)
  end if

  ! Read in initial diffusion thickness for Rayleigh-Taylor initialazing
  call parser_read('initial diffusion thickness', initialThicknessInverse)
  initialThicknessInverse = 1.0_WP / initialThicknessInverse

  ! Correct the velocity to prevent unphysical waves at the interface
  call parser_read('rayleigh taylor velocity correction', velocityCorrection, .false.)
  if (velocityCorrection .and. .not. useViscosity) call die('perturbation_controller_setup: &
       &velocity correction cannot be used with inviscid flow!')
  if (velocityCorrection.and.waveNumberSensitivity) call die('perturbation_controller_setup: &
       &not yet implemented for including both velocity correction and wave number &
       &sensitivity!')
  if (velocityCorrection.and.phaseSensitivity) call die('perturbation_controller_setup: not &
       &yet implemented for for including both velocity correction and phase sensitivity!')

  return
end subroutine perturbation_controller_setup


! ================================= !
! Cleanup the perturbation actuator !
! ================================= !
subroutine perturbation_controller_cleanup

  ! Internal modules
  use perturbation_controller

  implicit none

  if (allocated(waveNumberx)) deallocate(waveNumberx)
  if (allocated(waveNumberz)) deallocate(waveNumberz)
  if (allocated(amplitude1)) deallocate (amplitude1)
  if (allocated(amplitude2)) deallocate(amplitude2)
  if (allocated(amplitude3)) deallocate(amplitude3, amplitude4)
  if (allocated(phase1)) deallocate(phase1)
  if (allocated(phase2)) deallocate(phase2)
  if (allocated(equalityConstraint)) deallocate(equalityConstraint)

  return
end subroutine perturbation_controller_cleanup


! ===================================== !
! Update control parameters using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine perturbation_controller_update_parameters(actuationAmount)

  ! Internal modules
  use perturbation_controller

  ! External modules

  implicit none

  ! Arguments
  real(WP), intent(in) :: actuationAmount

  ! Local variables
  integer  :: i, j, m, n, nModesz2
  real(WP) :: energySum
  real(WP), parameter :: smallNumber = 1.0E-9_WP

  if (abs(actuationAmount) .le. 0.0_WP) return

  nModesz2 = 0 ! ... may be called in for loops later
  if (nDimensions.eq.3) nModesz2 = nModesz

  n = 0
  ! Update the wavenumber
  if (waveNumberSensitivity) then
     do i = 1, nModesx
        waveNumberx(i) = baselineValue(i)                                                    &
             + real(gradientDirection, WP) * actuationAmount * controlGradient(i)
     end do
     n = n + nModesx

     do i = 1, nModesz2
        waveNumberz(i) = baselineValue(n+i)                                                  &
             + real(gradientDirection, WP) * actuationAmount * controlGradient(n+i)
     end do
     n = n + nModesz2
  end if

  ! Update amplitude1
  if (amplitudeSensitivity) then
     do j = 1, nModesz
        do i = 1, nModesx
           m = (j-1)*nModesx + i
           if (constantEnergy .and. m.eq.iMode) then
              ! Skip updating amplitude1(iMode) for now
              n = n - 1
              cycle 
           end if
           amplitude1(m) = max(minAmplitude, baselineValue(n+m)                              &
                + real(gradientDirection, WP) * actuationAmount * controlGradient(n+m))
           amplitude1(m) = min(maxAmplitude, amplitude1(m))
        end do
     end do
     n = n + nModesx*nModesz
  end if

  ! Update amplitude2-4 or phase1-2
  if (phaseSensitivity) then
     if (sineTerms) then
        do j = 1, nModesz
           do i = 1, nModesx
              m = (j-1)*nModesx + i
              
              amplitude3(m) = max(minAmplitude, baselineValue(n + nModesx*nModesz2 + m)      &
                   + real(gradientDirection, WP) * actuationAmount                           &
                   * controlGradient(n + nModesx*nModesz2 + m))
              amplitude3(m) = min(maxAmplitude, amplitude3(m))

              amplitude4(m) = max(minAmplitude, baselineValue(n + 2*nModesx*nModesz2 + m)    &
                   + real(gradientDirection, WP) * actuationAmount                           &
                   * controlGradient(n + 2*nModesx*nModesz2 + m))
              amplitude4(m) = min(maxAmplitude, amplitude4(m))

              ! Updating amplitude"2", to prevent being ovewrited in
              ! some cases, it comes after updating amplitude3 and amplitude4
              amplitude2(m) = max(minAmplitude, baselineValue(n+m)                           &
                   + real(gradientDirection, WP) * actuationAmount                           &
                   * controlGradient(n+m))
              amplitude2(m) = min(maxAmplitude, amplitude2(m))
           end do
        end do
     else
        do i = 1, nModesx
           phase1(i) = baselineValue(n+i)                                                    &
                + real(gradientDirection, WP) * actuationAmount * controlGradient(n+i)
        end do
        n = n + nModesx

        do i = 1, nModesz2
           phase2(i) = baselineValue(n+i)                                                    &
                + real(gradientDirection, WP) * actuationAmount * controlGradient(n+i)
        end do
     end if

  end if

  if (constantEnergy) then
     ! Update amplitde1(iMode) for having constant perturbation energy
     energySum = sum(amplitude1**2) - amplitude1(iMode)**2
     if (sineTerms .and. phaseSensitivity) then
        energySum = energySum + sum(amplitude2**2) + sum(amplitude3**2) + sum(amplitude4**2)
     end if

     if (perturbationEnergy .ge. energySum) then
        amplitude1(iMode) = sign(1.0_WP, amplitude1(iMode)) *                                &
             sqrt(perturbationEnergy - energySum)
     else
        amplitude1(iMode) = smallNumber
     end if
  end if

  return
end subroutine perturbation_controller_update_parameters


! ===================================== !
! Update initial conditions using the   !
! adjoint gradient and acutation amount !
! ===================================== !
subroutine perturbation_controller_update_ic(controlIteration)

  ! Internal modules
  use perturbation_controller

  ! External modules
  use geometry
  use grid
  use grid_functions
  use state, only : conservedVariables
  use equation_of_state

  implicit none

  ! Arguments
  integer, intent(in)  :: controlIteration

  ! Local variables
  integer  :: ii, jj, kk, gridIndex, nx, nz, m, n, p
  real(WP) :: y, yn, eta, molarFraction1, massFraction1, pressure, molecularWeightOfMixture, &
       integralOfW, constant1, constant2, constant3, constant4, constant5, temperature(1),   &
       massDiffusivity(1,1), dEta_dX, dEta_dZ
  real(WP), allocatable, dimension(:,:) :: eta_2d, cosX, sinX, cosZ, sinZ

  ! Compute constants
  constant1 = 0.5_WP*( 1.0_WP/molecularWeightInverse(1) + 1.0_WP/molecularWeightInverse(2) )
  constant2 = 0.5_WP*( 1.0_WP/molecularWeightInverse(1) - 1.0_WP/molecularWeightInverse(2) )
  constant3 = constant2 / (initialThicknessInverse * sqrt(pi))
  constant4 = ratioOfSpecificHeats/( (ratioOfSpecificHeats-1.0_WP)*T0(1) ) ! used if T0=const
  constant5 = 2.0_WP * constant2 * initialThicknessInverse / sqrt(pi) ! if velocityCorrection

  nx = globalGridSize(1)
  nz = globalGridSize(3)

  ! (Re)allocate and calclate eta, cosX, sinX, cosZ, sinZ
  allocate(eta_2d(nx,nz))
  allocate(cosX(nx,nModesx), sinX(nx,nModesx), cosZ(nz,nModesz), sinZ(nz,nModesz))
  call calculate_eta_2d(eta_2d, cosX, sinX, cosZ, sinZ)

  ! Compute the mass fraction and density and pressure
  do kk = iStart(3), iEnd(3)
     do jj = iStart(2), iEnd(2)
        do ii = iStart(1), iEnd(1)

           gridIndex = grid_index(ii,jj,kk)

           ! Vertical coordinate
           y = coordinates(gridIndex, 2)

           ! Compute the perturbated fluid interface
           eta = eta_2d(ii,kk) * controlMollifier(gridIndex,1) * maxControlMollifierInverse

           ! Compute the normalized perturbed vertical position
           yn = (y - eta) * initialThicknessInverse

           ! Compute the molecular weight of mixture
           molecularWeightOfMixture = constant1 + constant2 * erf(yn)

           ! Compute the molar and mass fractions
           molarFraction1 = 0.5_WP * (1.0_WP + erf(yn))
           massFraction1 = molarFraction1 / (molecularWeightInverse(1) *                     &
                molecularWeightOfMixture)

           ! Compute the integral of molecular weight of mixture
           integralOfW = molecularWeightOfMixture * (y - eta) +                              &
                constant3 * ( exp( - yn**2 ) - 1.0_WP )

           ! Compute the pressure and overwrite the density
           if (isothermalInitialConditions) then           
              ! Compute the pressure
              pressure = P0 * exp( - constant4 * froudeNumberInverse * integralOfW )

              ! Overwrite the density
              conservedVariables(gridIndex,1) = constant4*molecularWeightOfMixture*pressure
           else
              ! Overwrite the density
              conservedVariables(gridIndex,1) = molecularWeightOfMixture

              ! Compute the pressure 
              pressure = P0 - froudeNumberInverse * integralOfW
           end if

           ! Overwrite the momentum
           conservedVariables(gridIndex, 2:nDimensions+1) = 0.0_WP
           if (velocityCorrection) then
              ! Compute the local temperature if necessary
              temperature = T0
              if (.not.isothermalInitialConditions) then ! ... molecularWeight = density
                 temperature = ratioOfSpecificHeats/(ratioOfSpecificHeats - 1.0_WP)*pressure
              end if

              call compute_transport_variables(temperature, massDiffusivity=massDiffusivity)

              ! Compute the vertical momentum
              conservedVariables(gridIndex, 3) = - massDiffusivity(1,1) * constant5 *        &
                   exp( - yn**2 ) / molecularWeightOfMixture

              ! Compute the horizontal momentum(s)
              dEta_dX = 0.0_WP; dEta_dZ = 0.0_WP
              if (sineTerms) then
                 do n = 1, nModesz
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m
                       dEta_dX = dEta_dX                                                     &
                            - amplitude1(p) * waveNumberx(m) * sinX(ii,m) * cosZ(kk,n)       &
                            + amplitude2(p) * waveNumberx(m) * cosX(ii,m) * cosZ(kk,n)       &
                            - amplitude3(p) * waveNumberx(m) * sinX(ii,m) * sinZ(kk,n)       &
                            + amplitude4(p) * waveNumberx(m) * cosX(ii,m) * sinZ(kk,n)

                       dEta_dZ = dEta_dZ                                                     &
                            - amplitude1(p) * waveNumberz(n) * cosX(ii,m) * sinZ(kk,n)       &
                            - amplitude2(p) * waveNumberz(n) * sinX(ii,m) * sinZ(kk,n)       &
                            + amplitude3(p) * waveNumberz(n) * cosX(ii,m) * cosZ(kk,n)       &
                            + amplitude4(p) * waveNumberz(n) * sinX(ii,m) * cosZ(kk,n)
                    end do
                 end do
              else
                 do n = 1, nModesz
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m
                       dEta_dX = dEta_dX - amplitude1(p) * waveNumberx(m) * sinX(ii,m)       &
                            * cosZ(kk,n)

                       dEta_dZ = dEta_dZ - amplitude1(p) * waveNumberz(n) * cosX(ii,m)       &
                            * sinZ(kk,n)
                    end do
                 end do
              end if
              ! Consider the control mollifier effect
              dEta_dX = dEta_dX * controlMollifier(gridIndex,1) * maxControlMollifierInverse
              conservedVariables(gridIndex,2) = - conservedVariables(gridIndex,3) * dEta_dX
              if (nDimensions.eq.3) then
                 dEta_dZ = dEta_dZ * controlMollifier(gridIndex,1)*maxControlMollifierInverse
                 conservedVariables(gridIndex,4) = - conservedVariables(gridIndex,3) *dEta_dZ
              end if
           end if ! (velocityCorrection)

           ! Overwrite other conserved variables
           conservedVariables(gridIndex, nDimensions+2) = pressure /                         &
                (ratioOfSpecificHeats - 1.0_WP) + 0.5_WP *                                   &
                sum(conservedVariables(gridIndex,2:nDimensions+1)**2) /                      &
                conservedVariables(gridIndex,1)
           conservedVariables(gridIndex,nDimensions+3) = conservedVariables(gridIndex,1) *   &
                massFraction1

        end do
     end do
  end do

  ! Clean up
  deallocate(eta_2d, cosX, sinX, cosZ, sinZ)

  call write_initial_perturbation(controlIteration)

  return
end subroutine perturbation_controller_update_ic


! =============================== !
! Compute the adjoint sensitivity !
! =============================== !
subroutine perturbation_controller_compute_sensitivity

  ! Internal modules
  use perturbation_controller

  implicit none

  instantaneousSensitivity = 0.0_WP
  return

  return
end subroutine perturbation_controller_compute_sensitivity


! ======================================================== !
! Compute the adjoint sensitivity of the initial condition !
! ======================================================== !
subroutine perturbation_controller_compute_sensitivity_ic

  ! Internal modules
  use perturbation_controller

  ! External modules

  use parallel
  use geometry
  use solver_options
  use grid, only : coordinates
  use grid_functions, only : inner_product, grid_index
  use state
  use time_info

  implicit none

  ! Local variables
  integer  :: ii, jj, kk, gridIndex, m, n, p, s, nx, nz, nModesz2
  real(WP) ::  x, y, z, yn, molecularWeightOfMixture, eta, dEta, dPressure, dLnP, dW, dLnW,  &
       innerProduct, constant1, constant2, constant5, temp, tempX, tempZ, ddEta_dX, ddEta_dZ
  real(WP), dimension(nUnknowns) :: dQ0
  real(WP), allocatable, dimension(:,:) :: eta_2d, cosX, sinX, cosZ, sinZ

  ! Only measure sensitivity of initial condition
  if (timestep .ne. 0) return

  ! Compute constants
  constant1 = 0.5_WP*( 1.0_WP/molecularWeightInverse(1) + 1.0_WP/molecularWeightInverse(2) )
  constant2 = 0.5_WP*( 1.0_WP/molecularWeightInverse(1) - 1.0_WP/molecularWeightInverse(2) )
  constant5 = 2.0_WP * constant2 * initialThicknessInverse / sqrt(pi)

  nx = globalGridSize(1)
  nz = globalGridSize(3)

  ! (Re)allocate and calclate eta, cosX, sinX, cosZ, sinZ
  allocate(eta_2d(nx,nz))
  allocate(cosX(nx,nModesx), sinX(nx,nModesx), cosZ(nz,nModesz), sinZ(nz,nModesz))
  call calculate_eta_2d(eta_2d, cosX, sinX, cosZ, sinZ)

  nModesz2 = 0 ! .. be called in for loops later
  if (nDimensions.eq.3) nModesz2 = nModesz

  ! Compute Sensitivities
  instantaneousSensitivity = 0.0_WP
  do kk = iStart(3), iEnd(3)
     do jj = iStart(2), iEnd(2)
        do ii = iStart(1), iEnd(1)

           gridIndex = grid_index(ii,jj,kk)

           ! Coordinates
           x = coordinates(gridIndex, 1)
           y = coordinates(gridIndex, 2)
           z = 0.0_WP
           if (nz .gt. 1) z = coordinates(gridIndex, 3)

           eta = eta_2d(ii,kk) * controlMollifier(gridIndex,1) * maxControlMollifierInverse

           ! Compute the normalized perturbed vertical position
           yn = (y - eta) * initialThicknessInverse

           ! Compute the molecular weight of mixture and its derivative w.r.t eta 
           molecularWeightOfMixture = constant1 + constant2 * erf(yn)
           dW = - constant5 * exp( - yn**2 )
           dLnW = dW / molecularWeightOfMixture ! ... derivative of Ln(W)

           ! Compute the derivative of pressure w.r.t eta
           dPressure = conservedVariables(gridIndex, 1) * froudeNumberInverse
           dLnP = dPressure / pressure(gridIndex, 1) ! ... derivative of Ln(pressure)

           ! Compute the gradient of density w.r.t eta
           if (isothermalInitialConditions) then
              dQ0(1) = conservedVariables(gridIndex, 1) * (dLnW + dLnP)
           else
              dQ0(1) = dW ! ... molecularWeightOfMixture = density
           end if

           ! Compute the gradient of momentum w.r.t eta
           dQ0(2:nDimensions+1) = 0.0_WP
           tempX = 0.0_WP; tempZ = 0.0_WP
           if (velocityCorrection) then
              temp = 2.0_WP * yn * initialThicknessInverse - dLnW
              if (.not.isothermalInitialConditions .and. powerLawExponent.gt.0.0_WP)         &
                   temp = temp + powerLawExponent * dlnP

              dQ0(2:nDimensions+1) = conservedVariables(gridIndex, 2:nDimensions+1) * temp

              ! Will be used later
              tempX = - gridNorm(gridIndex, 1) * ( adjointVariables(gridIndex, 2) +          &
                   velocity(gridIndex, 1) * adjointVariables(gridIndex, nDimensions+2) ) *   &
                   conservedVariables(gridIndex, 3) * controlMollifier(gridIndex, 1) *       &
                   maxControlMollifierInverse
              if (nDimensions.eq.3)                                                          &
                   tempZ = - gridNorm(gridIndex, 1) * ( adjointVariables(gridIndex, 4) +     &
                   velocity(gridIndex, 3) * adjointVariables(gridIndex, nDimensions+2) ) *   &
                   conservedVariables(gridIndex, 3) * controlMollifier(gridIndex, 1) *       &
                   maxControlMollifierInverse
           end if

           ! Compute the gradient of energy w.r.t eta
           dQ0(nDimensions+2) = dPressure / ( ratioOfSpecificHeats - 1.0_WP )                &
                + sum( velocity(gridIndex, 1:nDimensions) * dQ0(2:nDimensions+1) )           &
                - 0.5_WP * sum( velocity(gridIndex, 1:nDimensions)**2 ) * dQ0(1)

           ! Compute the gradient of mass fraction w.r.t eta
           dQ0(nDimensions+3) = massFraction(gridIndex, 1) * dQ0(1) +                        &
                conservedVariables(gridIndex, 1) * dLnW *                                    &
                ( 0.5_WP / (molecularWeightInverse(1)*constant2) - massFraction(gridIndex,1) )

           ! Compute inner product of dQ0 with adjoint variables
           innerProduct = gridNorm(gridIndex,1) * sum(adjointVariables(gridIndex,:)*dQ0(:))  &
                 * controlMollifier(gridIndex,1) * maxControlMollifierInverse

           s = 0
           ! Update the waveNumber sensitivities
           ddEta_dX = 0.0_WP; ddEta_dZ = 0.0_WP
           if (waveNumberSensitivity) then
              if (sineTerms) then
                 do m = 1, nModesx
                    dEta = 0.0_WP
                    do n = 1, nModesz
                       p = (n-1)*nModesx + m
                       dEta = dEta                                                           &
                            - amplitude1(p) * sinX(ii,m) * cosZ(kk,n)                        &
                            + amplitude2(p) * cosX(ii,m) * cosZ(kk,n)                        &
                            - amplitude3(p) * sinX(ii,m) * sinZ(kk,n)                        &
                            + amplitude4(p) * cosX(ii,m) * sinZ(kk,n)
                    end do
                    dEta = dEta * x
                    instantaneousSensitivity(m) = instantaneousSensitivity(m) +              &
                         dEta * innerProduct
                 end do
                 s = s + nModesx

                 do n = 1, nModesz2
                    dEta = 0.0_WP
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m
                       dEta = dEta                                                           &
                            - amplitude1(p) * cosX(ii,m) * sinZ(kk,n)                        &
                            - amplitude2(p) * sinX(ii,m) * sinZ(kk,n)                        &
                            + amplitude3(p) * cosX(ii,m) * cosZ(kk,n)                        &
                            + amplitude4(p) * sinX(ii,m) * cosZ(kk,n)
                    end do
                    dEta = dEta * z
                    instantaneousSensitivity(s+n) = instantaneousSensitivity(s+n) +          &
                         dEta * innerProduct
                 end do
                 s = s + nModesz2

              else !if not simeTerms

                 do m = 1, nModesx
                    dEta = 0.0_WP
                    do n = 1, nModesz
                       p = (n-1)*nModesx + m
                       dEta = dEta - amplitude1(p) * sinX(ii,m) * cosZ(kk,n)
                    end do
                    dEta = dEta * x
                    instantaneousSensitivity(m) = instantaneousSensitivity(m) +              &
                         dEta * innerProduct
                 end do
                 s = s + nModesx

                 do n = 1, nModesz2
                    dEta = 0.0_WP
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m
                       dEta = dEta - amplitude1(p) * cosX(ii,m) * sinZ(kk,n)
                    end do
                    dEta = dEta * z
                    instantaneousSensitivity(s+n) = instantaneousSensitivity(s+n) +          &
                         dEta * innerProduct
                 end do
                 s = s + nModesz2
              end if
           end if

           ! Update the amplitude1 sensitivities
           ddEta_dX = 0.0_WP; ddEta_dZ = 0.0_WP
           if (amplitudeSensitivity) then
              do n = 1, nModesz
                 do m = 1, nModesx
                    p = (n-1)*nModesx + m
                    if (constantEnergy .and. p.eq.iMode) then
                       ! Skip computing sensitivity of amplitude1(iMode)
                       s = s - 1
                       cycle 
                    end if
                    
                    dEta = cosX(ii,m) * cosZ(kk,n)
                    if (velocityCorrection) then
                       ddEta_dX = - waveNumberx(m) * sinX(ii,m) * cosZ(kk,n)
                       ddEta_dZ = - waveNumberz(n) * cosX(ii,m) * sinZ(kk,n)
                    end if
                    if (constantEnergy) then
                       dEta = dEta - amplitude1(p) / amplitude1(iMode)                       &
                            * cosX(ii,iModex) * cosZ(kk,iModez)
                       if (velocityCorrection) then
                          ddEta_dX = ddEta_dX + amplitude1(p) / amplitude1(iMode) *          &
                               waveNumberx(iModex) * sinX(ii,iModex) * cosZ(kk,iModez)
                          ddEta_dZ = ddEta_dZ + amplitude1(p) / amplitude1(iMode) *          &
                               waveNumberz(iModez) * cosX(ii,iModex) * sinZ(kk,iModez)
                       end if
                    end if
                    instantaneousSensitivity(s + p) = instantaneousSensitivity(s + p)        &
                         + dEta * innerProduct + tempX * ddEta_dX + tempZ * ddEta_dZ
                 end do
              end do
              s = s + nModesx * nModesz
           end if

           ! Update the amplitude2-4 or phase1-2 sensitivities
           ddEta_dX = 0.0_WP; ddEta_dZ = 0.0_WP
           if (phaseSensitivity) then
              if (sineTerms) then

                 do n = 1, nModesz
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m

                       ! Updating amplitude"3" sensitivity
                       dEta = cosX(ii,m) * sinZ(kk,n)
                       instantaneousSensitivity(s + nModesx*nModesz2 + p) =                  &
                            instantaneousSensitivity(s + nModesx*nModesz2 + p) + dEta *      &
                            innerProduct

                       ! Updating amplitude"4" sensitivity
                       dEta = sinX(ii,m) * sinZ(kk,n)
                       instantaneousSensitivity(s + 2*nModesx*nModesz2 + p) =                &
                            instantaneousSensitivity(s + 2*nModesx*nModesz2 + p) + dEta *    &
                            innerProduct

                       ! Updating amplitude"2" sensitivity, to prevent ovewriting it in some
                       ! cases, it comes after updating amplitude3 and amplitude4
                       dEta = sinX(ii,m) * cosZ(kk,n)
                       instantaneousSensitivity(s + p) = instantaneousSensitivity(s + p)     &
                            + dEta * innerProduct
                    end do
                 end do
              else ! if not sineTerms

                 do m = 1, nModesx
                    dEta = 0.0_WP
                    do n = 1, nModesz
                       p = (n-1)*nModesx + m
                       dEta = dEta - amplitude1(p) * sinX(ii,m) * cosZ(kk,n)
                    end do
                    instantaneousSensitivity(s+m) = instantaneousSensitivity(s+m) +          &
                         dEta * innerProduct
                 end do
                 s = s + nModesx

                 do n = 1, nModesz2
                    dEta = 0.0_WP
                    do m = 1, nModesx
                       p = (n-1)*nModesx + m
                       dEta = dEta - amplitude1(p) * cosX(ii,m) * sinZ(kk,n)  
                    end do
                    instantaneousSensitivity(s+n) = instantaneousSensitivity(s+n) +          &
                         dEta * innerProduct
                 end do

              end if
           end if

        end do
     end do
  end do

  ! Clean up
  deallocate(eta_2d, cosX, sinX, cosZ, sinZ)

  call parallel_sum(instantaneousSensitivity)

  return
end subroutine perturbation_controller_compute_sensitivity_ic


! ========================= !
! Set the parameters bounds !
! ========================= !
subroutine perturbation_controller_parameter_bound(parameterLowerBound, parameterUpperBound, &
     parameterScale)

  ! Internal modules
  use perturbation_controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(out) :: parameterLowerBound,               &
       parameterUpperBound
  real(WP), dimension(nControlParameters), intent(out), optional :: parameterScale

  ! Local variables and parameters
  integer  :: i, f

  i = 0
  ! Set lower and upper bounds for the wavenumbers
  if (waveNumberSensitivity) then
     i = 1
     f = nModesx + (nDimensions-2)*nModesz
     parameterLowerBound(i:f) = -huge(1.0_WP)
     parameterUpperBound(i:f) = huge(1.0_WP)
     if (present(parameterScale)) parameterScale(i:f) = 1.0_WP
     i = f
  end if

  ! Set lower and upper bounds for amplitude1
  if (amplitudeSensitivity) then
     f = i + nModesx*nModesz
     if (constantEnergy) f = f - 1
     i = i + 1
     parameterLowerBound(i:f) = minAmplitude
     parameterUpperBound(i:f) = maxAmplitude
     if (present(parameterScale)) parameterScale(i:f) = maxAmplitude
     i = f
  end if

  ! Set lower and upper bounds for amplitude2-4 or phase1-2
  if (phaseSensitivity) then
     i = i + 1
     f = nControlParameters
     if (sineTerms) then ! amplitude2-4
        parameterLowerBound(i:f) = minAmplitude
        parameterUpperBound(i:f) = maxAmplitude
        if (present(parameterScale)) parameterScale(i:f) = maxAmplitude
     else ! phase1-2
        parameterLowerBound(i:f) = 0.0_WP
        parameterUpperBound(i:f) = 2.0_WP * pi
        if (present(parameterScale)) parameterScale(i:f) = 2.0_WP * pi
     end if
  end if

  return
end subroutine perturbation_controller_parameter_bound


! ============================ !
! Compute constraint functions !
! ============================ !
subroutine perturbation_controller_constraint(controlParameters, functionValue)

  ! Internal modules
  use perturbation_controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(in) :: controlParameters
  real(WP), dimension(nConstraint), intent(out) :: functionValue

  ! Local variables
  integer :: startIndex, endIndex

  ! Set start and end indices for amplitudes in controlParameters
  startIndex = 1
  endIndex = nControlParameters
  if (waveNumberSensitivity) startIndex = 1 + nModesx + (nDimensions-2)*nModesz
  if (.not.sineTerms.and.phaseSensitivity) endIndex = nControlParameters - nModesx -         &
       (nDimensions-2)*nModesz

  ! Calculate the constraint
  functionValue(1) = sum(controlParameters(startIndex:endIndex)**2) - perturbationEnergy /   &
       maxAmplitude**2

  return
end subroutine perturbation_controller_constraint


! ============================ !
! Compute constraint functions !
! ============================ !
subroutine perturbation_controller_constraint_gradient(controlParameters, gradientValue)

  ! Internal modules
  use perturbation_controller

  implicit none

  ! Arguments
  real(WP), dimension(nControlParameters), intent(in) :: controlParameters
  real(WP), dimension(nConstraint,nControlParameters), intent(out) :: gradientValue

  ! Local variables
  integer :: startIndex, endIndex

  ! Set start and end indices for amplitudes in controlParameters
  startIndex = 1
  endIndex = nControlParameters
  if (waveNumberSensitivity) startIndex = 1 + nModesx + (nDimensions-2)*nModesz
  if (.not.sineTerms.and.phaseSensitivity) endIndex = nControlParameters - nModesx -         &
       (nDimensions-2)*nModesz

  ! Initialize gradientValue
  gradientValue = 0.0_WP
  ! Calculate gradientValue
  gradientValue(1,startIndex:endIndex) = 2.0_WP * controlParameters(startIndex:endIndex)

  return
end subroutine perturbation_controller_constraint_gradient
