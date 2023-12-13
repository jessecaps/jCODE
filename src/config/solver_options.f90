module solver_options

  ! External modules
  use precision
  use string

  implicit none

  ! Solver mode
  integer, parameter ::                                                                      &
       FORWARD   = +1,                                                                       &
       ADJOINT   = -1,                                                                       &
       UNDEFINED =  0

  ! I/O options
  integer, parameter ::                                                                      &
       IO_GRID           = 100,                                                              &
       IO_FORWARD_STATE  = 101,                                                              &
       IO_TARGET_STATE   = 102,                                                              &
       IO_ADJOINT_STATE  = 103,                                                              &
       IO_PARTICLE       = 104,                                                              &
       IO_IBM            = 105,                                                              &
       IO_LEVELSET       = 106

  ! Equation of state
  integer, parameter ::                                                                      &
       IDEAL_GAS         = 1,                                                                &
       IDEAL_GAS_MIXTURE = 2

  ! Viscosity
  integer, parameter ::                                                                      &
       POWER_LAW         = 1,                                                                &
       SUTHERLANDS_LAW   = 2

  ! Sensor types
  integer, parameter ::                                                                      &
       DUCROS_SENSOR  = 1,                                                                   &
       BOUNDED_SENSOR = 2

  ! General solver options
  real(WP) :: reynoldsNumberInverse, prandtlNumberInverse, froudeNumberInverse,              &
       ratioOfSpecificHeats, powerLawExponent, sutherlandConstant, bulkViscosityRatio,       &
       shockCoefficient, minSpecies, maxSpecies, minTemperature, maxTemperature, inputCFL,   &
       inputTimeStepSize,  targetDensity  
  real(WP), dimension(:), allocatable ::  schmidtNumberInverse, molecularWeightInverse
  integer :: nSpecies, nUnknowns, equationOfState, viscosityModel
  character(len = str_medium) :: simulationName
  character(len = str_medium), dimension(:), allocatable :: speciesName

contains

  ! Get the number of unknowns from the solution file
  ! -------------------------------------------------
  subroutine get_nUnknowns

    ! External modules
    use string
    use parser
    use fileio
    use parallel
    use simulation_flags
    use geometry

    implicit none

    integer :: ifile, ierror
    integer, dimension(4) :: dims
    character(len=str_medium) :: filename_, filename
    logical :: fileIsThere

    ! Get the solution filename
    call parser_read('solution file to read', filename_)

    ! Get the name of the header file in case serial i/o is used
    fileIsThere = .false.
    if (useSerialIO) then
       filename = trim(filename_) // '/data.header'
       inquire(file = filename, exist = fileIsThere)
    end if
    if (.not. fileIsThere) filename = trim(filename_)

    ! Root reads the header file
    if (irank.eq.iroot) then

       ! Open the file
       call BINARY_FILE_OPEN(ifile, trim(filename), "r", ierror)

       ! Check number of processors are correct if serial file is used
       if (fileIsThere) then
          call BINARY_FILE_READ(ifile, dims(1), 1, kind(dims), ierror)
          if (dims(1) .ne. nProcs) then
             print*, 'Expected ', dims(1), ' processor(s)'
             print*, 'Currently ', nProcs, ' processor(s)'
             call die('Number of processors incompatible with serial data file')
          end if
       end if

       ! Read dimensions from header
       call BINARY_FILE_READ(ifile, dims, 4, kind(dims), ierror)
       if ((dims(1).ne.globalGridSize(1)) .or. (dims(2).ne.globalGridSize(2)) .or.           &
            (dims(3).ne.globalGridSize(3))) then
          print*, 'grid = ',globalGridSize(1), globalGridSize(2), globalGridSize(3)
          print*, trim(filename_), ' = ',dims(1),dims(2),dims(3)
          call die('The size of the data file does not correspond to the grid file')
       end if
       ! Set the number of unknowns
       nUnknowns = dims(4)

       ! Close the file
       call BINARY_FILE_CLOSE(ifile, ierror)
    end if

    ! Broadcast information
    call parallel_bc(nUnknowns)

    ! Do we transport any scalars?
    nSpecies = nUnknowns - nDimensions - 2
    if (nSpecies .lt. 0) then
       call die('Number of species must be non-negative!')
    end if

    return
  end subroutine get_nUnknowns

end module solver_options


! ============================= !
! Initialize the solver options !
! ============================= !
subroutine solver_options_setup

  ! Internal modules
  use solver_options

  ! External modules
  use simulation_flags
  use geometry
  use parser
  use thermochem, only : get_molecular_weight

  ! Local variables
  integer :: i, nInput
  character(len = str_medium) :: inputString
  real(wp) :: referenceMolecularWeight, atwoodNumber
  logical :: present

  ! Clean slate
  call solver_options_cleanup

  ! Get the name of the simulation
  call parser_read('simulation name', simulationName, 'NULL')

  ! Reference ratio of specific heats
  call parser_read('ratio of specific heats', ratioOfSpecificHeats, 1.4_WP)

  ! Gravity
  if (useGravity) then
     call parser_read('Froude number', froudeNumberInverse, 0.0_WP)
     froudeNumberInverse = max(0.0_WP, froudeNumberInverse)
     if (froudeNumberInverse .gt. 0.0_WP) froudeNumberInverse = 1.0_WP / froudeNumberInverse
  end if

  ! Transport coefficients
  if (useViscosity) then

     call parser_read('Reynolds number', reynoldsNumberInverse, 0.0_WP)
     reynoldsNumberInverse = max(0.0_WP, reynoldsNumberInverse)
     call parser_read('Prandtl number', prandtlNumberInverse, 0.7_WP)
     prandtlNumberInverse = max(0.0_WP, prandtlNumberInverse)
     if (nSpecies .gt. 0) then
        allocate(schmidtNumberInverse(nSpecies + 1))
        call parser_is_defined('Schmidt number', present)
        if (present) then
           call parser_getsize('Schmidt number', nInput)
        else
           nInput = 0
        end if
        select case (nInput)
        case (0)
           ! Default: set Schmidt number equal to Prandtl number (unity Lewis number)
           do i = 1, nSpecies + 1
              schmidtNumberInverse(i) = prandtlNumberInverse
           end do
        case (1)
           ! All species have equal Schmidt numbers
           call parser_read('Schmidt number', schmidtNumberInverse(1))
           do i = 2, nSpecies + 1
              schmidtNumberInverse(i) = schmidtNumberInverse(1)
           end do
        case (2:)
           ! Assign the individual Schmidt numbers
           if (nInput .ne. nSpecies + 1) call die("Error reading value of 'Schmidt number'.&
                & Please input either 1 value or nSpecies + 1 values")
           call parser_read('Schmidt number', schmidtNumberInverse)
        end select

        do i = 1, nSpecies + 1
           schmidtNumberInverse(i) = max(0.0_WP,schmidtNumberInverse(i))
        end do
     end if

     if (reynoldsNumberInverse .le. 0.0_WP .or. prandtlNumberInverse .le. 0.0_WP) then
        viscosityModel = 0
        powerLawExponent = 0.0_WP
        bulkViscosityRatio = 0.0_WP
        if (nSpecies .gt. 0) schmidtNumberInverse = 0.0_WP
     else
        ! Setup transport coefficients
        call parser_read('bulk viscosity ratio', bulkViscosityRatio, 0.6_WP)
        reynoldsNumberInverse = 1.0_WP / reynoldsNumberInverse
        prandtlNumberInverse = 1.0_WP / prandtlNumberInverse
        if (nSpecies .gt. 0) schmidtNumberInverse = 1.0_WP / schmidtNumberInverse
        if (nSpecies .gt. 0 .and. usePTKE)  schmidtNumberInverse = 0.0_WP

        ! Viscosity model
        call parser_read('viscosity model', inputString, 'power law')
        select case (trim(inputString))

        case ('power law', 'POWER LAW')
           viscosityModel = POWER_LAW
           call parser_read('viscosity power law exponent', powerLawExponent, 2.0_WP / 3.0_WP)
           if (powerLawExponent.le.0.0_WP) viscosityModel = 0

        case ('sutherlands law', 'Sutherlands law', 'SUTHERLANDS LAW')
           viscosityModel = SUTHERLANDS_LAW
           call parser_read('sutherland constant', sutherlandConstant, 110.56_WP / 273.11_WP)
           if (.not.predictionOnly) then
              call die("solver_options_setup: adjoints not setup for Sutherland's law!")
           end if

        case ('', 'none', 'NONE')
           viscosityModel = 0

        case default
           call die("solver_options_setup: invalid viscosity model '" //                     &
                trim(inputString) // "'!")
        end select

        ! Adoints are only compatible with repeated first derivative
        if (useSplitViscosity .and. .not.predictionOnly)                                     &
             call die('adjoint requires non-split viscosity!')

     end if

  else

     ! Hack if subgrid-scale models used and inviscid
     if (useLES .or. useShockCapturing) then
        useViscosity = .true.
        viscosityModel = 0
        reynoldsNumberInverse = 0.0_WP
        prandtlNumberInverse = 0.0_WP
        bulkViscosityRatio = 0.0_WP
        powerLawExponent = 0.0_WP
        if (nSpecies .gt. 0) then
           allocate(schmidtNumberInverse(nSpecies + 1))
           schmidtNumberInverse = 0.0_WP
        end if
     end if

  end if

  if (useConstantCfl) then
     call parser_read('cfl', inputCFL, 0.5_WP)
  else
     call parser_read('time step size', inputTimeStepSize)
  end if

  call parser_read('equation of state', inputString, 'ideal gas')
  select case (trim(inputString))
  case ('ideal gas', 'IDEAL GAS')
     equationOfState = IDEAL_GAS

  case ('ideal gas mixture', 'IDEAL GAS MIXTURE')
     if (nSpecies .gt. 0) then

        equationOfState = IDEAL_GAS_MIXTURE

        allocate(speciesName(nSpecies+1))
        allocate(molecularWeightInverse(nSpecies+1))

        ! Check if this is a binary mixture
        call parser_read('Atwood number', atwoodnumber, -1.0_WP)
        if (atwoodNumber .gt. 0.0_WP) then ! ... get molecular weights from the Atwood number
           if (atwoodNumber .ge. 1.0_WP) call die('Atwood number must be < 1')
           if (nSpecies .ne. 1) call die('Using Atwood number requires nSpecies = 1')
           molecularWeightInverse(2) = 1.0_WP ! ... reference species
           molecularWeightInverse(1) = (1.0_WP + atwoodNumber) / (1.0_WP - atwoodNumber)
           molecularWeightInverse = 1.0_WP / molecularWeightInverse
           speciesName(1) = 'fluid1'
        else ! ... read in the molecular weights

           ! Get the molecular weight of the reference species
           call parser_read('reference species', inputString, 'air')
           call get_molecular_weight(trim(inputString), referenceMolecularWeight)
           call parser_getsize('active species', nInput)
           if (nInput .ne. nSpecies) call die("Error reading 'active species'.&
                & Please input nSpecies values")
           call parser_read('active species', speciesName(1:nSpecies))
           do i = 1, nSpecies
              call get_molecular_weight(trim(speciesName(i)), molecularWeightInverse(i))
              molecularWeightInverse(i) = referenceMolecularWeight /                         &
                   molecularWeightInverse(i)
           end do

           ! Get the molecular weight of the inert species
           call parser_read('inert species', speciesName(nSpecies + 1), 'N2')
           call get_molecular_weight(trim(speciesName(nSpecies+1)),                          &
                molecularWeightInverse(nSpecies+1))
           molecularWeightInverse(nSpecies+1) = referenceMolecularWeight /                   &
                molecularWeightInverse(nSpecies+1)

        end if

     else

        call die('solver_options_init: ideal gas mixture requires nSpecies > 0!')

     end if

  case default
     call die("solver_options_setup: invalid equation of state '" //                         &
          trim(inputString) // "'!")

  end select

  ! Shock capturing
  if (useShockCapturing) call parser_read('shock capturing coefficient',                     &
       shockCoefficient, 1.0_WP)

  ! Bound scalars
  if (nSpecies .gt. 0) then
     if (boundScalars) then
        call parser_read('minimum scalar value', minSpecies, -huge(1.0_WP))
        call parser_read('maximum scalar value', maxSpecies, +huge(1.0_WP))
     end if
  else
     boundScalars = .false.
  end if
  
  ! Mass correction
  if (useMassCorrection) call parser_read('target density', targetDensity)

  return
end subroutine solver_options_setup


! ========================== !
! Cleanup the solver options !
! ========================== !
subroutine solver_options_cleanup

  ! Internal modules
  use solver_options

  implicit none

  if (allocated(speciesName)) deallocate(speciesName)
  if (allocated(schmidtNumberInverse)) deallocate(schmidtNumberInverse)
  if (allocated(molecularWeightInverse)) deallocate(molecularWeightInverse)

  return
end subroutine solver_options_cleanup


! ==================== !
! Lookup species index !
! ==================== !
subroutine get_species_index(species, index)

  ! Internal modules
  use solver_options

  implicit none

  ! Arguments
  character(len = *), intent(in) :: species
  integer, intent(out) :: index

  ! Local variables
  integer :: i

  if (nSpecies .le. 0) call die('get_species_index: nSpecies must be > 0!')
  if (.not. allocated(speciesName)) call die('get_species_index: species names not specified')

  index = 0

  do i = 1, nSpecies + 1
     if (trim(species) .eq. speciesName(i)) index = i
  end do

  if (index .eq. 0)  call die('get_species_index: unknown species: ' // trim(species))

  return
end subroutine get_species_index
