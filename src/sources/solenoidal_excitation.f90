module solenoidal_excitation

  ! External modules
  use precision

  implicit none

  integer :: nModes
  real(WP) :: location(2), speed(2), amplitude, gaussianFactor, frequency
  real(WP), allocatable :: angularFrequencies(:), phases(:,:), strength(:),                  &
       spatialFunctionsCache(:,:,:)
  logical :: useExcitation

end module solenoidal_excitation


! =============================== !
! Setup the solenoidal excitation !
! =============================== !
subroutine solenoidal_excitation_setup

  ! Internal modules
  use solenoidal_excitation

  ! External modules
  use parser
  use string
  use math
  use parallel
  use geometry
  use grid

  implicit none

  ! Local variables
  integer :: i, j, n, seed
  integer, allocatable :: seed_(:)
  real(WP) :: distance(2)

  call parser_read('use solenoidal excitation', useExcitation, .false.)

  if (.not. useExcitation) return

  ! Read excitation parameters
  call parser_read('excitation number of modes', nModes)
  call parser_read('excitation x', location(1), 0.0_WP)
  call parser_read('excitation y', location(2), 0.0_WP)
  call parser_read('excitation u', speed(1), 0.0_WP)
  call parser_read('excitation v', speed(2), 0.0_WP)
  call parser_read('excitation amplitude', amplitude)
  call parser_read('excitation frequency', frequency)
  call parser_read('excitation radius', gaussianFactor)
  gaussianFactor = 9.0_WP / (2.0_WP * gaussianFactor ** 2)

  if (nModes .gt. 0) then
     allocate(angularFrequencies(nModes))
     allocate(phases(nModes, 3))
     if (irank.eq.iroot) then
        call parser_read('excitation random seed', seed, -1)
        call random_seed(size = n)
        allocate(seed_(n))
        seed_ = seed
        call random_seed(put = seed_)
        deallocate(seed_)

        do i = 1, 3
           call random_number(phases(:,i))
        end do
        call random_number(angularFrequencies)
        phases = 2.0_WP * pi * phases

        do i = 1, nModes
           angularFrequencies(i) = 2.0_WP * pi * frequency *                                 &
                (real(i, WP) + (angularFrequencies(i) - 0.5_WP)) / (0.5_WP * real(nModes, WP))
        end do
     end if
     call parallel_bc(angularFrequencies)
     call parallel_bc(phases)
  end if

  allocate(strength(nGridPoints))
  strength = 0.0_WP

  if (nModes .gt. 0) allocate(spatialFunctionsCache(nGridPoints, nModes, 4))
  
  do i = 1, nGridPoints

     ! Get minimal distance based on periodicity
     distance(1:2) = coordinates(i,1:2) - location(1:2)
     do j = 1, 2
        if (isPeriodic(j)) then
           if (abs(coordinates(i,j) - periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) - periodicLength(j)         &
                - location(j)
           if (abs(coordinates(i,j) + periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) + periodicLength(j)         &
                - location(j)
        end if
     end do

     strength(i) = amplitude * exp(- gaussianFactor * sum(distance ** 2))

     do j = 1, nModes
        spatialFunctionsCache(i,j,1) = sin(angularFrequencies(j) * distance(1) + phases(j,1))
        spatialFunctionsCache(i,j,2) = cos(angularFrequencies(j) * distance(1) + phases(j,1))
        spatialFunctionsCache(i,j,3) = sin(angularFrequencies(j) * distance(2) + phases(j,2))
        spatialFunctionsCache(i,j,4) = cos(angularFrequencies(j) * distance(2) + phases(j,2))
     end do
  end do

  return
end subroutine solenoidal_excitation_setup


! ============================ !
! Cleanup the ignition sources !
! ============================ !
subroutine solenoidal_excitation_cleanup

  ! Internal modules
  use solenoidal_excitation

  implicit none

  if (allocated(angularFrequencies)) deallocate(angularFrequencies)
  if (allocated(phases)) deallocate(phases)
  if (allocated(strength)) deallocate(strength)
  if (allocated(spatialFunctionsCache)) deallocate(spatialFunctionsCache)

  return
end subroutine solenoidal_excitation_cleanup


! =============================================== !
! Add the ignition sources during the forward run !
! =============================================== !
subroutine solenoidal_excitation_forward(time, source)

  ! Internal modules
  use solenoidal_excitation

  ! External modules
  use math
  use solver_options
  use geometry
  use grid

  implicit none

  ! Arguments
  real(WP), intent(in) :: time
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j
  real(WP) :: distance(2)
  real(WP), allocatable :: temporalFunctions(:,:), temp(:,:,:)

  if (.not. useExcitation) return

  if (nModes .gt. 0) then
     allocate(temporalFunctions(nModes, 4))
     do i = 1, nModes
        temporalFunctions(i,1) = sin(angularFrequencies(i) * speed(1) * time)
        temporalFunctions(i,2) = cos(angularFrequencies(i) * speed(1) * time)
        temporalFunctions(i,3) = sin(angularFrequencies(i) * speed(2) * time)
        temporalFunctions(i,4) = cos(angularFrequencies(i) * speed(2) * time)
     end do
  end if

  allocate(temp(nGridPoints, nModes, 4))
  do i = 1, nModes
     temp(:,i,1) = spatialFunctionsCache(:,i,1) * temporalFunctions(i,2) -                   &
          spatialFunctionsCache(:,i,2) * temporalFunctions(i,1)
     temp(:,i,2) = spatialFunctionsCache(:,i,2) * temporalFunctions(i,2) +                   &
          spatialFunctionsCache(:,i,1) * temporalFunctions(i,1)
     temp(:,i,3) = spatialFunctionsCache(:,i,3) * temporalFunctions(i,4) -                   &
          spatialFunctionsCache(:,i,4) * temporalFunctions(i,3)
     temp(:,i,4) = spatialFunctionsCache(:,i,4) * temporalFunctions(i,4) +                   &
          spatialFunctionsCache(:,i,3) * temporalFunctions(i,3)
  end do

  do i = 1, nGridPoints
     ! Get minimal distance based on periodicity
     distance(1:2) =  coordinates(i,1:2) - location(1:2)
     do j = 1, 2
        if (isPeriodic(j)) then
           if (abs(coordinates(i,j) - periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) - periodicLength(j)         &
                - location(j)
           if (abs(coordinates(i,j) + periodicLength(j) - location(j)) .lt.                  &
                abs(distance(j))) distance(j) = coordinates(i,j) + periodicLength(j)         &
                - location(j)
        end if
     end do

     ! Add the solenoidal excitation source to x-momentum
     source(i, 2) = source(i, 2) + strength(i) * sum(temp(i,:,1) *                           &
          (angularFrequencies * temp(i,:,4) - 2.0_WP * gaussianFactor *                      &
          distance(2) * temp(i,:,3)))

     ! Add the solenoidal excitation source to y-momentum
     source(i, 3) = source(i, 3) - strength(i) * sum(temp(i,:,3) *                           &
          (angularFrequencies * temp(i,:,2) - 2.0_WP * gaussianFactor *                      &
          distance(1) * temp(i,:,1)))
  end do

  return
end subroutine solenoidal_excitation_forward
