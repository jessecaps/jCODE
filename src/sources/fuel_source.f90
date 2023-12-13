module fuel_source

  ! External modules
  use precision

  implicit none

  type :: t_FuelSource
     integer :: fuelIndex
     real(WP) :: location(3), amplitude, gaussianFactor, angularFrequency, phase
  end type t_FuelSource

  type(t_FuelSource), allocatable, dimension(:) :: fuelSources

end module fuel_source


! ====================== !
! Setup the fuel sources !
! ====================== !
subroutine fuel_source_setup

  ! Internal modules
  use fuel_source

  ! External modules
  use parser
  use string
  use math
  use solver_options

  implicit none

  ! Local variables
  integer :: i, n
  real(WP) :: radius, frequency
  character(len = str_medium) :: key, fuel
  
  call parser_read('number of fuel sources', n, 0)

  if (n .le. 0) return

  if (nSpecies .le. 0) call die('fuel_source_setup: fuel sources require nSpecies > 0!')

  allocate(fuelSources(n))

  do i = 1, n
     write(key, '(A,I1.1)') 'fuel source ', i
     call parser_read(trim(key) // ' fuel', fuel)
     call get_species_index(trim(fuel), fuelSources(i)%fuelIndex)
     call parser_read(trim(key) // ' x', fuelSources(i)%location(1), 0.0_WP)
     call parser_read(trim(key) // ' y',  fuelSources(i)%location(2), 0.0_WP)
     call parser_read(trim(key) // ' z',  fuelSources(i)%location(3), 0.0_WP)
     call parser_read(trim(key) // ' amplitude',  fuelSources(i)%amplitude, 1.0_WP)
     call parser_read(trim(key) // ' phase',  fuelSources(i)%phase, 0.0_WP)
     call parser_read(trim(key) // ' frequency',  frequency, 1.0_WP)
     call parser_read(trim(key) // ' radius', radius, 1.0_WP)
     fuelSources(i)%angularFrequency = 2.0_WP * pi * frequency
     fuelSources(i)%gaussianFactor = 9.0_WP / (2.0_WP * radius ** 2)
  end do

  return
end subroutine fuel_source_setup


! ======================== !
! Cleanup the fuel sources !
! ======================== !
subroutine fuel_source_cleanup

  ! Internal modules
  use fuel_source

  implicit none

  if (allocated(fuelSources)) deallocate(fuelSources)

  return
end subroutine fuel_source_cleanup


! ======================================= !
! Add fuel sources during the forward run !
! ======================================= !
subroutine fuel_source_forward(time, source)

  ! Internal modules
  use fuel_source

  ! External modules
  use math
  use grid
  use solver_options

  implicit none

  ! Arguments
  real(WP), intent(in) :: time
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j
  real(WP) :: a, r

  if (.not. allocated(fuelSources)) return

  do i = 1, size(fuelSources)
     a = fuelSources(i)%amplitude *                                                          &
          abs( cos(fuelSources(i)%angularFrequency * time + fuelSources(i)%phase) )
     do j = 1, nGridPoints
        r = real(sum((coordinates(j,:) - fuelSources(i)%location(1:nDimensions)) ** 2), WP)
        source(j,nDimensions+2+fuelSources(i)%fuelIndex) =                                   &
             source(j,nDimensions+2+fuelSources(i)%fuelIndex) +                              &
             a * exp(-fuelSources(i)%gaussianFactor * r)
     end do
  end do

  return
end subroutine fuel_source_forward
