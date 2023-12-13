module acoustic_source

  ! External modules
  use precision

  implicit none

  type :: t_AcousticSource
     real(WP) :: location(3), amplitude, gaussianFactor, angularFrequency, phase
  end type t_AcousticSource

  type(t_AcousticSource), allocatable, dimension(:) :: acousticSources

end module acoustic_source


! ========================== !
! Setup the acoustic sources !
! ========================== !
subroutine acoustic_source_setup

  ! Internal modules
  use acoustic_source

  ! External modules
  use parser
  use string
  use math

  implicit none

  ! Local variables
  integer :: i, n
  real(WP) radius, frequency
  character(len = str_medium) :: key
  
  call parser_read('number of acoustic sources', n, 0)

  if (n .le. 0) return

  allocate(acousticSources(n))

  do i = 1, n
     write(key, '(A,I1.1)') 'acoustic source ', i
     call parser_read(trim(key) // ' x', acousticSources(i)%location(1), 0.0_WP)
     call parser_read(trim(key) // ' y',  acousticSources(i)%location(2), 0.0_WP)
     call parser_read(trim(key) // ' z',  acousticSources(i)%location(3), 0.0_WP)
     call parser_read(trim(key) // ' amplitude',  acousticSources(i)%amplitude)
     call parser_read(trim(key) // ' frequency',  frequency, 1.0_WP)
     call parser_read(trim(key) // ' radius', radius, 1.0_WP)
     call parser_read(trim(key) // ' phase',  acousticSources(i)%phase, 0.0_WP)

     acousticSources(i)%angularFrequency = 2.0_WP * pi * frequency
     acousticSources(i)%gaussianFactor = 9.0_WP / (2.0_WP * radius ** 2)
  end do

  return
end subroutine acoustic_source_setup


! ============================ !
! Cleanup the acoustic sources !
! ============================ !
subroutine acoustic_source_cleanup

  ! Internal modules
  use acoustic_source

  implicit none

  if (allocated(acousticSources)) deallocate(acousticSources)

  return
end subroutine acoustic_source_cleanup


! =========================================== !
! Add acoustic sources during the forward run !
! =========================================== !
subroutine acoustic_source_forward(time, source)

  ! Internal modules
  use acoustic_source
  use solver_options
  use geometry
  use grid

  implicit none

  ! Arguments
  real(WP), intent(in) :: time
  real(WP), dimension(nGridPoints, nUnknowns), intent(inout) :: source

  ! Local variables
  integer :: i, j
  real(WP) :: a, r

  if (.not. allocated(acousticSources)) return

  do i = 1, size(acousticSources)
     a = acousticSources(i)%amplitude * cos(acousticSources(i)%angularFrequency *            &
          time + acousticSources(i)%phase)

     do j = 1, nGridPoints
        r = sum((coordinates(j,:) - acousticSources(i)%location(1:nDimensions)) ** 2)
        source(j,nDimensions+2) = source(j,nDimensions+2) +                                  &
             a * exp(-acousticSources(i)%gaussianFactor * r)
     end do
  end do

  return
end subroutine acoustic_source_forward
