module stat_part

  ! External modules
  use stat
  
  implicit none

  ! Global variables
  integer, parameter :: nbin = 100
  real(WP), dimension(:), allocatable :: bin_Map, bin_Rep, bin_VFp
  real(WP), dimension(:,:), allocatable :: jPDF1, jPDF2

contains

  subroutine compute_part_data(part,Map,Rep,VFp)

    ! External modules
    use solver_options
    use particle
    use particle_exchange

    implicit none

    ! Arguments
    type(t_Particle), intent(in) :: part
    real(WP), intent(out) :: Map, Rep, VFp

    ! Local variables
    real(WP) :: density, temperature, viscosity, velocity(3), stress(3), vf, soundSpeed,       &
         relativeVelocity

    call interpolate_fluid_to_particle(part%gridIndex, part%position, density, temperature,    &
         viscosity, velocity, stress, vf)

    relativeVelocity = vf * sqrt(sum((part%velocity(1:nDimensions) - velocity(1:nDimensions))**2))
    soundSpeed = sqrt((ratioOfSpecificHeats-1.0_WP) * temperature)

    ! Output quantities
    Vfp = 1.0_WP - vf
    Rep = density * part%diameter * relativeVelocity / viscosity
    Map = relativeVelocity / soundSpeed

    return
  end subroutine compute_part_data
  
end module stat_part


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_part_setup

  ! Internal modules
  use stat_part

  ! External modules
  use simulation_flags

  implicit none
  
  ! Return if particles are not used
  if (.not. useParticles) return

  ! Allocate
  allocate(bin_Map(nbin+1))
  allocate(bin_Rep(nbin+1))
  allocate(bin_VFp(nbin+1))
  allocate(jPDF1(nbin,nbin)); jPDF1 = 0.0_WP
  allocate(jPDF2(nbin,nbin)); jPDF2 = 0.0_WP
 
  return
end subroutine stat_part_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_part_cleanup

  ! Internal modules
  use stat_part
  
  implicit none

  if (allocated(bin_Map)) deallocate(bin_Map)
  if (allocated(bin_Rep)) deallocate(bin_Rep)
  if (allocated(bin_VFp)) deallocate(bin_VFp)
  if (allocated(jPDF1)) deallocate(jPDF1)
  if (allocated(jPDF2)) deallocate(jPDF2)

  return
end subroutine stat_part_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_part_compute

  ! Internal modules
  use stat_part

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use particle

  implicit none

  ! Local variables
  integer :: i, iRep, iMap, iVFp
  real(WP) :: maxRep, maxMap, maxVFp, Map, Rep, VFp

  ! Return if particles are not used
  if (.not. useParticles) return

  ! Prepare the particle arrays for interpolation
  call prepare_cartesian_arrays

  ! Determine max values
  maxRep=0.0_WP; maxMap=0.0_WP; maxVFp=0.0_WP
  do i = 1, nParticles
     call compute_part_data(particles(i), Map, Rep, VFp)
     maxMap = max(maxMap, Map)
     maxRep = max(maxRep, Rep)
     maxVFp = max(maxVFp, VFp)
  end do
  call parallel_max(maxMap)
  call parallel_max(maxRep)
  call parallel_max(maxVFp)

  ! Compute bins
  do i = 1, nbin + 1
     bin_Rep(i) = real(i-1, WP) * maxRep / real(nbin, WP)
     bin_Map(i) = real(i-1, WP) * maxMap / real(nbin, WP)
     bin_VFp(i) = real(i-1, WP) * maxVFp / real(nbin, WP)
  end do

  ! Populate the PDF
  jPDF1 = 0.0_WP; jPDF2 = 0.0_WP
  do i = 1, nParticles
     call compute_part_data(particles(i), Map, Rep, VFp)
     
     iRep = floor( Rep * real(nbin, WP) / maxval(bin_Rep) ) + 1
     iRep = max(iRep, 1); iRep = min(iRep, nbin)
     iMap = floor( Map * real(nbin, WP) / maxval(bin_Map) ) + 1
     iMap = max(iMap, 1); iMap = min(iMap, nbin)
     iVFp = floor( VFp * real(nbin, WP) / maxval(bin_VFp) ) + 1
     iVFp = max(iVFp, 1); iVFp = min(iVFp, nbin)

     jPDF1(iVFp, iRep) = jPDF1(iVFp, iRep) + 1.0_WP
     jPDF2(iVFp, iMap) = jPDF2(iVFp, iMap) + 1.0_WP
  end do

  call parallel_sum(jPDF1)
  call parallel_sum(jPDF2)

  return
end subroutine stat_part_compute


! ======================== !
! Read particle statistics !
! ======================== !
subroutine stat_part_read

  ! Internal modules
  use stat_part

  implicit none
  
  ! Nothing to do
  
  return
end subroutine stat_part_read


! ======================== !
! Read particle statistics !
! ======================== !
subroutine stat_part_write

  ! Internal modules
  use stat_part

  ! External modules
  use time_info
  use parallel
  use fileio

  implicit none

  ! Local variables
  character(len=str_medium) :: filename
  integer  :: i, j, iunit, ierror

  ! Root process writes
  if (irank.eq.iroot) then

     ! ------- WRITE THE ASCI FILE ------
     iunit = iopen()
     write(filename,'(a18,I8.8,a4)') "stat/jPDF-VFp-Rep-", timestep, ".txt"
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(3a20)') 'bin-VFp','bin-Rep','jPDF'
     do i = 1, nbin
        do j = 1, nbin
           write(iunit,'(10000ES20.12)') 0.5_WP * (bin_VFp(i) + bin_VFp(i+1)), &
                0.5_WP * (bin_Rep(j) + bin_Rep(j+1)), jPDF1(i,j)
        end do
     end do
     close(iclose(iunit))

     iunit = iopen()
     write(filename,'(a18,I8.8,a4)') "stat/jPDF-VFp-Map-", timestep, ".txt"
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(3a20)') 'bin-VFp','bin-Map','jPDF'
     do i = 1, nbin
        do j = 1, nbin
           write(iunit,'(10000ES20.12)') 0.5_WP * (bin_VFp(i) + bin_VFp(i+1)), &
                0.5_WP * (bin_Map(j) + bin_Map(j+1)), jPDF2(i,j)
        end do
     end do
     close(iclose(iunit)) 
  end if

  ! Log
  call monitor_log("PART STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_part_write
