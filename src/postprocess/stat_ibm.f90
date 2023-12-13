module stat_ibm

  ! External modules
  use stat
  use string
  
  implicit none
 
end module stat_ibm


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_ibm_setup

  ! Internal modules
  use stat_ibm

  ! External modules
  use parser
  use simulation_flags
  use geometry

  implicit none
  
  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Nothing to do
 
  return
end subroutine stat_ibm_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_ibm_cleanup

  ! Internal modules
  use stat_ibm
  
  implicit none

  ! Nothing to do

  return
end subroutine stat_ibm_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_ibm_compute

  ! Internal modules
  use stat_ibm

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use ibm

  implicit none

  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Nothing to do either

  return
end subroutine stat_ibm_compute


! =================== !
! Read ibm statistics !
! =================== !
subroutine stat_ibm_read

  ! Internal modules
  use stat_ibm

  ! External modules
  use parallel
  use ibm

  implicit none

  ! Nothing to do
  
  return
end subroutine stat_ibm_read


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_ibm_write

  ! Internal modules
  use stat_ibm

  ! External modules
  use time_info
  use parallel
  use fileio
  use math
  use geometry
  use ibm

  implicit none

  ! Local variables
  character(len=str_medium) :: filename
  integer  :: i, iunit, ierror
  real(WP) :: D, F(3)

  ! Get forces acting on IBM objects computed on the grid
  call ibm_integrate_forces

  ! Root process writes
  if (iRank .eq. iRoot) then
     ! ------- WRITE THE BINARY FILE ------
    
     ! Open the file
!!$     filename = "stat/stat-ibm"
!!$     write(filename, '(2A,I8.8)') trim(filename), ".", timestep
!!$     call BINARY_FILE_OPEN(iunit,trim(filename),"w",ierror)
!!$
!!$     ! Write the current time
!!$     call BINARY_FILE_WRITE(iunit, time, 1, kind(time), ierror)
!!$     do i = 1, nObjects
!!$     call BINARY_FILE_WRITE(iunit, buf, (nbin+1) * 4, kind(buf), ierror)
!!$
!!$     ! Close the file
!!$     call BINARY_FILE_CLOSE(iunit,ierror)

     ! ------- WRITE THE ASCI FILE ------
     filename = "stat/stat-ibm"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES20.12)') "t=",time
     write(iunit,'(5a20)') 'X','D','Fx','Fy','Fz'
     F = 0.0_WP
     do i = 1, nObjects
        D = (2.0_WP * real(nDimensions, WP) * object(i)%volume / pi)                         &
             ** (1.0_WP / real(nDimensions, WP))
        F(1:nDimensions) =  object(i)%pForce(1:nDimensions) +  object(i)%vForce(1:nDimensions)
        write(iunit,'(10000ES20.12)') object(i)%position(1), D, F(1), F(2), F(3)
     end do
     close(iclose(iunit)) 

  end if

  ! Log
  call monitor_log("IBM STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_ibm_write
