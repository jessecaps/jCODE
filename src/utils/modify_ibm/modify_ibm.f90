program modify_ibm

  ! External modules
  use precision
  use string
  use parser
  use geometry
  use fileio
  use parallel
  use simulation_flags
  use solver_options
  use ibm
  
  implicit none

  type t_Object_old
     real(WP)               :: volume          !... Object volume
     real(WP), dimension(3) :: position        !... Center of mass
     real(WP), dimension(3) :: velocity        !... Translational velocity
     real(WP), dimension(3) :: angularVelocity !... Angular velocity
     real(WP), dimension(3) :: pForce          !... Force due to pressure
     real(WP), dimension(3) :: vForce          !... Force due to viscous stress
     real(WP), dimension(3) :: cForce          !... Force due to collisions
     real(WP), dimension(3) :: hTorque         !... Torque due to hydrodynamic stress
     real(WP), dimension(3) :: cTorque         !... Torque due to collisions
     real(WP), dimension(3) :: dudt            !... RHS for advancing velocity
     real(WP), dimension(3) :: dwdt            !... RHS for advancing ang. velocity
  end type t_Object_old
  type(t_Object_old), dimension(:), allocatable :: objectOld
  integer, parameter :: old_size=248

  character(len = str_medium) :: input, gridFile
  
  ! File data
  character(len=str_medium) :: ibmFile
  integer :: i,iunit,ierr
  integer :: size
  real(WP) :: buf

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  ! Set up the grid and stencil operators
  call simulation_flags_setup
  call get_dimensions(gridFile)
  call geometry_setup
  call get_nUnknowns
  call solver_options_setup

  write (*,"(a20)", advance = 'no') ' IBM file to read : '
  read "(A)", ibmFile
  call BINARY_FILE_OPEN(iunit,trim(ibmFile),"r",ierr)

    ! Read local number of particles
  call binary_file_read(iunit,nObjects,1,kind(nObjects),ierr)
  print *, 'Number of objects: ', nObjects
  call binary_file_read(iunit,size,1,kind(size),ierr)
  if (size.ne.old_size) then
     print *, trim(ibmFILE), ' is wrong type! size:', size,' should be:', old_size
     stop
  end if
  call binary_file_read(iunit,buf,1,kind(buf),ierr)

  allocate(objectOld(nObjects))
  if (nObjects.gt.0) call binary_file_read(iunit, objectOld(1:nObjects), nObjects,          &
       size, ierr)
  
  ! Close the files
  call BINARY_FILE_CLOSE(iunit,ierr)

  ! Set new object type
  allocate(object(nObjects))
  do i = 1, nObjects
     object(i)%volume = objectOld(i)%volume
     object(i)%position = objectOld(i)%position
     object(i)%velocity = objectOld(i)%velocity
     object(i)%angularVelocity= objectOld(i)%angularVelocity
     object(i)%pForce = objectOld(i)%pForce
     object(i)%vForce = objectOld(i)%vForce
     object(i)%cForce = objectOld(i)%cForce
     object(i)%hTorque = objectOld(i)%hTorque
     object(i)%dudt = objectOld(i)%dudt
     object(i)%dwdt = objectOld(i)%dwdt
     object(i)%remove = .false.
  end do

  ! Write the IBM file
  ibmFile = trim(ibmFile) // '_new'
  call simulation_write(IO_IBM, ibmFile)

  ! Finalize the parallel environment
  call parallel_finalize

end program modify_ibm

! -------------------------------------------
subroutine die(errorText)

  ! External modules
  use parallel

  implicit none

  ! Arguments
  character(len = *), intent(in) :: errorText

  call parallel_kill(errorText)

  return
end subroutine die
