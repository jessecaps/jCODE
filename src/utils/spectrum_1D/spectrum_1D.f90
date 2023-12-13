program spectrum_1D

  ! External modules
  use precision
  use string
  use parser
  use parallel
  use geometry
  use grid
  use grid_functions
  use grid_patch
  use state
  use equation_of_state
  use time_info
  use dump_viz

  implicit none

#ifdef USE_FFTW

  ! Global variables
  integer :: i, j, k, dim, nx, ny, nz, ierror
  real(WP), dimension(:), allocatable :: y
  character(len = str_medium) :: input, filename, gridFile

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '=============================================================='
     write (*,*) '| jCODE - Spectrum: compute 1D spectrum of data at a given y |'
     write (*,*) '=============================================================='
     write (*,*)
  end if

  ! Read file data
  call parser_read('spectrum solution file', filename)

  ! Set up the grid and stencil operators
  call simulation_flags_setup
  call get_dimensions(gridFile)
  call geometry_setup
  call get_nUnknowns
  call solver_options_setup
  call operator_setup
  call grid_setup
  call simulation_read(IO_GRID, trim(gridFile))
  call grid_patch_setup

  ! Setup the metrics and prepare the state
  call grid_metrics_setup
  call allocate_state
  call particle_solver_setup
  call particle_exchange_setup

  ! Only decompose in non-periodic directions
  do dim = 1, nDimensions
     if (nProcsDir(dim).gt.1 .and. isPeriodic(dim)) then
        call die('Cannot deompose procs in periodic direction!')
     end if
  end do

  ! Get local grid size
  nx = globalGridsize(1)
  ny = globalGridSize(2)
  nz = globalGridsize(3)

  ! Store y
  allocate(y(ny)); y = 0.0_WP
  ! Populate the y-vector
  i = iStart(1); k = iStart(3)
  do j = iStart(2), iEnd(2)
     y(j) = coordinates(grid_index(i,j,k), 2)
  end do
  call MPI_ALLREDUCE(MPI_IN_PLACE, y, ny, MPI_REAL_WP, MPI_SUM, commDir(2), ierror)

  ! Read in the solution file
  if (iRank .eq. iRoot) print *, 'Reading: ' // trim(filename)
  call simulation_read(IO_FORWARD_STATE, trim(filename))

  ! Update the state
  call update_state

  ! Compute various spectra
  call spectrum_compute(velocity(:,1), 1, 'spect_E11')
  call spectrum_compute(pressure(:,1), 1, 'spect_p1')
  if (nDimensions.gt.2) then
     call spectrum_compute(velocity(:,3), 3, 'spect_E33')
     call spectrum_compute(pressure(:,1), 3, 'spect_p3')
  end if

  ! Finalize the parallel environment
  call parallel_finalize

contains

  ! Compute spectrum on input data
  ! ------------------------------
  subroutine spectrum_compute(data, dir, name)

    implicit none

    ! Arguments
    integer, intent(in) :: dir
    real(WP), dimension(:), intent(in) :: data
    character(len = *), intent(in) :: name

    ! Local variables
    integer :: n
    real(WP), dimension(:,:), allocatable :: data2D
    real(WP), dimension(:,:), allocatable :: spect,data1D
    real(WP), dimension(:), allocatable :: meanVal, varVal, vol

    ! Determine number of grid points in direction `dir'
    n = globalGridSize(dir)

    ! Allocate data arrays
    allocate(data2D(1:nx,1:nz)); data2D = 0.0_WP
    allocate(data1D(1:n+1,1:ny)); data1D = 0.0_WP
    allocate(spect(1:n+1,2)); spect = 0.0_WP
    allocate(vol(ny)); vol=0.0_WP
    allocate(meanVal(ny)); meanVal=0.0_WP
    allocate(varVal(ny)); varVal=0.0_WP

    ! Compute the mean
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             vol(j) = vol(j) + gridNorm(grid_index(i,j,k), 1)
             meanVal(j) = meanVal(j) + data(grid_index(i,j,k)) *                             &
                  gridNorm(grid_index(i,j,k), 1)
          end do
       end do
    end do
    call parallel_sum(vol)
    call parallel_sum(meanVal)
    meanVal = meanVal / vol

    ! Compute autocorrelation using 1D Fourier transforms
    do j = iStart(2), iEnd(2)
       do k = iStart(3), iend(3)
          do i = iStart(1), iEnd(1)
             data2D(i,k) = data(grid_index(i,j,k)) - meanVal(j)
             varVal(j) = varVal(j) + (data(grid_index(i,j,k)) - meanVal(j))**2 *             &
                  gridNorm(grid_index(i,j,k), 1)
          end do
       end do
       call spectrum_fftw(data2D,spect,n,periodicLength(dir))
       data1D(:,j) = spect(:,2)
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,data1D,(n+1)*ny,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierror)

    ! Normalize the variance
    call parallel_sum(varVal)
    varVal = varVal / vol

    ! Write spectrum to file
    call spectrum_write(spect(:,1), data1D, varVal, trim(name), n)

    ! Cleanup
    deallocate(data2D)
    deallocate(data1D)
    deallocate(spect)
    deallocate(vol)
    deallocate(meanVal)

    return
  end subroutine spectrum_compute


  ! Compute the spectrum via FFTW
  ! -----------------------------
  subroutine spectrum_fftw(A,S,n,L)

    implicit none

    include 'fftw3.f'  

    ! Arguments
    integer, intent(in) :: n
    real(WP), intent(in) :: L
    real(WP), dimension(nx,nz), intent(in) :: A
    real(WP), dimension(n+1,2), intent(out) :: S

    ! Local variables
    complex(WP), dimension(n) :: in,out
    real(WP), dimension(n) :: B
    real(WP) :: pi,dk,kc,eps,kx,kz
    integer(KIND=8) :: plan
    integer :: i,k,ik
    complex(WP), parameter :: ii =(0.0_WP,1.0_WP)

    ! Reset
    S = 0.0_WP

    ! Some parameters
    pi=acos(-1.0_WP)
    dk=2.0_WP*pi/L
    kc=pi*real(n,WP)/L
    eps=kc/1000000.0_WP

    ! Create the plan
    call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)

    if (n.eq.nx) then
       ! FFT in x
       B = 0.0_WP
       do k=1,nz
          in = real(A(:,k),WP)
          call dfftw_execute(plan)
          out = out / real(nx,WP)
          do i=1,nx
             B(i) = B(i) + real(out(i)*conjg(out(i)))
          end do
       end do
       B = B / real(nz,WP)

       do i=1,nx
          ! Wavenumbers
          kx=real(i-1,WP)*dk
          if (i.gt.(nx/2+1)) kx=-real(nx-i+1,WP)*dk
          ik=1+idint(kx/dk+0.5_WP)
          ! Spectrum
          if ((kx.gt.eps).and.(kx.le.kc)) then
             S(ik,2)=S(ik,2)+0.5_WP*B(i)/dk
          end if
       end do
    elseif (n.eq.nz) then
       ! FFT in z
       B = 0.0_WP
       do i=1,nx
          ! FFT
          in = real(A(i,:),WP)
          call dfftw_execute(plan)
          out = out / real(nz,WP)
          do k=1,nz
             B(k) = B(k) + real(out(k)*conjg(out(k)))
          end do
       end do
       B = B / real(nx,WP)

       do k=1,nz
          ! Wavenumbers
          kz=real(k-1,WP)*dk
          if (k.gt.(nz/2+1)) kz=-real(nz-k+1,WP)*dk
          ik=1+idint(kz/dk+0.5_WP)
          ! Spectrum
          if ((kz.gt.eps).and.(kz.le.kc)) then
             S(ik,2)=S(ik,2)+0.5_WP*B(k)/dk
          end if
       end do
    else
       call die('WARNING: Attempting to FFT on inhomogeneous direction!')
    end if

    ! Return the wave number
    do ik=1,n+1
       S(ik,1)=dk*real(ik-1,WP)
    end do

    ! Clean up
    call dfftw_destroy_plan(plan)

    return
  end subroutine spectrum_fftw


  ! Write spectrum to file
  ! ----------------------
  subroutine spectrum_write(kk,data,variance,name,n)
    use fileio
    implicit none

    ! Arguments
    real(WP), dimension(n+1), intent(in) :: kk
    real(WP), dimension(n+1,ny), intent(in) :: data
    real(WP), dimension(ny), intent(in) :: variance
    character(len=*), intent(in) :: name

    ! Local variables
    integer :: i,j,iunit
    integer, intent(in) :: n

    ! Only the root process writes
    if (irank.ne.iroot) return

    ! Write spectra
    print *, 'Writing ',trim(adjustl(name))
    iunit=iopen()
    open(iunit,file=trim(name),form='formatted')
    do i=1,n+1
       if (abs(kk(i)).gt.epsilon(1.0_WP) .and. maxval(abs(data(i,:))).gt.epsilon(1.0_WP)) then
          write(iunit,*) kk(i),',',(data(i,j), j=1,ny)
       end if
    end do
    close(iClose(iunit))

    ! Write the corresponding variance
    iunit=iopen()
    open(iunit,file=trim(name)//'_var',form='formatted')
    do j=1,ny
       write(iunit,*) y(j), variance(j)
    end do
    close(iClose(iunit))

    return
  end subroutine spectrum_write

#endif
end program spectrum_1D


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
