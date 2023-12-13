program analyzer

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
  use ibm

  implicit none

  ! Local variables
  integer :: nx, ny, nz
  real(WP), dimension(:), allocatable :: y, meanRho, meanRhoU, norm
  character(len = str_medium) :: input, solutionFile, gridFile

  ! Initialize parallel environment and parse the input file
  call parallel_init
  call parallel_get_inputname(input)
  call parser_init
  call parser_parsefile(input)
  disableManualDecomp = .true.

  if (iRank .eq. iRoot) then
     write (*,*)
     write (*,*) '============================================='
     write (*,*) '| jCODE - Analyzer: postprocess output data |'
     write (*,*) '| This version computes mean stats vs. y    |'
     write (*,*) '============================================='
     write (*,*)
  end if

  ! Read file data
  call parser_read('analyzer solution file', solutionFile)

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

  ! Setup the metrics and levelset
  call grid_metrics_setup
  call grid_levelset_setup
  
  ! Prepare the state
  call allocate_state

  ! Store the total grid points in y
  nx = globalGridsize(1)
  ny = globalGridSize(2)
  nz = globalGridsize(3)

  ! Read in the solution file
  if (iRank .eq. iRoot) print *, 'Reading: ' // trim(solutionFile)
  call simulation_read(IO_FORWARD_STATE, trim(solutionFile))
   
  ! Update the state
  call update_state

  ! We are all set! Let's analyze the data...
  call analyze_data

  ! Finalize the parallel environment
  call parallel_finalize

contains


  ! Compute statistics at a given timestap
  ! --------------------------------------
  subroutine analyze_data

    ! External modules
    use fileio

    implicit none
    
    ! Local variables
    integer :: i, j, k, gridIndex, fileUnit, iostat, ierror
    character(len = str_medium) :: filename

    ! Allocate 1D arrays
    allocate(y(ny)); y = 0.0_WP
    allocate(norm(ny)); norm = 0.0_WP
    allocate(meanRho(ny)); meanRho = 0.0_WP
    allocate(meanRhoU(ny)); meanRhoU = 0.0_WP

    ! Populate the y-vector
    i = iStart(1); k = iStart(3)
    do j = iStart(2), iEnd(2)
       y(j) = coordinates(grid_index(i,j,k), 2)
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, y, ny, MPI_REAL_WP, MPI_SUM, commDir(2), ierror)

    !norm = 0.0_WP
    do k = iStart(3), iEnd(3)
       do j = iStart(2), iEnd(2)
          do i = iStart(1), iEnd(1)
             gridIndex = i - gridOffset(1) + localGridSize(1) *                              &
                  (j - 1 - gridOffset(2) + localGridSize(2) *                                &
                  (k - 1 - gridOffset(3)))
             norm(j) = norm(j) + gridNorm(gridIndex, 1)
             meanRho(j) = meanRho(j) + conservedVariables(gridIndex, 1) *                    &
                  gridNorm(gridIndex, 1)
             meanRhoU(j) = meanRhoU(j) + conservedVariables(gridIndex, 2) *                  &
                  gridNorm(gridIndex, 1)
          end do
       end do
    end do

    ! Average
    call parallel_sum(norm)
    call parallel_sum(meanRho); meanRho = meanRho / norm
    call parallel_sum(meanRhoU); meanRhoU = meanRhoU / norm

    ! Write to file
    filename = 'mean_stats.txt'
    if (iRank .eq. iRoot) then
       fileUnit = iOpen()
       open(unit = fileUnit, file = trim(filename), action = 'write', status = 'unknown',    &
            iostat = iostat)
       write(fileUnit, '(3A20)') "y","<rho>","<rhoU>"
       do j = 1, ny
          write(fileUnit, '(3(ES20.5))') y(j), meanRho(j), meanRhoU(j)                     
       end do
       close(iClose(fileUnit))
    end if
    
    return
  end subroutine analyze_data

end program analyzer


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
