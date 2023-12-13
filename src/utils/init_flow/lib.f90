module lib

  ! External modules
  use precision
  use string
  use parser
  use simulation_flags
  use solver_options
  use geometry
  use grid
  use grid_functions

  implicit none

  ! Global variables
  integer :: nx, ny, nz
  real(WP) :: Lx, Ly, Lz, dx, dy, dz, v0

end module lib

subroutine lib_grid

  ! Internal modules
  use lib

  ! External modules
  use parallel
  use math, only : pi

  implicit none

  ! Local variables
  integer :: i, j, k

  ! Simplify
  nx = globalGridSize(1)
  ny = globalGridSize(2)
  nz = globalGridSize(3)

  ! Read in the grid size
  call parser_read('Lx', Lx, 0.0_WP)
  call parser_read('Ly', Ly, 0.0_WP)
  call parser_read('Lz', Lz, 0.0_WP)

  ! Compute the grid spacing
  dx = Lx / real(nx-1, WP)
  dy = Ly / real(ny-1, WP)
  dz = Lz / real(nz-1, WP)

  ! Generate the grid
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X
           if (nx .gt. 1) then
              coordinates(grid_index(i,j,k), 1) = dx * real(i - 1, WP) - Lx * 0.5
           end if

           ! Create X
           if (ny .gt. 1) then
              coordinates(grid_index(i,j,k), 2) = dy * real(j - 1, WP) - Ly * 0.5
           end if


           ! Create Z
           if (nz .gt. 1) then
              coordinates(grid_index(i,j,k), 3) = dz * real(k - 1, WP) - Lz * 0.5
           end if

        end do
     end do
  end do

  ! Setup the grid metrics
  call grid_metrics_setup

end subroutine lib_grid

subroutine lib_data

  ! Internal modules
  use lib

  ! External modules
  use math, only : pi
  use parallel
  use parser
  use state, only : conservedVariables

  implicit none

  ! Local variables
  integer :: i, j, k, l

  ! LIB data reader
  real(WP) :: dummyVal
  character(len = 1000) :: libFilename, dummyline
  integer :: nskip, s

  ! Set the number of conserved variables and allocate the array
  nUnknowns = nDimensions + 2
  allocate(conservedVariables(nGridPoints, nUnknowns))

  !probably add asserts if ND .eq. 3
  !do later

  !add get LIB filename

  ! Zero-out the conserved variables
  conservedVariables = 0.0_WP

  libFilename = 'lib.dat'
  print *, "Reading lib data file: ", trim(libFilename)
  !27 lines of header information (skip)
  nskip=27
  OPEN(1,FILE=trim(libFilename))
  do s=1,nskip
     READ(1,*) dummyline
  end do

  ! Read in the x,y,z,rho,rhou,rhov,rhow,rhoE fields
  do k=1,nz
     do j=1,ny
        do i=1,nx
           l = (k-1)*nx*ny + &
                (j-1)*nx + &
                (i-1+1)

           read(1,*),dummyVal,dummyVal,dummyVal,& !skip x,y,z
                conservedVariables(l,1), & 
                conservedVariables(l,2), &
                conservedVariables(l,3), &
                conservedVariables(l,4), &
                conservedVariables(l,5) 
        end do
     end do
  end do

  close(1)


end subroutine lib_data
