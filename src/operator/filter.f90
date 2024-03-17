module filter

  ! External modules
  use precision
  use operator

  implicit none

  type(t_StencilOperator), allocatable, dimension(:) :: filterStencil, gaussianFilter
  logical :: useFilter

contains

  subroutine filter_apply(direction, x)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('filter_apply: size of `x` inconsistent with local grid size')

    ! Explicit step
    select case (direction)
    case (1)
       call operator_apply_1(filterStencil(1), x)
    case (2)
       call operator_apply_2(filterStencil(2), x)
    case (3)
       call operator_apply_3(filterStencil(3), x)
    end select

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(filterStencil(1), x)
    case (2)
       call operator_apply_implicit_2(filterStencil(2), x)
    case (3)
       call operator_apply_implicit_3(filterStencil(3), x)
    end select

    return
  end subroutine filter_apply

  subroutine gaussian_filter_apply(direction, x)

    ! External modules
    use geometry

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('filter_apply: size of `x` inconsistent with local grid size')

    select case (direction)
    case (1)
       call operator_apply_1(gaussianFilter(1), x)
    case (2)
       call operator_apply_2(gaussianFilter(2), x)
    case (3)
       call operator_apply_3(gaussianFilter(3), x)
    end select

    return
  end subroutine gaussian_filter_apply

end module filter


! ========================= !
! Setup the filter operator !
! ========================= !
subroutine filter_setup

  ! Internal modules
  use filter

  ! External modules
  use string
  use parser
  use simulation_flags
  use operator, only : discretizationType

  implicit none

  ! Local variables
  integer :: i, order
  real(WP) :: alpha, a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
  character(len = str_medium) :: stencilScheme

  ! Decide if filter is used
  call parser_read('filter solution', useFilter, .false.)
  if (useFilter) allocate(filterStencil(nDimensions))
  if (useShockCapturing .or. useIBM .or. useParticles) allocate(gaussianFilter(nDimensions))
  if (.not.allocated(filterStencil) .and. .not.allocated(gaussianFilter)) return

  ! Setup the filter to be used on the solution
  if (allocated(filterStencil)) then
     do i = 1, nDimensions

        if (globalGridSize(i) .gt. 1) then
           call parser_read('filtering scheme', stencilScheme)
           if (trim(stencilScheme) .eq. 'implicit') then
              ! Pade-type implicit filter (tri-diagonal), details can be found in
              ! Gaitonde & Visbal (1999) AIAA.
              call parser_read('filter parameter', alpha, 0.49_WP)
              if (abs(alpha) .gt. 0.5_WP)                                                    &
                   call die('filter_setup: filter parameter must be <= 0.5')
              call parser_read('filter order', order, 10)
              stencilScheme = trim(stencilScheme) // ' order ' // num2str(order)
           end if
        else
           stencilScheme = 'null matrix'
        end if

        select case (trim(stencilScheme))

        case ('implicit order 10')

           ! Determine if the operator is implicit
           if (abs(alpha) .gt. 0.0_WP) then
              filterStencil(i)%implicit = .true.
           else
              filterStencil(i)%implicit = .false.
           end if

           filterStencil(i)%symmetryType = SYMMETRIC
           filterStencil(i)%interiorWidth = 11
           filterStencil(i)%boundaryWidth = 11
           filterStencil(i)%boundaryDepth = 5
           filterStencil(i)%nDiagonals    = 3
           call allocate_operator(filterStencil(i))

           a0 = (193.0_WP + 126.0_WP * alpha) / 256.0_WP
           a1 = (105.0_WP + 302.0_WP * alpha) / 256.0_WP
           a2 = 15.0_WP * (-1.0_WP + 2.0_WP * alpha) / 64.0_WP
           a3 = 45.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a4 = 5.0_WP * (-1.0_WP + 2.0_WP * alpha) / 256.0_WP
           a5 = (1.0_WP - 2.0_WP * alpha) / 512.0_WP

           filterStencil(i)%rhsInterior(0:5) = 0.5_WP * (/ 2.0_WP * a0, a1, a2, a3, a4, a5 /)
           filterStencil(i)%rhsInterior(-1:-5:-1) = filterStencil(i)%rhsInterior(1:5)

           filterStencil(i)%normBoundary = 1.0_WP

           ! Left-hand side coefficients
           filterStencil(i)%lhsInterior(-1:+1) = (/ alpha, 1.0_WP, alpha /)

           ! (i = 1) on boundary filter: no filtering
           filterStencil(i)%rhsBoundary1(1:1,1) = (/ 1.0_WP /)
           filterStencil(i)%lhsBoundary1(0:0,1) = (/ 1.0_WP /)

           ! (i = 2) boundary filter: from Gaitonde & Visbal (1999)
           a1 = (1.0_WP + 1022.0_WP * alpha) / 1024.0_WP
           a2 = (507.0_WP + 10.0_WP * alpha) / 512.0_WP
           a3 = (45.0_WP + 934.0_WP * alpha) / 1024.0_WP
           a4 = 15.0_WP * (-1.0_WP + 2.0_WP * alpha) / 128.0_WP
           a5 = 105.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a6 = 63.0_WP * (-1.0_WP + 2.0_WP * alpha) / 256.0_WP
           a7 = 105.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a8 = 15.0_WP * (-1.0_WP + 2.0_WP * alpha) / 128.0_WP
           a9 = 45.0_WP * (1.0_WP - 2.0_WP * alpha) / 1024.0_WP
           a10= 5.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a11= (1.0_WP - 2.0_WP * alpha) / 1024.0_WP
           filterStencil(i)%rhsBoundary1(1:11,2) =                                           &
                (/ a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 /)
           filterStencil(i)%lhsBoundary1(-1:+1,2) = (/ alpha, 1.0_WP, alpha /)

           ! (i = 3) boundary filter: from Gaitonde & Visbal (1999)
           a1 = (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           a2 = (5.0_WP + 502.0_WP * alpha) / 512.0_WP
           a3 = (979.0_WP + 90.0_WP * alpha) / 1024.0_WP
           a4 = (15.0_WP + 98.0_WP * alpha) / 128.0_WP
           a5 = 105.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a6 = 63.0_WP * (1.0_WP - 2.0_WP * alpha) / 256.0_WP
           a7 = 105.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a8 = 15.0_WP * (1.0_WP - 2.0_WP * alpha) / 128.0_WP
           a9 = 45.0_WP * (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           a10= 5.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a11= (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           filterStencil(i)%rhsBoundary1(1:11,3) =                                           &
                (/ a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 /)
           filterStencil(i)%lhsBoundary1(-1:+1,3) = (/ alpha, 1.0_WP, alpha /)

           ! (i = 4) boundary filter: from Gaitonde & Visbal (1999)
           a1 = (1.0_WP - 2.0_WP * alpha) / 1024.0_WP
           a2 = 5.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a3 = (45.0_WP + 934.0_WP * alpha) / 1024.0_WP
           a4 = (113.0_WP + 30.0_WP * alpha) / 128.0_WP
           a5 = (105.0_WP + 302.0_WP * alpha) / 512.0_WP
           a6 = 63.0_WP * (-1.0_WP + 2.0_WP * alpha) / 256.0_WP
           a7 = 105.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a8 = 15.0_WP * (-1.0_WP + 2.0_WP * alpha) / 128.0_WP
           a9 = 45.0_WP * (1.0_WP - 2.0_WP * alpha) / 1024.0_WP
           a10= 5.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a11= (1.0_WP - 2.0_WP * alpha) / 1024.0_WP
           filterStencil(i)%rhsBoundary1(1:11,4) =                                           &
                (/ a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 /)
           filterStencil(i)%lhsBoundary1(-1:+1,4) = (/ alpha, 1.0_WP, alpha /)

           ! (i = 5) boundary filter: from Gaitonde & Visbal (1999)
           a1 = (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           a2 = 5.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a3 = 45.0_WP * (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           a4 = (15.0_WP + 98.0_WP * alpha) / 128.0_WP
           a5 = (407.0_WP + 210.0_WP * alpha) / 512.0_WP
           a6 = (63.0_WP + 130.0_WP * alpha) / 256.0_WP
           a7 = 105.0_WP * (-1.0_WP + 2.0_WP * alpha) / 512.0_WP
           a8 = 15.0_WP * (1.0_WP - 2.0_WP * alpha) / 128.0_WP
           a9 = 45.0_WP * (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           a10= 5.0_WP * (1.0_WP - 2.0_WP * alpha) / 512.0_WP
           a11= (-1.0_WP + 2.0_WP * alpha) / 1024.0_WP
           filterStencil(i)%rhsBoundary1(1:11,5) =                                           &
                (/ a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 /)
           filterStencil(i)%lhsBoundary1(-1:+1,5) = (/ alpha, 1.0_WP, alpha /)

        case ('DRP 9-point')

           filterStencil(i)%implicit = .false.
           filterStencil(i)%symmetryType = SYMMETRIC
           filterStencil(i)%interiorWidth = 9
           filterStencil(i)%boundaryWidth = 7
           filterStencil(i)%boundaryDepth = 4
           call allocate_operator(filterStencil(i))

           filterStencil(i)%rhsInterior(0:4) = (/ 0.75647250688_WP, 0.204788880640_WP,       &
                -0.120007591680_WP, 0.045211119360_WP,                                       &
                -0.008228661760_WP /)
           filterStencil(i)%rhsInterior(-1:-4:-1) = filterStencil(i)%rhsInterior(1:4)

           filterStencil(i)%normBoundary = 1.0_WP

           filterStencil(i)%rhsBoundary1(1:1,1) = (/ 1.0_WP /)
           filterStencil(i)%rhsBoundary1(1:3,2) = (/ 1.0_WP / 4.0_WP, &
                1.0_WP / 2.0_WP, &
                1.0_WP / 4.0_WP /)
           filterStencil(i)%rhsBoundary1(3:5,3) = (/ 5.0_WP / 8.0_WP, &
                1.0_WP / 4.0_WP, &
                -1.0_WP / 16.0_WP /)
           filterStencil(i)%rhsBoundary1(2:1:-1,3) = filterStencil(i)%rhsBoundary1(4:5,3)
           filterStencil(i)%rhsBoundary1(4:7,4) = (/ 11.0_WP / 16.0_WP, 15.0_WP / 64.0_WP,   &
                -3.0_WP / 32.0_WP, 1.0_WP / 64.0_WP /)
           filterStencil(i)%rhsBoundary1(3:1:-1,4) = filterStencil(i)%rhsBoundary1(5:7,4)

        case ('DRP 13-point')

           filterStencil(i)%implicit = .false.
           filterStencil(i)%symmetryType = SYMMETRIC
           filterStencil(i)%interiorWidth = 13
           filterStencil(i)%boundaryWidth = 11
           filterStencil(i)%boundaryDepth = 6
           call allocate_operator(filterStencil(i))

           filterStencil(i)%rhsInterior(0:6) = (/ 0.809100488494_WP, 0.171503832236_WP,      &
                -0.123632891797_WP, 0.069975429105_WP, &
                -0.029662754736_WP, 0.008520738659_WP, &
                -0.001254597714_WP /)
           filterStencil(i)%rhsInterior(-1:-6:-1) = filterStencil(i)%rhsInterior(1:6)

           filterStencil(i)%normBoundary = 1.0_WP

           filterStencil(i)%rhsBoundary1(1:4,1) = (/ 0.679117647059_WP, 0.465_WP,            &
                -0.179117647059_WP, 0.035_WP /)

           filterStencil(i)%rhsBoundary1(1:7,2) = (/ 0.085777408970_WP, 0.722371828476_WP,   &
                0.356848072173_WP, -0.223119093072_WP, &
                0.057347064865_WP, 0.000747264596_WP,  &
                0.000027453993_WP /)

           filterStencil(i)%rhsBoundary1(1:7,3) = (/ -0.032649010764_WP, 0.143339502575_WP,  &
                0.72667882202_WP, 0.294622121167_WP,  &
                -0.186711738069_WP, 0.062038376258_WP, &
                -0.007318073189_WP /)

           filterStencil(i)%rhsBoundary1(1:11,4) = (/ 0.000054596010_WP, -0.042124772446_WP, &
                0.173103107841_WP, 0.700384128648_WP,  &
                0.276543612935_WP, -0.131223506571_WP, &
                0.023424966418_WP, -0.013937561779_WP, &
                0.024565095706_WP, -0.013098287852_WP, &
                0.002308621090_WP /)

           filterStencil(i)%rhsBoundary1(1:11,5) = (/ -0.008391235145_WP, 0.047402506444_WP, &
                -0.121438547725_WP, 0.200063042812_WP, &
                0.759930952164_WP, 0.207269200140_WP, &
                -0.122263107844_WP, 0.047121062819_WP, &
                -0.009014891495_WP, -0.001855812216_WP,&
                0.001176830044_WP /)

           filterStencil(i)%rhsBoundary1(6:11,6) = (/ 0.784955115888_WP, 0.187772883589_WP,  &
                -0.123755948787_WP, 0.059227575576_WP, &
                -0.018721609157_WP, 0.002999540835_WP /)
           filterStencil(i)%rhsBoundary1(5:1:-1,6) = filterStencil(i)%rhsBoundary1(7:11,6)

        case ('null matrix')

           filterStencil(i)%implicit = .false.
           filterStencil(i)%symmetryType = SYMMETRIC
           filterStencil(i)%interiorWidth = 0
           filterStencil(i)%boundaryWidth = 1
           filterStencil(i)%boundaryDepth = 1
           call allocate_operator(filterStencil(i))

           filterStencil(i)%rhsInterior(0:0) = 0.0_WP
           filterStencil(i)%rhsBoundary1 = 0.0_WP

        case default

           call die('filter_setup: unknown stencil scheme: ' // trim(stencilScheme))

        end select

        ! Fill the right-boundary coefficients
        if (allocated(filterStencil(i)%rhsBoundary1) .and.                                   &
             allocated(filterStencil(i)%rhsBoundary2)) then
           select case (filterStencil(i)%symmetryType)
           case (SYMMETRIC)
              filterStencil(i)%rhsBoundary2(1:filterStencil(i)%boundaryWidth,:) =            &
                   +filterStencil(i)%rhsBoundary1(filterStencil(i)%boundaryWidth:1:-1,:)
           case (SKEW_SYMMETRIC)
              filterStencil(i)%rhsBoundary2(1:filterStencil(i)%boundaryWidth,:) =            &
                   -filterStencil(i)%rhsBoundary1(filterStencil(i)%boundaryWidth:1:-1,:)
           end select
        end if
        if (allocated(filterStencil(i)%lhsBoundary1) .and.                                   &
             allocated(filterStencil(i)%lhsBoundary2)) then
           filterStencil(i)%lhsBoundary2(                                                    &
                -filterStencil(i)%nDiagonals/2:filterStencil(i)%ndiagonals/2,:) =            &
                filterStencil(i)%lhsBoundary1(                                               &
                filterStencil(i)%nDiagonals/2:-filterStencil(i)%nDiagonals/2:-1,:)
        end if

        call operator_update(filterStencil(i), i)
     end do
  end if

  ! Setup Gaussian filter used during shock capturing
  if (allocated(gaussianFilter)) then
     do i = 1, nDimensions
        if (globalGridSize(i) .gt. 1) then
           ! Truncated 9-point Gaussian filter of Cook & Cabot (JCP 2004)
           ! Boundary treatment is same as Kawai & Lele (JCP 2008)
           gaussianFilter(i)%implicit = .false.
           gaussianFilter(i)%symmetryType = SYMMETRIC
           gaussianFilter(i)%interiorWidth = 9
           gaussianFilter(i)%boundaryWidth = 8
           gaussianFilter(i)%boundaryDepth = 4
           call allocate_operator(gaussianFilter(i))

           a0 = 3565.0_WP / 10368.0_WP
           a1 = 3091.0_WP / 12960.0_WP
           a2 = 1997.0_WP / 25920.0_WP
           a3 = 149.0_WP  / 12960.0_WP
           a4 = 107.0_WP  / 103680.0_WP
           gaussianFilter(i)%rhsInterior(0:4) = (/ a0, a1, a2, a3, a4 /)
           gaussianFilter(i)%rhsInterior(-1:-4:-1) = gaussianFilter(i)%rhsInterior(1:4)

           gaussianFilter(i)%normBoundary = 1.0_WP

           ! (i = 1): all lefthand side points are mirrored
           gaussianFilter(i)%rhsBoundary1(1:5,1) = (/ a0, 2.0_WP*a1, 2.0_WP*a2, 2.0_WP*a3,   &
                2.0_WP*a4 /)

           ! (i = 2): left three points are mirrored
           gaussianFilter(i)%rhsBoundary1(1:6,2) = (/ a1, a0, a1, 2.0_WP*a2, 2.0_WP*a3,      &
                2.0_WP*a4/)

           ! (i = 3): left two points are mirrored
           gaussianFilter(i)%rhsBoundary1(1:7,3) = (/ a2, a1, a0, a1, a2, 2.0_WP*a3, 2.0_WP*a4 /)

           ! (i = 4): left one points are mirrored
           gaussianFilter(i)%rhsBoundary1(1:8,4) = (/ a3, a2, a1, a0, a1, a2, a3, 2.0_WP*a4 /)

        else
           ! Null matrix
           gaussianFilter(i)%implicit = .false.
           gaussianFilter(i)%symmetryType = SYMMETRIC
           gaussianFilter(i)%interiorWidth = 0
           gaussianFilter(i)%boundaryWidth = 1
           gaussianFilter(i)%boundaryDepth = 1
           call allocate_operator(gaussianFilter(i))

           gaussianFilter(i)%rhsInterior(0:0) = 0.0_WP
           gaussianFilter(i)%rhsBoundary1 = 0.0_WP
        end if

        ! Fill the right-boundary coefficients
        if (allocated(gaussianFilter(i)%rhsBoundary1) .and.                                  &
             allocated(gaussianFilter(i)%rhsBoundary2)) then
           gaussianFilter(i)%rhsBoundary2(1:gaussianFilter(i)%boundaryWidth,:) =             &
                +gaussianFilter(i)%rhsBoundary1(gaussianFilter(i)%boundaryWidth:1:-1,:)
        end if
        call operator_update(gaussianFilter(i), i)

     end do
  end if

  return
end subroutine filter_setup


! =========================== !
! Cleanup the filter operator !
! =========================== !
subroutine filter_cleanup

  ! Internal modules
  use filter

  implicit none

  ! Local variables
  integer :: i

  if (allocated(filterStencil)) then
     do i = 1, size(filterStencil)
        call deallocate_operator(filterStencil(i))
     end do
     deallocate(filterStencil)
  end if

  if (allocated(gaussianFilter)) then
     do i = 1, size(gaussianFilter)
        call deallocate_operator(gaussianFilter(i))
     end do
     deallocate(gaussianFilter)
  end if

  return
end subroutine filter_cleanup
