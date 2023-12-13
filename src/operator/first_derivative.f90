module first_derivative

  ! External modules
  use precision
  use operator

  implicit none

  type(t_StencilOperator), allocatable, dimension(:) :: firstDerivative,                     &
       adjointFirstDerivative, upwindLeft, upwindRight

contains
  
  subroutine first_derivative_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('first_derivative_apply: size of `x` inconsistent with local grid size')

    ! Explicit step
    select case (direction)
    case (1)
       call operator_apply_1(firstDerivative(1), x)
    case (2)
       call operator_apply_2(firstDerivative(2), x)
    case (3)
       call operator_apply_3(firstDerivative(3), x)
    end select

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(firstDerivative(1), x)
    case (2)
       call operator_apply_implicit_2(firstDerivative(2), x)
    case (3)
       call operator_apply_implicit_3(firstDerivative(3), x)
    end select

    return
  end subroutine first_derivative_apply


  subroutine first_derivative_apply_interior(direction, xWithGhostPoints, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(in) :: xWithGhostPoints(:,:,:,:)
    real(WP), intent(out) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints) call die                                                &
         ('first_derivative_apply_interior: size of `x` inconsistent with local grid size')

    ! Explicit step
    select case (direction)
    case (1)
       call operator_apply_interior_1(firstDerivative(1), xWithGhostPoints, x)
    case (2)
       call operator_apply_interior_2(firstDerivative(2), xWithGhostPoints, x)
    case (3)
       call operator_apply_interior_3(firstDerivative(3), xWithGhostPoints, x)
    end select

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(firstDerivative(1), x)
    case (2)
       call operator_apply_implicit_2(firstDerivative(2), x)
    case (3)
       call operator_apply_implicit_3(firstDerivative(3), x)
    end select

    return
  end subroutine first_derivative_apply_interior


  subroutine first_derivative_apply_norm(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints) call die                                                &
         ('first_derivative_apply_norm: size of `x` inconsistent with local grid size')

    select case (direction)
    case (1)
       call operator_apply_norm_1(firstDerivative(1), x)
    case (2)
       call operator_apply_norm_2(firstDerivative(2), x)
    case (3)
       call operator_apply_norm_3(firstDerivative(3), x)
    end select

    return
  end subroutine first_derivative_apply_norm


  subroutine first_derivative_apply_norm_inverse(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints) call die                                                &
         ('first_derivative_apply_norm: size of `x` inconsistent with local grid size')

    select case (direction)
    case (1)
       call operator_apply_norm_inverse_1(firstDerivative(1), x)
    case (2)
       call operator_apply_norm_inverse_2(firstDerivative(2), x)
    case (3)
       call operator_apply_norm_inverse_3(firstDerivative(3), x)
    end select

    return
  end subroutine first_derivative_apply_norm_inverse


  subroutine adjoint_first_derivative_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    select case (direction)
    case (1)
       call operator_apply_1(adjointFirstDerivative(1), x)
    case (2)
       call operator_apply_2(adjointFirstDerivative(2), x)
    case (3)
       call operator_apply_3(adjointFirstDerivative(3), x)
    end select

    return
  end subroutine adjoint_first_derivative_apply


  subroutine adjoint_first_derivative_apply_interior(direction, xWithGhostPoints, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(in) :: xWithGhostPoints(:,:,:,:)
    real(WP), intent(out) :: x(:,:)

    select case (direction)
    case (1)
       call operator_apply_interior_1(adjointFirstDerivative(1), xWithGhostPoints, x)
    case (2)
       call operator_apply_interior_2(adjointFirstDerivative(2), xWithGhostPoints, x)
    case (3)
       call operator_apply_interior_3(adjointFirstDerivative(3), xWithGhostPoints, x)
    end select

    return
  end subroutine adjoint_first_derivative_apply_interior


  subroutine adjoint_first_derivative_apply_norm(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints) call die                                                &
         ('first_derivative_apply_norm: size of `x` inconsistent with local grid size')

    select case (direction)
    case (1)
       call operator_apply_norm_1(adjointFirstDerivative(1), x)
    case (2)
       call operator_apply_norm_2(adjointFirstDerivative(2), x)
    case (3)
       call operator_apply_norm_3(adjointFirstDerivative(3), x)
    end select

    return
  end subroutine adjoint_first_derivative_apply_norm


  subroutine adjoint_first_derivative_apply_norm_inverse(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    select case (direction)
    case (1)
       call operator_apply_norm_inverse_1(adjointFirstDerivative(1), x)
    case (2)
       call operator_apply_norm_inverse_2(adjointFirstDerivative(2), x)
    case (3)
       call operator_apply_norm_inverse_3(adjointFirstDerivative(3), x)
    end select

    return
  end subroutine adjoint_first_derivative_apply_norm_inverse

  
  subroutine adjoint_first_derivative_project_boundary_and_apply(direction, x,               &
       faceOrientation)

    implicit none

    ! Arguments
    integer, intent(in) :: direction, faceOrientation
    real(WP), intent(inout) :: x(:,:)

    select case (direction)
    case (1)
       call operator_project_boundary_and_apply_1(adjointFirstDerivative(1), x,              &
            faceOrientation)
    case (2)
       call operator_project_boundary_and_apply_2(adjointFirstDerivative(2), x,              &
            faceOrientation)
    case (3)
       call operator_project_boundary_and_apply_3(adjointFirstDerivative(3), x,              &
            faceOrientation)
    end select

    return
  end subroutine adjoint_first_derivative_project_boundary_and_apply


  subroutine upwind_apply(direction, x, y)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:), y(:,:)

    if ((size(x, 1) .ne. nGridPoints) .or. (size(y, 1) .ne. nGridPoints))                    &
         call die ('upwind_apply: size of `x` or `y` inconsistent with local grid size')

    ! Left operator
    select case (direction)
    case (1)
       call operator_apply_1(upwindLeft(1), x)
    case (2)
       call operator_apply_2(upwindLeft(2), x)
    case (3)
       call operator_apply_3(upwindLeft(3), x)
    end select

    ! Right operator
    select case (direction)
    case (1)
       call operator_apply_1(upwindRight(1), y)
    case (2)
       call operator_apply_2(upwindRight(2), y)
    case (3)
       call operator_apply_3(upwindRight(3), y)
    end select

    return
  end subroutine upwind_apply


  subroutine upwind_apply_central(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    ! Local variables
    real(WP), allocatable :: left(:,:), right(:,:)
    integer :: i, j

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('upwind_apply: size of `x` inconsistent with local grid size')

    allocate(left(size(x, 1), size(x,2)))
    allocate(right(size(x, 1), size(x,2)))

    ! Left operator
    left = x
    select case (direction)
    case (1)
       call operator_apply_1(upwindLeft(1), left)
    case (2)
       call operator_apply_2(upwindLeft(2), left)
    case (3)
       call operator_apply_3(upwindLeft(3), left)
    end select

    ! Right operator
    right = x
    select case (direction)
    case (1)
       call operator_apply_1(upwindRight(1), right)
    case (2)
       call operator_apply_2(upwindRight(2), right)
    case (3)
       call operator_apply_3(upwindRight(3), right)
    end select

    do j = 1, size(x,2)
       do i = 1, size(x,1)
          x(i,j) = 0.5_WP * (right(i,j) + left(i,j))
       end do
    end do

    ! Cleanup
    deallocate(left, right)

    return
  end subroutine upwind_apply_central

end module first_derivative


! ==================================== !
! Setup the first derivative operators !
! ==================================== !
subroutine first_derivative_setup

  ! Internal modules
  use first_derivative

  ! External modules
  use string
  use parser
  use simulation_flags
  use operator, only : discretizationType

  implicit none

  ! Local variables
  integer :: i, lb, ub
  real(WP) :: x1, x2, x3, alpha, a, b, c
  character(len = str_medium) :: stencilScheme

  allocate(firstDerivative(nDimensions))
  if (.not. predictionOnly) allocate(adjointFirstDerivative(nDimensions))
  if (useUpwinding) then
     allocate(upwindLeft(nDimensions))
     allocate(upwindRight(nDimensions))
  end if

  do i = 1, nDimensions

     if (globalGridSize(i) .gt. 1) then
        call parser_read('first derivative scheme', stencilScheme, trim(discretizationType))
     else
        stencilScheme = 'null matrix'
     end if

     select case (trim(stencilScheme))

     case ('SBP 1-2')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 3
        firstDerivative(i)%boundaryWidth = 2
        firstDerivative(i)%boundaryDepth = 1
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:1) = (/ 1.0_WP / 2.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-1:-1) = - firstDerivative(i)%rhsInterior(1:1)

        firstDerivative(i)%normBoundary = (/ 1.0_WP / 2.0_WP /)

        firstDerivative(i)%rhsBoundary1(1:2,1) = (/ -1.0_WP, 1.0_WP /)

     case ('SBP 2-4')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 5
        firstDerivative(i)%boundaryWidth = 6
        firstDerivative(i)%boundaryDepth = 4
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:2) = (/ 2.0_WP / 3.0_WP, -1.0_WP / 12.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-2:-1) = - firstDerivative(i)%rhsInterior(1:2)

        firstDerivative(i)%normBoundary = (/ 17.0_WP / 48.0_WP, &
                                             59.0_WP / 48.0_WP, &
                                             43.0_WP / 48.0_WP, &
                                             49.0_WP / 48.0_WP /)

        firstDerivative(i)%rhsBoundary1(1:4,1) = (/ -24.0_WP / 17.0_WP, &
                                                     59.0_WP / 34.0_WP, &
                                                     -4.0_WP / 17.0_WP, &
                                                     -3.0_WP / 34.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:3,2) = (/  -1.0_WP / 2.0_WP,  &
                                                               0.0_WP,  &
                                                      1.0_WP / 2.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:5,3) = (/   4.0_WP / 43.0_WP, &
                                                    -59.0_WP / 86.0_WP, &
                                                                0.0_WP, &
                                                     59.0_WP / 86.0_WP, &
                                                     -4.0_WP / 43.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:6,4) = (/   3.0_WP / 98.0_WP, &
                                                                0.0_WP, &
                                                    -59.0_WP / 98.0_WP, &
                                                                0.0_WP, &
                                                     32.0_WP / 49.0_WP, &
                                                     -4.0_WP / 49.0_WP /)

     case ('SBP 3-6')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 7
        firstDerivative(i)%boundaryWidth = 9
        firstDerivative(i)%boundaryDepth = 6
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:3) = (/ 3.0_WP / 4.0_WP, &
                                                -3.0_WP / 20.0_WP,&
                                                 1.0_WP / 60.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-3:-1) = - firstDerivative(i)%rhsInterior(1:3)

        firstDerivative(i)%normBoundary = (/ 13649.0_WP / 43200.0_WP, &
                                             12013.0_WP /  8640.0_WP, &
                                              2711.0_WP /  4320.0_WP, &
                                              5359.0_WP /  4320.0_WP, &
                                              7877.0_WP /  8640.0_WP, &
                                             43801.0_WP / 43200.0_WP /)

        firstDerivative(i)%rhsBoundary1(1:6,1) = (/  -21600.0_WP /  13649.0_WP, &
                                                     104009.0_WP /  54596.0_WP, &
                                                      30443.0_WP /  81894.0_WP, &
                                                     -33311.0_WP /  27298.0_WP, &
                                                      16863.0_WP /  27298.0_WP, &
                                                     -15025.0_WP / 163788.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:6,2) = (/ -104009.0_WP / 240260.0_WP, &
                                                                        0.0_WP, &
                                                       -311.0_WP /  72078.0_WP, &
                                                      20229.0_WP /  24026.0_WP, &
                                                     -24337.0_WP /  48052.0_WP, &
                                                      36661.0_WP / 360390.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:6,3) = (/  -30443.0_WP / 162660.0_WP, &
                                                        311.0_WP /  32532.0_WP, &
                                                                        0.0_WP, &
                                                     -11155.0_WP /  16266.0_WP, &
                                                      41287.0_WP /  32532.0_WP, &
                                                     -21999.0_WP /  54220.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:7,4) = (/   33311.0_WP / 107180.0_WP, &
                                                     -20229.0_WP /  21436.0_WP, &
                                                        485.0_WP /   1398.0_WP, &
                                                                        0.0_WP, &
                                                       4147.0_WP /  21436.0_WP, &
                                                      25427.0_WP / 321540.0_WP, &
                                                         72.0_WP /   5359.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:8,5) = (/  -16863.0_WP /  78770.0_WP, &
                                                      24337.0_WP /  31508.0_WP, &
                                                     -41287.0_WP /  47262.0_WP, &
                                                      -4147.0_WP /  15754.0_WP, &
                                                                        0.0_WP, &
                                                     342523.0_WP / 472620.0_WP, &
                                                      -1296.0_WP /   7877.0_WP, &
                                                        144.0_WP /   7877.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:9,6) = (/   15025.0_WP / 525612.0_WP, &
                                                     -36661.0_WP / 262806.0_WP, &
                                                      21999.0_WP /  87602.0_WP, &
                                                     -25427.0_WP / 262806.0_WP, &
                                                    -342523.0_WP / 525612.0_WP, &
                                                                        0.0_WP, &
                                                      32400.0_WP /  43801.0_WP, &
                                                      -6480.0_WP /  43801.0_WP, &
                                                        720.0_WP /  43801.0_WP /)


     case ('SBP 4-8')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 9
        firstDerivative(i)%boundaryWidth = 12
        firstDerivative(i)%boundaryDepth = 8
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:4) = (/  4.0_WP /   5.0_WP, &
                                                 -1.0_WP /   5.0_WP, &
                                                  4.0_WP / 105.0_WP, &
                                                 -1.0_WP / 280.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-4:-1) = - firstDerivative(i)%rhsInterior(1:4)

        firstDerivative(i)%normBoundary = (/ 1498139.0_WP / 5080320.0_WP, &
                                             1107307.0_WP /  725760.0_WP, &
                                               20761.0_WP /   80640.0_WP, &
                                             1304999.0_WP /  725760.0_WP, &
                                              299527.0_WP /  725760.0_WP, &
                                              103097.0_WP /   80640.0_WP, &
                                              670091.0_WP /  725760.0_WP, &
                                             5127739.0_WP / 5080320.0_WP /)

        x1 =  541.0_WP / 1000.0_WP
        x2 = - 27.0_WP /  400.0_WP
        x3 =  187.0_WP /  250.0_WP

        firstDerivative(i)%rhsBoundary1(1,1) = -2540160.0_WP / 1498139.0_WP
        firstDerivative(i)%rhsBoundary1(2,1) = 9.0_WP * (2257920.0_WP * x1 +                 &
             11289600.0_WP * x2 + 22579200.0_WP * x3 - 15849163.0_WP) / 5992556.0_WP
        firstDerivative(i)%rhsBoundary1(3,1) = 3.0_WP * (-33868800.0_WP * x1 -               &
             162570240.0_WP * x2 - 304819200.0_WP * x3 + 235236677.0_WP) / 5992556.0_WP
        firstDerivative(i)%rhsBoundary1(4,1) = (609638400.0_WP * x1 +                        &
             2743372800.0_WP * x2 + 4572288000.0_WP * x3 - 3577778591.0_WP) / 17977668.0_WP
        firstDerivative(i)%rhsBoundary1(5,1) = 3.0_WP * (-16934400 * x1 -                    &
             67737600.0_WP * x2 - 84672000.0_WP * x3 + 67906303.0_WP) / 1498139.0_WP
        firstDerivative(i)%rhsBoundary1(6,1) = 105.0_WP * (967680.0_WP * x1 +                &
             2903040.0_WP * x2 - 305821.0_WP) / 5992556.0_WP
        firstDerivative(i)%rhsBoundary1(7,1) = 49.0_WP * (-1244160.0_WP * x1 +               &
             18662400.0_WP * x3 - 13322233.0_WP) / 17977668.0_WP
        firstDerivative(i)%rhsBoundary1(8,1) = 3.0_WP * (-6773760.0_WP * x2 -                &
             33868800.0_WP * x3 + 24839327.0_WP) / 5992556.0_WP

        firstDerivative(i)%rhsBoundary1(1,2) = 9.0_WP * (-2257920.0_WP * x1 -                &
             11289600.0_WP * x2 - 22579200.0_WP * x3 + 15849163.0_WP) / 31004596.0_WP
        firstDerivative(i)%rhsBoundary1(2,2) = 0.0_WP
        firstDerivative(i)%rhsBoundary1(3,2) = 3.0_WP * (7257600.0_WP * x1 +                 &
             33868800.0_WP * x2 + 60963840.0_WP * x3 - 47167457.0_WP) / 2214614.0_WP
        firstDerivative(i)%rhsBoundary1(4,2) = 3.0_WP * (-9676800.0_WP * x1 -                &
             42336000.0_WP * x2 - 67737600.0_WP * x3 + 53224573.0_WP) / 1107307.0_WP
        firstDerivative(i)%rhsBoundary1(5,2) = 7.0_WP * (55987200.0_WP * x1 +                &
             217728000.0_WP * x2 + 261273600.0_WP * x3 - 211102099.0_WP) / 13287684.0_WP
        firstDerivative(i)%rhsBoundary1(6,2) = 3.0_WP * (-11612160.0_WP * x1 -               &
             33868800.0_WP * x2  + 3884117.0_WP) / 2214614.0_WP
        firstDerivative(i)%rhsBoundary1(7,2) = 150.0_WP * (24192.0_WP * x1 -                 &
             338688.0_WP * x3 + 240463.0_WP) / 1107307.0_WP
        firstDerivative(i)%rhsBoundary1(8,2) = (152409600.0_WP * x2 +                        &
             731566080.0_WP * x3 - 536324953.0_WP) / 46506894.0_WP

        firstDerivative(i)%rhsBoundary1(1,3) = (33868800.0_WP * x1 +                         &
             162570240.0_WP * x2 + 304819200.0_WP * x3 - 235236677.0_WP) / 1743924.0_WP
        firstDerivative(i)%rhsBoundary1(2,3) = (-7257600.0_WP * x1 -                         &
             33868800.0_WP * x2 - 60963840.0_WP * x3 + 47167457.0_WP) / 124566.0_WP
        firstDerivative(i)%rhsBoundary1(3,3) = 0.0_WP
        firstDerivative(i)%rhsBoundary1(4,3) = (24192000.0_WP * x1 +                         &
             101606400.0_WP * x2 + 152409600.0_WP * x3 - 120219461.0_WP) / 124566.0_WP
        firstDerivative(i)%rhsBoundary1(5,3) = (-72576000.0_WP * x1 -                        &
             270950400.0_WP * x2 - 304819200.0_WP * x3 + 249289259.0_WP) / 249132.0_WP
        firstDerivative(i)%rhsBoundary1(6,3) = 9.0_WP * (806400.0_WP * x1 +                  &
             2257920.0_WP * x2 - 290167.0_WP) / 41522.0_WP
        firstDerivative(i)%rhsBoundary1(7,3) = 6.0_WP * (-134400.0_WP * x1 +                 &
             1693440.0_WP * x3 - 1191611.0_WP) / 20761.0_WP
        firstDerivative(i)%rhsBoundary1(8,3) = 5.0_WP * (-2257920.0_WP * x2 -                &
             10160640.0_WP * x3 + 7439833.0_WP) / 290654.0_WP

        firstDerivative(i)%rhsBoundary1(1,4) = (-609638400.0_WP * x1 -                       &
             2743372800.0_WP * x2 - 4572288000.0_WP * x3 + 3577778591.0_WP) / 109619916.0_WP
        firstDerivative(i)%rhsBoundary1(2,4) = 3.0_WP * (9676800.0_WP * x1 +                 &
             42336000.0_WP * x2 + 67737600.0_WP * x3 - 53224573.0_WP) / 1304999.0_WP
        firstDerivative(i)%rhsBoundary1(3,4) = 3.0_WP * (-24192000.0_WP * x1 -               &
             101606400.0_WP * x2 - 152409600.0_WP * x3 + 120219461.0_WP) / 2609998.0_WP
        firstDerivative(i)%rhsBoundary1(4,4) = 0.0_WP
        firstDerivative(i)%rhsBoundary1(5,4) = 9.0_WP * (16128000.0_WP * x1 +                &
             56448000.0_WP * x2 + 56448000.0_WP * x3 - 47206049.0_WP) / 5219996.0_WP
        firstDerivative(i)%rhsBoundary1(6,4) = 3.0_WP * (-19353600.0_WP * x1 -               &
             50803200.0_WP * x2 + 7628371.0_WP) / 2609998.0_WP
        firstDerivative(i)%rhsBoundary1(7,4) = 2.0_WP * (10886400.0_WP * x1 -                &
             114307200.0_WP * x3  + 79048289.0_WP) / 3914997.0_WP
        firstDerivative(i)%rhsBoundary1(8,4) = 75.0_WP * (1354752.0_WP * x2 +                &
             5419008.0_WP * x3  - 3952831.0_WP) / 18269986.0_WP

        firstDerivative(i)%rhsBoundary1(1,5) = 3.0_WP * (16934400.0_WP * x1 +                &
             67737600.0_WP * x2 + 84672000.0_WP * x3 - 67906303.0_WP) / 2096689.0_WP
        firstDerivative(i)%rhsBoundary1(2,5) = 7.0_WP * (-55987200.0_WP * x1 -               &
             217728000.0_WP * x2 - 261273600.0_WP * x3 + 211102099.0_WP) / 3594324.0_WP
        firstDerivative(i)%rhsBoundary1(3,5) = 3.0_WP * (72576000.0_WP * x1 +                &
             270950400.0_WP * x2 + 304819200.0_WP * x3 - 249289259.0_WP) / 1198108.0_WP
        firstDerivative(i)%rhsBoundary1(4,5) = 9.0_WP * (-16128000.0_WP * x1 -               &
             56448000.0_WP * x2 - 56448000.0_WP * x3 + 47206049.0_WP) / 1198108.0_WP
        firstDerivative(i)%rhsBoundary1(5,5) = 0.0_WP
        firstDerivative(i)%rhsBoundary1(6,5) = 105.0_WP * (414720.0_WP * x1 +                &
             967680.0_WP * x2 - 165527.0_WP) / 1198108.0_WP
        firstDerivative(i)%rhsBoundary1(7,5) = 15.0_WP * (-967680.0_WP * x1 +                &
             6773760.0_WP * x3 - 4472029.0_WP) / 1198108.0_WP
        firstDerivative(i)%rhsBoundary1(8,5) = (-304819200.0_WP * x2 -                       &
             914457600.0_WP * x3 + 657798011.0_WP) / 25160268.0_WP
        firstDerivative(i)%rhsBoundary1(9,5) = -2592.0_WP / 299527.0_WP

        firstDerivative(i)%rhsBoundary1(1,6)  = 5.0_WP * (-967680.0_WP * x1 -                &
             2903040.0_WP * x2 + 305821.0_WP) / 1237164.0_WP
        firstDerivative(i)%rhsBoundary1(2,6)  = (11612160.0_WP * x1 +                        &
             33868800.0_WP * x2 - 3884117.0_WP) / 618582.0_WP
        firstDerivative(i)%rhsBoundary1(3,6)  = 9.0_WP * (-806400.0_WP * x1 -                &
             2257920.0_WP * x2 + 290167.0_WP) / 206194.0_WP
        firstDerivative(i)%rhsBoundary1(4,6)  = (19353600.0_WP * x1 +                        &
             50803200.0_WP * x2 - 7628371.0_WP) / 618582.0_WP
        firstDerivative(i)%rhsBoundary1(5,6)  = 35.0_WP * (-414720.0_WP * x1 -               &
             967680.0_WP * x2 + 165527.0_WP) / 1237164.0_WP
        firstDerivative(i)%rhsBoundary1(6,6)  = 0.0_WP
        firstDerivative(i)%rhsBoundary1(7,6)  = 80640.0_WP * x1 / 103097.0_WP
        firstDerivative(i)%rhsBoundary1(8,6)  = 80640.0_WP * x2 / 103097.0_WP
        firstDerivative(i)%rhsBoundary1(9,6)  = 3072.0_WP / 103097.0_WP
        firstDerivative(i)%rhsBoundary1(10,6) = -288.0_WP / 103097.0_WP

        firstDerivative(i)%rhsBoundary1(1,7)  = 7.0_WP * (1244160.0_WP * x1 -                &
             18662400.0_WP * x3 + 13322233.0_WP) / 8041092.0_WP
        firstDerivative(i)%rhsBoundary1(2,7)  = 150.0_WP * (-24192.0_WP * x1 +               &
             338688.0_WP * x3 - 240463.0_WP) / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(3,7)  = 54.0_WP * (134400.0_WP * x1 -                &
             1693440.0_WP * x3 + 1191611.0_WP) / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(4,7)  = 2.0_WP * (-10886400.0_WP * x1 +              &
             114307200.0_WP * x3 - 79048289.0_WP) / 2010273.0_WP
        firstDerivative(i)%rhsBoundary1(5,7)  = 15.0_WP * (967680.0_WP * x1 -                &
             6773760.0_WP * x3 + 4472029.0_WP) / 2680364.0_WP
        firstDerivative(i)%rhsBoundary1(6,7)  = -725760.0_WP * x1 / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(7,7)  = 0.0_WP
        firstDerivative(i)%rhsBoundary1(8,7)  = 725760.0_WP * x3 / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(9,7)  = -145152.0_WP / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(10,7) = 27648.0_WP / 670091.0_WP
        firstDerivative(i)%rhsBoundary1(11,7) = -2592.0_WP / 670091.0_WP

        firstDerivative(i)%rhsBoundary1(1,8)  = 3.0_WP * (6773760.0_WP * x2 +                &
             33868800.0_WP * x3 - 24839327.0_WP) / 20510956.0_WP
        firstDerivative(i)%rhsBoundary1(2,8)  = (-152409600.0_WP * x2 -                      &
             731566080.0_WP * x3 + 536324953.0_WP) / 30766434.0_WP
        firstDerivative(i)%rhsBoundary1(3,8)  = 45.0_WP * (2257920.0_WP * x2 +               &
             10160640.0_WP * x3 - 7439833.0_WP) / 10255478.0_WP
        firstDerivative(i)%rhsBoundary1(4,8)  = 75.0_WP * (-1354752.0_WP * x2 -              &
             5419008.0_WP * x3 + 3952831.0_WP) / 10255478.0_WP
        firstDerivative(i)%rhsBoundary1(5,8)  = (304819200.0_WP * x2 +                       &
             914457600.0_WP * x3 - 657798011.0_WP) / 61532868.0_WP
        firstDerivative(i)%rhsBoundary1(6,8)  = -5080320.0_WP * x2 / 5127739.0_WP
        firstDerivative(i)%rhsBoundary1(7,8)  = -5080320.0_WP * x3 / 5127739.0_WP
        firstDerivative(i)%rhsBoundary1(8,8)  = 0.0_WP
        firstDerivative(i)%rhsBoundary1(9,8)  = 4064256.0_WP / 5127739.0_WP
        firstDerivative(i)%rhsBoundary1(10,8) = -1016064.0_WP / 5127739.0_WP
        firstDerivative(i)%rhsBoundary1(11,8) = 193536.0_WP / 5127739.0_WP
        firstDerivative(i)%rhsBoundary1(12,8) = -18144.0_WP / 5127739.0_WP

     case ('SBP 2-2-8')
        ! Summation-by-parts Operators with Minimal Dispersion Error for Accurate and
        ! Efficient Flow Calculations AIAA SciTech 2016

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 9
        firstDerivative(i)%boundaryWidth = 12
        firstDerivative(i)%boundaryDepth = 8
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:4) = (/  0.820909522170212_WP , &
             -0.221666309292960_WP, &
             0.047921361451327_WP, &
             -0.005335246984568_WP /)

        firstDerivative(i)%rhsInterior(-1:-4:-1) = - firstDerivative(i)%rhsInterior(1:4)

        firstDerivative(i)%normBoundary = (/ 0.302148955111481_WP, &
             1.517433450301689_WP, &
             0.285932266776763_WP, &
             1.555010218521753_WP, &
             1.012284063361710_WP, &
             0.625467352483611_WP, &
             1.259534176052987_WP, &
             1.259534176052987_WP /)

        firstDerivative(i)%rhsBoundary1(2,1) = 0.666853094480032_WP
        firstDerivative(i)%rhsBoundary1(3,1) = 0.076479756121075_WP
        firstDerivative(i)%rhsBoundary1(4,1) = -0.606965769217888_WP
        firstDerivative(i)%rhsBoundary1(5,1) = 0.507594157334931_WP
        firstDerivative(i)%rhsBoundary1(6,1) = -0.112550043333316_WP
        firstDerivative(i)%rhsBoundary1(7,1) = -0.055485611063658_WP
        firstDerivative(i)%rhsBoundary1(8,1) = 0.024074415678823_WP

        firstDerivative(i)%rhsBoundary1(3,2) = -0.124693288718265_WP
        firstDerivative(i)%rhsBoundary1(4,2) = 1.849263324563050_WP
        firstDerivative(i)%rhsBoundary1(5,2) = -1.613862462196115_WP
        firstDerivative(i)%rhsBoundary1(6,2) = 0.650971333814076_WP
        firstDerivative(i)%rhsBoundary1(7,2) = -0.083403924642150_WP
        firstDerivative(i)%rhsBoundary1(8,2) = -0.011421888340565_WP

        firstDerivative(i)%rhsBoundary1(4,3) = -0.976617502001401_WP
        firstDerivative(i)%rhsBoundary1(5,3) = 2.003068861125611_WP
        firstDerivative(i)%rhsBoundary1(6,3) = -1.700521224095258_WP
        firstDerivative(i)%rhsBoundary1(7,3) = 0.799572166580460_WP
        firstDerivative(i)%rhsBoundary1(8,3) = -0.173715834206602_WP

        firstDerivative(i)%rhsBoundary1(5,4) = -0.828851116783346_WP
        firstDerivative(i)%rhsBoundary1(6,4) = 2.270511280380005_WP
        firstDerivative(i)%rhsBoundary1(7,4) = -1.645747376085641_WP
        firstDerivative(i)%rhsBoundary1(8,4) = 0.469767265832744_WP

        firstDerivative(i)%rhsBoundary1(6,5) = -1.016252629771842_WP
        firstDerivative(i)%rhsBoundary1(7,5) = 1.584810115859904_WP
        firstDerivative(i)%rhsBoundary1(8,5) = -0.495272799622414_WP
        firstDerivative(i)%rhsBoundary1(9,5) = -0.005335246984568_WP

        firstDerivative(i)%rhsBoundary1(7,6)  = 0.060442326278035_WP
        firstDerivative(i)%rhsBoundary1(8,6)  = -0.010869723751127_WP
        firstDerivative(i)%rhsBoundary1(9,6)  = 0.047921361451327_WP
        firstDerivative(i)%rhsBoundary1(10,6) = -0.005335246984568_WP

        firstDerivative(i)%rhsBoundary1(8,7)  = 0.839267891753151_WP
        firstDerivative(i)%rhsBoundary1(9,7)  = -0.221666309292960_WP
        firstDerivative(i)%rhsBoundary1(10,7) = 0.047921361451327_WP
        firstDerivative(i)%rhsBoundary1(11,7) = -0.005335246984568_WP

        firstDerivative(i)%rhsBoundary1(9,8)  = 0.820909522170212_WP
        firstDerivative(i)%rhsBoundary1(10,8) = -0.221666309292960_WP
        firstDerivative(i)%rhsBoundary1(11,8) = 0.047921361451327_WP
        firstDerivative(i)%rhsBoundary1(12,8) = -0.005335246984568_WP

     case ('STANDARD 2-4')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 5
        firstDerivative(i)%boundaryWidth = 3
        firstDerivative(i)%boundaryDepth = 2
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:2) = (/ 2.0_WP / 3.0_WP, -1.0_WP / 12.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-2:-1) = - firstDerivative(i)%rhsInterior(1:2)

        firstDerivative(i)%normBoundary = 1.0_WP

        firstDerivative(i)%rhsBoundary1(1:3,1) = (/  -3.0_WP / 2.0_WP, &
             2.0_WP, &
             -1.0_WP / 2.0_WP /)
        firstDerivative(i)%rhsBoundary1(1:3,2) = (/  -1.0_WP / 2.0_WP,  &
             0.0_WP,  &
             1.0_WP / 2.0_WP /)

     case ('STANDARD 3-6')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 7
        firstDerivative(i)%boundaryWidth = 4
        firstDerivative(i)%boundaryDepth = 3
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:3) = (/ 3.0_WP / 4.0_WP, &
             -3.0_WP / 20.0_WP,&
             1.0_WP / 60.0_WP /)
        firstDerivative(i)%rhsInterior(-1:-3:-1) = - firstDerivative(i)%rhsInterior(1:3)

        firstDerivative(i)%normBoundary = 1.0_WP

        firstDerivative(i)%rhsBoundary1(1:4,1) = (/ -11.0_WP / 6.0_WP,  &
             3.0_WP,  &
             -3.0_WP / 2.0_WP,  &
             1.0_WP / 3.0_WP /)

        firstDerivative(i)%rhsBoundary1(1:4,2) = (/  -1.0_WP / 3.0_WP,  &
             -1.0_WP / 2.0_WP,  &
             1.0_WP,  &
             -1.0_WP / 6.0_WP /)

        firstDerivative(i)%rhsBoundary1(1:4,3) = (/   1.0_WP / 6.0_WP,  &
             -1.0_WP,  &
             1.0_WP / 2.0_WP,  &
             1.0_WP / 3.0_WP /)

     case ('Pade 3-6')

        firstDerivative(i)%implicit = .true.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 5
        firstDerivative(i)%boundaryWidth = 3
        firstDerivative(i)%boundaryDepth = 2
        firstDerivative(i)%nDiagonals    = 3
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%normBoundary = 1.0_WP

        ! Interior coefficients
        alpha = 1.0_WP / 3.0_WP
        firstDerivative(i)%lhsInterior(-1:+1) = (/ alpha, 1.0_WP, alpha /)

        a = 1.0_WP / 6.0_WP * (alpha + 9.0_WP)
        b = 1.0_WP / 15.0_WP * (32.0_WP * alpha - 9.0_WP)
        firstDerivative(i)%rhsInterior(1:2) = (/  a / 2.0_WP, b / 4.0_WP  /)
        firstDerivative(i)%rhsInterior(-1:-2:-1) = - firstDerivative(i)%rhsInterior(1:2)

        ! (i = 1) on boundary first derivative
        alpha = 2.0_WP
        a = -(3.0_WP + alpha) / 2.0_WP
        b = 2.0_WP
        c = -(1.0_WP - alpha) / 2.0_WP
        firstDerivative(i)%lhsBoundary1(0:1,1) = (/ 1.0_WP, alpha /)
        firstDerivative(i)%rhsBoundary1(1:3,1) = (/ a, b, c /)

        ! (i = 2) boundary first derivative
        alpha = 0.25_WP
        a = 0.75_WP
        firstDerivative(i)%lhsBoundary1(-1:1,2) = (/ alpha, 1.0_WP, alpha /)
        firstDerivative(i)%rhsBoundary1(1:3,2) = (/ -a, 0.0_WP, a /)

     case ('DRP 13 point')
        ! Bogey, C., & Bailly, C. (2004). A family of low dispersive and low dissipative
        ! explicit schemes for flow and noise computations. Journal of Computational physics,
        ! 194 (1), 194-214.


        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = 13
        firstDerivative(i)%boundaryWidth = 11
        firstDerivative(i)%boundaryDepth = 6
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(1:6) = (/ +9.1083989961046158884449223488881E-001_WP, &
                                                 -3.4195287358140815340035015906213E-001_WP, &
                                                 +1.3803907328137065312584595523507E-001_WP, &
                                                 -4.8269781409397925994599952280221E-002_WP, &
                                                 +1.2472818178268679717084650313541E-002_WP, &
                                                 -1.7227229242514893380588708194318E-003_WP /)
        firstDerivative(i)%rhsInterior(-1:-6:-1) = - firstDerivative(i)%rhsInterior(1:6)

        firstDerivative(i)%normBoundary = 1.0_WP

        firstDerivative(i)%rhsBoundary1(1:7,1) =                                             &
                                              (/ -2.2054462468058574498358102539027E+000_WP, &
                                                 +4.7907207504743750847052990415757E+000_WP, &
                                                 -5.1219100502840136759893413993387E+000_WP, &
                                                 +4.3560242991761195202872469263174E+000_WP, &
                                                 -2.6621263984801656044415289901374E+000_WP, &
                                                 +1.0228938290312966274670491142147E+000_WP, &
                                                 -1.8015618311175450219291443872908E-001_WP /)

        firstDerivative(i)%rhsBoundary1(1:7,2) =                                             &
                                              (/ -2.0839945703453454625946405976346E-001_WP, &
                                                 -1.0890523690038639260460742477278E+000_WP, &
                                                 +2.1545870338706711574556654684244E+000_WP, &
                                                 -1.3931686380866877775059715344834E+000_WP, &
                                                 +7.6849925659269384213662516663408E-001_WP, &
                                                 -2.8018214718148207379084200629556E-001_WP, &
                                                 +4.7716320843203324010061213211805E-002_WP /)

        firstDerivative(i)%rhsBoundary1(1:7,3) =                                             &
                                              (/ +4.8653402919582049864719144773788E-002_WP, &
                                                 -4.6718942991066221653035665769951E-001_WP, &
                                                 -4.7718722757375299865233721644263E-001_WP, &
                                                 +1.2742418176766598292851992139565E+000_WP, &
                                                 -5.1750883227457007860879493782542E-001_WP, &
                                                 +1.6506785384999144743480951940672E-001_WP, &
                                                 -2.6077584687248032793239066169436E-002_WP /)

        firstDerivative(i)%rhsBoundary1(1:7,4) =                                             &
                                              (/ -2.6250993375175630031627008311165E-002_WP, &
                                                 +1.8833730683403585345984136657799E-001_WP, &
                                                 -7.9792163354254481682480170822249E-001_WP, &
                                                 +0.0000000000000000000000000000000E+000_WP, &
                                                 +7.9792163354254481682480170822249E-001_WP, &
                                                 -1.8833730683403585345984136657799E-001_WP, &
                                                 +2.6250993375175630031627008311165E-002_WP /)

        firstDerivative(i)%rhsBoundary1(1:9,5) =                                             &
                                              (/ +8.3010422318690440302147315416164E-003_WP, &
                                                 -6.2400883586085285768298993582152E-002_WP, &
                                                 +2.4992644535898403610438199224578E-001_WP, &
                                                 -8.4585440888718839102472592991157E-001_WP, &
                                                 +0.0000000000000000000000000000000E+000_WP, &
                                                 +8.4585440888718839102472592991157E-001_WP, &
                                                 -2.4992644535898403610438199224578E-001_WP, &
                                                 +6.2400883586085285768298993582152E-002_WP, &
                                                 -8.3010422318690440302147315416164E-003_WP /)

        firstDerivative(i)%rhsBoundary1(1:11,6) =                                            &
                                              (/ -2.6494900234774859217489035169923E-003_WP, &
                                                 +2.1641176020156560356254918102549E-002_WP, &
                                                 -9.2327847148325566766361264237684E-002_WP, &
                                                 +2.8922276219461971527120727959842E-001_WP, &
                                                 -8.7477923690750154205960592130903E-001_WP, &
                                                 +0.0000000000000000000000000000000E+000_WP, &
                                                 +8.7477923690750154205960592130903E-001_WP, &
                                                 -2.8922276219461971527120727959842E-001_WP, &
                                                 +9.2327847148325566766361264237684E-002_WP, &
                                                 -2.1641176020156560356254918102549E-002_WP, &
                                                 +2.6494900234774859217489035169923E-003_WP /)

     case ('null matrix')

        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SYMMETRIC
        firstDerivative(i)%interiorWidth = 0
        firstDerivative(i)%boundaryWidth = 1
        firstDerivative(i)%boundaryDepth = 1
        call allocate_operator(firstDerivative(i))

        firstDerivative(i)%rhsInterior(0:0) = 0.0_WP
        firstDerivative(i)%rhsBoundary1 = 0.0_WP

     case default

        call die('first_derivative_setup: unknown stencil scheme: ' // trim(stencilScheme))

     end select

     ! Fill the right-boundary coefficients
     if (allocated(firstDerivative(i)%rhsBoundary1) .and.                                    &
          allocated(firstDerivative(i)%rhsBoundary2)) then
        select case (firstDerivative(i)%symmetryType)
        case (SYMMETRIC)
           firstDerivative(i)%rhsBoundary2(1:firstDerivative(i)%boundaryWidth,:) =           &
                +firstDerivative(i)%rhsBoundary1(firstDerivative(i)%boundaryWidth:1:-1,:)
        case (SKEW_SYMMETRIC)
           firstDerivative(i)%rhsBoundary2(1:firstDerivative(i)%boundaryWidth,:) =           &
                -firstDerivative(i)%rhsBoundary1(firstDerivative(i)%boundaryWidth:1:-1,:)
        end select
     end if
     if (allocated(firstDerivative(i)%lhsBoundary1) .and.                                    &
          allocated(firstDerivative(i)%lhsBoundary2)) then
        firstDerivative(i)%lhsBoundary2(                                                     &
             -firstDerivative(i)%nDiagonals/2:firstDerivative(i)%ndiagonals/2,:) =           &
             firstDerivative(i)%lhsBoundary1(                                                &
             firstDerivative(i)%nDiagonals/2:-firstDerivative(i)%nDiagonals/2:-1,:)
     end if
     
     call operator_update(firstDerivative(i), i)

     ! Adjoint first derivative operator
     if (allocated(adjointFirstDerivative)) then
        if (useContinuousAdjoint .or. trim(stencilScheme) .eq. 'null matrix') then
           adjointFirstDerivative(i)%implicit = .false.
           adjointFirstDerivative(i)%symmetryType = firstDerivative(i)%symmetryType
           adjointFirstDerivative(i)%interiorWidth = firstDerivative(i)%interiorWidth
           adjointFirstDerivative(i)%boundaryWidth = firstDerivative(i)%boundaryWidth
           adjointFirstDerivative(i)%boundaryDepth = firstDerivative(i)%boundaryDepth
           call allocate_operator(adjointFirstDerivative(i))
           adjointFirstDerivative(i)%normBoundary = firstDerivative(i)%normBoundary
           adjointFirstDerivative(i)%rhsInterior = -firstDerivative(i)%rhsInterior
           adjointFirstDerivative(i)%rhsBoundary1 = -firstDerivative(i)%rhsBoundary1
           adjointFirstDerivative(i)%rhsBoundary2 = -firstDerivative(i)%rhsBoundary2
        else
           call operator_adjoint(firstDerivative(i), adjointFirstDerivative(i))
        end if
        call operator_update(adjointFirstDerivative(i), i)
     end if

  end do

  ! Upwind operators
  if (useUpwinding) then
     do i = 1, nDimensions

        if (globalGridSize(i) .gt. 1) then
           call parser_read('upwind scheme', stencilScheme, trim(stencilScheme))
        else
           stencilScheme = 'null matrix'
        end if

        select case (trim(stencilScheme))

        case ('SBP 1-2')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 3
           upwindRight(i)%boundaryDepth = 2
           upwindRight(i)%boundaryWidth = 4
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(0:2))

           upwindRight(i)%rhsInterior(0:2) = (/ -3.0_WP / 2.0_WP, &
                2.0_WP, &
                -1.0_WP / 2.0_WP /)

           upwindRight(i)%normBoundary = (/ 1.0_WP/4.0_WP, 5.0_WP/4.0_WP /)

           upwindRight(i)%rhsBoundary1(1:2,1) = (/ -1.0_WP/4.0_WP,  5.0_WP/4.0_WP /)
           upwindRight(i)%rhsBoundary1(1:2,2) = (/ -1.0_WP/4.0_WP, -5.0_WP/4.0_WP /)

        case ('SBP 1-3')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 4
           upwindRight(i)%boundaryDepth = 2
           upwindRight(i)%boundaryWidth = 4
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-1:2))

           upwindRight(i)%rhsInterior(-1:2) = (/ -1.0_WP / 3.0_WP, &
                -1.0_WP / 2.0_WP, &
                1.0_WP, &
                -1.0_WP / 6.0_WP /)

           upwindRight(i)%normBoundary = (/ 5.0_WP/12.0_WP, 13.0_WP/12.0_WP /)

           upwindRight(i)%rhsBoundary1(1:2,1) = (/ -1.0_WP/12.0_WP,  3.0_WP/4.0_WP /)
           upwindRight(i)%rhsBoundary1(1:2,2) = (/ -5.0_WP/12.0_WP, -5.0_WP/12.0_WP /)

        case ('SBP 2-4')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 5
           upwindRight(i)%boundaryDepth = 4
           upwindRight(i)%boundaryWidth = 7
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-1:3))

           upwindRight(i)%rhsInterior(-1:3) = (/ -1.0_WP /  4.0_WP, &
                -5.0_WP /  6.0_WP, &
                3.0_WP /  2.0_WP, &
                - 1.0_WP /  2.0_WP, &
                1.0_WP / 12.0_WP /)

           upwindRight(i)%normBoundary = (/  49.0_WP / 144.0_WP, &
                61.0_WP /  48.0_WP, &
                41.0_WP /  48.0_WP, &
                149.0_WP / 144.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,1) = (/ -1.0_WP /  48.0_WP, &
                205.0_WP / 288.0_WP, &
                -29.0_WP / 144.0_WP, &
                1.0_WP /  96.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,2) = (/ -169.0_WP / 288.0_WP, &
                -11.0_WP /  48.0_WP, &
                33.0_WP /  32.0_WP, &
                -43.0_WP / 144.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,3) = (/ 11.0_WP / 144.0_WP, &
                -13.0_WP /  32.0_WP, &
                -29.0_WP /  48.0_WP, &
                389.0_WP / 288.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,4) = (/  1.0_WP /  32.0_WP, &
                -11.0_WP / 144.0_WP, &
                -65.0_WP / 288.0_WP, &
                -13.0_WP /  16.0_WP /)


        case ('SBP 2-5')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 6
           upwindRight(i)%boundaryDepth = 4
           upwindRight(i)%boundaryWidth = 7
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-2:3))

           upwindRight(i)%rhsInterior(-2:3) = (/ 1.0_WP / 20.0_WP, &
                -1.0_WP /  2.0_WP, &
                -1.0_WP /  3.0_WP, &
                1.0_WP, &
                -1.0_WP /  4.0_WP, &
                1.0_WP / 30.0_WP /)

           upwindRight(i)%normBoundary = (/ 251.0_WP / 720.0_WP, &
                299.0_WP / 240.0_WP, &
                211.0_WP / 240.0_WP, &
                739.0_WP / 720.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,1) = (/ -1.0_WP /  120.0_WP, &
                941.0_WP / 1440.0_WP, &
                -47.0_WP /  360.0_WP, &
                -7.0_WP /  480.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,2) = (/ -869.0_WP / 1440.0_WP, &
                -11.0_WP /  120.0_WP, &
                25.0_WP /   32.0_WP, &
                -43.0_WP /  360.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,3) = (/   29.0_WP /  360.0_WP, &
                -17.0_WP /   32.0_WP, &
                -29.0_WP /  120.0_WP, &
                1309.0_WP / 1440.0_WP /)

           upwindRight(i)%rhsBoundary1(1:4,4) = (/    1.0_WP /   32.0_WP, &
                -11.0_WP /  360.0_WP, &
                -661.0_WP / 1440.0_WP, &
                -13.0_WP /   40.0_WP /)

        case ('SBP 3-6')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 7
           upwindRight(i)%boundaryDepth = 6
           upwindRight(i)%boundaryWidth = 10
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-2:4))

           upwindRight(i)%rhsInterior(-2:4) = (/  1.0_WP / 30.0_WP, &
                -2.0_WP /  5.0_WP, &
                -7.0_WP / 12.0_WP, &
                4.0_WP /  3.0_WP, &
                -1.0_WP /  2.0_WP, &
                2.0_WP / 15.0_WP, &
                -1.0_WP / 60.0_WP /)

           upwindRight(i)%normBoundary = (/ 13613.0_WP / 43200.0_WP, &
                12049.0_WP /  8640.0_WP, &
                535.0_WP /   864.0_WP, &
                1079.0_WP /   864.0_WP, &
                7841.0_WP /  8640.0_WP, &
                43837.0_WP / 43200.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,1) = (/      -265.0_WP /     128688.0_WP, &
                1146190567.0_WP / 1737288000.0_WP, &
                -1596619.0_WP /   18384000.0_WP, &
                -55265831.0_WP /  579096000.0_WP, &
                26269819.0_WP / 3474576000.0_WP, &
                2464501.0_WP /  144774000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,2) = (/ -1116490567.0_WP / 1737288000.0_WP, &
                -8839.0_WP /     214480.0_WP, &
                190538869.0_WP /  347457600.0_WP, &
                102705469.0_WP /  694915200.0_WP, &
                413741.0_WP /    9651600.0_WP, &
                -191689861.0_WP / 3474576000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,3) = (/   1096619.0_WP /  18384000.0_WP, &
                -135385429.0_WP / 347457600.0_WP, &
                -61067.0_WP /    321720.0_WP, &
                45137333.0_WP /  57909600.0_WP, &
                -253641811.0_WP / 694915200.0_WP, &
                70665929.0_WP / 579096000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,4) = (/    66965831.0_WP /  579096000.0_WP, &
                -208765789.0_WP /  694915200.0_WP, &
                -17623253.0_WP /   57909600.0_WP, &
                -18269.0_WP /      45960.0_WP, &
                410905829.0_WP /  347457600.0_WP, &
                -477953317.0_WP / 1158192000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,5) = (/  -49219819.0_WP / 3474576000.0_WP, &
                293299.0_WP /    9651600.0_WP, &
                26422771.0_WP /  694915200.0_WP, &
                -141938309.0_WP /  347457600.0_WP, &
                -346583.0_WP /     643400.0_WP, &
                2217185207.0_WP / 1737288000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,6) = (/   -2374501.0_WP /  144774000.0_WP, &
                142906261.0_WP / 3474576000.0_WP, &
                -3137129.0_WP /  579096000.0_WP, &
                -29884283.0_WP / 1158192000.0_WP, &
                -630168407.0_WP / 1737288000.0_WP, &
                -3559.0_WP /       6128.0_WP /)

        case ('SBP 3-7')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 8
           upwindRight(i)%boundaryDepth = 6
           upwindRight(i)%boundaryWidth = 10
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-3:4))

           upwindRight(i)%rhsInterior(-3:4) = (/ -1.0_WP / 105.0_WP, &
                1.0_WP /  10.0_WP, &
                -3.0_WP /   5.0_WP, &
                -1.0_WP /   4.0_WP, &
                1.0_WP, &
                -3.0_WP /  10.0_WP, &
                1.0_WP /  15.0_WP, &
                -1.0_WP / 140.0_WP /)

           upwindRight(i)%normBoundary = (/ 19087.0_WP / 60480.0_WP, &
                84199.0_WP / 60480.0_WP, &
                18869.0_WP / 30240.0_WP, &
                37621.0_WP / 30240.0_WP, &
                55031.0_WP / 60480.0_WP, &
                61343.0_WP / 60480.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,1) = (/       -265.0_WP /      300272.0_WP, &
                1587945773.0_WP /  2432203200.0_WP, &
                -1926361.0_WP /    25737600.0_WP, &
                -84398989.0_WP /   810734400.0_WP, &
                48781961.0_WP /  4864406400.0_WP, &
                3429119.0_WP /   202683600.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,2) = (/ -1570125773.0_WP / 2432203200.0_WP, &
                -26517.0_WP /    1501360.0_WP, &
                240029831.0_WP /  486440640.0_WP, &
                202934303.0_WP /  972881280.0_WP, &
                118207.0_WP /   13512240.0_WP, &
                -231357719.0_WP / 4864406400.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,3) = (/    1626361.0_WP /  25737600.0_WP, &
                -206937767.0_WP / 486440640.0_WP, &
                -61067.0_WP /    750680.0_WP, &
                49602727.0_WP /  81073440.0_WP, &
                -43783933.0_WP / 194576256.0_WP, &
                51815011.0_WP / 810734400.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,4) = (/   91418989.0_WP /  810734400.0_WP, &
                -53314099.0_WP /  194576256.0_WP, &
                -33094279.0_WP /   81073440.0_WP, &
                -18269.0_WP /     107240.0_WP, &
                440626231.0_WP /  486440640.0_WP, &
                -365711063.0_WP / 1621468800.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,5) = (/  -62551961.0_WP / 4864406400.0_WP, &
                799.0_WP /      35280.0_WP, &
                82588241.0_WP /  972881280.0_WP, &
                -279245719.0_WP /  486440640.0_WP, &
                -346583.0_WP /    1501360.0_WP, &
                2312302333.0_WP / 2432203200.0_WP /)

           upwindRight(i)%rhsBoundary1(1:6,6) = (/   -3375119.0_WP /  202683600.0_WP, &
                202087559.0_WP / 4864406400.0_WP, &
                -11297731.0_WP /  810734400.0_WP, &
                61008503.0_WP / 1621468800.0_WP, &
                -1360092253.0_WP / 2432203200.0_WP, &
                -10677.0_WP /      42896.0_WP /)

        case ('SBP 4-8')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 9
           upwindRight(i)%boundaryDepth = 8
           upwindRight(i)%boundaryWidth = 13
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-3:5))

           upwindRight(i)%rhsInterior(-3:5) = (/ -1.0_WP / 168.0_WP, &
                1.0_WP /  14.0_WP, &
                -1.0_WP /   2.0_WP, &
                -9.0_WP /  20.0_WP, &
                5.0_WP /   4.0_WP, &
                -1.0_WP /   2.0_WP, &
                1.0_WP /   6.0_WP, &
                -1.0_WP /  28.0_WP, &
                1.0_WP / 280.0_WP /)

           upwindRight(i)%normBoundary = (/ 7489399.0_WP / 25401600.0_WP, &
                5537831.0_WP /  3628800.0_WP, &
                103373.0_WP /   403200.0_WP, &
                261259.0_WP /   145152.0_WP, &
                298231.0_WP /   725760.0_WP, &
                515917.0_WP /   403200.0_WP, &
                3349159.0_WP /  3628800.0_WP, &
                25639991.0_WP / 25401600.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,1) = (/ -16683.0_WP / 63018560.0_WP, &
                29345822969987.0_WP /  43354248537600.0_WP, &
                -2734625426591.0_WP /  40644608004000.0_WP, &
                -113480208109603.0_WP / 780376473676800.0_WP, &
                -830250230261.0_WP /  26012549122560.0_WP, &
                2500519492033.0_WP /  32515686403200.0_WP, &
                2235718279643.0_WP / 390188236838400.0_WP, &
                -388481888477.0_WP /  26543417472000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,2) = (/ -29227665839987.0_WP /   43354248537600.0_WP, &
                -493793.0_WP /         63018560.0_WP, &
                8302717120817.0_WP /   26543417472000.0_WP, &
                3739408501537.0_WP /    9290196115200.0_WP, &
                2684481534461.0_WP /   13935294172800.0_WP, &
                -4450185662513.0_WP /   18580392230400.0_WP, &
                -1221838279381.0_WP /   37160784460800.0_WP, &
                90595000956023.0_WP / 1950941184192000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,3) = (/ 2505689537591.0_WP / 40644608004000.0_WP, &
                -7312922392817.0_WP / 26543417472000.0_WP, &
                -69332623.0_WP /     1323389760.0_WP, &
                10994933811709.0_WP / 18580392230400.0_WP, &
                -9270952411151.0_WP / 18580392230400.0_WP, &
                3191238635141.0_WP / 20644880256000.0_WP, &
                4442211176987.0_WP / 92901961152000.0_WP, &
                -940661365031.0_WP / 32515686403200.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,4) = (/ 118016946570403.0_WP / 780376473676800.0_WP, &
                -4173878828737.0_WP /   9290196115200.0_WP, &
                -7990503962509.0_WP /  18580392230400.0_WP, &
                -207799621.0_WP /      1323389760.0_WP, &
                2044021560341.0_WP /   2477385630720.0_WP, &
                511197701761.0_WP /  18580392230400.0_WP, &
                1237681717213.0_WP /  13935294172800.0_WP, &
                -7784834666617.0_WP / 130062745612800.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,5) = (/ 68609076271.0_WP / 2364777192960.0_WP, &
                -2235651762161.0_WP /  13935294172800.0_WP, &
                6527681584751.0_WP /  18580392230400.0_WP, &
                -1115980068821.0_WP /   2477385630720.0_WP, &
                -55386253.0_WP /       189055680.0_WP, &
                3208334350649.0_WP /   3716078446080.0_WP, &
                -407569013461.0_WP /    844563283200.0_WP, &
                136474842626653.0_WP / 780376473676800.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,6) = (/ -2487637785013.0_WP / 32515686403200.0_WP, &
                4244231077313.0_WP / 18580392230400.0_WP, &
                -1550378843141.0_WP / 20644880256000.0_WP, &
                -5726967564961.0_WP / 18580392230400.0_WP, &
                -1017898941929.0_WP /  3716078446080.0_WP, &
                -526653889.0_WP /     1323389760.0_WP, &
                45241297077547.0_WP / 37160784460800.0_WP, &
                -2279608411897.0_WP /  5080576000500.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,7) = (/ -2164019088443.0_WP / 390188236838400.0_WP, &
                1263196075861.0_WP /  37160784460800.0_WP, &
                -6600697610987.0_WP /  92901961152000.0_WP, &
                556610591687.0_WP /  13935294172800.0_WP, &
                926842346471.0_WP /   9290196115200.0_WP, &
                -18757693936747.0_WP /  37160784460800.0_WP, &
                -584765899.0_WP /      1323389760.0_WP, &
                204646287449.0_WP /    168431424000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,8) = (/ 387091928477.0_WP /   26543417472000.0_WP, &
                -90231551688023.0_WP / 1950941184192000.0_WP, &
                1032404418251.0_WP /   32515686403200.0_WP, &
                3502353445417.0_WP /  130062745612800.0_WP, &
                -15385068876253.0_WP /  780376473676800.0_WP, &
                262499068919.0_WP /   10161152001000.0_WP, &
                -867004691939.0_WP /    1852745664000.0_WP, &
                -85017967.0_WP /        189055680.0_WP /)

        case ('SBP 4-9')

           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 10
           upwindRight(i)%boundaryDepth = 8
           upwindRight(i)%boundaryWidth = 13
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-4:5))

           upwindRight(i)%rhsInterior(-4:5) = (/ 1.0_WP / 504.0_WP, &
                -1.0_WP /  42.0_WP, &
                1.0_WP /   7.0_WP, &
                -2.0_WP /   3.0_WP, &
                -1.0_WP /   5.0_WP, &
                1.0_WP           , &
                -1.0_WP /   3.0_WP, &
                2.0_WP /  21.0_WP, &
                -1.0_WP /  56.0_WP, &
                1.0_WP / 630.0_WP /)

           upwindRight(i)%normBoundary = (/ 1070017.0_WP / 3628800.0_WP, &
                5537111.0_WP / 3628800.0_WP, &
                103613.0_WP /  403200.0_WP, &
                261115.0_WP /  145152.0_WP, &
                298951.0_WP /  725760.0_WP, &
                515677.0_WP /  403200.0_WP, &
                3349879.0_WP / 3628800.0_WP, &
                3662753.0_WP / 3628800.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,1) = (/ -5561.0_WP / 47263920.0_WP, &
                4186300102421.0_WP /   6193464076800.0_WP, &
                -377895002003.0_WP /   5806372572000.0_WP, &
                -16485548951749.0_WP / 111482353382400.0_WP, &
                -113245973003.0_WP /   3716078446080.0_WP, &
                355360297339.0_WP /   4645098057600.0_WP, &
                321012170669.0_WP /  55741176691200.0_WP, &
                -388397049437.0_WP /  26543417472000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,2) = (/ -4178798062421.0_WP /   6193464076800.0_WP, &
                -493793.0_WP /       141791760.0_WP, &
                725405227507.0_WP /   2413037952000.0_WP, &
                3904159533697.0_WP /   9290196115200.0_WP, &
                2483046570341.0_WP /  13935294172800.0_WP, &
                -4336328670953.0_WP /  18580392230400.0_WP, &
                -1258688487061.0_WP /  37160784460800.0_WP, &
                12931584852209.0_WP / 278705883456000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,3) = (/ 363359390003.0_WP /  5806372572000.0_WP, &
                -7539548734577.0_WP / 26543417472000.0_WP, &
                -69332623.0_WP /     2977626960.0_WP, &
                9994352248429.0_WP / 18580392230400.0_WP, &
                -8195655811631.0_WP / 18580392230400.0_WP, &
                7361486640463.0_WP / 61934640768000.0_WP, &
                5539855071347.0_WP / 92901961152000.0_WP, &
                -12898722943.0_WP /   422281641600.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,4) = (/ 16773595838149.0_WP / 111482353382400.0_WP, &
                -372477950627.0_WP /    844563283200.0_WP, &
                -8659050093229.0_WP /  18580392230400.0_WP, &
                -207799621.0_WP /      2977626960.0_WP, &
                1734921317461.0_WP /   2477385630720.0_WP, &
                2530020015841.0_WP /  18580392230400.0_WP, &
                441856623253.0_WP /  13935294172800.0_WP, &
                -115132773073.0_WP /   2654341747200.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,5) = (/ 108449122763.0_WP /   3716078446080.0_WP, &
                -2283566671541.0_WP /  13935294172800.0_WP, &
                6976424333231.0_WP /  18580392230400.0_WP, &
                -440819477447.0_WP /    825795210240.0_WP, &
                -55386253.0_WP /       425375280.0_WP, &
                2479572560009.0_WP /   3716078446080.0_WP, &
                -40258468963.0_WP /    120651897600.0_WP, &
                11808221047099.0_WP / 111482353382400.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,6) = (/ -32231128289.0_WP /   422281641600.0_WP, &
                4244793299753.0_WP / 18580392230400.0_WP, &
                -5173673584463.0_WP / 61934640768000.0_WP, &
                -4848139955041.0_WP / 18580392230400.0_WP, &
                -1506045711689.0_WP /  3716078446080.0_WP, &
                -526653889.0_WP /     2977626960.0_WP, &
                36411368691307.0_WP / 37160784460800.0_WP, &
                -825434105779.0_WP /  2903186286000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,7) = (/ -316459841069.0_WP / 55741176691200.0_WP, &
                1277069729941.0_WP / 37160784460800.0_WP, &
                -6499182375347.0_WP / 92901961152000.0_WP, &
                355606625147.0_WP / 13935294172800.0_WP, &
                1519272420551.0_WP /  9290196115200.0_WP, &
                -2240079855137.0_WP /  3378253132800.0_WP, &
                -584765899.0_WP /     2977626960.0_WP, &
                2301241355533.0_WP /  2382101568000.0_WP /)

           upwindRight(i)%rhsBoundary1(1:8,8) = (/ 387779289437.0_WP /  26543417472000.0_WP, &
                -12908508708209.0_WP / 278705883456000.0_WP, &
                147710908133.0_WP /   4645098057600.0_WP, &
                534025841911.0_WP /  18580392230400.0_WP, &
                -4119981443899.0_WP / 111482353382400.0_WP, &
                279819152779.0_WP /   2903186286000.0_WP, &
                -1510324515533.0_WP /   2382101568000.0_WP, &
                -85017967.0_WP /       425375280.0_WP /)

        case ('SBP DRP 1-3')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 6
           upwindRight(i)%boundaryDepth = 4
           upwindRight(i)%boundaryWidth = 7
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-2:3))

           upwindRight(i)%rhsInterior(-2:3) = (/ -5.0_WP/48.0_WP, &
                29.0_WP / 360.0_WP, &
                -401.0_WP / 360.0_WP, &
                7.0_WP / 5.0_WP, &
                -187.0_WP / 720.0_WP, &
                -1.0_WP / 360.0_WP /)

           upwindRight(i)%normBoundary = (/ 0.407206_WP, 1.05763_WP, 1.07979_WP, 0.955377_WP /)

           upwindRight(i)%rhsBoundary1(1:4,1) = (/ -0.0371647_WP, 0.690309_WP, -0.176330_WP, 0.0231858_WP /)
           upwindRight(i)%rhsBoundary1(1:4,2) = (/ -0.523249_WP, -0.290858_WP, 1.09105_WP, -0.274169_WP /)
           upwindRight(i)%rhsBoundary1(1:4,3) = (/ 0.0651977_WP, -0.431581_WP, -0.677496_WP, 1.30638_WP /)
           upwindRight(i)%rhsBoundary1(1:4,4) = (/ -0.00478431_WP, 0.0321293_WP, 0.133060_WP, -1.03178_WP /)

        case ('SBP DRP 3-6 draft')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 10
           upwindRight(i)%boundaryDepth = 8
           upwindRight(i)%boundaryWidth = 13
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-4:5))

           upwindRight(i)%rhsInterior(-4:5) = (/  -1.0_WP / 168.0_WP, &
                149.0_WP /  3150.0_WP, &
                -199.0_WP / 1575.0_WP, &
                -8.0_WP / 75_WP, &
                -8.0_WP / 9.0_WP, &
                67.0_WP / 45.0_WP, &
                -37.0_WP / 75.0_WP, &
                124.0_WP / 1575.0_WP, &
                139.0_WP / 12600_WP, &
                -1.0_WP / 210.0_WP /)

           upwindRight(i)%normBoundary = (/ 0.294425_WP, &
                1.52829_WP, &
                0.251092_WP, &
                1.80762_WP, &
                0.403131_WP, &
                1.28497_WP, &
                0.920654_WP, &
                1.00981_WP /)

           upwindRight(i)%rhsBoundary1(1:8,1) = (/ -0.00903009_WP, 0.710170_WP, -0.103201_WP, -0.149256_WP, -0.00889415_WP, 0.0775492_WP, -0.00762020_WP, -0.00971753_WP /)
           upwindRight(i)%rhsBoundary1(1:8,2) = (/ -0.678377_WP, -0.0306467_WP, 0.383435_WP, 0.382131_WP, 0.144309_WP, -0.260057_WP, 0.0457140_WP, 0.0134917_WP /)
           upwindRight(i)%rhsBoundary1(1:8,3) = (/ 0.109170_WP, -0.384446_WP, -0.0114785_WP, 0.597174_WP, -0.372833_WP, -0.00555940_WP, 0.0897052_WP, -0.0217330_WP /)
           upwindRight(i)%rhsBoundary1(1:8,4) = (/ 0.145873_WP, -0.369601_WP, -0.597113_WP, -0.0223665_WP, 0.509298_WP, 0.736223_WP, -0.541190_WP, 0.143639_WP /)
           upwindRight(i)%rhsBoundary1(1:8,5) = (/ -0.0410135_WP, -0.0712000_WP, 0.413473_WP, -0.467929_WP, -0.150629_WP, 0.0494434_WP, 0.470024_WP, -0.208439_WP /)
           upwindRight(i)%rhsBoundary1(1:8,6) = (/ -0.0648459_WP, 0.227413_WP, 0.00266513_WP, -0.704684_WP, 0.160351_WP, -0.339735_WP, 0.894812_WP, -0.260975_WP /)
           upwindRight(i)%rhsBoundary1(1:8,7) = (/ 0.0512884_WP, -0.0996057_WP, -0.127011_WP, 0.476373_WP, -0.448256_WP, -0.0480874_WP, -0.833357_WP, 1.43699_WP /)
           upwindRight(i)%rhsBoundary1(1:8,8) = (/ -0.0130649_WP, 0.0179173_WP, 0.0392310_WP, -0.111443_WP, 0.172608_WP, -0.251125_WP, -0.0330888_WP, -0.901590_WP /)

        case ('SBP DRP 2-4')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 8
           upwindRight(i)%boundaryDepth = 6
           upwindRight(i)%boundaryWidth = 10
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-3:4))

           upwindRight(i)%rhsInterior(-3:4) = (/ 2.0_WP/65.0_WP, -443.0_WP/2860.0_WP, 57.0_WP/1430.0_WP, -889.0_WP/858.0_WP, 205.0_WP/143.0_WP, -873.0_WP/2860.0_WP, -133.0_WP/4290.0_WP, 3.0_WP/130.0_WP /)

           upwindRight(i)%normBoundary = &
                (/ 0.25292350710583155211867593031500623221984793e44_WP / 0.80144000202690995277057891212902866905988720e44_WP, &
                0.12366522918918608201850133537536530111621851e44_WP / 0.8904888911410110586339765690322540767332080e43_WP, &
                0.431945550921990834855601327346624756668871e42_WP / 0.684991454723854660487674283870964674410160e42_WP, &
                0.19955473972660125574814490349612056132376685e44_WP / 0.16028800040538199055411578242580573381197744e44_WP, &
                0.159295736114548212938119642e27_WP / 0.176585123967931181608910805e27_WP, &
                0.3004824804395338925078707877e28_WP / 0.2951613866135020453161224430e28_WP /)

           upwindRight(i)%rhsBoundary1(1:6,1) = &
                (/ -0.299337951148667946616390345482513051882664999823e48_WP / 0.13889756675128376391466903126108195863476905063200e50_WP, &
                0.334203102386451445389405940943896029825001960741e48_WP / 0.462991889170945879715563437536939862115896835440e48_WP, &
                -0.154761287365716184513371737750091076641062676551e48_WP / 0.925983778341891759431126875073879724231793670880e48_WP, &
                -0.47551546840572027745239855308325563115118731109e47_WP / 0.1388975667512837639146690312610819586347690506320e49_WP, &
                -0.3249244246682129529368477e25_WP / 0.128425544703949950261026040e27_WP, &
                0.389387454207582946160914171e27_WP / 0.14758069330675102265806122150e29_WP /)

           upwindRight(i)%rhsBoundary1(1:6,2) = &
                (/ -0.163867315630929559163701592975418537505497389e45_WP / 0.243423706188720231185890345708170274508883720e45_WP, &
                -0.30656114770770165599996417750517210097119869413e47_WP / 0.925983778341891759431126875073879724231793670880e48_WP, &
                0.2031111734644347800339460058542491781926664601e46_WP / 0.2967896725454781280227970753441922193050620740e46_WP, &
                0.5440031102252193970520342577039764722606917e43_WP / 0.324564941584960308247853794277560366011844960e45_WP, &
                0.24341883735526123464906853e26_WP / 0.1059510743807587089653464830e28_WP, &
                -0.52575461023873067818412747e26_WP / 0.2951613866135020453161224430e28_WP /)

           upwindRight(i)%rhsBoundary1(1:6,3) = &
                (/ 0.24402651831679141629765081741317534139059529921e47_WP / 0.102887086482432417714569652785986636025754852320e48_WP, &
                -0.115008374618384762928474480582274121580026720269e48_WP / 0.154330629723648626571854479178979954038632278480e48_WP, &
                -0.25337437961186700206692845540531216985830930791e47_WP / 0.308661259447297253143708958357959908077264556960e48_WP, &
                0.68971e5_WP / 0.92432e5_WP, &
                -0.5410e4_WP / 0.92211e5_WP, &
                -0.11219e5_WP / 0.93126e5_WP /)

           upwindRight(i)%rhsBoundary1(1:6,4) = &
                (/ -0.365347013678543672843654601844986483417906215e45_WP / 0.9921254767948840279619216518648711331054932188e46_WP, &
                0.70833766889864192195958289226074523178045259801e47_WP / 0.925983778341891759431126875073879724231793670880e48_WP, &
                -0.141702981263929354814860763599550983301777836687e48_WP / 0.231495944585472939857781718768469931057948417720e48_WP, &
                -0.88930074506280496946385186439556417902700092547e47_WP / 0.555590267005135055658676125044327834539076202528e48_WP, &
                0.8316e4_WP / 0.10333e5_WP, &
                -0.6562e4_WP / 0.101953e6_WP /)

           upwindRight(i)%rhsBoundary1(1:6,5) = &
                (/ -0.2056582668234943118496241e25_WP / 0.94178732782896630191419096e26_WP, &
                -0.2847973112561150368136731e25_WP / 0.529755371903793544826732415e27_WP, &
                0.10755e5_WP / 0.53963e5_WP, &
                -0.3235e4_WP / 0.6007e4_WP, &
                -0.51051360223797094328698381e26_WP / 0.86490672963884660379874680e26_WP, &
                0.41187e5_WP / 0.32434e5_WP /)

           upwindRight(i)%rhsBoundary1(1:6,6) = &
                (/ 0.39873929407274324907444906e26_WP / 0.2459678221779183710967687025e28_WP, &
                -0.43208769313625564285244989e26_WP / 0.2951613866135020453161224430e28_WP, &
                -0.4629e4_WP / 0.207343e6_WP, &
                -0.3184e4_WP / 0.52305e5_WP, &
                -0.3961e4_WP / 0.134610e6_WP, &
                -0.29791922049731222907256317667e29_WP / 0.29516138661350204531612244300e29_WP /)

        case ('SBP DRP 2-5')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 8
           upwindRight(i)%boundaryDepth = 6
           upwindRight(i)%boundaryWidth = 10
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-3:4))

           upwindRight(i)%rhsInterior(-3:4) = (/ 13.0_WP / 525.0_WP, &
                -109.0_WP / 1050.0_WP, &
                -17.0_WP / 175.0_WP, &
                -127.0_WP / 140.0_WP, &
                31.0_WP / 21.0_WP, &
                -167.0_WP / 350.0_WP, &
                47.0_WP / 525.0_WP, &
                -11.0_WP / 2100.0_WP /)

           upwindRight(i)%normBoundary = &
                (/ 0.906210599613069573511257867906670294043316871e45_WP / 0.2849010651464160240006827153805046624772811600e46_WP, &
                0.1314264380885890713019873837069561991169202483e46_WP / 0.949670217154720080002275717935015541590937200e45_WP, &
                0.1801014065190369285712759682017676542405561049e46_WP / 0.2849010651464160240006827153805046624772811600e46_WP, &
                0.3544873730015507083491260332026958478606873099e46_WP / 0.2849010651464160240006827153805046624772811600e46_WP, &
                0.16373485002812522035248089567e29_WP / 0.18079306480429549386044404300e29_WP, &
                0.1331657974094248518055091423e28_WP / 0.1310745007841312569668932475e28_WP /)

           upwindRight(i)%rhsBoundary1(1:6,1) = &
                (/ -0.25041807479373797745150180079771434881201393483e47_WP / 0.1385210082524478146915319412751523335976576354080e49_WP, &
                0.33575067810193006193264690074734087621366726055151e50_WP / 0.46750840285201137458392030180363912589209451950200e50_WP, &
                -0.3792076605224842645592745442883314164325629945541e49_WP / 0.20778151237867172203729791191272850039648645311200e50_WP, &
                -0.5091541262828815498082963563648745549700549757e46_WP / 0.421178741308118355481009280904179392695580648200e48_WP, &
                -0.6222986245399120158716767907e28_WP / 0.216951677765154592632532851600e30_WP, &
                0.8994644988289303196582066593e28_WP / 0.387980522321028520622004012600e30_WP /)

           upwindRight(i)%rhsBoundary1(1:6,2) = &
                (/ -0.844192226267689525316450178585450773685541811e45_WP / 0.1253424139984212165593576958332477514892273200e46_WP, &
                -0.2014745257661048877151542919979509115504239523621e49_WP / 0.62334453713601516611189373573818550118945935933600e50_WP, &
                0.9413289869368022003302464265727179181229704347487e49_WP / 0.13357382938628896416683437194389689311202700557200e50_WP, &
                -0.279406459093804970571346265420012025077411e42_WP / 0.9678950887908974251687852960096351466349600e43_WP, &
                0.2297181358702433625294300242e28_WP / 0.40678439580966486118599909675e29_WP, &
                -0.146875746940652049694861873e27_WP / 0.5542578890300407437457200180e28_WP /)

           upwindRight(i)%rhsBoundary1(1:6,3) = &
                (/ 0.1586236945888525038433939537009336304249401707393e49_WP / 0.6926050412622390734576597063757616679882881770400e49_WP, &
                -0.34908961647816784350959167524764879510543037765649e50_WP / 0.46750840285201137458392030180363912589209451950200e50_WP, &
                -0.366290752645263039490807984567188757445239510981e48_WP / 0.6926050412622390734576597063757616679882881770400e49_WP, &
                0.18524e5_WP / 0.23949e5_WP, &
                -0.7389e4_WP / 0.40157e5_WP, &
                -0.3009e4_WP / 0.220076e6_WP /)

           upwindRight(i)%rhsBoundary1(1:6,4) = &
                (/ -0.1056122861806108514511401272102581256909299815773e49_WP / 0.31167226856800758305594686786909275059472967966800e50_WP, &
                0.1381533772781156753248785269380866777134377268311e49_WP / 0.14384873933908042294889855440111973104372139061600e50_WP, &
                -0.62111250282499028018134580845475927427441764682551e50_WP / 0.93501680570402274916784060360727825178418903900400e50_WP, &
                -0.987988638532518851912453211991365609048454507e45_WP / 0.5348301476928487053727101979735611335816897120e46_WP, &
                0.19662e5_WP / 0.20015e5_WP, &
                -0.54609e5_WP / 0.195175e6_WP /)

           upwindRight(i)%rhsBoundary1(1:6,5) = &
                (/ -0.3644478718032016407646992871e28_WP / 0.216951677765154592632532851600e30_WP, &
                -0.9640303890867238936682396837e28_WP / 0.325427516647731888948799277400e30_WP, &
                0.21820e5_WP / 0.102449e6_WP, &
                -0.20872e5_WP / 0.41245e5_WP, &
                -0.419036161755782772084170213183e30_WP / 0.650855033295463777897598554800e30_WP, &
                0.42397e5_WP / 0.30808e5_WP /)

           upwindRight(i)%rhsBoundary1(1:6,6) = &
                (/ 0.1120165022724575268042491e25_WP / 0.84564194054278230301221450e26_WP, &
                -0.14585214695327843265197651e26_WP / 0.2621490015682625139337864950e28_WP, &
                -0.11247e5_WP / 0.623410e6_WP, &
                -0.4961e4_WP / 0.74597e5_WP, &
                -0.15296e5_WP / 0.148107e6_WP, &
                -0.35077251164449944159758501e26_WP / 0.38836889121224076138338740e26_WP /)

        case ('SBP DRP 3-6')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 10
           upwindRight(i)%boundaryDepth = 8
           upwindRight(i)%boundaryWidth = 13
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-4:5))

           upwindRight(i)%rhsInterior(-4:5) = &
                (/ -1.0_WP / 168.0_WP, &
                149.0_WP / 3150.0_WP, &
                -199.0_WP / 1575.0_WP, &
                -8.0_WP / 75.0_WP, &
                -8.0_WP / 9.0_WP, &
                67.0_WP / 45.0_WP, &
                -37.0_WP / 75.0_WP, &
                124.0_WP / 1575.0_WP, &
                139.0_WP / 12600.0_WP, &
                -1.0_WP / 210.0_WP /)

           upwindRight(i)%normBoundary = &
                (/ 0.9332474760507988611613909984697980130625228656933410369e55_WP / 0.31697237956293195112204875285069437638526538665203032000e56_WP, &
                0.461358579388732672194183754182312353634531103631318081e54_WP / 0.301878456726601858211475002714947025128824177763838400e54_WP, &
                0.113698985717781971959457551849125386390652850420495609e54_WP / 0.452817685089902787317212504072420537693236266645757600e54_WP, &
                0.818522599690756243349430192726111796225380289930283513e54_WP / 0.452817685089902787317212504072420537693236266645757600e54_WP, &
                0.365089939389748562274173538583206470381976103350108947e54_WP / 0.905635370179805574634425008144841075386472533291515200e54_WP, &
                0.831225138134787260258030069078690949041125229074787969e54_WP / 0.646882407271289696167446434389172196704623238065368000e54_WP, &
                0.160865935583274169227438367745404037e36_WP / 0.174730090913054890018556018162864970e36_WP, &
                0.46962085995528305694995629e26_WP / 0.46505724599538673682179050e26_WP /)

           upwindRight(i)%rhsBoundary1(1:8,1) = &
                (/ -0.180775940537992658045526482429536261001256306443696241791373542184316062920711e78_WP / 0.20018921269138658207947234508273516066105750435795494156871262962114773057440000e80_WP, &
                0.96713189095217748845257723954559643286493509970860915494069923414784294527e74_WP / 0.136183137885296994611885949035874258953100343100649620114770496340916823520e75_WP, &
                -0.101352954845097745898821172233578958699359566208599338773962876583566129e72_WP / 0.982089936672814864989562132470247059757935166591223221981518002458534785e72_WP, &
                -0.8254031724170300416127954338199887869685823614030483112183e58_WP / 0.55301265426974557706689211484852503006861865536646438415200e59_WP, &
                -0.90327572941632962609565963221732791904875901597488414854074500229e65_WP / 0.10155430091210340998241044224596294880764593679933668288077757318400e68_WP, &
                0.39159451033029416705758645043110485680470361280292102747960469231596187e71_WP / 0.504960514831672562930935500311001711663858135372314034281635613422230000e72_WP, &
                -0.1597776579136356343460698268105575e34_WP / 0.209676109095665868022267221795437964e36_WP, &
                -0.4197892495607803036908028427072026436452829482058405104525727e61_WP / 0.431989325088324009224011719758603591437500867681347892518476200e63_WP /)

           upwindRight(i)%rhsBoundary1(1:8,2) = &
                (/ -0.5389031545800344393765837033581973558132191660503280048753102124819207346419e76_WP / 0.7944016376642324685693347027092665105597520014204561173361612286553481372000e76_WP, &
                -0.2504257501903645414594428064616987136329388735336980697085512957205781619077e76_WP / 0.81709882731178196767131569421524555371860205860389772068862297804550094112000e77_WP, &
                0.5221743803551880059495482465697247098785243319061663042232094049679136157e73_WP / 0.13618313788529699461188594903587425895310034310064962011477049634091682352e74_WP, &
                0.204177116807698308109736293822872169160274596227394051221e57_WP / 0.534311743255792828083953734153164280259534932721221627200e57_WP, &
                0.1282344035671760522407967033031686584307603380741241404834403650901e67_WP / 0.8886001329809048373460913696521758020669019469941959752068037653600e67_WP, &
                -0.220116276091032799580704465544553390002870626829362271639628264818025223e72_WP / 0.846410005813089248341377600521298107169895541195497809843503504402976000e72_WP, &
                0.14263609738251635139141484822112252e35_WP / 0.312018019487598017890278603862258875e36_WP, &
                0.416311564576878185060284177554555771229294243616194211326933e60_WP / 0.30856380363451714944572265697043113674107204834381992322748300e62_WP /)

           upwindRight(i)%rhsBoundary1(1:8,3) = &
                (/ 0.1734483714914131982638212569774527237762864059789300854717888971713604215169e76_WP / 0.15888032753284649371386694054185330211195040028409122346723224573106962744000e77_WP, &
                -0.21814496814965829048420466963809752220988809335008179845982022759715892247e74_WP / 0.56742974118873747754952478764947607897125142958604008381154373475382009800e74_WP, &
                -0.46899447741653031874756366546060913088167964891146847363109386493631626051e74_WP / 0.4085494136558909838356578471076227768593010293019488603443114890227504705600e76_WP, &
                0.524198819865168380389604893701341357773650248502131495731e57_WP / 0.877797863920231074709352563251627031854950246613435530400e57_WP, &
                -0.110433657234063890651203615263301163549417271507741041770308528623e66_WP / 0.296200044326968279115363789884058600688967315664731991735601255120e66_WP, &
                -0.949866069749029570714748402376965151520827339660415625500277540238959e69_WP / 0.170909712712258405915085861643723656255459676587552442372245899927524000e72_WP, &
                0.940451378570533970096091476278074353e36_WP / 0.10483805454783293401113361089771898200e38_WP, &
                -0.894146945530470316900987178091374331801872823027903069498311e60_WP / 0.41141840484602286592763020929390818232142939779175989763664400e62_WP /)

           upwindRight(i)%rhsBoundary1(1:8,4) = &
                (/ 0.5214678338745675810795555491310271345896269039814316578457915408911593663053e76_WP / 0.35748073694890461085620061621916992975188840063920525280127255289490666174000e77_WP, &
                -0.4339100416857417971551659877432649666620853113775736194721182764169734801e73_WP / 0.11739925679766982294128099054816746461474167508676691389204353132837657200e74_WP, &
                -0.1219750844773001309344516354781628097986388790164928195024119194257418211531e76_WP / 0.2042747068279454919178289235538113884296505146509744301721557445113752352800e76_WP, &
                -0.224893090002137414756567713661079266431060316909514028717e57_WP / 0.10054775532177192310307129360882273273974884643026625166400e59_WP, &
                0.4525628116417329186104717270483447143619279124211502949994529099377e67_WP / 0.8886001329809048373460913696521758020669019469941959752068037653600e67_WP, &
                0.3271516476790879957448755516432433094390179267678901374712403119431581597e73_WP / 0.4443652530518718553792232402736815062641951591276363501678393398115624000e73_WP, &
                -0.9171e4_WP / 0.16946e5_WP, &
                0.4432192463292494848012456350107723030385464750031624230777293e61_WP / 0.30856380363451714944572265697043113674107204834381992322748300e62_WP /)

           upwindRight(i)%rhsBoundary1(1:8,5) = &
                (/ -0.11201237309749528205897173874662178391589788187071829274100311e62_WP / 0.273111851205669033219564828218534260630883962670476578021139200e63_WP, &
                -0.15625977221154819697766923854048395267458649495808504864449e59_WP / 0.219464880433126901694293165532750745149817470003061535909844e60_WP, &
                0.28231093723191470411864587057203096030108769000203523326816173e62_WP / 0.68277962801417258304891207054633565157720990667619144505284800e62_WP, &
                -0.7835e4_WP / 0.16744e5_WP, &
                -0.16097758434557993144330470080539549010981374789748275377787267e62_WP / 0.106869854819609621694612324085513406333824159305838660964793600e63_WP, &
                0.7298e4_WP / 0.147603e6_WP, &
                0.15202e5_WP / 0.32343e5_WP, &
                -0.24913e5_WP / 0.119522e6_WP /)

           upwindRight(i)%rhsBoundary1(1:8,6) = &
                (/ -0.542281927897990944194381154029631202734854652778936771323727275623e66_WP / 0.8362620848790935734128093607584035597443351405979334020007108780000e67_WP, &
                0.1825697755524986590952938367139708686350713461381264060081038344517e67_WP / 0.8028116014839298304762969863280674173545617349740160659206824428800e67_WP, &
                0.24069044410446302657114849457614991389564550594917237667684548363e65_WP / 0.9031630516694210592858341096190758445238819518457680741607677482400e67_WP, &
                -0.28687e5_WP / 0.40709e5_WP, &
                0.9275e4_WP / 0.57842e5_WP, &
                -0.55788436217580687600435149585922190876767848325211433022737920572287e68_WP / 0.164211463939894738051969838112559244458887627608321468029230499680000e69_WP, &
                0.25546e5_WP / 0.28549e5_WP, &
                -0.16211e5_WP / 0.62117e5_WP /)

           upwindRight(i)%rhsBoundary1(1:8,7) = &
                (/ 0.3763886388709655470705217623826631133e37_WP / 0.73386638183483053807793527628403287400e38_WP, &
                -0.6961648168952882024293793425163113747e37_WP / 0.69892036365221956007422407265145988000e38_WP, &
                -0.16069e5_WP / 0.126517e6_WP, &
                0.11079e5_WP / 0.23257e5_WP, &
                -0.9577e4_WP / 0.21365e5_WP, &
                -0.7397e4_WP / 0.153824e6_WP, &
                -0.173348352280286186280998547303748683e36_WP / 0.208012012991732011926852402574839250e36_WP, &
                0.42182341844716490888368602504391516733e38_WP / 0.29354655273393221523117411051361314960e38_WP /)

           upwindRight(i)%rhsBoundary1(1:8,8) = &
                (/ -0.28354285528443696161146e23_WP / 0.2170267147978471438501689e25_WP, &
                0.357109270498310227962559e24_WP / 0.19931024828373717292362450e26_WP, &
                0.5767e4_WP / 0.147001e6_WP, &
                -0.3621e4_WP / 0.32492e5_WP, &
                0.8987e4_WP / 0.52066e5_WP, &
                -0.7531e4_WP / 0.29989e5_WP, &
                -0.6463e4_WP / 0.195323e6_WP, &
                -0.3522043255326282016490092819e28_WP / 0.3906480866361248589303040200e28_WP /)


        case ('SBP DRP 3-7')
           upwindRight(i)%implicit = .false.
           upwindRight(i)%symmetryType = ASYMMETRIC
           upwindRight(i)%interiorWidth = 10
           upwindRight(i)%boundaryDepth = 8
           upwindRight(i)%boundaryWidth = 13
           call allocate_operator(upwindRight(i))

           allocate(upwindRight(i)%rhsInterior(-4:5))

           upwindRight(i)%rhsInterior(-4:5) = &
                (/ -0.43e2_WP / 0.7056e4_WP, &
                +0.4859e4_WP / 0.117600e6_WP, &
                -0.107e3_WP / 0.1225e4_WP, &
                -0.841e3_WP / 0.4200e4_WP, &
                -0.1111e4_WP / 0.1400e4_WP, &
                +0.119e3_WP / 0.80e2_WP, &
                -0.617e3_WP / 0.1050e4_WP, &
                +0.5113e4_WP / 0.29400e5_WP, &
                -0.587e3_WP / 0.19600e5_WP, &
                +0.737e3_WP / 0.352800e6_WP /)

           upwindRight(i)%normBoundary = &
                (/ 0.32612589119845772001865474176308615053989456881054591315219e59_WP / 0.110588548303924429111746555990950439992763334355264722747200e60_WP, &
                0.4820609787469451380456598357197370327847417072166890177949e58_WP / 0.3159672808683555117478473028312869714078952410150420649920e58_WP, &
                0.2034603860470626768046523146406119413070350909537066994339e58_WP / 0.7899182021708887793696182570782174285197381025376051624800e58_WP, &
                0.2840521456436738681256835049188074812009893965236437952833e58_WP / 0.1579836404341777558739236514156434857039476205075210324960e58_WP, &
                0.186316294181664167494102878452926008162536969135702585471e57_WP / 0.451381829811936445354067575473267102011278915735774378560e57_WP, &
                0.20198039880886397420119082053804443554406154945206382773553e59_WP / 0.15798364043417775587392365141564348570394762050752103249600e59_WP, &
                0.106739398921576344786855185165609254577e39_WP / 0.115608875908367172736668570302045884680e39_WP, &
                0.131059199350960203297730429021e30_WP / 0.129846691063143520265326690800e30_WP /)

           upwindRight(i)%rhsBoundary1(1:8,1) = &
                (/ -0.52261358142498941499949464975602780739602300028629113102569797578605941788630379e80_WP / 0.5547082942311126114222813721552630983159480953277788863314811583369275683571488000e82_WP, &
                0.91157186794483056826076122312771547107546769500083270111292901645516074349699e77_WP / 0.128642925378272869068247071464578640611305216912750205549972439317469287652400e78_WP, &
                -0.27225729524148140056459092345613438474804978386469892068744099001474388174961e77_WP / 0.283014435832200311950143557222073009344871477208050452209939366498432432835280e78_WP, &
                -0.2149454530683064334530890818723014283823218000361118312814403541e64_WP / 0.13492638626536677151961652887011942433156500674254780961420303840e65_WP, &
                -0.57533797620871627209083635859617442585536009288685516740869951931723e68_WP / 0.25951610770911040127027989948392250077888649728363405931750712193964800e71_WP, &
                0.8989899972139712619891237695552685461986690664313677002112081560260422701e73_WP / 0.117700032368852384219078724921221110626095418337000196382650888111003532000e75_WP, &
                -0.29405475871990759059321419563628367437e38_WP / 0.3468266277251015182100057109061376540400e40_WP, &
                -0.1192379384157571322250302840179951193878136027913692721352515821e64_WP / 0.127330641940743179871210930138577044467701519437371611378508425950e66_WP /)

           upwindRight(i)%rhsBoundary1(1:8,2) = &
                (/ -0.4456849697544281990037602771332172357172000055509857801309595741741079452369319e79_WP / 0.6603670169418007278836683001848370218047001134854510551565251884963423432823200e79_WP, &
                -0.3625911323575924556320286941909191809627511202149651257923960578055629030910609e79_WP / 0.113205774332880124780057422888829203737948590883220180883975746599372973134112000e81_WP, &
                0.346619171489406086207467675483522630741630418021465025683339713130980226216651e78_WP / 0.943381452774001039833811857406910031149571590693501507366464554994774776117600e78_WP, &
                0.336253495753891266248708478214387692179940175030119462432975373e63_WP / 0.817735674335556191027978962849208632312515192379077634025472960e63_WP, &
                0.2876443529250687529805457958743475865047231221979852652520463948768483e70_WP / 0.22707659424547160111149491204843218818152568512317980190281873169719200e71_WP, &
                -0.16268936218819694290268092004100895487335023888368224772310739127881009343e74_WP / 0.62773350596721271583508653291317925667250889779733438070747140325868550400e74_WP, &
                0.78037278977406893193465057435409608473e38_WP / 0.1651555370119531039095265290029226924000e40_WP, &
                0.3945767210261829923162869611560011163274275358760457560096215697e64_WP / 0.291041467293127268277053554602461815926174901571135111722304973600e66_WP /)

           upwindRight(i)%rhsBoundary1(1:8,3) = &
                (/ 0.40996921736863733989842971400575398623394478761135606145482187727338020075343e77_WP / 0.400222434510182259323435333445355770790727341506333972822136477876571117140800e78_WP, &
                -0.350576161720426820159228612672170342491193858261290762514072270415103108481939e78_WP / 0.943381452774001039833811857406910031149571590693501507366464554994774776117600e78_WP, &
                -0.48899970355913152469421244969734422337777596958271859949009534652803713277551e77_WP / 0.5660288716644006239002871144441460186897429544161009044198787329968648656705600e79_WP, &
                0.94425415917615061578696656160409311801070656769277245907599411e62_WP / 0.166575785512798483357551270210023980656253094743886184708892640e63_WP, &
                -0.8640808500511739680553095116332756794443989460372707372253893493343e67_WP / 0.24575388987605151635443172299613873179818797091253225314157871395800e68_WP, &
                0.23944004261466051019846442017994377972431693571528776347724801662732261e71_WP / 0.5885001618442619210953936246061055531304770916850009819132544405550176600e73_WP, &
                0.2986748932789792703157318508467071101e37_WP / 0.38536291969455724245556190100681961560e38_WP, &
                -0.307381857998921420452838571171768433875552195605971097066501747e63_WP / 0.16168970405173737126502975255692323107009716753951950651239165200e65_WP /)

           upwindRight(i)%rhsBoundary1(1:8,4) = &
                (/ 0.427952120942894743815617286397805998548874928084254607389010858562531278133083e78_WP / 0.2830144358322003119501435572220730093448714772080504522099393664984324328352800e79_WP, &
                -0.1111912812783381464882252678283331672976876100870930629098604124508767100663793e79_WP / 0.2830144358322003119501435572220730093448714772080504522099393664984324328352800e79_WP, &
                -0.1599696759745610934745700490934071301517063017138530390239775891279468609161069e79_WP / 0.2830144358322003119501435572220730093448714772080504522099393664984324328352800e79_WP, &
                -0.395815202366413293129524140713233825835788407966186586263205209e63_WP / 0.26985277253073354303923305774023884866313001348509561922840607680e65_WP, &
                0.121538486776787523189893798810514400404088463266354249275628155260933e69_WP / 0.249534718951067693529115287965310096902775478157340441651449155711200e69_WP, &
                0.703111769387557124919890815741790867136274052156521275044037336226316917e72_WP / 0.1001702403139169227396414680180605196817833347548937841554475643497902400e73_WP, &
                -0.16966e5_WP / 0.35549e5_WP, &
                0.22628840955829608944685243738416306699578101193273282288777477e62_WP / 0.209986628638619962681856821502497702688437879921453904561547600e63_WP /)

           upwindRight(i)%rhsBoundary1(1:8,5) = &
                (/ -0.1019249422313673686509620616193453584613953685274622269162522248483e67_WP / 0.24096203129908115252579377853660399329515923610365279416667327942400e68_WP, &
                -0.3565831596484310882429109554637417306934245316265552009783606045551e67_WP / 0.63252533216008802538020866865858548239979299477208858468751735848800e68_WP, &
                0.163349564987719227850726389247121861527043470090458694571157109709e66_WP / 0.425942984619587895878928396403087866935887538567063019991594180800e66_WP, &
                -0.27815e5_WP / 0.59562e5_WP, &
                -0.764116593479884538611420087225553116841412769956341580793574982559e66_WP / 0.5560662260748026596749087196998553691426751602391987557692460294400e67_WP, &
                0.9311e4_WP / 0.88804e5_WP, &
                0.16345e5_WP / 0.46859e5_WP, &
                -0.14519e5_WP / 0.136686e6_WP /)

           upwindRight(i)%rhsBoundary1(1:8,6) = &
                (/ -0.39145054214004491098464056761028939716099959050728114120233505288171e68_WP / 0.589062716738747424886159907317593855262252542863435562875800830348000e69_WP, &
                0.660109568193310770426644319759422203113432643741474485262183853770301e69_WP / 0.2827501040345987639453567555124450505258812205744490701803843985670400e70_WP, &
                -0.506859236681854929589687434721939753221451455062634409339924298299e66_WP / 0.106031289012974536479508783317166893947205457715418401317644149462640e69_WP, &
                -0.25629e5_WP / 0.36134e5_WP, &
                0.12919e5_WP / 0.67851e5_WP, &
                -0.350551717829392183514568406785916141815918474476420360806039476098269e69_WP / 0.902393949046591799825606666529079948486854959280156606958673612448000e69_WP, &
                0.39152e5_WP / 0.40241e5_WP, &
                -0.18054e5_WP / 0.48289e5_WP /)

           upwindRight(i)%rhsBoundary1(1:8,7) = &
                (/ 0.58552039474615468102759060262347727833e38_WP / 0.1103539270034413921577290898337710717400e40_WP, &
                -0.1917722481913915714354280142018472252441e40_WP / 0.17341331386255075910500285545306882702000e41_WP, &
                -0.7238e4_WP / 0.65577e5_WP, &
                0.12791e5_WP / 0.26114e5_WP, &
                -0.30171e5_WP / 0.60752e5_WP, &
                -0.9709e4_WP / 0.257928e6_WP, &
                -0.25903743628751519829385853082982169837e38_WP / 0.32596487568148638929511814934787373500e38_WP, &
                0.70337969501417118523361802238856954351173e41_WP / 0.48555727881514212549400799526859271565600e41_WP /)

           upwindRight(i)%rhsBoundary1(1:8,8) = &
                (/ -0.24772283621689535123847576463e29_WP / 0.1817853674884009283714573671200e31_WP, &
                0.2384271497669209460992645529e28_WP / 0.111297163768408731655994306400e30_WP, &
                0.18854e5_WP / 0.549607e6_WP, &
                -0.10883e5_WP / 0.92516e5_WP, &
                0.8014e4_WP / 0.42219e5_WP, &
                -0.4529e4_WP / 0.19097e5_WP, &
                -0.10524e5_WP / 0.92341e5_WP, &
                -0.4412432084832700153493001602867e31_WP / 0.5453561024652027851143721013600e31_WP /)

        case default

           call die('first_derivative_setup: unknown upwind scheme: ' // trim(stencilScheme))

        end select

        ! Fill the right-boundary coefficients of Q+
        call operator_fill_rhsBoundary(upwindRight(i), useUpwinding)

        ! Compute transpose(Q+)
        call operator_transpose(upwindRight(i), upwindLeft(i), useUpwinding)

        ! Compute (Q+ + B/2)
        upwindRight(i)%rhsBoundary1(1,1) = upwindRight(i)%rhsBoundary1(1,1) - 1.0_WP / 2.0_WP
        upwindRight(i)%rhsBoundary2(upwindRight(i)%boundaryWidth,1) =                        &
             upwindRight(i)%rhsBoundary2(upwindRight(i)%boundaryWidth,1) + 1.0_WP / 2.0_WP

        ! Compute -(Q- + B/2) = transpose(Q+) - B/2
        upwindLeft(i)%rhsBoundary1(1,1) = upwindLeft(i)%rhsBoundary1(1,1) + 1.0_WP / 2.0_WP
        upwindLeft(i)%rhsBoundary2(upwindLeft(i)%boundaryWidth,1) =                          &
             upwindLeft(i)%rhsBoundary2(upwindLeft(i)%boundaryWidth,1) - 1.0_WP / 2.0_WP

        ! Compute D+ = H^(-1)(Q+ + B/2)
        call operator_multiply_inverse_diagonal_norm(upwindRight(i), +1.0_WP)

        ! Compute D- = H^(-1)(Q- + B/2)
        call operator_multiply_inverse_diagonal_norm(upwindLeft(i), -1.0_WP)

        ! Update
        call operator_update(upwindRight(i), i)
        call operator_update(upwindLeft(i), i)

        ! Overwrite the first derivative operator using upwind coefficients
        call deallocate_operator(firstDerivative(i))
        firstDerivative(i)%implicit = .false.
        firstDerivative(i)%symmetryType = SKEW_SYMMETRIC
        firstDerivative(i)%interiorWidth = (upwindRight(i)%boundaryDepth/2+1) * 2 + 1
        firstDerivative(i)%boundaryDepth = upwindRight(i)%boundaryDepth
        firstDerivative(i)%boundaryWidth = upwindRight(i)%boundaryWidth
        call allocate_operator(firstDerivative(i))
        firstDerivative(i)%rhsInterior = 0.0_WP

        lb = lbound(upwindRight(i)%rhsInterior,1)
        ub = ubound(upwindRight(i)%rhsInterior,1)
        firstDerivative(i)%rhsInterior(lb:ub) = 0.5_WP * upwindRight(i)%rhsInterior(lb:ub)

        lb = lbound(upwindLeft(i)%rhsInterior,1)
        ub = ubound(upwindleft(i)%rhsInterior,1)
        firstDerivative(i)%rhsInterior(lb:ub) = firstDerivative(i)%rhsInterior(lb:ub) +      &
             0.5_WP * upwindLeft(i)%rhsInterior(lb:ub)

        firstDerivative(i)%normBoundary = upwindRight(i)%normBoundary
        firstDerivative(i)%rhsBoundary1 =                                                    &
             0.5_WP * (upwindRight(i)%rhsBoundary1 + upwindLeft(i)%rhsBoundary1)
        firstDerivative(i)%rhsBoundary2 =                                                    &
             0.5_WP * (upwindRight(i)%rhsBoundary2 + upwindLeft(i)%rhsBoundary2)
        
        call operator_update(firstDerivative(i), i)

     end do

  end if

  return
end subroutine first_derivative_setup


! ====================================== !
! Cleanup the first derivative operators !
! ====================================== !
subroutine first_derivative_cleanup

  ! Internal modules
  use first_derivative

  implicit none

  ! Local variables
  integer :: i

  if (allocated(firstDerivative)) then
     do i = 1, size(firstDerivative)
        call deallocate_operator(firstDerivative(i))
     end do
     deallocate(firstDerivative)
  end if

  if (allocated(adjointFirstDerivative)) then
     do i = 1, size(adjointFirstDerivative)
        call deallocate_operator(adjointFirstDerivative(i))
     end do
     deallocate(adjointFirstDerivative)
  end if

  if (allocated(upwindLeft)) then
     do i = 1, size(upwindLeft)
        call deallocate_operator(upwindLeft(i))
     end do
     deallocate(upwindLeft)
  end if

  if (allocated(upwindRight)) then
     do i = 1, size(upwindRight)
        call deallocate_operator(upwindRight(i))
     end do
     deallocate(upwindRight)
  end if

  return
end subroutine first_derivative_cleanup
