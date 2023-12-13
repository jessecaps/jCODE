module fourth_derivative

  ! External modules
  use precision
  use operator

  implicit none

  type(t_StencilOperator), allocatable, dimension(:) :: fourthDerivative

contains

  subroutine fourth_derivative_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('fourth_derivative_apply: size of `x` inconsistent with local grid size')

    ! Explicit step
    select case (direction)
    case (1)
       call operator_apply_1(fourthDerivative(1), x)
    case (2)
       call operator_apply_2(fourthDerivative(2), x)
    case (3)
       call operator_apply_3(fourthDerivative(3), x)
    end select

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(fourthDerivative(1), x)
    case (2)
       call operator_apply_implicit_2(fourthDerivative(2), x)
    case (3)
       call operator_apply_implicit_3(fourthDerivative(3), x)
    end select

    return
  end subroutine fourth_derivative_apply

end module fourth_derivative


! ==================================== !
! Setup the fourth derivative operator !
! ==================================== !
subroutine fourth_derivative_setup

  ! Internal modules
  use fourth_derivative

  ! External modules
  use string
  use parser
  use simulation_flags
  use operator, only : discretizationType

  implicit none

  ! Local variables
  integer :: i
  real(WP) :: alpha, a, b
  character(len = str_medium) :: stencilScheme

  ! Decide if fourth derivative is used
  if (.not. useShockCapturing) return

  allocate(fourthDerivative(nDimensions))

  do i = 1, nDimensions

     if (globalGridSize(i) .gt. 1) then
        call parser_read('fourth derivative scheme', stencilScheme, 'Pade 3-6')
     else
        stencilScheme = 'null matrix'
     end if

     select case (trim(stencilScheme))

     case ('explicit 2nd-order')

        fourthDerivative(i)%implicit = .false.
        fourthDerivative(i)%symmetryType = SYMMETRIC
        fourthDerivative(i)%interiorWidth = 5
        fourthDerivative(i)%boundaryWidth = 5
        fourthDerivative(i)%boundaryDepth = 3
        call allocate_operator(fourthDerivative(i))

        fourthDerivative(i)%rhsInterior(0:2) = (/ 6.0_WP, -4.0_WP, 1.0_WP /)
        fourthDerivative(i)%rhsInterior(-1:-2:-1) = fourthDerivative(i)%rhsInterior(1:2)

        fourthDerivative(i)%normBoundary = 1.0_WP

        fourthDerivative(i)%rhsBoundary1(1:5,1) = (/ 1.0_WP, -4.0_WP, 6.0_WP, -4.0_WP, 1.0_WP /)
        fourthDerivative(i)%rhsBoundary1(1:5,2) = (/ 1.0_WP, -4.0_WP, 6.0_WP, -4.0_WP, 1.0_WP /)
        fourthDerivative(i)%rhsBoundary1(1:5,3) = (/ 1.0_WP, -4.0_WP, 6.0_WP, -4.0_WP, 1.0_WP /)

     case ('Pade 3-6')

        fourthDerivative(i)%implicit = .true.
        fourthDerivative(i)%symmetryType = SYMMETRIC
        fourthDerivative(i)%interiorWidth = 7
        fourthDerivative(i)%boundaryWidth = 6
        fourthDerivative(i)%boundaryDepth = 3
        fourthDerivative(i)%nDiagonals = 3
        call allocate_operator(fourthDerivative(i))

        ! Right-hand side interior
        alpha = 7.0_WP / 26.0_WP
        a = 2.0_WP * (1.0_WP - alpha)
        b = 4.0_WP * alpha - 1.0_WP
        fourthDerivative(i)%rhsInterior(0:3) = (/                                            &
             16.0_WP * b / 6.0_WP + 6.0_WP * a,                                              &
             -9.0_WP * b / 6.0_WP - 4.0_WP * a,                                              &
             a, b / 6.0_WP /)
        fourthDerivative(i)%rhsInterior(-1:-3:-1) = fourthDerivative(i)%rhsInterior(1:3)

        ! Left-hand side coefficients
        fourthDerivative(i)%lhsInterior(-1:1) = (/ alpha, 1.0_WP, alpha /)

        fourthDerivative(i)%normBoundary = 1.0_WP

        ! (i = 1) on boundary (one-sided explicit)
        fourthDerivative(i)%lhsBoundary1(-1:1,1) = (/ 0.0_WP, 1.0_WP, 0.0_WP /)
        fourthDerivative(i)%rhsBoundary1(1:6,1) = (/ 3.0_WP, -14.0_WP, 26.0_WP, -24.0_WP,    &
             11.0_WP, -2.0_WP /)

        ! (i = 2) on boundary (one-sided explicit)
        fourthDerivative(i)%lhsBoundary1(-1:1,2) = (/ 0.0_WP, 1.0_WP, 0.0_WP /)
        fourthDerivative(i)%rhsBoundary1(1:6,2) = (/ 2.0_WP, -9.0_WP, 16.0_WP, -14.0_WP,     &
             6.0_WP, -1.0_WP /)
        
        ! (i = 3) on boundary (one-sided explicit)
        fourthDerivative(i)%lhsBoundary1(-1:1,3) = (/ 0.0_WP, 1.0_WP, 0.0_WP /)
        fourthDerivative(i)%rhsBoundary1(1:6,3) = (/ 1.0_WP, -4.0_WP, 6.0_WP, -4.0_WP,       &
             1.0_WP, 0.0_WP /)
        
     case ('DRP 13 point')

        fourthDerivative(i)%implicit = .false.
        fourthDerivative(i)%symmetryType = SYMMETRIC
        fourthDerivative(i)%interiorWidth = 13
        fourthDerivative(i)%boundaryWidth = 11
        fourthDerivative(i)%boundaryDepth = 6
        call allocate_operator(fourthDerivative(i))

        fourthDerivative(i)%rhsInterior(0:6) = (/                                            &
             15.53434242252894_WP, -11.70502896612137_WP, 5.06622328028830_WP,               &
             -1.41994874382408_WP, 0.34825585767883_WP, -0.06252689371137_WP,                &
             0.00585425442522_WP /)
        fourthDerivative(i)%rhsInterior(-1:-6:-1) = fourthDerivative(i)%rhsInterior(1:6)

        fourthDerivative(i)%normBoundary = 1.0_WP

        ! (i = 1) on boundary fourth derivative: 2nd order one-sided (Kawai & Lele, JCP 2008)
        fourthDerivative(i)%rhsBoundary1(1:6,1) = (/                                         &
             3.0_WP, -14.0_WP, 26.0_WP, -24.0_WP, 11.0_WP, -2.0_WP /)

        ! (i = 2) on boundary fourth derivative: 2nd order one-sided (Kawai & Lele, JCP 2008)
        fourthDerivative(i)%rhsBoundary1(1:6,2) = (/                                         &
             2.0_WP, -9.0_WP, 16.0_WP, -14.0_WP, 6.0_WP, -1.0_WP /)

        ! (i = 3) on boundary fourth derivative: 2nd order non-optimized centered
        fourthDerivative(i)%rhsBoundary1(1:5,3) = (/                                         &
             1.0_WP, -4.0_WP, 6.0_WP, -4.0_WP, 1.0_WP /)

        ! (i = 4) on boundary fourth derivative: 4th order non-optimized centered
        fourthDerivative(i)%rhsBoundary1(1:7,4) = (/                                         &
             -1.0_WP/6.0_WP, 2.0_WP, -6.5_WP, 28.0_WP/3.0_WP, -6.5_WP, 2.0_WP,               &
             -1.0_WP/6.0_WP /)

        ! (i = 5) on boundary fourth derivative: 4th order 9pt DRP
        fourthDerivative(i)%rhsBoundary1(1:9,5) = (/                                         &
             0.04427486965170_WP, -0.52086562388028_WP, 3.23969635024767_WP,                 &
             -8.97939270049536_WP, 12.43257420895253_WP, -8.97939270049536_WP,               &
             3.23969635024767_WP, -0.52086562388028_WP, 0.04427486965170_WP/)

        ! (i = 6) on boundary fourth derivative: 4th order 11pt DRP
        fourthDerivative(i)%rhsBoundary1(1:11,5) = (/                                        &
             -0.01214058678882_WP, 0.14402579928738_WP, -0.89395252335694_WP,                &
             4.09022849383507_WP, -10.31623938400544_WP, 13.97615640205749_WP,               &
             -10.31623938400544_WP, 4.09022849383507_WP, -0.89395252335694_WP,               &
             0.14402579928738_WP, -0.01214058678882_WP /)

     case ('null matrix')

        fourthDerivative(i)%symmetryType = SYMMETRIC
        fourthDerivative(i)%interiorWidth = 0
        fourthDerivative(i)%boundaryWidth = 1
        fourthDerivative(i)%boundaryDepth = 1
        call allocate_operator(fourthDerivative(i))

        fourthDerivative(i)%rhsInterior(0:0) = 0.0_WP
        fourthDerivative(i)%rhsBoundary1 = 0.0_WP

     case default

        call die('fourth_derivative_setup: unknown stencil scheme: ' // trim(stencilScheme))

     end select
     
     ! Fill the right-boundary coefficients
     if (allocated(fourthDerivative(i)%rhsBoundary1) .and.                                   &
          allocated(fourthDerivative(i)%rhsBoundary2)) then
        select case (fourthDerivative(i)%symmetryType)
        case (SYMMETRIC)
           fourthDerivative(i)%rhsBoundary2(1:fourthDerivative(i)%boundaryWidth,:) =         &
                +fourthDerivative(i)%rhsBoundary1(fourthDerivative(i)%boundaryWidth:1:-1,:)
        case (SKEW_SYMMETRIC)
           fourthDerivative(i)%rhsBoundary2(1:fourthDerivative(i)%boundaryWidth,:) =         &
                -fourthDerivative(i)%rhsBoundary1(fourthDerivative(i)%boundaryWidth:1:-1,:)
        end select
     end if
     if (allocated(fourthDerivative(i)%lhsBoundary1) .and.                                   &
          allocated(fourthDerivative(i)%lhsBoundary2)) then
        fourthDerivative(i)%lhsBoundary2(                                                    &
             -fourthDerivative(i)%nDiagonals/2:fourthDerivative(i)%ndiagonals/2,:) =         &
             fourthDerivative(i)%lhsBoundary1(                                               &
             fourthDerivative(i)%nDiagonals/2:-fourthDerivative(i)%nDiagonals/2:-1,:)
     end if

     call operator_update(fourthDerivative(i), i)

  end do

  return
end subroutine fourth_derivative_setup


! ====================================== !
! Cleanup the fourth derivative operator !
! ====================================== !
subroutine fourth_derivative_cleanup

  ! Internal modules
  use fourth_derivative

  implicit none

  ! Local variables
  integer :: i

  if (allocated(fourthDerivative)) then
     do i = 1, size(fourthDerivative)
        call deallocate_operator(fourthDerivative(i))
     end do
     deallocate(fourthDerivative)
  end if

  return
end subroutine fourth_derivative_cleanup
