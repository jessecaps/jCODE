module dissipation

  ! External modules
  use precision
  use operator

  implicit none

  integer :: dissSensorType
  type(t_StencilOperator), allocatable, dimension(:) :: artificialDissipation,               &
       dissipationTranspose, artificialDissipation2
  real(WP) :: dissipationAmount
  real(WP), dimension(:,:), allocatable :: dissipationSource, dissipationSensor
  logical :: useDissipation, compositeDissipation, hybridDissipation

contains

  subroutine dissipation_apply(direction, x, phi)

    ! External modules
    use simulation_flags

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)
    real(WP), intent(in), optional :: phi(:)

    ! Explicit step
    if (.not.present(phi) .or. .not.hybridDissipation) then
       select case (direction)
       case (1)
          call operator_apply_1(artificialDissipation(1), x)
       case (2)
          call operator_apply_2(artificialDissipation(2), x)
       case (3)
          call operator_apply_3(artificialDissipation(3), x)
       end select
    else
       select case (direction)
       case (1)
          call operator_hybrid_apply_1(artificialDissipation(1), artificialDissipation2(1),  &
               phi, x)
       case (2)
          call operator_hybrid_apply_2(artificialDissipation(2), artificialDissipation2(2),  &
               phi, x)
       case (3)
          call operator_hybrid_apply_3(artificialDissipation(3), artificialDissipation2(3),  &
               phi, x)
       end select
    end if

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(artificialDissipation(1), x)
    case (2)
       call operator_apply_implicit_2(artificialDissipation(2), x)
    case (3)
       call operator_apply_implicit_3(artificialDissipation(3), x)
    end select

    return
  end subroutine dissipation_apply


  subroutine dissipation_transpose_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    select case (direction)
    case (1)
       call operator_apply_1(dissipationTranspose(1), x)
    case (2)
       call operator_apply_2(dissipationTranspose(2), x)
    case (3)
       call operator_apply_3(dissipationTranspose(3), x)
    end select

    return
  end subroutine dissipation_transpose_apply

end module dissipation


! ============================== !
! Setup the dissipation operator !
! ============================== !
subroutine dissipation_setup

  ! Internal modules
  use dissipation

  ! External modules
  use string
  use parser
  use simulation_flags
  use solver_options
  use operator, only : discretizationType

  implicit none

  ! Local variables
  integer :: i
  character(len = str_medium) :: stencilScheme, sensorType

  ! Decide if dissipation is used
  call parser_read('add dissipation', useDissipation, .false.)
  hybridDissipation = .false.
  if (.not. useDissipation) return

  ! Read from input
  call parser_read('dissipation amount', dissipationAmount)
  call parser_read('composite dissipation', compositeDissipation, .true.)
  call parser_read('use hybrid dissipation', hybridDissipation, .false.)

  ! Hybrid dissipation only works for composite true
  if (hybridDissipation .and. .not.compositeDissipation)                                     &
       call die('dissipation_setup: hybrid dissipation requires composite dissipation: true')

  ! Sensor type for hybrid dissipation
  if (hybridDissipation) then
     allocate(dissipationSensor(nGridPoints, 1))
     call parser_read('dissipation sensor', sensorType)
     select case (trim(sensorType))
     case ('ducros', 'Ducros', 'DUCROS')
        dissSensorType = DUCROS_SENSOR
     case ('bounded', 'Bounded', 'BOUNDED')
        if (nSpecies .le. 0) then
           call die('dissipation_setup: bounded sensor requires nSpecies > 0')
        end if
        dissSensorType = BOUNDED_SENSOR
     case default
        call die("dissipation_setup: invalid dissipation sensor '" //                        &
             trim(sensorType) // "'!")
     end select
  end if

  ! Monitor dissipation if homogeneous
  if (all(isPeriodic(1:nDimensions))) allocate(dissipationSource(nGridPoints, nUnknowns))

  ! Allocate arrays
  allocate(artificialDissipation(nDimensions))
  if (.not. compositeDissipation) allocate(dissipationTranspose(nDimensions))
  if (hybridDissipation) allocate(artificialDissipation2(nDimensions))

  do i = 1, nDimensions

     if (globalGridSize(i) .gt. 1) then
        call parser_read('artificial dissipation scheme', stencilScheme,                     &
             trim(discretizationType))
        if (compositeDissipation) stencilScheme = trim(stencilScheme) //                     &
             ' composite dissipation'
     else
        stencilScheme = 'null matrix'
     end if

     select case (trim(stencilScheme))

     case ('SBP 1-2')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = ASYMMETRIC
        artificialDissipation(i)%interiorWidth = 2
        artificialDissipation(i)%boundaryWidth = 2
        artificialDissipation(i)%boundaryDepth = 1
        call allocate_operator(artificialDissipation(i))

        allocate(artificialDissipation(i)%rhsInterior(-1:0))
        artificialDissipation(i)%rhsInterior(-1:0) = (/ -1.0_WP, 1.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:2,1) = (/ 0.0_WP, 0.0_WP /)
        !     artificialDissipation(i)%rhsInterior(-1:0)

        artificialDissipation(i)%rhsBoundary2(1:2,1) = (/ -1.0_WP, 1.0_WP /)

        ! Transpose dissipation operator
        dissipationTranspose(i)%implicit = .false.
        dissipationTranspose(i)%symmetryType = ASYMMETRIC
        dissipationTranspose(i)%interiorWidth = 2
        dissipationTranspose(i)%boundaryWidth = 2
        dissipationTranspose(i)%boundaryDepth = 1
        call allocate_operator(dissipationTranspose(i))

        allocate(dissipationTranspose(i)%rhsInterior(0:1))
        dissipationTranspose(i)%rhsInterior(0:1) = (/ 1.0_WP, -1.0_WP /)

        dissipationTranspose(i)%rhsBoundary1(1:2,1) = (/ 0.0_WP, -1.0_WP /)

        dissipationTranspose(i)%rhsBoundary2(1:2,1) = (/ 0.0_WP,  1.0_WP /)

     case ('SBP 2-4')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 3
        artificialDissipation(i)%boundaryWidth = 3
        artificialDissipation(i)%boundaryDepth = 1
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(-1:1) = (/ 1.0_WP, -2.0_WP, 1.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:3,1) = artificialDissipation(i)%rhsInterior

        ! Transpose dissipation operator
        dissipationTranspose(i)%implicit = .false.
        dissipationTranspose(i)%symmetryType = SYMMETRIC
        dissipationTranspose(i)%interiorWidth = 3
        dissipationTranspose(i)%boundaryWidth = 4
        dissipationTranspose(i)%boundaryDepth = 3
        call allocate_operator(dissipationTranspose(i))

        dissipationTranspose(i)%rhsInterior(-1:1) = (/ 1.0_WP, -2.0_WP, 1.0_WP /)

        dissipationTranspose(i)%rhsBoundary1(1,1:3) = dissipationTranspose(i)%rhsInterior(-1:1)
        dissipationTranspose(i)%rhsBoundary1(2,1:3) = dissipationTranspose(i)%rhsInterior(-1:1)
        dissipationTranspose(i)%rhsBoundary1(3,2:3) = dissipationTranspose(i)%rhsInterior(-1:0)
        dissipationTranspose(i)%rhsBoundary1(4,3:3) = dissipationTranspose(i)%rhsInterior(-1:-1)

     case ('SBP 3-6')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = ASYMMETRIC
        artificialDissipation(i)%interiorWidth = 4
        artificialDissipation(i)%boundaryWidth = 4
        artificialDissipation(i)%boundaryDepth = 2
        call allocate_operator(artificialDissipation(i))

        allocate(artificialDissipation(i)%rhsInterior(-2:1))
        artificialDissipation(i)%rhsInterior(-2:1) = (/ -1.0_WP, 3.0_WP, -3.0_WP, 1.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:4,1) =                                       &
             artificialDissipation(i)%rhsInterior(-2:1)
        artificialDissipation(i)%rhsBoundary1(1:4,2) =                                       &
             artificialDissipation(i)%rhsInterior(-2:1)

        artificialDissipation(i)%rhsBoundary2(1:4,:) =                                       &
             artificialDissipation(i)%rhsBoundary1(1:4,:)

        ! Transpose dissipation operator
        dissipationTranspose(i)%implicit = .false.
        dissipationTranspose(i)%symmetryType = ASYMMETRIC
        dissipationTranspose(i)%interiorWidth = 4
        dissipationTranspose(i)%boundaryWidth = 6
        dissipationTranspose(i)%boundaryDepth = 4
        call allocate_operator(dissipationTranspose(i))

        allocate(dissipationTranspose(i)%rhsInterior(-1:2))
        dissipationTranspose(i)%rhsInterior(-1:2) = (/ 1.0_WP, -3.0_WP, 3.0_WP, -1.0_WP /)

        dissipationTranspose(i)%rhsBoundary1(1,1:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:-1:-1)
        dissipationTranspose(i)%rhsBoundary1(2,1:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:-1:-1)
        dissipationTranspose(i)%rhsBoundary1(3,1:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:-1:-1)
        dissipationTranspose(i)%rhsBoundary1(4,2:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:0:-1)
        dissipationTranspose(i)%rhsBoundary1(5,3:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:1:-1)
        dissipationTranspose(i)%rhsBoundary1(6,4:4) =                                        &
             dissipationTranspose(i)%rhsInterior(2:2:-1)

!!$        dissipationTranspose(i)%rhsBoundary2(6,1:4) =                                        &
!!$             dissipationTranspose(i)%rhsInterior(2:-1:-1)
!!$        dissipationTranspose(i)%rhsBoundary2(5,1:4) =                                        &
!!$             dissipationTranspose(i)%rhsInterior(2:-1:-1)
!!$        dissipationTranspose(i)%rhsBoundary2(4,2:4) =                                        &
!!$             dissipationTranspose(i)%rhsInterior(-1:1)
!!$        dissipationTranspose(i)%rhsBoundary2(3,3:4) =                                        &
!!$             dissipationTranspose(i)%rhsInterior(-1:0)
!!$        dissipationTranspose(i)%rhsBoundary2(2,4:4) =                                        &
!!$             dissipationTranspose(i)%rhsInterior(-1:-1)
!!$        dissipationTranspose(i)%rhsBoundary2(1,1:4) = 0.0_WP ! ... the first column is zero

        dissipationTranspose(i)%rhsBoundary2(6,1:4) =                                        &
             dissipationTranspose(i)%rhsInterior(-1:2)
        dissipationTranspose(i)%rhsBoundary2(5,1:4) =                                        &
             dissipationTranspose(i)%rhsInterior(-1:2)
        dissipationTranspose(i)%rhsBoundary2(4,2:4) =                                        &
             dissipationTranspose(i)%rhsInterior(-1:1)
        dissipationTranspose(i)%rhsBoundary2(3,3:4) =                                        &
             dissipationTranspose(i)%rhsInterior(-1:0)
        dissipationTranspose(i)%rhsBoundary2(2,4:4) =                                        &
             dissipationTranspose(i)%rhsInterior(-1:-1)
        dissipationTranspose(i)%rhsBoundary2(1,1:4) = 0.0_WP ! ... the first column is zero

     case ('SBP 4-8')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 5
        artificialDissipation(i)%boundaryWidth = 5
        artificialDissipation(i)%boundaryDepth = 2
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:2) = (/ 6.0_WP, -4.0_WP, 1.0_WP /)
        artificialDissipation(i)%rhsInterior(-1:-2:-1) =                                     &
             artificialDissipation(i)%rhsInterior(1:2)

        artificialDissipation(i)%rhsBoundary1(1:5,1) =                                       &
             artificialDissipation(i)%rhsInterior(-2:2)
        artificialDissipation(i)%rhsBoundary1(1:5,2) =                                       &
             artificialDissipation(i)%rhsInterior(-2:2)

        ! Transpose dissipation operator
        dissipationTranspose(i)%implicit = .false.
        dissipationTranspose(i)%symmetryType = SYMMETRIC
        dissipationTranspose(i)%interiorWidth = 5
        dissipationTranspose(i)%boundaryWidth = 7
        dissipationTranspose(i)%boundaryDepth = 5
        call allocate_operator(dissipationTranspose(i))

        dissipationTranspose(i)%rhsInterior(0:2) = (/ 6.0_WP, -4.0_WP, 1.0_WP /)
        dissipationTranspose(i)%rhsInterior(-1:-2:-1) =                                      &
             dissipationTranspose(i)%rhsInterior(1:2)

        dissipationTranspose(i)%rhsBoundary1(1,1:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:2)
        dissipationTranspose(i)%rhsBoundary1(2,1:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:2)
        dissipationTranspose(i)%rhsBoundary1(3,1:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:2)
        dissipationTranspose(i)%rhsBoundary1(4,2:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:1)
        dissipationTranspose(i)%rhsBoundary1(5,3:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:0)
        dissipationTranspose(i)%rhsBoundary1(6,4:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:-1)
        dissipationTranspose(i)%rhsBoundary1(7,5:5) =                                        &
             dissipationTranspose(i)%rhsInterior(-2:-2)

     case ('SBP 1-2 composite dissipation')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 3
        artificialDissipation(i)%boundaryWidth = 2
        artificialDissipation(i)%boundaryDepth = 1
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:1) = (/ -2.0_WP, 1.0_WP /)
        artificialDissipation(i)%rhsInterior(-1:-1:-1) =                                     &
             artificialDissipation(i)%rhsInterior(1:1)
        artificialDissipation(i)%rhsInterior = artificialDissipation(i)%rhsInterior / 4.0_WP

        artificialDissipation(i)%normBoundary = (/ 1.0_WP / 2.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:2,1) = (/ -2.0_WP, 2.0_WP /)

        artificialDissipation(i)%rhsBoundary1 = artificialDissipation(i)%rhsBoundary1 / 4.0_WP

     case ('SBP 2-4 composite dissipation')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 5
        artificialDissipation(i)%boundaryWidth = 6
        artificialDissipation(i)%boundaryDepth = 4
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:2) = (/ -6.0_WP, 4.0_WP, -1.0_WP /)
        artificialDissipation(i)%rhsInterior(-1:-2:-1) =                                     &
             artificialDissipation(i)%rhsInterior(1:2)
        artificialDissipation(i)%rhsInterior =                                               &
             artificialDissipation(i)%rhsInterior / 16.0_WP

        artificialDissipation(i)%normBoundary = (/ 17.0_WP / 48.0_WP, &
             59.0_WP / 48.0_WP, &
             43.0_WP / 48.0_WP, &
             49.0_WP / 48.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:3,1) = (/  -96.0_WP / 17.0_WP, &
             192.0_WP / 17.0_WP, &
             -96.0_WP / 17.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:4,2) = (/  192.0_WP / 59.0_WP, &
             -432.0_WP / 59.0_WP, &
             288.0_WP / 59.0_WP, &
             -48.0_WP / 59.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:5,3) = (/  -96.0_WP / 43.0_WP, &
             288.0_WP / 43.0_WP, &
             -336.0_WP / 43.0_WP, &
             192.0_WP / 43.0_WP, &
             -48.0_WP / 43.0_WP /)
        artificialDissipation(i)%rhsBoundary1(2:6,4) = (/  -48.0_WP / 49.0_WP, &
             192.0_WP / 49.0_WP, &
             -288.0_WP / 49.0_WP, &
             192.0_WP / 49.0_WP, &
             -48.0_WP / 49.0_WP /)

        artificialDissipation(i)%rhsBoundary1 =                                              &
             artificialDissipation(i)%rhsBoundary1 / 16.0_WP

     case ('SBP 3-6 composite dissipation')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 7
        artificialDissipation(i)%boundaryWidth = 9
        artificialDissipation(i)%boundaryDepth = 6
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:3) = (/ -20.0_WP, 15.0_WP, -6.0_WP, 1.0_WP /)
        artificialDissipation(i)%rhsInterior(-1:-3:-1) =                                    &
             artificialDissipation(i)%rhsInterior(1:3)
        artificialDissipation(i)%rhsInterior = artificialDissipation(i)%rhsInterior / 64.0_WP

        artificialDissipation(i)%normBoundary = (/ 13649.0_WP / 43200.0_WP, &
             12013.0_WP /  8640.0_WP, &
             2711.0_WP /  4320.0_WP, &
             5359.0_WP /  4320.0_WP, &
             7877.0_WP /  8640.0_WP, &
             43801.0_WP / 43200.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:4,1) = (/ -129600.0_WP / 13649.0_WP, &
             388800.0_WP / 13649.0_WP, &
             -388800.0_WP / 13649.0_WP, &
             129600.0_WP / 13649.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:5,2) = (/   77760.0_WP / 12013.0_WP, &
             -241920.0_WP / 12013.0_WP, &
             259200.0_WP / 12013.0_WP, &
             -103680.0_WP / 12013.0_WP, &
             8640.0_WP / 12013.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:6,3) = (/  -38880.0_WP /  2711.0_WP, &
             129600.0_WP /  2711.0_WP, &
             -159840.0_WP /  2711.0_WP, &
             90720.0_WP /  2711.0_WP, &
             -25920.0_WP /  2711.0_WP, &
             4320.0_WP /  2711.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:7,4) = (/   12960.0_WP /  5359.0_WP, &
             -51840.0_WP /  5359.0_WP, &
             90720.0_WP /  5359.0_WP, &
             -95040.0_WP /  5359.0_WP, &
             64800.0_WP /  5359.0_WP, &
             -25920.0_WP /  5359.0_WP, &
             4320.0_WP /  5359.0_WP /)
        artificialDissipation(i)%rhsBoundary1(2:8,5) = (/    8640.0_WP /  7877.0_WP, &
             -51840.0_WP /  7877.0_WP, &
             129600.0_WP /  7877.0_WP, &
             -172800.0_WP /  7877.0_WP, &
             129600.0_WP /  7877.0_WP, &
             -51840.0_WP /  7877.0_WP, &
             8640.0_WP /  7877.0_WP /)
        artificialDissipation(i)%rhsBoundary1(3:9,6) = (/   43200.0_WP / 43801.0_WP, &
             -259200.0_WP / 43801.0_WP, &
             648000.0_WP / 43801.0_WP, &
             -864000.0_WP / 43801.0_WP, &
             648000.0_WP / 43801.0_WP, &
             -259200.0_WP / 43801.0_WP, &
             43200.0_WP / 43801.0_WP /)

        artificialDissipation(i)%rhsBoundary1 =                                             &
             artificialDissipation(i)%rhsBoundary1 / 64.0_WP

     case ('SBP 4-8 composite dissipation')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 9
        artificialDissipation(i)%boundaryWidth = 12
        artificialDissipation(i)%boundaryDepth = 8
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:4) = (/ -70.0_WP, 56.0_WP, -28.0_WP, 8.0_WP,  &
             -1.0_WP /)
        artificialDissipation(i)%rhsInterior(-1:-4:-1) =                                     &
             artificialDissipation(i)%rhsInterior(1:4)
        artificialDissipation(i)%rhsInterior =                                               &
             artificialDissipation(i)%rhsInterior / 256.0_WP

        artificialDissipation(i)%normBoundary = (/ 1498139.0_WP / 5080320.0_WP, &
             1107307.0_WP /  725760.0_WP, &
             20761.0_WP /   80640.0_WP, &
             1304999.0_WP /  725760.0_WP, &
             299527.0_WP /  725760.0_WP, &
             103097.0_WP /   80640.0_WP, &
             670091.0_WP /  725760.0_WP, &
             5127739.0_WP / 5080320.0_WP /)

        artificialDissipation(i)%rhsBoundary1(1:5,1)  = (/  -15240960.0_WP / 1498139.0_WP, &
             60963840.0_WP / 1498139.0_WP, &
             -91445760.0_WP / 1498139.0_WP, &
             60963840.0_WP / 1498139.0_WP, &
             -15240960.0_WP / 1498139.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:6,2)  = (/    8709120.0_WP / 1107307.0_WP, &
             -35562240.0_WP / 1107307.0_WP, &
             55157760.0_WP / 1107307.0_WP, &
             -39191040.0_WP / 1107307.0_WP, &
             11612160.0_WP / 1107307.0_WP, &
             -725760.0_WP / 1107307.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:7,3)  = (/   -1451520.0_WP /   20761.0_WP, &
             6128640.0_WP /   20761.0_WP, &
             -10080000.0_WP /   20761.0_WP, &
             8064000.0_WP /   20761.0_WP, &
             -3225600.0_WP /   20761.0_WP, &
             645120.0_WP /   20761.0_WP, &
             -80640.0_WP /   20761.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:8,4)  = (/    8709120.0_WP / 1304999.0_WP, &
             -39191040.0_WP / 1304999.0_WP, &
             72576000.0_WP / 1304999.0_WP, &
             -73301760.0_WP / 1304999.0_WP, &
             46448640.0_WP / 1304999.0_WP, &
             -20321280.0_WP / 1304999.0_WP, &
             5806080.0_WP / 1304999.0_WP, &
             -725760.0_WP / 1304999.0_WP /)
        artificialDissipation(i)%rhsBoundary1(1:9,5)  = (/   -2177280.0_WP /  299527.0_WP, &
             11612160.0_WP /  299527.0_WP, &
             -29030400.0_WP /  299527.0_WP, &
             46448640.0_WP /  299527.0_WP, &
             -52254720.0_WP /  299527.0_WP, &
             40642560.0_WP /  299527.0_WP, &
             -20321280.0_WP /  299527.0_WP, &
             5806080.0_WP /  299527.0_WP, &
             -725760.0_WP /  299527.0_WP /)
        artificialDissipation(i)%rhsBoundary1(2:10,6) = (/     -80640.0_WP /  103097.0_WP, &
             645120.0_WP /  103097.0_WP, &
             -2257920.0_WP /  103097.0_WP, &
             4515840.0_WP /  103097.0_WP, &
             -5644800.0_WP /  103097.0_WP, &
             4515840.0_WP /  103097.0_WP, &
             -2257920.0_WP /  103097.0_WP, &
             645120.0_WP /  103097.0_WP, &
             -80640.0_WP /  103097.0_WP /)
        artificialDissipation(i)%rhsBoundary1(3:11,7) = (/    -725760.0_WP /  670091.0_WP, &
             5806080.0_WP /  670091.0_WP, &
             -20321280.0_WP /  670091.0_WP, &
             40642560.0_WP /  670091.0_WP, &
             -50803200.0_WP /  670091.0_WP, &
             40642560.0_WP /  670091.0_WP, &
             -20321280.0_WP /  670091.0_WP, &
             5806080.0_WP /  670091.0_WP, &
             -725760.0_WP /  670091.0_WP /)
        artificialDissipation(i)%rhsBoundary1(4:12,8) = (/   -5080320.0_WP / 5127739.0_WP, &
             40642560.0_WP / 5127739.0_WP, &
             -142248960.0_WP / 5127739.0_WP, &
             284497920.0_WP / 5127739.0_WP, &
             -355622400.0_WP / 5127739.0_WP, &
             284497920.0_WP / 5127739.0_WP, &
             -142248960.0_WP / 5127739.0_WP, &
             40642560.0_WP / 5127739.0_WP, &
             -5080320.0_WP / 5127739.0_WP /)

        artificialDissipation(i)%rhsBoundary1 = artificialDissipation(i)%rhsBoundary1 / 256.0_WP

     case ('null matrix')

        artificialDissipation(i)%implicit = .false.
        artificialDissipation(i)%symmetryType = SYMMETRIC
        artificialDissipation(i)%interiorWidth = 0
        artificialDissipation(i)%boundaryWidth = 1
        artificialDissipation(i)%boundaryDepth = 1
        call allocate_operator(artificialDissipation(i))

        artificialDissipation(i)%rhsInterior(0:0) = 0.0_WP
        artificialDissipation(i)%rhsBoundary1 = 0.0_WP

        if (.not. compositeDissipation) then
           dissipationTranspose(i)%implicit = .false.
           dissipationTranspose(i)%symmetryType = SYMMETRIC
           dissipationTranspose(i)%interiorWidth = 0
           dissipationTranspose(i)%boundaryWidth = 1
           dissipationTranspose(i)%boundaryDepth = 1
           call allocate_operator(dissipationTranspose(i))

           dissipationTranspose(i)%rhsInterior(0:0) = 0.0_WP
           dissipationTranspose(i)%rhsBoundary1 = 0.0_WP
        end if

     case default

        call die('dissipation_setup: unknown stencil scheme: ' // trim(stencilScheme))

     end select

     ! Fill the right-boundary coefficients
     if (allocated(artificialDissipation(i)%rhsBoundary1) .and.                              &
          allocated(artificialDissipation(i)%rhsBoundary2)) then
        select case (artificialDissipation(i)%symmetryType)
        case (SYMMETRIC)
           artificialDissipation(i)%rhsBoundary2(1:artificialDissipation(i)%boundaryWidth,:) =&
                +artificialDissipation(i)%rhsBoundary1(artificialDissipation(i)%boundaryWidth:1:-1,:)
        case (SKEW_SYMMETRIC)
           artificialDissipation(i)%rhsBoundary2(1:artificialDissipation(i)%boundaryWidth,:) =&
                -artificialDissipation(i)%rhsBoundary1(artificialDissipation(i)%boundaryWidth:1:-1,:)
        end select
     end if
     if (allocated(artificialDissipation(i)%lhsBoundary1) .and.                              &
          allocated(artificialDissipation(i)%lhsBoundary2)) then
        artificialDissipation(i)%lhsBoundary2(                                               &
             -artificialDissipation(i)%nDiagonals/2:artificialDissipation(i)%ndiagonals/2,:) =&
             artificialDissipation(i)%lhsBoundary1(                                          &
             artificialDissipation(i)%nDiagonals/2:-artificialDissipation(i)%nDiagonals/2:-1,:)
     end if

     ! Ali's correction to right-boundary for SBP 3-6 composite dissipation
     if (trim(stencilScheme) .eq. 'SBP 3-6 composite dissipation') then
        artificialDissipation(i)%rhsBoundary2(6:9,4) = (/ -90720.0_WP /  5359.0_WP, &
             77760.0_WP /  5359.0_WP, &
             -38880.0_WP /  5359.0_WP, &
             8640.0_WP /  5359.0_WP /) / 64.0_WP
        artificialDissipation(i)%rhsBoundary2(6:9,3) = (/ 77760.0_WP /  2711.0_WP, &
             -120960.0_WP /  2711.0_WP, &
             90720.0_WP  /  2711.0_WP, &
             -25920.0_WP /  2711.0_WP /) / 64.0_WP
        artificialDissipation(i)%rhsBoundary2(6:9,2) = (/ -77760.0_WP /  12013.0_WP, &
             181440.0_WP /  12013.0_WP, &
             -164160.0_WP /  12013.0_WP, &
             51840.0_WP /  12013.0_WP /) / 64.0_WP
        artificialDissipation(i)%rhsBoundary2(6:9,1) = (/ 86400.0_WP /  13649.0_WP, &
             -259200.0_WP /  13649.0_WP, &
             259200.0_WP /  13649.0_WP, &
             -86400.0_WP /  13649.0_WP /) / 64.0_WP
     end if

     ! Fill the right-boundary coefficients for dissipation transpose
     if (.not. compositeDissipation) then
        if (allocated(dissipationTranspose(i)%rhsBoundary1) .and.                            &
             allocated(dissipationTranspose(i)%rhsBoundary2)) then
           select case (dissipationTranspose(i)%symmetryType)
           case (SYMMETRIC)
              dissipationTranspose(i)%rhsBoundary2(1:dissipationTranspose(i)%boundaryWidth,:) =&
                   +dissipationTranspose(i)%rhsBoundary1(dissipationTranspose(i)%boundaryWidth:1:-1,:)
           case (SKEW_SYMMETRIC)
              dissipationTranspose(i)%rhsBoundary2(1:dissipationTranspose(i)%boundaryWidth,:) =&
                   -dissipationTranspose(i)%rhsBoundary1(dissipationTranspose(i)%boundaryWidth:1:-1,:)
           end select
        end if
        if (allocated(dissipationTranspose(i)%lhsBoundary1) .and.                               &
             allocated(dissipationTranspose(i)%lhsBoundary2)) then
           dissipationTranspose(i)%lhsBoundary2(                                                &
                -dissipationTranspose(i)%nDiagonals/2:dissipationTranspose(i)%ndiagonals/2,:) = &
                dissipationTranspose(i)%lhsBoundary1(                                           &
                dissipationTranspose(i)%nDiagonals/2:-dissipationTranspose(i)%nDiagonals/2:-1,:)
        end if
     end if

     ! Update the operators
     call operator_update(artificialDissipation(i), i)
     if (.not. compositeDissipation) call operator_update(dissipationTranspose(i), i)

     ! Store 2nd-order stencil if hybrid dissipation is used
     if (hybridDissipation) then
        artificialDissipation2(i)%implicit = .false.
        artificialDissipation2(i)%symmetryType = SYMMETRIC
        artificialDissipation2(i)%interiorWidth = 3
        artificialDissipation2(i)%boundaryWidth = 2
        artificialDissipation2(i)%boundaryDepth = 1
        call allocate_operator(artificialDissipation2(i))

        artificialDissipation2(i)%rhsInterior(0:1) = (/ -2.0_WP, 1.0_WP /)
        artificialDissipation2(i)%rhsInterior(-1:-1:-1) =                                    &
             artificialDissipation2(i)%rhsInterior(1:1)
        artificialDissipation2(i)%rhsInterior = artificialDissipation2(i)%rhsInterior / 4.0_WP

        artificialDissipation2(i)%normBoundary = (/ 1.0_WP / 2.0_WP /)

        artificialDissipation2(i)%rhsBoundary1(1:2,1) = (/ -2.0_WP, 2.0_WP /)

        artificialDissipation2(i)%rhsBoundary1 = artificialDissipation2(i)%rhsBoundary1 / 4.0_WP

        ! Fill the right-boundary coefficients for low order derivative
        artificialDissipation2(i)%rhsBoundary2(1:artificialDissipation2(i)%boundaryWidth,:) =&
             + artificialDissipation2(i)%rhsBoundary1(                                       &
             artificialDissipation2(i)%boundaryWidth:1:-1,:) !... SYMMETRIC

        ! Update operator based on higher-order stencil
        artificialDissipation2(i)%hasDomainBoundary = artificialDissipation(i)%hasDomainBoundary
        artificialDissipation2(i)%nGhost = artificialDissipation(i)%nGhost
        artificialDissipation2(i)%periodicOffset = artificialDissipation(i)%periodicOffset
     end if
     
  end do

  return
end subroutine dissipation_setup


! ================================ !
! Cleanup the dissipation operator !
! ================================ !
subroutine dissipation_cleanup

  ! Internal modules
  use dissipation

  implicit none

  ! Local variables
  integer :: i

  if (allocated(artificialDissipation)) then
     do i = 1, size(artificialDissipation)
        call deallocate_operator(artificialDissipation(i))
     end do
     deallocate(artificialDissipation)
  end if

  if (allocated(dissipationTranspose)) then
     do i = 1, size(dissipationTranspose)
        call deallocate_operator(dissipationTranspose(i))
     end do
     deallocate(dissipationTranspose)
  end if

  if (allocated(artificialDissipation2)) then
     do i = 1, size(artificialDissipation2)
        call deallocate_operator(artificialDissipation2(i))
     end do
     deallocate(artificialDissipation2)
  end if

  if (allocated(dissipationSource)) deallocate(dissipationSource)
  if (allocated(dissipationSensor)) deallocate(dissipationSensor)

  return
end subroutine dissipation_cleanup
