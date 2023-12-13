module second_derivative

  ! External modules
  use precision
  use operator

  implicit none

  type(t_StencilOperator), allocatable, dimension(:) :: secondDerivative,                    &
       adjointSecondDerivative

contains

  subroutine second_derivative_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    if (size(x, 1) .ne. nGridPoints)                                                         &
         call die ('second_derivative_apply: size of `x` inconsistent with local grid size')

    ! Explicit step
    select case (direction)
    case (1)
       call operator_apply_1(secondDerivative(1), x)
    case (2)
       call operator_apply_2(secondDerivative(2), x)
    case (3)
       call operator_apply_3(secondDerivative(3), x)
    end select

    ! Implicit step
    select case (direction)
    case (1)
       call operator_apply_implicit_1(secondDerivative(1), x)
    case (2)
       call operator_apply_implicit_2(secondDerivative(2), x)
    case (3)
       call operator_apply_implicit_3(secondDerivative(3), x)
    end select

    return
  end subroutine second_derivative_apply

  subroutine adjoint_second_derivative_apply(direction, x)

    implicit none

    ! Arguments
    integer, intent(in) :: direction
    real(WP), intent(inout) :: x(:,:)

    select case (direction)
    case (1)
       call operator_apply_1(adjointSecondDerivative(1), x)
    case (2)
       call operator_apply_2(adjointSecondDerivative(2), x)
    case (3)
       call operator_apply_3(adjointSecondDerivative(3), x)
    end select

    return
  end subroutine adjoint_second_derivative_apply

end module second_derivative


! ==================================== !
! Setup the second derivative operator !
! ==================================== !
subroutine second_derivative_setup

  ! Internal modules
  use second_derivative

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

  allocate(secondDerivative(nDimensions))
  if (.not. predictionOnly) allocate(adjointSecondDerivative(nDimensions))

  do i = 1, nDimensions

     if (globalGridSize(i) .gt. 1) then
        call parser_read('second derivative scheme', stencilScheme, trim(discretizationType))
     else
        stencilScheme = 'null matrix'
     end if

     select case (trim(stencilScheme))

     case ('SBP 1-2')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 3
        secondDerivative(i)%boundaryWidth = 3
        secondDerivative(i)%boundaryDepth = 1
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:1) = (/ -2.0_WP, 1.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-1:-1) = secondDerivative(i)%rhsInterior(1:1)

        secondDerivative(i)%normBoundary = (/ 1.0_WP / 2.0_WP /)

        secondDerivative(i)%rhsBoundary1(1:3,1) = (/ 1.0_WP, -2.0_WP, 1.0_WP /)

     case ('SBP 2-4')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 5
        secondDerivative(i)%boundaryWidth = 6
        secondDerivative(i)%boundaryDepth = 4
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:2) = (/ -5.0_WP / 2.0_WP, &
             4.0_WP / 3.0_WP, &
             -1.0_WP / 12.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-2:-1) = secondDerivative(i)%rhsInterior(1:2)

        secondDerivative(i)%normBoundary = (/ 17.0_WP / 48.0_WP, &
             59.0_WP / 48.0_WP, &
             43.0_WP / 48.0_WP, &
             49.0_WP / 48.0_WP /)

        secondDerivative(i)%rhsBoundary1(1:4,1) = (/ 2.0_WP, &
             -5.0_WP, &
             4.0_WP, &
             -1.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:3,2) = (/ 1.0_WP, &
             -2.0_WP, &
             1.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:5,3) = (/ -4.0_WP / 43.0_WP, &
             59.0_WP / 43.0_WP, &
             -110.0_WP / 43.0_WP, &
             59.0_WP / 43.0_WP, &
             -4.0_WP / 43.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:6,4) = (/ -1.0_WP / 49.0_WP, &
             0.0_WP, &
             59.0_WP / 49.0_WP, &
             -118.0_WP / 49.0_WP, &
             64.0_WP / 49.0_WP, &
             -4.0_WP / 49.0_WP /)

     case ('SBP 3-6')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 7
        secondDerivative(i)%boundaryWidth = 9
        secondDerivative(i)%boundaryDepth = 6
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:3) = (/ -49.0_WP / 18.0_WP, &
             3.0_WP /  2.0_WP, &
             -3.0_WP / 20.0_WP, &
             1.0_WP / 90.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-3:-1) = secondDerivative(i)%rhsInterior(1:3)

        secondDerivative(i)%normBoundary = (/ 13649.0_WP / 43200.0_WP, &
             12013.0_WP /  8640.0_WP, &
             2711.0_WP /  4320.0_WP, &
             5359.0_WP /  4320.0_WP, &
             7877.0_WP /  8640.0_WP, &
             43801.0_WP / 43200.0_WP /)

        secondDerivative(i)%rhsBoundary1(1:6,1) = (/  114170.0_WP /  40947.0_WP, &
             -438107.0_WP /  54596.0_WP, &
             336409.0_WP /  40947.0_WP, &
             -276997.0_WP /  81894.0_WP, &
             3747.0_WP /  13649.0_WP, &
             21035.0_WP / 163788.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:6,2) = (/    6173.0_WP /   5860.0_WP, &
             -2066.0_WP /    879.0_WP, &
             3283.0_WP /   1758.0_WP, &
             -303.0_WP /    293.0_WP, &
             2111.0_WP /   3516.0_WP, &
             -601.0_WP /   4395.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:6,3) = (/  -52391.0_WP /  81330.0_WP, &
             134603.0_WP /  32532.0_WP, &
             -21982.0_WP /   2711.0_WP, &
             112915.0_WP /  16266.0_WP, &
             -46969.0_WP /  16266.0_WP, &
             30409.0_WP /  54220.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:7,4) = (/   68603.0_WP / 321540.0_WP, &
             -12423.0_WP /  10718.0_WP, &
             112915.0_WP /  32154.0_WP, &
             -75934.0_WP /  16077.0_WP, &
             53369.0_WP /  21436.0_WP, &
             -54899.0_WP / 160770.0_WP, &
             48.0_WP /   5359.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,5) = (/   -7053.0_WP /  39385.0_WP, &
             86551.0_WP /  94524.0_WP, &
             -46969.0_WP /  23631.0_WP, &
             53369.0_WP /  15754.0_WP, &
             -87904.0_WP /  23631.0_WP, &
             820271.0_WP / 472620.0_WP, &
             -1296.0_WP /   7877.0_WP, &
             96.0_WP /   7877.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:9,6) = (/   21035.0_WP / 525612.0_WP, &
             -24641.0_WP / 131403.0_WP, &
             30409.0_WP /  87602.0_WP, &
             -54899.0_WP / 131403.0_WP, &
             820271.0_WP / 525612.0_WP, &
             -117600.0_WP /  43801.0_WP, &
             64800.0_WP /  43801.0_WP, &
             -6480.0_WP /  43801.0_WP, &
             480.0_WP /  43801.0_WP /)

     case ('SBP 4-8')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 9
        secondDerivative(i)%boundaryWidth = 12
        secondDerivative(i)%boundaryDepth = 8
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:4) = (/ -205.0_WP /  72.0_WP, &
             8.0_WP /   5.0_WP, &
             -1.0_WP /   5.0_WP, &
             8.0_WP / 315.0_WP, &
             -1.0_WP / 560.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-4:-1) = secondDerivative(i)%rhsInterior(1:4)

        secondDerivative(i)%normBoundary = (/ 1498139.0_WP / 5080320.0_WP, &
             1107307.0_WP /  725760.0_WP, &
             20761.0_WP /   80640.0_WP, &
             1304999.0_WP /  725760.0_WP, &
             299527.0_WP /  725760.0_WP, &
             103097.0_WP /   80640.0_WP, &
             670091.0_WP /  725760.0_WP, &
             5127739.0_WP / 5080320.0_WP /)

        secondDerivative(i)%rhsBoundary1(1:7,1)  = (/  4870382994799.0_WP / 1358976868290.0_WP, &
             -893640087518.0_WP /   75498714905.0_WP, &
             926594825119.0_WP /   60398971924.0_WP, &
             -1315109406200.0_WP /  135897686829.0_WP, &
             39126983272.0_WP /   15099742981.0_WP, &
             12344491342.0_WP /   75498714905.0_WP, &
             -451560522577.0_WP / 2717953736580.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,2)  = (/   333806012194.0_WP /  390619153855.0_WP, &
             -154646272029.0_WP /  111605472530.0_WP, &
             1168338040.0_WP /   33481641759.0_WP, &
             82699112501.0_WP /  133926567036.0_WP, &
             -171562838.0_WP /   11160547253.0_WP, &
             -28244698346.0_WP /  167408208795.0_WP, &
             11904122576.0_WP /  167408208795.0_WP, &
             -2598164715.0_WP /  312495323084.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,3)  =  (/    7838984095.0_WP /   52731029988.0_WP, &
             1168338040.0_WP /    5649753213.0_WP, &
             -88747895.0_WP /     144865467.0_WP, &
             423587231.0_WP /     627750357.0_WP, &
             -43205598281.0_WP /   22599012852.0_WP, &
             4876378562.0_WP /    1883251071.0_WP, &
             -5124426509.0_WP /    3766502142.0_WP, &
             10496900965.0_WP /   39548272491.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,4)  = (/   -94978241528.0_WP /  828644350023.0_WP, &
             82699112501.0_WP /  157837019052.0_WP, &
             1270761693.0_WP /   13153084921.0_WP, &
             -167389605005.0_WP /  118377764289.0_WP, &
             48242560214.0_WP /   39459254763.0_WP, &
             -31673996013.0_WP /   52612339684.0_WP, &
             43556319241.0_WP /  118377764289.0_WP, &
             -44430275135.0_WP /  552429566682.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:9,5)  = (/     1455067816.0_WP /   21132528431.0_WP, &
             -171562838.0_WP /    3018932633.0_WP, &
             -43205598281.0_WP /   36227191596.0_WP, &
             48242560214.0_WP /    9056797899.0_WP, &
             -52276055645.0_WP /    6037865266.0_WP, &
             57521587238.0_WP /    9056797899.0_WP, &
             -80321706377.0_WP /   36227191596.0_WP, &
             8078087158.0_WP /   21132528431.0_WP, &
             -1296.0_WP /        299527.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:10,6) = (/    10881504334.0_WP /  327321118845.0_WP, &
             -28244698346.0_WP /  140280479505.0_WP, &
             4876378562.0_WP /    9352031967.0_WP, &
             -10557998671.0_WP /   12469375956.0_WP, &
             57521587238.0_WP /   28056095901.0_WP, &
             -278531401019.0_WP /   93520319670.0_WP, &
             73790130002.0_WP /   46760159835.0_WP, &
             -137529995233.0_WP /  785570685228.0_WP, &
             2048.0_WP /        103097.0_WP, &
             -144.0_WP /        103097.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:11,7) = (/  -135555328849.0_WP / 8509847458140.0_WP, &
             11904122576.0_WP /  101307707835.0_WP, &
             -5124426509.0_WP /   13507694378.0_WP, &
             43556319241.0_WP /   60784624701.0_WP, &
             -80321706377.0_WP /   81046166268.0_WP, &
             73790130002.0_WP /   33769235945.0_WP, &
             -950494905688.0_WP /  303923123505.0_WP, &
             239073018673.0_WP /  141830790969.0_WP, &
             -145152.0_WP /        670091.0_WP, &
             18432.0_WP /        670091.0_WP, &
             -1296.0_WP /        670091.0_WP /)
        secondDerivative(i)%rhsBoundary1(2:12,8) = (/    -2598164715.0_WP /  206729925524.0_WP, &
             10496900965.0_WP /  155047444143.0_WP, &
             -44430275135.0_WP /  310094888286.0_WP, &
             425162482.0_WP /    2720130599.0_WP, &
             -137529995233.0_WP /  620189776572.0_WP, &
             239073018673.0_WP /  155047444143.0_WP, &
             -144648000000.0_WP /   51682481381.0_WP, &
             8128512.0_WP /       5127739.0_WP, &
             -1016064.0_WP /       5127739.0_WP, &
             129024.0_WP /       5127739.0_WP, &
             -9072.0_WP /       5127739.0_WP /)

     case ('SBP 2-2-8')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 9
        secondDerivative(i)%boundaryWidth = 12
        secondDerivative(i)%boundaryDepth = 8
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:4) = (/ -205.0_WP /  72.0_WP, &
             8.0_WP /   5.0_WP, &
             -1.0_WP /   5.0_WP, &
             8.0_WP / 315.0_WP, &
             -1.0_WP / 560.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-4:-1) = secondDerivative(i)%rhsInterior(1:4)

        secondDerivative(i)%normBoundary = (/ 1498139.0_WP / 5080320.0_WP, &
             1107307.0_WP /  725760.0_WP, &
             20761.0_WP /   80640.0_WP, &
             1304999.0_WP /  725760.0_WP, &
             299527.0_WP /  725760.0_WP, &
             103097.0_WP /   80640.0_WP, &
             670091.0_WP /  725760.0_WP, &
             5127739.0_WP / 5080320.0_WP /)

        secondDerivative(i)%rhsBoundary1(1:7,1)  = (/  4870382994799.0_WP / 1358976868290.0_WP, &
             -893640087518.0_WP /   75498714905.0_WP, &
             926594825119.0_WP /   60398971924.0_WP, &
             -1315109406200.0_WP /  135897686829.0_WP, &
             39126983272.0_WP /   15099742981.0_WP, &
             12344491342.0_WP /   75498714905.0_WP, &
             -451560522577.0_WP / 2717953736580.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,2)  = (/   333806012194.0_WP /  390619153855.0_WP, &
             -154646272029.0_WP /  111605472530.0_WP, &
             1168338040.0_WP /   33481641759.0_WP, &
             82699112501.0_WP /  133926567036.0_WP, &
             -171562838.0_WP /   11160547253.0_WP, &
             -28244698346.0_WP /  167408208795.0_WP, &
             11904122576.0_WP /  167408208795.0_WP, &
             -2598164715.0_WP /  312495323084.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,3)  =  (/    7838984095.0_WP /   52731029988.0_WP, &
             1168338040.0_WP /    5649753213.0_WP, &
             -88747895.0_WP /     144865467.0_WP, &
             423587231.0_WP /     627750357.0_WP, &
             -43205598281.0_WP /   22599012852.0_WP, &
             4876378562.0_WP /    1883251071.0_WP, &
             -5124426509.0_WP /    3766502142.0_WP, &
             10496900965.0_WP /   39548272491.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:8,4)  = (/   -94978241528.0_WP /  828644350023.0_WP, &
             82699112501.0_WP /  157837019052.0_WP, &
             1270761693.0_WP /   13153084921.0_WP, &
             -167389605005.0_WP /  118377764289.0_WP, &
             48242560214.0_WP /   39459254763.0_WP, &
             -31673996013.0_WP /   52612339684.0_WP, &
             43556319241.0_WP /  118377764289.0_WP, &
             -44430275135.0_WP /  552429566682.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:9,5)  = (/     1455067816.0_WP /   21132528431.0_WP, &
             -171562838.0_WP /    3018932633.0_WP, &
             -43205598281.0_WP /   36227191596.0_WP, &
             48242560214.0_WP /    9056797899.0_WP, &
             -52276055645.0_WP /    6037865266.0_WP, &
             57521587238.0_WP /    9056797899.0_WP, &
             -80321706377.0_WP /   36227191596.0_WP, &
             8078087158.0_WP /   21132528431.0_WP, &
             -1296.0_WP /        299527.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:10,6) = (/    10881504334.0_WP /  327321118845.0_WP, &
             -28244698346.0_WP /  140280479505.0_WP, &
             4876378562.0_WP /    9352031967.0_WP, &
             -10557998671.0_WP /   12469375956.0_WP, &
             57521587238.0_WP /   28056095901.0_WP, &
             -278531401019.0_WP /   93520319670.0_WP, &
             73790130002.0_WP /   46760159835.0_WP, &
             -137529995233.0_WP /  785570685228.0_WP, &
             2048.0_WP /        103097.0_WP, &
             -144.0_WP /        103097.0_WP /)
        secondDerivative(i)%rhsBoundary1(1:11,7) = (/  -135555328849.0_WP / 8509847458140.0_WP, &
             11904122576.0_WP /  101307707835.0_WP, &
             -5124426509.0_WP /   13507694378.0_WP, &
             43556319241.0_WP /   60784624701.0_WP, &
             -80321706377.0_WP /   81046166268.0_WP, &
             73790130002.0_WP /   33769235945.0_WP, &
             -950494905688.0_WP /  303923123505.0_WP, &
             239073018673.0_WP /  141830790969.0_WP, &
             -145152.0_WP /        670091.0_WP, &
             18432.0_WP /        670091.0_WP, &
             -1296.0_WP /        670091.0_WP /)
        secondDerivative(i)%rhsBoundary1(2:12,8) = (/    -2598164715.0_WP /  206729925524.0_WP, &
             10496900965.0_WP /  155047444143.0_WP, &
             -44430275135.0_WP /  310094888286.0_WP, &
             425162482.0_WP /    2720130599.0_WP, &
             -137529995233.0_WP /  620189776572.0_WP, &
             239073018673.0_WP /  155047444143.0_WP, &
             -144648000000.0_WP /   51682481381.0_WP, &
             8128512.0_WP /       5127739.0_WP, &
             -1016064.0_WP /       5127739.0_WP, &
             129024.0_WP /       5127739.0_WP, &
             -9072.0_WP /       5127739.0_WP /)

     case ('Pade 3-6')

        secondDerivative(i)%implicit = .true.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 5
        secondDerivative(i)%boundaryWidth = 4
        secondDerivative(i)%boundaryDepth = 2
        secondDerivative(i)%nDiagonals    = 3
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%normBoundary = 1.0_WP

        ! Interior coefficients
        alpha = 2.0_WP / 11.0_WP
        secondDerivative(i)%lhsInterior(-1:+1) = (/ alpha, 1.0_WP, alpha /)

        a = 12.0_WP / 11.0_WP
        b = 3.0_WP / 11.0_WP
        secondDerivative(i)%rhsInterior(0:2) = (/ -2.0_WP * (b / 4.0_WP + a), a, b / 4.0_WP /)
        secondDerivative(i)%rhsInterior(-1:-2:-1) = secondDerivative(i)%rhsInterior(1:2)

        ! (i = 1) on boundary second derivative
        alpha = 11.0_WP
        secondDerivative(i)%lhsBoundary1(0:1,1) = (/ 1.0_WP, alpha /)
        secondDerivative(i)%rhsBoundary1(1:4,1) = (/ 13.0_WP, -27.0_WP, 15.0_WP, -1.0_WP /)

        ! (i = 2) boundary second derivative
        alpha = 0.1_WP
        secondDerivative(i)%lhsBoundary1(-1:1,2) = (/ alpha, 1.0_WP, alpha /)
        secondDerivative(i)%rhsBoundary1(1:3,2) = (/ 12.0_WP / 1.0_WP, -24.0_WP / 10.0_WP,   &
             12.0_WP / 1.0_WP /)

     case ('DRP 13 point')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 13
        secondDerivative(i)%boundaryWidth = 11
        secondDerivative(i)%boundaryDepth = 6
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:6) =                                              &
             (/ -3.0932689374140908764950319934624E+000_WP,&
             +1.8117592073153699737816408065692E+000_WP,&
             -3.3438269916183987705583208818595E-001_WP,&
             +8.7349401467551123473691980072727E-002_WP,&
             -2.1869332809975397294053625502616E-002_WP,&
             +4.2243461063809429482875929492585E-003_WP,&
             -4.4645421044132760621866917137460E-004_WP /)
        secondDerivative(i)%rhsInterior(-1:-6:-1) = secondDerivative(i)%rhsInterior(1:6)

        secondDerivative(i)%normBoundary = 1.0_WP

        secondDerivative(i)%rhsBoundary1(1:7,1) =                                              &
             (/ +4.7331413118533270375404984391079E+000_WP,&
             -1.8732181204453295558576323967981E+001_WP,&
             +3.2580453011133238896440809919952E+001_WP,&
             -3.2662826237066540750809968782159E+001_WP,&
             +1.9830453011133238896440809919952E+001_WP,&
             -6.7321812044532955585763239679810E+000_WP,&
             +9.8314131185332703754049843910794E-001_WP /)

        secondDerivative(i)%rhsBoundary1(1:7,2) =                                              &
             (/ +8.2022895972330285853254027778687E-001_WP,&
             -1.1713737583398171511952416667212E+000_WP,&
             -5.2989893748379045534522916653024E-001_WP,&
             +1.4287541388672761626825277775959E+000_WP,&
             -6.9656560415045712201189583319691E-001_WP,&
             +1.6195957499351618213809166661210E-001_WP,&
             -1.3104373610030474800793055546461E-002_WP /)

        secondDerivative(i)%rhsBoundary1(1:7,3) =                                              &
             (/ -7.7535624543688550227579420261105E-002_WP,&
             +1.2985470805954646346988098549000E+000_WP,&
             -2.4130343681553282534136913039166E+000_WP,&
             +1.2173791575404376712182550718888E+000_WP,&
             +3.6322985113384132529753627500979E-003_WP,&
             -3.4786252737868698634523478433373E-002_WP,&
             +5.7977087896447831057539130722288E-003_WP /)

        secondDerivative(i)%rhsBoundary1(1:7,4) =                                              &
             (/ +1.5660656422279830804838779522893E-002_WP,&
             -1.7729727186701231816236601047069E-001_WP,&
             +1.5682431796675307954059150261767E+000_WP,&
             -2.8132131284455966160967755904579E+000_WP,&
             +1.5682431796675307954059150261767E+000_WP,&
             -1.7729727186701231816236601047069E-001_WP,&
             +1.5660656422279830804838779522893E-002_WP /)

        secondDerivative(i)%rhsBoundary1(1:9,5) =                                              &
             (/ -3.5741930674022842916259048164782E-003_WP,&
             +3.8095462267220224407907676108308E-002_WP,&
             -2.4042224558860899394826129365362E-001_WP,&
             +1.6760169110278905047878905667033E+000_WP,&
             -2.9402318692781989019118220886831E+000_WP,&
             +1.6760169110278905047878905667033E+000_WP,&
             -2.4042224558860899394826129365362E-001_WP,&
             +3.8095462267220224407907676108308E-002_WP,&
             -3.5741930674022842916259048164782E-003_WP /)

        secondDerivative(i)%rhsBoundary1(1:11,6) =                                             &
             (/ +8.9801406812550448657719127266081E-004_WP,&
             -9.7460853822832689972936919112537E-003_WP,&
             +5.8160065952544240616240402171360E-002_WP,&
             -2.8227272480920862141376147176946E-001_WP,&
             +1.7391373200773310119011515562991E+000_WP,&
             -3.0123531798130177331858279721249E+000_WP,&
             +1.7391373200773310119011515562991E+000_WP,&
             -2.8227272480920862141376147176946E-001_WP,&
             +5.8160065952544240616240402171360E-002_WP,&
             -9.7460853822832689972936919112537E-003_WP,&
             +8.9801406812550448657719127266081E-004_WP /)

     case ('null matrix')

        secondDerivative(i)%implicit = .false.
        secondDerivative(i)%symmetryType = SYMMETRIC
        secondDerivative(i)%interiorWidth = 0
        secondDerivative(i)%boundaryWidth = 1
        secondDerivative(i)%boundaryDepth = 1
        call allocate_operator(secondDerivative(i))

        secondDerivative(i)%rhsInterior(0:0) = 0.0_WP
        secondDerivative(i)%rhsBoundary1 = 0.0_WP

     case default

        call die('second_derivative_setup: unknown stencil scheme: ' // trim(stencilScheme))

     end select

     ! Fill the right-boundary coefficients
     if (allocated(secondDerivative(i)%rhsBoundary1) .and.                                   &
          allocated(secondDerivative(i)%rhsBoundary2)) then
        select case (secondDerivative(i)%symmetryType)
        case (SYMMETRIC)
           secondDerivative(i)%rhsBoundary2(1:secondDerivative(i)%boundaryWidth,:) =         &
                +secondDerivative(i)%rhsBoundary1(secondDerivative(i)%boundaryWidth:1:-1,:)
        case (SKEW_SYMMETRIC)
           secondDerivative(i)%rhsBoundary2(1:secondDerivative(i)%boundaryWidth,:) =         &
                -secondDerivative(i)%rhsBoundary1(secondDerivative(i)%boundaryWidth:1:-1,:)
        end select
     end if
     if (allocated(secondDerivative(i)%lhsBoundary1) .and.                                   &
          allocated(secondDerivative(i)%lhsBoundary2)) then
        secondDerivative(i)%lhsBoundary2(                                                    &
             -secondDerivative(i)%nDiagonals/2:secondDerivative(i)%ndiagonals/2,:) =         &
             secondDerivative(i)%lhsBoundary1(                                               &
             secondDerivative(i)%nDiagonals/2:-secondDerivative(i)%nDiagonals/2:-1,:)
     end if

     call operator_update(secondDerivative(i), i)

     ! Adjoint second derivative operator
     if (allocated(adjointSecondDerivative)) then
        if (useContinuousAdjoint .or. trim(stencilScheme) .eq. 'null matrix') then
           adjointSecondDerivative(i)%implicit = .false.
           adjointSecondDerivative(i)%symmetryType = secondDerivative(i)%symmetryType
           adjointSecondDerivative(i)%interiorWidth = secondDerivative(i)%interiorWidth
           adjointSecondDerivative(i)%boundaryWidth = secondDerivative(i)%boundaryWidth
           adjointSecondDerivative(i)%boundaryDepth = secondDerivative(i)%boundaryDepth
           call allocate_operator(adjointSecondDerivative(i))
           adjointSecondDerivative(i)%normBoundary = secondDerivative(i)%normBoundary
           adjointSecondDerivative(i)%rhsInterior = -secondDerivative(i)%rhsInterior
           adjointSecondDerivative(i)%rhsBoundary1 = -secondDerivative(i)%rhsBoundary1
           adjointSecondDerivative(i)%rhsBoundary2 = -secondDerivative(i)%rhsBoundary2
        else
           call operator_adjoint(secondDerivative(i), adjointSecondDerivative(i))
        end if
        call operator_update(adjointSecondDerivative(i), i)
     end if

  end do

  return
end subroutine second_derivative_setup


! ====================================== !
! Cleanup the second derivative operator !
! ====================================== !
subroutine second_derivative_cleanup

  ! Internal modules
  use second_derivative

  implicit none
  
  ! Local variables
  integer :: i

  if (allocated(secondDerivative)) then
     do i = 1, size(secondDerivative)
        call deallocate_operator(secondDerivative(i))
     end do
     deallocate(secondDerivative)
  end if

  if (allocated(adjointSecondDerivative)) then
     do i = 1, size(adjointSecondDerivative)
        call deallocate_operator(adjointSecondDerivative(i))
     end do
     deallocate(adjointSecondDerivative)
  end if

  return
end subroutine second_derivative_cleanup
