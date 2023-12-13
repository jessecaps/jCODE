module thermochem

  ! External modules
  use precision
  use string

  implicit none
  public

  ! Thermodynamic constants
  real(WP), parameter :: universalGasConstant = 8.314472_WP ! [J/(mol.K)]

contains

  ! Standard values of density at STP [kg/m^3]
  ! ------------------------------------------
  subroutine get_density(species, density)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(out) :: density

    select case (trim(species))
    case('AIR', 'air')
       density = 1.205_WP
    case ('Ar', 'ARGON')
       density = 1.661_WP
    case ('CH4', 'METHANE')
       density = 0.668_WP
    case ('CO', 'CARBON_MONOXIDE')
       density = 1.165_WP
    case ('CO2', 'CARBON_DIOXIDE')
       density = 1.842_WP
    case ('H')
       density = 1.00794_WP
    case ('H2', 'HYDROGEN')
       density = 0.0899_WP
    case('H2O', 'WATER')
       density = 18.01528_WP
    case ('N')
       density = 14.0067_WP
    case ('N2', 'NITROGEN')
       density = 28.0134_WP
    case ('O')
       density = 15.9994_WP
    case ('O2', 'OXYGEN')
       density = 31.9988_WP
    case default
       call die('Unknown species: ' // trim(species) // '!')
    end select

    return
  end subroutine get_density

  ! Return molecular weight in [g/mol]
  ! -----------------------------------
  subroutine get_molecular_weight(species, molecularWeight)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(out) :: molecularWeight

    select case (trim(species))
    case('air', 'AIR')
       molecularWeight = 28.97_WP
    case ('Ar', 'ARGON', 'argon')
       molecularWeight = 39.944_WP
    case ('C', 'CARBON', 'carbon')
       molecularWeight = 12.011_WP
    case ('CH')
       molecularWeight = 13.019_WP
    case ('CHO')
       molecularWeight = 29.019_WP
    case ('CH2')
       molecularWeight = 14.027_WP
    case ('CH2O')
       molecularWeight = 30.027_WP
    case ('CHCO')
       molecularWeight = 41.030_WP
    case ('CH3')
       molecularWeight = 15.035_WP
    case ('CH4', 'METHANE', 'methane')
       molecularWeight = 16.043_WP
    case ('C2H')
       molecularWeight = 25.038_WP
    case ('C2H2')
       molecularWeight = 26.046_WP
    case ('C2H3')
       molecularWeight = 27.054_WP
    case ('C2H4')
       molecularWeight = 28.062_WP
    case ('C2H5')
       molecularWeight = 29.070_WP
    case ('C2H6')
       molecularWeight = 30.078_WP
    case ('C3H3')
       molecularWeight = 39.065_WP
    case ('C3H4')
       molecularWeight = 40.073_WP
    case ('C3H5')
       molecularWeight = 41.081_WP
    case ('C3H6')
       molecularWeight = 42.050_WP
    case ('CO', 'CARBON MONOXIDE', 'carbon monoxide')
       molecularWeight = 28.011_WP
    case ('CO2', 'CARBON DIOXIDE', 'carbon dioxide')
       molecularWeight = 44.011_WP
    case ('H')
       molecularWeight = 1.008_WP !1.00794_WP
    case ('H2', 'HYDROGEN', 'hydrogen')
       molecularWeight = 2.016_WP !2.01588_WP
    case('H2O', 'WATER', 'water')
       molecularWeight = 18.016_WP !18.01528_WP
    case ('HO2')
       molecularWeight = 33.008_WP !33.00674_WP
    case ('H2O2', 'HYDROGEN PEROXIDE', 'hydrogen peroxide')
       molecularWeight = 34.016_WP !34.0147_WP
    case ('N')
       molecularWeight = 14.010_WP !14.0067_WP
    case ('N2', 'NITROGEN', 'nitrogen')
       molecularWeight = 28.014_WP !28.0134_WP
    case ('O')
       molecularWeight = 16.000_WP !15.9994_WP
    case ('O2' ,'OXYGEN', 'oxygen')
       molecularWeight = 32.000_WP !31.9988_WP
    case ('OH', 'HYDROXIDE', 'hydroxide')
       molecularWeight = 17.008_WP !17.0074_WP
    case default
       call die('Unknown species: ' // trim(species) // '!')
    end select

    return
  end subroutine get_molecular_weight


  ! Return mass formation enthalpy at 298 K [J/kg]
  ! ----------------------------------------------
  subroutine get_enthalpy_formation(species, enthalpyOfFormation)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(out) :: enthalpyOfFormation

    select case (trim(species))
    case('C')
       enthalpyOfFormation = 5.967E+07_WP !5.95307478394E7_WP
    case('CO')
       enthalpyOfFormation = -3.947E+06_WP !-3.946090681899E6_WP
    case('CO2')
       enthalpyOfFormation = -8.942E+06_WP !-8.9417057680728E6_WP
    case ('H')
       enthalpyOfFormation = 2.162E+08_WP !2.1799E8_WP
    case ('H2', 'HYDROGEN', 'hydrogen')
       enthalpyOfFormation = 0.0_WP
    case('H2O', 'WATER', 'water')
       enthalpyOfFormation = -1.342E+07_WP !-1.3435E7_WP
    case ('HO2')
       enthalpyOfFormation = 3.166E+05_WP !3.80364E5_WP
    case ('H2O2', 'HYDROGEN PEROXIDE', 'hydrogen peroxide')
       enthalpyOfFormation = -4.001E+06_WP !-5.6202E6_WP
    case ('N')
       enthalpyOfFormation = 3.374E+07_WP !3.375456031756E7_WP
    case ('N2', 'NITROGEN', 'nitrogen')
       enthalpyOfFormation = 0.0_WP
    case('NO')
       enthalpyOfFormation = 3.009E+06_WP !3.0090548255E6_WP
    case('NO2')
       enthalpyOfFormation = 7.520539479424E5_WP
    case ('O')
       enthalpyOfFormation = 1.557E+07_WP !1.5574375E7_WP
    case ('O2' ,'OXYGEN', 'oxygen')
       enthalpyOfFormation = 0.0_WP
    case ('OH', 'HYDROXIDE', 'hydroxide')
       enthalpyOfFormation = 2.292E+06_WP !2.19291E6_WP
    case('CH')
       enthalpyOfFormation = 4.564E+07_WP
    case('CHO')
       enthalpyOfFormation = 1.499E+06_WP
    case('CH2')
       enthalpyOfFormation = 2.759E+07_WP
    case('CH2O')
       enthalpyOfFormation = -3.861E+06_WP
    case('CHCO')
       enthalpyOfFormation = 4.328E+06_WP
    case('CH3')
       enthalpyOfFormation = 9.690E+06_WP
    case('CH4')
       enthalpyOfFormation = -4.669E+06_WP
    case('C2H')
       enthalpyOfFormation = 2.239E+07_WP
    case('C2H2')
       enthalpyOfFormation = 8.706E+06_WP
    case('C2H3')
       enthalpyOfFormation = 1.058E+07_WP
    case('C2H4')
       enthalpyOfFormation = 1.869E+06_WP
    case('C2H5')
       enthalpyOfFormation = 4.032E+06_WP
    case('C2H6')
       enthalpyOfFormation = -2.788E+06_WP
    case('C3H4')
       enthalpyOfFormation = 4.723E+06_WP
    case('C3H5')
       enthalpyOfFormation = 3.936E+06_WP
    case('C3H6')
       enthalpyOfFormation = 4.862E+05_WP

    case default
       call die('Unknown species: ' // trim(species) // '!')
    end select

    return
  end subroutine get_enthalpy_formation


  ! Return specific heat capacity at constant pressure [J/(kg K)]
  ! -------------------------------------------------------------
  subroutine get_Cp(species, T, Cp)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(in) :: T
    real(WP), intent(out) :: Cp

    select case (trim(species))
    case ('air')
       Cp = 1005.0_WP
    case ('water', 'H20')
       Cp = 4184.0_WP
    case ('water vapor', 'H20 vapor')
       Cp = 1864.0_WP
    case ('hexane')
       Cp = 78.848_WP+8.8729E-1_WP*T-2.9482E-3_WP*T**2+4.1999E-6_WP*T**3 ! J/(mol.K)
       Cp = Cp/0.086177_WP ! J/(kg.K)
    case ('heptane')
       Cp = 101.121_WP+9.7739E-1_WP*T-3.0712E-3_WP*T**2+4.1844E-6_WP*T**3 ! J/(mol.K)
       Cp = Cp/0.100204_WP ! J/(kg.K)
    case ('dodecane')
       Cp = 84.485_WP+2.0358_WP*T-5.0981E-3_WP*T**2+5.2186E-6_WP*T**3 ! J/(mol.K)
       Cp = Cp/0.170338_WP ! J/(kg.K)
    case ('acetone')
       Cp = 46.878_WP+6.2652E-1_WP*T-2.0761E-3_WP*T**2+2.9583E-6_WP*T**3 ! J/(mol.K)
       Cp = Cp/0.05808_WP ! J/(kg.K)
    case default
       call die('Unknown species: ' // trim(species) // '!')
    end select

    return
  end subroutine get_Cp


  ! Return latent heat of vaporization [J/kg]
  ! -----------------------------------------
  subroutine get_latent_heat(species, T, Lv)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(in) :: T
    real(WP), intent(out) :: Lv

    select case (trim(species))
    case ('water', 'H20')
       Lv = 2.257E6_WP + 2.595E3_WP*(373.15_WP-T)
    case ('hexane')
       Lv = 45.61_WP*(1.0_WP-T/507.43_WP)**(0.401_WP) ! kJ/mol
       Lv = 1.0E3_WP*Lv/0.086177_WP ! J/kg
    case ('heptane')
       Lv = 49.73_WP*(1.0_WP-T/540.26_WP)**(0.386_WP) ! kJ/mol
       Lv = 1.0E3_WP*Lv/0.100204_WP ! J/kg
    case ('dodecane')
       Lv = 77.166_WP*(1.0_WP-T/658.2_WP)**(0.407_WP) ! kJ/mol
       Lv = 1.0E3_WP*Lv/0.17033_WP ! J/kg
    case ('acetone')
       Lv = 49.244_WP*(1.0_WP-T/508.2_WP)**(0.481_WP) ! kJ/mol
       Lv = 1.0E3_WP*Lv/0.05808_WP ! J/kg
    case default
       call die('Unknown species: ' // trim(species) // '!')
    end select

    return
  end subroutine get_latent_heat


  ! Compute the gas constnat [J/kg]
  ! -------------------------------
  subroutine get_gas_constant(species, gasConstant)

    implicit none

    ! Arguments
    character(len = *), intent(in) :: species
    real(WP), intent(out) :: gasConstant

    call get_molecular_weight(trim(species), gasConstant)
    gasConstant = universalGasConstant / gasConstant * 1000.0_WP

    return
  end subroutine get_gas_constant

end module thermochem
