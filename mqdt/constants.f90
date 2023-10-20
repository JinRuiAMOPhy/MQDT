module constants
!----------------------------------------------------------------------------
! CONSTANTS : basic constants for precision, numbers, and std input/output
!             physical constants from NIST
!----------------------------------------------------------------------------
!  use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN, IEEE_IS_NAN
  implicit none
  ! for precision : double precision
  !integer, parameter :: long = selected_real_kind(15)
  !integer, parameter :: long = kind(1.q0)
  integer, parameter :: long = selected_real_kind(24)
  integer, parameter :: sp = kind( 1.0 )
  integer, parameter :: dp = kind( 1.0d0 )
  !integer, parameter :: wp = dp
  integer, parameter :: wp = long
  integer, parameter :: short_str = 32
  integer, parameter :: max_col = 100
  integer, parameter :: longlong_str = short_str * max_col
  integer, parameter :: long_str = longlong_str / 2 
  ! print control
  integer, parameter :: print_col_l_fmt = 12
  integer, parameter :: print_col_s_fmt = 5
  !
  !
  !
  ! iteration control
  integer, parameter :: MAXITER = 2000
  !real(long), parameter :: THRESH = 1.0e-15_long
  real(long), parameter :: THRESH = 1.0e-26_long
  real(long), parameter :: JUMP_THRESH = 1.0e-5_long

  ! for numbers
  !
  real(long), parameter :: EPS = 1.0e-15_long
  real(long) :: nan 
  real(long), parameter :: ZERO  = 0.0_long
  real(long), parameter :: ONE   = 1.0_long
  real(long), parameter :: TWO   = 2.0_long
  real(long), parameter :: THREE = 3.0_long
  real(long), parameter :: FOUR  = 4.0_long
  real(long), parameter :: HALF  = 0.5_long
  !real(long), parameter :: PI    = 3.141592653589793_long
  real(long) :: PI ! internally evaluate in envset.f90
  complex(long), parameter :: IMAG_ONE = ( ZERO, ONE )
  complex(long), parameter :: CMPLX_ZERO = ( ZERO, ZERO )
  complex(long), parameter :: REAL_ONE = ( ONE, ZERO )
  !
  ! standard input/output parameters
  !
  integer, parameter :: STDIN  = 5      ! unit number of standard input
  integer, parameter :: STDOUT = 6      ! unit number of standard output
  integer, parameter :: STDERR = 0      ! unit number of standard error
  integer, parameter :: funit_dbg = 888      ! unit number of debug file
  !
  ! physical constants from NIST
  !   ref) Mohr & Taylor, J. Phys. Chem. Ref. Data 28, 1713 (1999)
  !        Mohr & Taylor, Rev. Mod. Phys. 72, 351 (2000)
  !   E: electron, P: proton, and N: neutron
  !
  real(long), parameter :: RYDBERG=109737.3177_long                  ! 1 Ryd.-> 1cm-1
  real(long), parameter :: SPEED_C_SI = 2.99792458e8_long       ! m/s (exact)
  real(long), parameter :: CHARGE_P_SI = 1.602176462e-19_long   ! C
  real(long), parameter :: CHARGE_E_SI = CHARGE_P_SI
  real(long), parameter :: MASS_E_SI = 9.10938188e-31_long      ! kg
  real(long), parameter :: MASS_P_SI = 1.67262158e-27_long      ! kg
  real(long), parameter :: MASS_N_SI = 1.67492715e-27_long      ! kg
  real(long), parameter :: AMU_SI = 1.66053873e-27_long         ! kg
  real(long), parameter :: PLANCK_H_SI = 6.62606876e-34_long    ! J s
  real(long), parameter :: BOLTZMANN_K_SI = 1.3806503e-23_long  ! J/K, k = R/NA
  real(long), parameter :: AVOGADRO_NA_SI = 6.02214199e23_long  ! 1/mol
  real(long), parameter :: GRAVITY_G_SI = 6.673e-11_long        ! m^3 / kg s^2
  !
  ! atomic units from NIST inspired by Burnus's diploma (2004)
  !   Atomic units defined by h_bar = 1, m_e = 1, e = 1, and a0 = 1.
  !   corollary) epsilon0 = 1/4PI, and E_h = e^2 / 4PI epsilon0 a0 = 1
  !
  ! a.u. of charge: e = 1 a.u. = 1.602176462e-19 C   (CHARGE_E_SI)
  ! a.u. of mass:   m_e = 1 a.u. = 9.10938188e-31 kg (MASS_E_SI)
  ! a.u. of action: h_bar = h / 2PI = 1 a.u.         (PLANCK_H_SI / TWO / PI)
  ! a.u. of length: a0 = h_bar^2 / 4PI epsilon0 m_e e^2 = 1 a.u. (bohr)
  ! a.u. of energy: E_h = e^2 / 4PI epsilon0 a0 = 1 a.u. (hartree)
  ! a.u. of time                        := h_bar / E_h
  ! a.u. of force                       := E_h / a0
  ! a.u. of velocity                    := a0 E_h / h_bar
  ! a.u. of momentum                    := h_bar / a0
  ! a.u. of current                     := e E_h / h_bar
  ! a.u. of charge density              := e / a0^3
  ! a.u. of electric potential          := E_h / e
  ! a.u. of electric field              := E_h / e a0
  ! a.u. of electric field gradient     := E_h / e a0^2
  ! a.u. of electric dipole moment      := e a0
  ! a.u. of electric quadrupole moment  := e a0^2
  ! a.u. of electric polarizability     := e^2 a0^2 / E_h
  ! a.u. of magnetic field              := h_bar / e a0^2
  !
  real(long), parameter :: AU_CHARGE = CHARGE_E_SI              ! 1.602e-19 C
  real(long), parameter :: AU_MASS   = MASS_E_SI                ! 9.109e-31 kg
  !real(long), parameter :: AU_ACTION = PLANCK_H_SI / TWO / PI   ! 1.055e-34 J s 
  ! later evatluate in envset after PI
  real(long) :: AU_ACTION ! = PLANCK_H_SI / TWO / PI   ! 1.055e-34 J s
  real(long), parameter :: AU_LENGTH = 0.5291772083e-10_long    ! m
  real(long), parameter :: AU_ENERGY = 4.35974381e-18_long      ! J
  real(long), parameter :: AU_TIME   = 2.418884326500e-17_long  ! sec
  real(long), parameter :: AU_FORCE  = 8.23872181e-8_long       ! N
  real(long), parameter :: AU_VELOCITY = 2.1876912529e6_long    ! m/s
  real(long), parameter :: AU_MOMENTUM = 1.99285151e24_long     ! kg m/s
  real(long), parameter :: AU_CURRENT  = 6.62361753e-3_long     ! A
  real(long), parameter :: AU_CHARGE_DENSITY   = 1.081202285e12_long    ! C/m^3
  real(long), parameter :: AU_E_POTENTIAL      = 27.2113834_long        ! V
  real(long), parameter :: AU_E_FIELD          = 5.14220624e11_long     ! V/m
  real(long), parameter :: AU_E_FIELD_GRAD     = 9.71736153e21_long     ! V/m^2
  real(long), parameter :: AU_E_DIPOLE_MOMENT  = 8.47835267e-30_long    ! C m
  real(long), parameter :: AU_E_QUAD_MOMENT    = 4.48655100e-40_long    ! C m^2
  real(long), parameter :: AU_E_POLARIZABILITY = 1.648777251e-41_long
                                                                  ! C^2 m^2 / J
  real(long), parameter :: AU_MAGNETIC_FIELD   = 2.35051735e5_long      ! T
  !
  ! a.u. conversion
  !
  real(long), parameter :: AU2J   = AU_ENERGY        ! 1 a.u. = 4.3597e-18 J
  real(long), parameter :: AU2EV  = AU_E_POTENTIAL   ! 1 a.u. = 27.2114 eV
  real(long), parameter :: AU2M   = AU_LENGTH        ! 1 a.u. = 5.29e-11 m
  real(long), parameter :: AU2ANG = AU_LENGTH * 1.0e10_long
                                                     ! 1 a.u. = 0.529 angstrom
  real(long), parameter :: AU2KG  = AU_MASS          ! 1 a.u. = 9.109e-31 kg
  real(long), parameter :: AU2SEC = AU_TIME          ! 1 a.u. = 2.4189e-17 sec
  real(long), parameter :: AU2ATTOSEC = AU_TIME * 1.0e18_long
                                                     ! 1 a.u. = 24.189 attosec.
  !
  ! a.u. reverse conversion
  !
  real(long), parameter :: J2AU   = ONE / AU2J           ! 1 J = 2.2937e17 a.u.
  real(long), parameter :: EV2AU  = ONE / AU2EV          ! 1 eV = 0.036749 a.u.
  real(long), parameter :: M2AU   = ONE / AU2M           ! 1 m = 1.8897e10 a.u.
  real(long), parameter :: ANG2AU = ONE / AU2ANG         ! 1 angs. = 1.8897 a.u.
  real(long), parameter :: KG2AU  = ONE / AU2KG          ! 1 kg = 1.0978e30 a.u.
  real(long), parameter :: SEC2AU = ONE / AU2SEC         ! 1 sec = 4.134e16 a.u.
  real(long), parameter :: ATTOSEC2AU = ONE / AU2ATTOSEC ! 1 asec = 0.04134 a.u.
  real(long), parameter :: AMU2AU = AMU_SI * KG2AU       ! 1 amu = 1822.9 a.u.
  !
  ! speed of light: c = 1 / sqrt( epsilon0 mu0 ) = 2.9979e8 m/s = 137.04 a.u.
  !   electric constant: epsilon0 = 8.854187817e-12 F / m (exact)
  !   magnetic constant: mu0 = 4PI * e-7 = 12.566370614e-7 N / A^2 (exact)
  !
  real(long), parameter :: SPEED_C_AU = SPEED_C_SI * M2AU / SEC2AU
  !
  ! physical constant conversion
  !
  real(long), parameter :: CAL2J = 4.184_long           ! 1 cal = 4.184 J
  real(long), parameter :: KCAL2J = CAL2J * 1.0e3_long  ! 1 kcal = 4184 J
  real(long), parameter :: EV2J = CHARGE_E_SI           ! 1 eV = 1.602e-19 J
  real(long), parameter :: EV2KCALMOL = EV2J / KCAL2J * AVOGADRO_NA_SI
                                                        ! 1 eV = 23.06 kcal/mol
  real(long), parameter :: HARTREE2J = AU2J             ! 1 a.u = 4.3597e-18 J
  real(long), parameter :: HARTREE2EV = AU2EV           ! 1 a.u = 27.2114 eV
  real(long), parameter :: HARTREE2KCALMOL = HARTREE2EV * EV2KCALMOL
                                                        ! 1 a.u = 627.5 kcal/mol
  real(long), parameter :: J2CAL     = ONE / CAL2J      ! 1 J = 0.239 cal
  real(long), parameter :: J2KCAL    = ONE / KCAL2J     ! 1 J = 2.39e-4 kcal
  real(long), parameter :: J2EV      = ONE / EV2J       ! 1 J = 6.24e18 eV
  real(long), parameter :: J2HARTREE = ONE / HARTREE2J  ! 1 J = 2.29e17 hartree
  !
  ! energy relations: E = h nu = h c sigma = h c / lambda = k T
  !   1 eV -> 2.419884e14 Hz -> 8065.5410 cm^-1 (1239.84244 nm) -> 11604.45 K
  !
  real(long), parameter :: EV2HZ = EV2J / PLANCK_H_SI   ! 1 eV = 2.419884e14 Hz
  real(long), parameter :: EV2WAVENUMBER = EV2HZ / SPEED_C_SI / 100.0_long
  real(long), parameter :: EV2KAYSER = EV2WAVENUMBER    ! 1 eV = 8065.5410 cm^-1
                                                        !     -> 1239.84244 nm
  real(long), parameter :: EV2K = EV2J / BOLTZMANN_K_SI ! 1 eV = 11604.45 K
contains
  function check_NaN(f) result(yes_or_no)
     real(long) :: f
     logical :: yes_or_no
!     yes_or_no = ieee_is_nan(f)
     yes_or_no = isnan(f)
  end function check_NaN
end module constants
