MODULE mod_param_iceplume

    implicit none

    ! local parameters
    integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
    real(r8), parameter :: pi            = 4.0d0*atan(1.0d0)
    real(r8), parameter :: E_0           = 1.d-1
    real(r8), parameter :: iceTemp       = -1.d1
    real(r8), parameter :: rho_ref       = 1.020d3
    real(r8), parameter :: g             = 9.81d0
    real(r8), parameter :: c_w           = 3.974d3
    real(r8), parameter :: c_i           = 2.009d3
    real(r8), parameter :: L             = 3.35d5
    real(r8), parameter :: lambda1       = -5.73d-2
    real(r8), parameter :: lambda2       = 8.32d-2
    real(r8), parameter :: lambda3       = 7.61d-4

    real(r8), parameter :: GamT          = 2.20d-2
    real(r8), parameter :: GamS          = 6.20d-4
    real(r8), parameter :: Cd            = 2.50d-3
    real(r8), parameter :: backgroundVel = 1.d-2

    real(r8), parameter :: T_sg          = 4.0d0
    real(r8), parameter :: S_sg          = 1.0d1
    ! real(r8), parameter :: T_sg          = 1.0d-3
    ! real(r8), parameter :: S_sg          = 1.0d-3
    real(r8), parameter :: w_sg          = 1.0d0
    real(r8), parameter :: r_sg          = 1.0d0
    ! real(r8), parameter :: theta_sg      = 2.0d0*atan(1.0d0)

    real(r8), parameter :: RTOL          = 1.0d-5
    real(r8), parameter :: ATOL          = 1.0d-5

    logical, parameter :: applyIcefrontTendT  = .true.
    logical, parameter :: applyIcefrontTendS  = .true.
    logical, parameter :: usePlumeDiagnostics = .false.
    logical, parameter :: conserveMass        = .false.
    logical, parameter :: useSheetPlume       = .false.
    logical, parameter :: useConePlume        = .true.

    ! profiles
    integer, parameter :: Nr = 40  ! total number of layers

    real(r8), dimension(Nr+1) :: zprofw, sprofw, tprofw
    real(r8), dimension(Nr) :: zprofr, sprofr, tprofr
    real(r8), dimension(Nr+1) :: zProf, zProfAbs, rProfPlume, wProfPlume, &
                               & tProfPlume, sProfPlume, aProfPlume, mIntProfPlume

END
