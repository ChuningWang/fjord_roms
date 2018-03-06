! This module stores all global variables.

MODULE mod_param_iceplume

    implicit none

    ! ====================================================================
    ! User defined variables and logical switches
    ! These parameters are read when coupling with GCMs.
    integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
    ! logical :: applyIcefrontTendT  = .true.
    ! logical :: applyIcefrontTendS  = .true.
    ! logical :: usePlumeDiagnostics = .false.
    logical :: conserveMass        = .false.
    logical :: useConePlume        = .true.
    logical :: useSheetPlume       = .false.
    integer, parameter :: Nr = 40            ! total number of layers
    real(r8), parameter :: delta_y = 200.d0  ! grid cell resolution [m]
    real(r8), parameter :: delta_x = 200.d0  ! grid cell resolution [m]
    ! ====================================================================

    ! Model parameters
    real(r8), parameter :: pi = 4.0d0*atan(1.0d0)  ! Pi

    ! For plume model
    real(r8), parameter :: E_0           = 1.d-1     ! entrainment rate
    real(r8), parameter :: iceTemp       = -1.d1     ! ice temperature [degC]
    real(r8), parameter :: rho_ref       = 1.020d3   ! reference density [kg m^-3]
    real(r8), parameter :: rho_ice       = 0.9167d3  ! ice density [kg m^-3]
    real(r8), parameter :: g             = 9.81d0    ! gravity acceleration [m s^-2]
    real(r8), parameter :: c_w           = 3.974d3   ! heat capacity of water [J kg^-1 degC^-1]
    real(r8), parameter :: c_i           = 2.009d3   ! heat capacity of ice [J kg^-1 degC^-1]
    real(r8), parameter :: L             = 3.35d5    ! latent heat of melting [J kg^-1]
    real(r8), parameter :: lambda1       = -5.73d-2  ! freezing point slope [degC PSU^-1]
    real(r8), parameter :: lambda2       = 8.32d-2   ! freezing point offset [degC]
    real(r8), parameter :: lambda3       = 7.61d-4   ! freezing point depth slope [degC m^-1]

    real(r8), parameter :: GamT          = 2.20d-2   ! thermal turbulent transfer coefficient
    real(r8), parameter :: GamS          = 6.20d-4   ! salt turbulent transfer coefficient
    real(r8), parameter :: Cd            = 2.50d-3   ! ice-plume drag coefficient
    real(r8), parameter :: backgroundVel = 1.d-2     ! background velocity [m s^-1]

    real(r8), parameter :: T_sg          = 1.0d-3    ! initial temp [degC]
    real(r8), parameter :: S_sg          = 1.0d-3    ! initial salt [PSU]
    real(r8), parameter :: w_sg          = 1.0d0     ! initial vertical velocity [m s^-1]
    real(r8), parameter :: r_sg          = 1.0d2     ! initial radius [m]
    real(r8) :: Q_sg                                 ! subglacial discharge [m^3 s^-1]

    integer :: iceDepthK  ! ice bottom layer index
    real(r8) :: iceDepth  ! ice bottom depth

    ! For iceplume_calc
    real(r8), dimension(Nr) :: delta_z       ! vertical resolution [m]
    real(r8), dimension(Nr+1) :: volFlux     ! volume flux [m^3 s^-1]
    real(r8), dimension(Nr) :: volFluxDiff   ! entrainment [m^3 s^-1]

    ! Taken from pkg/icefront
    real(r8), parameter :: mass2rUnit      = 1.d0/rho_ref  ! reciprocal density [kg^-1 m^3]

    ! Profiles
    ! ====================================================================
    ! AT RHO POINTS
    ! sProf - ambient salinity [PSU]
    ! tProf - ambient temperature [degC]
    ! vProf - horizontal velocity parallel to the glacier wall [m s^-1]
    ! wProf - vertical velocity parallel to the glacier wall [m s^-1]
    ! zProfAbs - depth (absolute value) [m]
    ! prProf - pressure [dbar]
    ! mProfAv - total melt rate [m s^-1]
    ! mProfPlume - plume melt rate [m s^-1]
    ! mProf - background melt rate [m s^-1]
    ! fwFlux - freshwater flux [kg m^-2 s^-1]
    ! heatFlux - heat flux [W m^-2]
    ! icefront_TendT - a tendency term calculated in MITgcm.
    ! icefront_TendS - a tendency term calculated in MITgcm.
    !
    ! AT OMEGA POINTS
    ! zProf - depth [m]
    ! rProfPlume - radius of plume [m]
    ! wProfPlume - vertical velocity of plume [m s^-1]
    ! tProfPlume - temperature of plume [degC]
    ! sProfPlume - salinity of plume [PSU]
    ! aProfPlume - area of plume [m^2]
    ! mIntProfPlume - area integrated melt of plume [m^3 s^-1]
    ! ====================================================================
    real(r8), dimension(Nr) :: sProf, tProf, vProf, wProf, zProfAbs, &
                             & prProf, mProfAv, mProfPlume, mProf
    real(r8), dimension(Nr) :: fwFlux, heatFlux, icefront_TendT, icefront_TendS
    real(r8), dimension(Nr+1) :: zProf, rProfPlume, wProfPlume, &
                               & tProfPlume, sProfPlume, aProfPlume, mIntProfPlume

END
