! This module stores all global variables.

MODULE mod_iceplume

    implicit none

    ! ====================================================================
    ! User defined variables and logical switches
    ! These parameters are read when coupling with GCMs.
    integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
    logical :: usePlumeDiagnostics = .false.
    logical :: conserveMass        = .false.
    logical :: correctMeltEnt      = .true.
    logical :: useConePlume        = .true.
    logical :: useSheetPlume       = .false.
    logical :: useTracers          = .true.
    logical :: useInputTracers     = .true.

    integer :: Dt           ! 3D time stepping interval [s]
    integer :: Nr           ! total number of layers
    real(r8) :: dy          ! grid cell resolution [m]
    real(r8) :: dx          ! grid cell resolution [m]
    real(r8) :: iceDepth    ! ice bottom depth [m], if iceDepth is positive,
                            ! then the ice bottom depth is always equal to
                            ! the water depth (iceDepthK = 1)

    ! Initial (discharge) conditions
    real(r8) :: wIni        ! initial vertical velocity [m s^-1]
    real(r8) :: rIni        ! initial radius [m]
    real(r8) :: TIni        ! initial temp [degC]
    real(r8) :: SIni        ! initial salt [PSU]

    ! Meltwater entrainment ratio
    ! This is a parameter to quantify how much background meltwater is
    ! entrained into the plume. This ratio varies between 0 and 1, where
    ! 0 means no entrainment into the plume and 1 means all entrainment
    real(r8) :: meltEnt     ! [ratio]

    ! ====================================================================

    ! Model parameters
    real(r8), parameter :: pi            = 4.0d0*atan(1.0d0)    ! Pi

    ! For plume model
    real(r8), parameter :: E_0           = 1.d-1     ! entrainment rate
    real(r8), parameter :: iceTemp       = -1.d1     ! ice temperature [degC]
    real(r8), parameter :: rho_ref       = 1.020d3   ! reference density [kg m^-3]
    real(r8), parameter :: rho_fresh     = 1.000d3   ! reference density [kg m^-3]
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

    integer :: iceDepthK                             ! ice bottom layer index [sigma points]
    integer :: plumeDepthK                           ! neutral buoyancy plume layer index [sigma points]
    real(r8) :: QIni                                 ! subglacial discharge [m^3 s^-1]

    ! Taken from pkg/icefront
    ! real(r8), parameter :: mass2rUnit = 1.d0/rho_ref  ! reciprocal density [kg^-1 m^3]

    ! Looping vars
    integer :: K  ! vertical layers

    ! ====================================================================
    ! Profiles
    ! ====================================================================
    ! AT RHO POINTS
    ! sProf         - ambient salinity [PSU]
    ! tProf         - ambient temperature [degC]
    ! vProf         - horizontal velocity parallel to the glacier wall [m s^-1]
    ! wProf         - vertical velocity parallel to the glacier wall [m s^-1]
    ! rhoProf       - ambient density [kg m^-3]
    ! zProfAbs      - depth (absolute value) [m]
    ! prProf        - pressure [dbar]
    ! mProfAv       - total melt rate [m s^-1]
    ! mProfPlume    - plume melt rate [m s^-1]
    ! mProf         - background melt rate [m s^-1]
    ! fwFlux        - meltwater freshwater flux into grid cell [kg m^-2 s^-1]
    ! heatFlux      - meltwater heat flux into grid cell [W m^-2]
    ! tendT         - a tendency term calculated as in MITgcm.
    ! tendS         - a tendency term calculated as in MITgcm.
    ! ptrProf       - passive tracer concentration
    ! delta_z       - veritcal resolution [m]
    ! volFluxDiff   - entrainment rate [m^3 s^-1]
    !
    ! AT OMEGA POINTS
    ! zProf - depth [m]
    ! rProfPlume    - radius of plume [m]
    ! wProfPlume    - vertical velocity of plume [m s^-1]
    ! tProfPlume    - temperature of plume [degC]
    ! sProfPlume    - salinity of plume [PSU]
    ! aProfPlume    - area of plume [m^2]
    ! mIntProfPlume - area integrated melt of plume [m^3 s^-1]
    ! rhoProfPlume  - density of plume [kg m^-3]
    ! volFLux       - vertical volume flux [m^3 s^-1]
    ! ====================================================================
    ! Input arguments
    real(r8), dimension(:), allocatable :: zProf, sProf, tProf, vProf, wProf
    ! Output arguments - from iceplume_plume_model
    real(r8), dimension(:), allocatable :: rProfPlume, tProfPlume, sProfPlume,    &
                                         & wProfPlume, aProfPlume, mIntProfPlume
    ! Output arguments - from iceplume_calc
    real(r8), dimension(:), allocatable :: mProfAv, mProfPlume, mProf,    &
                                         & volFluxDiff, fwFlux, heatFlux, &
                                         & tendT, tendS

    ! Internal use only
    real(r8), dimension(:), allocatable :: zProfAbs, prProf, dz, volFLux
    real(r8), dimension(:), allocatable :: TbProf, SbProf

    ! ===================================================================
    CONTAINS

        SUBROUTINE allocate_iceplume

            ! Allocate profiles
            allocate( &
                  &  zProf(Nr+1), sProf(Nr), tProf(Nr), vProf(Nr), wProf(Nr), &
                  & )
            allocate( &
                  &  rProfPlume(Nr+1), tProfPlume(Nr+1), sProfPlume(Nr+1),    &
                  &  wProfPlume(Nr+1), aProfPlume(Nr+1), mIntProfPlume(Nr+1), &
                  & )
            allocate( &
                  &  mProfAv(Nr), mProfPlume(Nr), mProf(Nr),                  &
                  &  volFLuxDiff(Nr), fwFlux(Nr), heatFlux(Nr),               &
                  &  tendT(Nr), tendS(Nr),                                    &
                  & )
            allocate( &
                  &  zProfAbs(Nr), prProf(Nr), dz(Nr), volFlux(Nr+1),         &
                  &  TbProf(Nr), SbProf(Nr),                                  &
                  & )

        END

END
