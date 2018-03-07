! This is the wrapper to combine iceplume_calc, iceplume_plume_model and
! icePlume_write. Use this subroutine to pass in input profiles.

SUBROUTINE iceplume_wrapper( &
                         &  NIn, dxIn, dyIn, iceDepthIn,                    &
                         &  wIniIn, rIniIn, TIniIn, SIniIn,                 &
                         &  zProfIn, sProfIn, tProfIn, vProfIn, wProfIn,    &
                         &  rProfPlumeOut, tProfPlumeOut, sProfPlumeOut,    &
                         &  wProfPlumeOut, aProfPlumeOut, mIntProfPlumeOut, &
                         &  mProfAvOut, mProfPlumeOut, mProfOut,            &
                         &  volFluxDiffOut, fwFluxOut, heatFluxOut)

    USE mod_param_iceplume
    USE mod_param_iceplume_tracers
    implicit none
    
    integer, intent(in) :: NIn
    real(r8), intent(in) :: &
                         &  dxIn, dyIn, iceDepthIn,                         &
                         &  wIniIn, rIniIn, TIniIn, SIniIn
    real(r8), intent(in), dimension(:), pointer :: &
                         &  zProfIn, sProfIn, tProfIn, vProfIn, wProfIn
    real(r8), intent(out), dimension(:), pointer :: &
                         &  rProfPlumeOut, tProfPlumeOut, sProfPlumeOut,    &
                         &  wProfPlumeOut, aProfPlumeOut, mIntProfPlumeOut, &
                         &  mProfAvOut, mProfPlumeOut, mProfOut,            &
                         &  volFluxDiffOut, fwFluxOut, heatFluxOut

    ! Allocate variables

END SUBROUTINE iceplume_wrapper
