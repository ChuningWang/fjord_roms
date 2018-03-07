! This module contains variable for tracers.

MODULE mod_param_iceplume_tracers

    USE mod_param_iceplume
    implicit none

    integer :: pTracersNum = 1    ! number of passive tracers
    integer :: iTracer  ! looping vars

    ! Profiles
    real(r8), allocatable :: ptrPlume(:), ptrPlumeCum(:)
    real(r8), allocatable :: ptrIni(:)
    real(r8), allocatable :: ptrProf(:, :)

    ! ===================================================================
    CONTAINS

        SUBROUTINE allocate_param_iceplume_tracers

            ! Allocate the variables
            allocate(ptrPlume(pTracersNum),    &
                  &  ptrPlumeCum(pTracersNum), &
                  &  ptrIni(pTracersNum),      &
                  &  ptrProf(Nr, pTracersNum), &
                  & )

        END SUBROUTINE allocate_param_iceplume_tracers

    ! ===================================================================

END
