! This module contains variable for tracers.

MODULE mod_param_iceplume_tracers

    USE mod_param_iceplume
    implicit none

    integer :: pTracersNum = 1    ! number of passive tracers
    integer :: iTracer  ! looping vars

    ! Profiles
    ! Input argument
    real(r8), allocatable :: ptrProf(:, :)
    real(r8), allocatable :: ptrIni(:)
    ! Output argument
    real(r8), allocatable :: ptrPlume(:), ptrPlumeCum(:)

    ! ===================================================================
    CONTAINS

        SUBROUTINE allocate_param_iceplume_tracers

            ! Allocate the variables
            allocate(ptrProf(Nr, pTracersNum), &
                  &  ptrIni(pTracersNum),      &
                  &  ptrPlume(pTracersNum),    &
                  &  ptrPlumeCum(pTracersNum), &
                  & )

        END SUBROUTINE allocate_param_iceplume_tracers

    ! ===================================================================

END
