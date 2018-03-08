! This module contains variable for tracers.

MODULE mod_iceplume_tracers

    USE mod_iceplume
    implicit none

    integer :: pTracersNum    ! number of passive tracers
    integer :: iTracer        ! looping vars

    ! Profiles
    ! Input argument
    real(r8), allocatable :: ptrProf(:, :)
    real(r8), allocatable :: ptrIni(:)
    ! Output argument
    real(r8), allocatable :: ptrPlume(:), ptrPlumeCum(:)

    ! ===================================================================
    CONTAINS

        SUBROUTINE allocate_iceplume_tracers

            ! Allocate the variables
            allocate(ptrProf(Nr, pTracersNum), &
                  &  ptrIni(pTracersNum),      &
                  &  ptrPlume(pTracersNum),    &
                  &  ptrPlumeCum(pTracersNum), &
                  & )

        END

    ! ===================================================================

END
