! Test the icefront model.

PROGRAM iceplume

    USE mod_param_iceplume
    USE mod_param_iceplume_tracers

    implicit none

    character(100) :: fileName
    real(r8), dimension(:), allocatable :: zprofw, sprofw, tprofw, vprofw, wprofw
    real(r8), allocatable :: ptrprofw(:, :)
    allocate(zprofw(Nr+1), sprofw(Nr+1), tprofw(Nr+1), vprofw(Nr+1), wprofw(Nr+1))
    allocate(ptrprofw(Nr+1, pTracersNum))
    CALL allocate_param_iceplume
    CALL allocate_param_iceplume_tracers

    fileName = 'iceplume_output.txt'

    print *, "IcePlume test run."

    ! ==================================================================
    ! read in profiles for testing
    open(unit=15, file='iceplume_test_input.txt', action='read')
    100 format(5 F12.6)
    DO K = 1, Nr+1
        read(15, 100)  zprofw(K), sprofw(K), tprofw(K), vprofw(K), wprofw(K)
    ENDDO
    close(unit=15)

    DO K = 1, Nr+1
        zProf(K) = zprofw(K)
    ENDDO
    ! get T and S at Rho points
    DO K = 1, Nr
        ! zProf(K) = 0.5*(zprofw(K) + zprofw(K+1))
        sProf(K) = 0.5*(sprofw(K) + sprofw(K+1))
        tProf(K) = 0.5*(tprofw(K) + tprofw(K+1))
        vProf(K) = 0.5*(vprofw(K) + vprofw(K+1))
        wProf(K) = 0.5*(wprofw(K) + wprofw(K+1))
    ENDDO

    IF (useTracers) THEN
        ! Read in tracer concentration
        open(unit=15, file='iceplume_test_input_tracers.txt', action='read')
        200 format(999 F12.6)
        DO K = 1, Nr+1
            read(15, 200)  ptrprofw(K, :)
            ! write(*, '(10 E20.10)')  ptrprofw(K, :)
        ENDDO
        close(unit=15)
    ENDIF

    DO K = 1, Nr
        DO iTracer = 1, pTracersNum
            ptrProf(K, iTracer) = 0.5*(ptrprofw(K, iTracer) + ptrprofw(K+1, iTracer))
        ENDDO
    ENDDO

    CALL iceplume_calc()
    CALL iceplume_write(fileName)
    ! CALL iceplume_plume_model()

END PROGRAM iceplume
