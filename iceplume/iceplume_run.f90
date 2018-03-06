! Test the icefront model.

PROGRAM iceplume

    USE mod_param_iceplume
    implicit none

    real(r8), dimension(Nr+1) :: zprofw, sprofw, tprofw, vprofw, wprofw
    ! looping integer
    integer :: K
    character(100) :: fileName

    fileName = 'iceplume_output.txt'

    print *, "IcePlume test run."

    ! ==================================================================
    ! read in profiles for testing
    open(unit=15, file='iceplume_test_input.txt', action='read')
    100 format(5 F12.6)
    DO K = 1, Nr+1
        read(15, 100)  zprofw(K), sprofw(K), tprofw(K), vprofw(K), wprofw(K)
    ENDDO
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

    CALL iceplume_calc()
    CALL iceplume_write(fileName)
    ! CALL iceplume_plume_model()

END PROGRAM iceplume
