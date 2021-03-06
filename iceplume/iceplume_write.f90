! This program write the computed plume info to a file.
SUBROUTINE ICEPLUME_WRITE(fileName)

    USE mod_iceplume

    implicit none
    character(*), intent(in) :: fileName
    
    open(unit=15, file='./zw_' // trim(fileName), action='write')
    200 format(7 E12.4)
    write(15, '(99 A12)') 'zProf', 'rProfP', 'wProfP', 'tProfP', 'sProfP', 'aProfP', 'mIntProfP'
    DO K = 1, Nr+1
        write(15, 200) abs(zProf(K)), rProfPlume(K), wProfPlume(K), tProfPlume(K), &
                     & sProfPlume(K), aProfPlume(K), mIntProfPlume(K)
    ENDDO
    close(unit=15)

    open(unit=15, file='./zr_' // trim(fileName), action='write')
    210 format(9 E12.4)
    write(15, '(99 A12)') 'zProf', 'mProfAv', 'mProfPlume', 'mProf', &
                        & 'fwFlux', 'heatFlux', 'Entrainment', 'tendT', 'tendS'
    DO K = 1, Nr
        write(15, 210) zProfAbs(K), mProfAv(K), mProfPlume(K), mProf(K), &
                     & fwFlux(K), heatFlux(K), volFluxDiff(K), tendT(K), tendS(K)
                     ! & fwFlux(K), heatFlux(K), volFluxDiff(K), icefront_TendT(K), icefront_TendS(K)
    ENDDO
    close(unit=15)

END
