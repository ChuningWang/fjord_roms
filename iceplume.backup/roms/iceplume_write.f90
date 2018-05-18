! This program write the computed plume info to a file.
SUBROUTINE ICEPLUME_WRITE(ng)

    USE mod_iceplume

    implicit none
    integer :: ng
    
    open(unit=15, file='./zw_out.txt', action='write')
    200 format(7 E12.4)
    write(15, '(99 A12)') 'zProf', 'rProfP', 'wProfP', 'tProfP', 'sProfP', 'aProfP', 'mIntProfP'
    DO K = 1, Nr+1
        write(15, 200) abs(PLUME(ng) % zW(K)), PLUME(ng) % r(K), &
                     & PLUME(ng) % w(K), PLUME(ng) % t(K), &
                     & PLUME(ng) % s(K), PLUME(ng) % a(K), PLUME(ng) % mInt(K)
    ENDDO
    close(unit=15)

    open(unit=15, file='./zr_out2.txt', action='write')
    210 format(9 E12.4)
    write(15, '(99 A12)') 'zProf', 'mProfAv', 'mProfPlume', 'mProf', &
                        & 'fwFlux', 'heatFlux', 'Entrainment', 'tendT', 'tendS'
    DO K = 1, Nr
        write(15, 210) PLUME(ng) % zRhoAbs(K), PLUME(ng) % mAv(K), &
                     & PLUME(ng) % m(K), PLUME(ng) % mAm(K), &
                     & PLUME(ng) % fwFlux(K), PLUME(ng) % heatFlux(K), &
                     & PLUME(ng) % ent(K), &
                     & PLUME(ng) % tendT(K), PLUME(ng) % tendS(K)
    ENDDO
    close(unit=15)

END
