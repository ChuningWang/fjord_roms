SUBROUTINE ICEPLUME_CALC()

    USE mod_param_iceplume
    implicit none

    ! ==================================================================
    ! LOCAL VARIABLES:
    ! K :: loop indices
    ! sProf, tProf, uProf, vProf :: salt, pot. temperature and
    !                        uVel and vVel profiles
    ! ptrProf  :: ambient ptracer profile
    ! ptrPlumeCum :: cumulative quantity of ptracer in plume
    ! ptrPlume :: quantity of ptracer in plume outflow
    ! eps5     :: for thermodynamics (see pkg icefront)
    ! maxDepth :: vertical extent of domain (m)
    ! plumeAreaInCell :: surface area of plume in contact with ice in that cell (m^2)
    ! negSum, posSum :: sum of negative and positive contributions to the plume volume
    ! posNegRatio    :: ratio of the above
    ! meanVel :: ice tangental velocity
    ! rho_0 :: average density of seawater
    ! delta_y :: grid cell resolution
    ! ==================================================================

    integer :: K
    real(r8) :: plumeAreaInCell, maxDepth
    real(r8) :: negSum, posSum, posNegRatio
    real(r8) :: meanVel, depth
    real(r8) :: sw_temp, sw_ptmp
    ! real(r8), parameter :: rho_0 = 1027.d0

    ! Read in the subglacial discharge for this cell
    IF (useSheetPlume) THEN
        Q_sg = w_sg*r_sg*delta_y
    ELSEIF (useConePlume) THEN
        Q_sg = w_sg*((pi*r_sg**2.)/2.)
    ELSE
        Q_sg = 0.d0
    ENDIF

    ! Create variables with temperature, salinity
    ! and velocity profiles for that column

    DO K = 1,Nr
        ! Tracers
        prProf(K) = zProfAbs(K)*rho_ref*g*1.0E-6  ! Pressure (dbar)
        ! sProf(K)  = salt(I,J,K,bi,bj)   ! Salinity
        ! ptProf(K) = theta(I,J,K,bi,bj)  ! Potential Temperature
        ! tProf(K)  = SW_TEMP(sProf(k),ptProf(k),prProf(k),0. _d 0) ! Temperature
        delta_z(K) = zProf(K+1) - zProf(K)
    ENDDO

    ! ==================================================================
    ! Set iceDepth to the bottom layer
    iceDepth = zProf(1)
    DO K = 1,Nr+1
        IF (zProf(K) .EQ. iceDepth) iceDepthK = K
    ENDDO

    ! ==================================================================
    IF (Q_sg .GT. 0) THEN
        ! Run the plume model
        CALL ICEPLUME_PLUME_MODEL()

        ! Calculate vertical plume volume flux
        DO k = 1, Nr
            ! After checking to see if we are above the base of the ice face...
            IF (K .GT. iceDepthK) THEN
                ! assuming specified plume horizontal extent (for sheet flow)...
                IF (useSheetPlume) THEN
                    volFlux(K) = rProfPlume(K)*wProfPlume(K)*delta_y
                ELSEIF (useConePlume) THEN
                    volFlux(K) = pi*(rProfPlume(K)**2)*wProfPlume(K)/2.
                ELSE
                    volFlux(K) = 0
                ENDIF
            ELSE
                volFlux(K) = 0.d0
            ENDIF
        ENDDO

        ! A couple of corrections:
        ! Even if plume is still buoyant, it cannot flow through the fjord surface
        volFlux(Nr+1) = 0.d0
        ! The initial volume flux is equal to runoff
        volFlux(iceDepthK) = Q_sg

        ! Calculate volume flux differential to give entrainment / extrainment
        ! First clear volfluxdiff

        DO K = 1,Nr
            volFluxDiff(K) = 0.d0
        ENDDO

        DO K = iceDepthK, Nr
            volFluxDiff(K) = volFlux(K+1) - volFlux(K)
            ! write(*, '(E20.10)') volFluxDiff(K)
        ENDDO

        IF (conserveMass) THEN
            ! Scale output to compensate for entrainment lost in expanding of output layer
            ! i.e. so that there is no net flow over boundary

            negSum = 0.D0
            posSum = 0.D0

            DO K = 1,Nr
                IF ( volFluxDiff(K) .LT. 0 ) THEN
                    negSum = volFluxDiff(K) + negSum
                ELSE
                    posSum = volFluxDiff(K) + posSum
                ENDIF
            ENDDO

            IF ( negSum .NE. 0 ) THEN
                posNegRatio = -posSum / negSum
                DO K = 1,Nr
                    IF ( volFluxDiff(K) .LT. 0 ) &
                        & volFluxDiff(K) = volFluxDiff(K) * posNegRatio
                ENDDO
            ENDIF
        ENDIF
    ELSE  ! (Q_sg .EQ. 0)
        ! If no subglacial output, then there is no plume
        DO k = 1,Nr
            rProfPlume(K) = 0.d0
            wProfPlume(K) = 0.d0
            tProfPlume(K) = 0.d0
            sProfPlume(K) = 0.d0
            aProfPlume(K) = 0.d0
            mIntProfPlume(K) = 0.d0
        ENDDO
    ENDIF

    ! ==================================================================
    ! Claculate melt rates
    DO K = 1, Nr

        ! Check if we are above sea bed
        IF (K .LE. iceDepthK) THEN
            ! If not then there is no melting
            mProfAv(K) = 0.d0
            mProfPlume(K) = 0.d0
            mProf(K) = 0.d0
        ELSE
            ! If there is a plume in that cell, then need to calculate plume melt rate
            ! distinct to background melt rate. Plume melt rate is already encorporated in 
            ! the plrume model, and taken into account in the temperature and salinity of the
            ! plume outflow. It is useful though to have it available as a diagnostic.
            plumeAreaInCell = 0.d0
            IF ((Q_sg .NE. 0) .AND. (useConePlume .OR. useSheetPlume)) THEN
                plumeAreaInCell = aProfPlume(K+1) - aProfPlume(K)

                IF (plumeAreaInCell .GT. 0.0) THEN
                    mProfPlume(K) =(mIntProfPlume(K)-mIntProfPlume(K+1))/ &
                        & plumeAreaInCell
                ELSE
                    mProfPlume (K) = 0.0
                ENDIF               
            ELSE
                ! If no plume in that cell set plume melt rate to zero
                mProfPlume(K) = 0.d0
            ENDIF

            ! ==================================================================
            ! Calculate the background melt rate (i.e. not generated by plumes). This will then be 
            ! used to update the temperature and salinity in the adjacent cells.
            ! Velocities are calculated at cell faces - find averages for cell centres.
            ! Does not include velocity perpendicular to ice - this differs depending on 
            ! orientation of ice front
            meanVel = ((vProf(K))**2. + (wProf(K))**2.)**0.5
            depth = 0.5d0*(zProf(K)+zProf(K+1))
            CALL ICEPLUME_MELTRATE(tProf(K), sProf(K), meanVel, depth, &
                                 & mProf(K))
            ! write(*, '(E20.10)') mProf(K)

            ! Get average melt rate. This is useful for visualizing melt patterns and
            ! assessing overall melt rate of glacier.
            ! The following should apply to both conical and sheet plume models
            IF (Q_sg .NE. 0) THEN
                plumeAreaInCell = aProfPlume(K+1) - aProfPlume(K)
                IF (plumeAreaInCell .LE. delta_y*delta_z(K)) THEN
                    IF (plumeAreaInCell .LE. 0) THEN
                        ! If there is no plume in cell, then the melt rate is
                        ! equal to the background melt rate.
                        mProfAv(K) = mProf(K)
                    ELSE
                        ! If there is a plume in cell, calculate average melt rate
                        mProfAv(K) = (mProfPlume(K)*plumeAreaInCell &
                            & +mProf(K)*(delta_y*delta_z(K)-plumeAreaInCell)) / &
                            & (delta_y * delta_z(K))

                        ! Scale down background melt rate to account for area occupied by plume
                        ! (necessary so that tendency terms aren't over estimated)

                        mProf(K) = mProf(K)*(1-plumeAreaInCell / &
                            & (delta_y*delta_z(K)))
                    ENDIF
                ELSE  ! plumeAreaInCell .GE. delta_y*delta_z(K)
                    ! If the plume contact area is larger than the cell area, we
                    ! assume there is no background melting
                    mProfAv(K) = mProfPlume(K)*plumeAreaInCell / &
                        & (delta_y*delta_z(K))
                    mProf(K) = 0.d0
                ENDIF
            ELSE  ! not coneplume or sheet plume
                ! If it is not a plume cell, then no plume melting.
                mProfPlume(K) = 0.d0
                mProfAv(K) = mProf(K)
            ENDIF  ! plume type
        ENDIF  ! above or below sea bed

        ! ==================================================================
        ! Tendencies
        ! These are applied in cells where there is no plume.
        ! The idea is that in these areas there are likely to be local subgrid convection cells.
        ! As such, it is most realistic to apply changes in T and S to ice edge cell.
        ! Where there is a plume, products of melt are transported upwards so no local changes applied.
        ! The thermodynamics in this section are taken from pkg/icefront in MITgcm.
        ! Tendencies are applied in S/R EXTERNAL_FORCING

        ! To convert from melt rate (m s^-1) to freshwater flux (kg m^-2 s^-1)
        fwFlux(K) = -mProf(K)*rho_ice

        ! Heat required to melt that much ice (W m^-2)
        heatFlux(K) = -fwFlux(K)*L

        ! Compute tendencies (as for pkg/icefront in MITgcm)
        icefront_TendT(K) = -heatFlux(K)*mass2rUnit/c_i
        icefront_TendS(K) = fwFlux(K)*mass2rUnit*sProf(K)

        ! Scale by icefrontlength, which is the ratio of the horizontal length
        ! of the ice front in each model grid cell divided by the grid cell area.
        ! (icefrontlength = dy / dxdy = 1 / dx)
        icefront_TendT(K) = icefront_TendT(K)/delta_x
        icefront_TendS(K) = icefront_TendS(K)/delta_x

    ENDDO  ! K = 1, Nr

END SUBROUTINE ICEPLUME_CALC
