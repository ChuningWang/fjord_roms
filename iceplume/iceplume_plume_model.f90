! Program to calculate the plume shape, S, T, V etc.
SUBROUTINE iceplume_plume_model

    USE mod_param_iceplume
    implicit none

    real(r8), dimension(Nr) :: proftemp

    ! ==================================================================
    ! Local variables for ODEPACK
    integer :: IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW, MF, NEQ
    real(r8), parameter :: RTOL = 1.0d-5
    real(r8), parameter :: ATOL = 1.0d-5
    real(r8) :: T, TOUT
    real(r8) :: RWORK(116), Y(6)


    ! Y is input/output vector for DLSODE
    !   Y(1) = plume thickness/radius
    !   Y(2) = plume velocity
    !   Y(3) = plume temperature
    !   Y(4) = plume salinity
    !   Y(5) = plume area
    !   Y(6) = area integrated melt

    ! Other local variables
    real(r8) :: RHO, temperature, salinity, depth
    real(r8) :: tambient, sambient
    real(r8) :: rho_plume, rho_ambient

    ! Plume model
    external JENKINS, HALFCONE, JEX

    ! ==================================================================
    ! For ODEPACK solver. See ODEPACK documentation and source code in Cowton et al. 2015.
    NEQ = 6
    LRW = 116
    LIW = 116

    ITOL = 1
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    MF = 10
    IWORK(7) = 2  ! To limit number of times repeat error messages are printed

    ! ==================================================================
    ! Initial conditions
    Y(1) = rIni  ! initial pume thickness (sheet model) or radius (halfcone model)
    Y(2) = wIni  ! initial vertical velocity
    Y(3) = TIni  ! initial temperature
    Y(4) = SIni  ! initial salinity
    Y(5) = 0.0   ! integrated contact area
    Y(6) = 0.0   ! integrated melt rate

    ! Prepare profiles
    DO K = 1, Nr+1
        rProfPlume(K) = 0.0
        wProfPlume(K) = 0.0
        tProfPlume(K) = 0.0
        sProfPlume(K) = 0.0
        aProfPlume(K) = 0.0
        mIntProfPlume(K) = 0.0
    ENDDO

    DO K = 1, Nr
        zProfAbs(K) = abs(0.5d0*(zProf(K) + zProf(K+1)))
    ENDDO

    ! Start at bottom of ice face
    T = zProf(iceDepthK)

    ! Next point at which to retrieve values
    TOUT = zProf(iceDepthK+1)

    ! Set initial conditions
    rProfPlume(iceDepthK) = Y(1)
    wProfPlume(iceDepthK) = Y(2)
    tProfPlume(iceDepthK) = Y(3)
    sProfPlume(iceDepthK) = Y(4)
    aProfPlume(iceDepthK) = Y(5)
    mIntProfPlume(iceDepthK) = Y(6)

    ! clean up some variables
    plumeDepthK = 0

    ! Move up through water column from lowest layer
    DO IOUT = iceDepthK+1, Nr+1
        ! Check to make sure plume hasn't reached neutral buoyancy in a lower layer
        IF (ISTATE .GT. -1) THEN

            IF (useSheetPlume) THEN
                CALL DLSODE (JENKINS, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                           & ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
            ELSEIF (useConePlume) THEN
                CALL DLSODE (HALFCONE, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                           & ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
            ENDIF

            ! Test to see if neutral buoyancy has now been reached.
            ! If solver returns ISTATE = -1, then it has been unable to meet required tolerances
            ! at this level. This generally occurs because plume has reached neutral buoyancy and
            ! run out of momentum, and so is no longer rising. At this point, we therefore end
            ! the call to the plume model.
            ! Our aim is to catch the plume at the point of neutral buoyancy. We therefore perform
            ! a manual comparrison of ambient and plume density. If plume density >= ambient density
            ! we assign ISTATE = -1, again ending the call to the plume model.

            temperature = Y(3)
            salinity = Y(4)
            depth = T

            ! Calculate plume density
            rho_plume = RHO(temperature, salinity, depth)

            ! Calculate ambient density
            IF (IOUT .EQ. Nr+1) THEN
                tambient = tProf(Nr)
                sambient = sProf(Nr)
            ELSE
                tambient = .5*(tProf(IOUT-1)+tProf(IOUT))
                sambient = .5*(sProf(IOUT-1)+sProf(IOUT))
            ENDIF
            rho_ambient = RHO(tambient, sambient, depth)

            ! write(*, '(i3, f12.6, f12.6)') ISTATE, rho_plume, rho_ambient
            IF (rho_plume .GT. rho_ambient) THEN
                ISTATE = -1
                plumeDepthK = K
                ! write(*, '(i20)') 'Reached ISTATE=-1!'
            ENDIF
            ! If ISTATE is now < 0, then plume has reached neutral buoyancy 
            IF (ISTATE .LT. 0) THEN
                ! If we have reached neutral buoyancy then there is no volume flux out of this cell,
                ! so plume area and velocity equal zero.
                ! Other values are kept for use in determining plume outflow properties.
                Y(1) = 0.d0
                Y(2) = 0.d0
            ELSE
                ! If the plume has not reached neutral buoyancy,
                ! then we assign a depth at which to calculate the next 
                ! value and loop round to call the plume model again.
                ! Make sure we're not at the surface
                IF (IOUT .NE. Nr+1) THEN
                    ! define present depth
                    T = TOUT
                    ! define next depth
                    TOUT = zProf(IOUT+1)
                ENDIF
            ENDIF

        ELSE
            ! This section is entered once the plume has reached neutral buoyancy
            ! once plume has reached neutral buoyancy, no plume values
            Y(1) = 0.0
            Y(2) = 0.0
            Y(3) = 0.0
            Y(4) = 0.0
            Y(5) = 0.0
            Y(6) = 0.0
        ENDIF

        ! Save results
        rProfPlume(IOUT) = Y(1)
        wProfPlume(IOUT) = Y(2)
        tProfPlume(IOUT) = Y(3)
        sProfPlume(IOUT) = Y(4)
        aProfPlume(IOUT) = Y(5)
        mIntProfPlume(IOUT) = Y(6)

    ENDDO

    ! For diagnostic purpose
    ! Write the plume profiles to file
    ! open(unit=16, file='iceplume_test_output.txt', action='write')
    ! 200 format(6 E20.10)
    ! DO K = 1, Nr+1
    !     write(16, 200)  rProfPlume(K), wProfPlume(K), tProfPlume(K), &
    !                   & sProfPlume(K), aProfPlume(K), mIntProfPlume(K)
END

! =========================================================================
! These subroutines are borrowed from Dr. Tom Cowton.
! =========================================================================

SUBROUTINE  HALFCONE (NEQ, T, Y, YDOT)

    USE mod_param_iceplume
    implicit none

    integer :: NEQ
    real(r8) :: T, Y(6), YDOT(6)
    real(r8) :: Tambient, Sambient, rho_0, rho_1
    real(r8) :: mdot, Sb, Tb
    real(r8) :: a, b, c
    real(r8) :: RHO

    ! Interpolate from imposed ambient profiles
    IF (abs(T) .LE. zProfAbs(1)) THEN
        Tambient = tProf(1)
        Sambient = sProf(1)
    ELSEIF (abs(T) .GE. zProfAbs(Nr)) THEN
        Tambient = tProf(Nr)
        Sambient = sProf(Nr)
    ELSE
        CALL linint(Nr, zProfAbs, tProf, T, Tambient)
        CALL linint(Nr, zProfAbs, sProf, T, Sambient)
    ENDIF

    rho_1 = RHO(Y(3), Y(4), T)
    rho_0 = RHO(Tambient, Sambient, T)

    ! Equations for Sb, Tb and mdot
    a = lambda1*(GamT*c_w-GamS*c_i)

    b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*T+ &
        &         iceTemp-(L/c_i)) &
        &        -GamT*c_w*(Y(3)-lambda2-lambda3*T)

    c = GamS*Y(4)*(c_i*(lambda2+lambda3*T-iceTemp)+L)

    Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
    Tb   = lambda1*Sb+lambda2+lambda3*T !Tb
    mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb ! mdot

    ! Differential equations
    ! Plume radius
    YDOT(1) = 2.*E_0+4.*mdot/(pi*Y(2))- &
            & Y(1)*g*(rho_0-rho_1)/(2.*Y(2)*Y(2)*rho_ref)+2.*Cd/pi

    ! Plume vertical velocity
    YDOT(2) = -2.*E_0*Y(2)/Y(1)-4.*mdot/(pi*Y(1))+g* &
        &          (rho_0-rho_1)/(Y(2)*rho_ref)-4.*Cd*Y(2)/(pi*Y(1))

    ! Plume temperature
    YDOT(3) = 2.*E_0*(TAMBIENT-Y(3))/Y(1)+4.*mdot* &
        &           (Tb-Y(3))/(pi*Y(1)*Y(2))-4.* &
        &           GamT*(Cd**0.5)*(Y(3)-Tb)/(pi*Y(1))

    ! Plume salinity
    YDOT(4) = 2.*E_0*(Sambient-Y(4))/Y(1)+4.*mdot*(Sb-Y(4))/ &
        &          (pi*Y(1)*Y(2))-4.*GamS*(Cd**0.5)*(Y(4)-Sb)/(pi*Y(1))

    ! Along-plume integrated contact area and melt rate
    YDOT(5) = 2.*Y(1)
    YDOT(6) = 2.*Y(1)*mdot

END SUBROUTINE HALFCONE

! =========================================================================

SUBROUTINE  JENKINS (NEQ, T, Y, YDOT)

    USE mod_param_iceplume
    implicit none

    integer ::  NEQ
    real(r8) :: T, Y(6), YDOT(6)
    real(r8) :: Tambient, Sambient, rho_0, rho_1
    real(r8) :: mdot, Sb, Tb
    real(r8) :: a,b,c
    real(r8) :: RHO

    ! Interpolate from imposed ambient profiles
    IF (abs(T) .LE. zProfAbs(1)) THEN
        Tambient = tProf(1)
        Sambient = sProf(1)
    ELSEIF (abs(T) .GE. zProfAbs(Nr)) THEN
        Tambient = tProf(Nr)
        Sambient = sProf(Nr)
    ELSE
        CALL linint(Nr, zProfAbs, tProf, T, Tambient)
        CALL linint(Nr, zProfAbs, sProf, T, Sambient)
    ENDIF

    rho_1 = RHO(Y(3), Y(4), T)
    rho_0 = RHO(Tambient, Sambient, T)

    ! Equations for Sb, Tb and mdot
    a = lambda1*(GamT*c_w-GamS*c_i)

    b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*T+ &
        &         iceTemp-(L/c_i)) &
        &        -GamT*c_w*(Y(3)-lambda2-lambda3*T)

    c = GamS*Y(4)*(c_i*(lambda2+lambda3*T-iceTemp)+L)

    Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
    Tb   = lambda1*Sb+lambda2+lambda3*T !Tb
    mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb ! mdot

    ! Differential equations
    ! Plume thickness
    YDOT(1)=2*E_0+Cd-(g*Y(1)/(Y(2)**2))*(rho_0-rho_1) &
    &  /rho_ref+2*mdot/Y(2)

    ! Plume vertical velocity
    YDOT(2)=-(Y(2)/Y(1))*(E_0+Cd+mdot/Y(2)) &
    &  +(g/Y(2))*(rho_0-rho_1)/rho_ref

    ! Plume temperature
    YDOT(3)=E_0*Tambient/Y(1)-(Y(3)/Y(1)) &
    &  *(E_0+mdot/Y(2))+(mdot/(Y(1)*Y(2))) &
    &  *(Tb-(L/c_w)-(c_i/c_w)*(Tb-iceTemp))

    ! Plume salinity
    YDOT(4)=E_0*Sambient/Y(1)-(Y(4)/Y(1)) &
    &  *(E_0+mdot/Y(2))

    ! along-plume integrated contact area and melt rate
    YDOT(5) = delta_y  ! This is constant in sheet model
    YDOT(6) = delta_y * mdot  ! This is constant in sheet model

END

! =========================================================================

DOUBLE PRECISION FUNCTION RHO(t,S,z)

    ! Equation of state (UNESCO 1983)
    !     T = temperature (deg C)
    !     S = salinity (PSU)
    !     z = depth (m)

    DOUBLE PRECISION T,S,z
    DOUBLE PRECISION rho_0, g, P
    DOUBLE PRECISION kw, Aw, Bw, k0
    DOUBLE PRECISION bulk_modulus
    DOUBLE PRECISION A, B, C, rho_w,rho_zero

    PARAMETER(rho_0=1027)
    PARAMETER(g=9.81)

    P= rho_0*g*abs(z)*1.0E-5
   
    ! RHO_1 (in situ)
    kw= 19652.21+ 148.4206*T- 2.327105*T**2+ &
     &    1.360477e-2*(T**3)-5.155288e-5*(T**4)
    Aw= 3.239908+ 1.43713e-3*T+ 1.16092e-4*T**2- &
     &    5.77905e-7*T**3
    Bw= 8.50935e-5- 6.12293e-6*T + 5.2787e-8*(T**2)
    k0= kw + (54.6746- 0.603459*T+ 1.09987e-2*(T**2) &
     &    -6.1670e-5*(T**3))*S +(7.944e-2 + 1.6483e-2* &
     &    T- 5.3009e-4*(T**2))*(S**1.5)
    A=  Aw+ (2.2838e-3- 1.0981e-5*T- 1.6078e-6*(T**2)) &
     &    *S+ 1.91075e-4*(S**1.5)
    B= Bw+ (-9.9348e-7+ 2.0816e-8*T+ 9.1697e-10*T**2)*S
    bulk_modulus= k0+ A*P+ B*P**2

    A= 8.24493e-1- 4.0899e-3*T+ 7.6438e-5*T**2- &
     &   8.2467e-7*T**3+5.3875e-9*T**4
    B= -5.72466e-3 + 1.0227e-4*T- 1.6546e-6*T**2
    C= 4.8314e-4
    rho_w= 999.842594 + 6.793952e-2*T- 9.095290e-3*T**2+ &
     &       1.001685e-4*T**3-1.120083e-6*T**4+ &
     &       6.536336e-9*T**5
    rho_zero= rho_w+ A*S + B*(S**1.5)+ C*(S**2)

    RHO= rho_zero/(1- (P/bulk_modulus))

END

! =========================================================================

subroutine linint(nx,xtab,ytab,x,y)

    ! Given a value of x return a value of y based on interpolation
    ! within a table of y values (ytab) corresponding to the x values
    ! contained in the array xtab.  The subroutine assumes that the
    ! values in xtab increase monotonically
    !
    ! John Mahaffy 2/12/95
    ! Modified slightly TRC 2014

    integer nx
    double precision xtab(nx), ytab(nx), x, y

    ! local variables
    integer i, i1
    double precision  wx

    if (x.lt.(xtab(1)).or.x.GT.(xtab(nx))) then
        write(6,*) 'x = ', x, '  is out of table range'
        stop
    endif
    do 100 i=2,nx
         if (x.le.xtab(i)) go to 200
    100 continue
    200 i1=i-1

    wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
    y=(1-wx)*ytab(i1)+wx*ytab(i1+1)

end

! =========================================================================
! Dummy routine for ODEPACK. Necessary for Jacobian matrix if stiff ODEs.

SUBROUTINE jex()
    RETURN
END

! =========================================================================
! Obsolete subroutines
! =========================================================================
!
! SUBROUTINE RHO_TO_W(profn, profr, profw)
! 
!     ! A quick program to extrapolate profiles from rho-points to w-points
!     implicit none
!     integer, parameter :: r8 = SELECTED_REAL_KIND(12,300)  ! 64-bit
!     integer, intent(in) :: profn
!     real(r8) :: profr(profn), profw(profn+1)
! 
!     integer :: K
! 
!     profw(1) = profr(1)
!     profw(profn+1) = profr(profn)
!     DO K = 2, profn
!         profw(K) = 0.5d0*(profr(K-1)+profr(K))
!     ENDDO
! 
! END SUBROUTINE RHO_TO_W
