SUBROUTINE ICEPLUME_MELTRATE(temperature, salinity, velocity, depth, &
                           & mdot)

    USE mod_iceplume
    implicit none

    real(r8) :: temperature, salinity, velocity, &
              & depth, absVelocity
    real(r8) :: a, b, c, tb, sb, mdot

    ! Routine can't cope with zero velocity.
    ! Unlikely to occur anyway with currents, waves, convection etc.
    ! This isn't very physical, but will do for now.
    IF ( velocity .LT. backgroundVel ) velocity = backgroundVel

    absVelocity = abs(velocity)

    ! Calculate melt rate from 3 equation formualtion (as for plume models)

    ! Equations for Sb, Tb and mdot

    a = lambda1*(GamT*c_w-GamS*c_i)

    b = GamS*c_i*(lambda1*salinity-lambda2-lambda3*depth+ &
     &         iceTemp-(L/c_i)) &
     &        -GamT*c_w*(temperature-lambda2-lambda3*depth)

    c = GamS*salinity*(c_i*(lambda2+lambda3*depth-iceTemp)+L)

    Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
    Tb   = lambda1*Sb+lambda2+lambda3*depth !Tb
    mdot = GamS*(Cd**0.5)*absVelocity*(salinity-Sb)/Sb ! mdot

END
