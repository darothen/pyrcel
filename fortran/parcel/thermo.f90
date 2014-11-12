module thermo

    use constants

    implicit none

    public

contains

    real (kind=8) function dv(Te, r, P, accom)
        real(kind=dp), intent(in)  :: Te, r, P, accom

        real(kind=dp) :: dv_t, denom, P_atm

        P_atm = P*1.01325d-5 ! convert pressure from Pa -> atm
        dv_t = (1.d-4)*(0.211d0/P_atm)*((Te/273.d0)**1.94d0)
        denom = 1.d0 + (dv_t/accom/r)*sqrt(2.d0*pi*Mw/R_uni/Te)
        dv = dv_t/denom

    end function dv

    real (kind=dp) function es(T_c)
        real(kind=dp), intent(in)  :: T_c

        es = 611.2*exp(17.67*T_c/(T_c + 243.5))

    end function es

    real (kind=dp) function ka(Te, rho, r)
        real(kind=dp), intent(in)  :: Te, rho, r

        real(kind=dp) :: ka_t, denom

        ka_t = (1d-3)*(4.39d0 + 0.071d0*Te)
        denom = 1d0 + (ka_t/at/r/rho/Cp)*sqrt(2.d0*pi*Ma/R_uni/Te)

        ka = ka_t/denom

    end function ka

    real (kind=dp) function kohler_crit(T, r_dry, kappa)
        real (kind=dp), intent(in) :: T, r_dry, kappa
        real (kind=dp) :: A

        A = 2d0*Mw*sigma_w(T)/R_uni/T/rho_w
        kohler_crit = sqrt( 3d0*kappa*(r_dry**3d0)/A )

    end function kohler_crit

    real (kind=dp) function rho_air(T, P, RH)
        real(kind=dp), intent(in)  :: T, P, RH
        real(kind=dp) :: qsat, Tv

        qsat    = RH*epsilon*(es(T - 273.15)/P)
        Tv      = T*(1.0 + 0.61*qsat)
        rho_air = P/Rd/Tv

    end function rho_air

    real (kind=dp) function Seq(r, r_dry, T, kappa)
        real (kind=dp), intent(in)  :: r, r_dry, T, kappa
        real (kind=dp) :: A, B

        A = (2d0*Mw*sigma_w(T))/R_uni/T/rho_w/r

        !Seq = A - kappa*(r_dry**3d0)/(r**3d0)
        B = (r**3d0 - r_dry**3d0)/(r**3d0 - (r_dry**3d0)*(1d0-kappa))
        Seq = exp(A)*B - 1d0

    end function Seq

    real (kind=dp) function sigma_w(T)
        real(kind=dp), intent(in)  :: T

        sigma_w = 0.0761d0 - (1.55d-4)*(T - 273.15d0)

    end function sigma_w

end module thermo
