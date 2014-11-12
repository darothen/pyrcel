module model
    
    use aerosol
    use constants
    use thermo
    use dvode_f90_m

    implicit none

    integer :: nbins, nmodes
    real (kind=dp) :: V
    real (kind=dp), dimension(:), allocatable :: kappas, Nis, r_drys, r0s
    real (kind=dp), dimension(:), allocatable :: aer_mus, aer_Ns, aer_sigs, aer_kaps

contains

    subroutine run_model( T0, S0, P0, Smax )
        real (kind=dp), intent(in) :: T0, S0, P0
        real (kind=dp), intent(out) :: Smax

        real (kind=dp) :: wv0, wc0, N
        integer :: i, j, k

        !---- DVODE solver vars
        real(kind=dp) atol, rtol, t, tout, y, rstats
        integer neq, itask, istate, istats, iout, ierror
        dimension y(5+nbins), atol(5+nbins), rstats(22), istats(31)    
        type (vode_opts) :: options

        integer :: ntout
        real(kind=dp) :: tfinal
        real(kind=dp) :: deltat
        logical :: Smax_found

        !-------------------------------------------
        ! SETUP METEOROLOGY
        !-------------------------------------------
        wv0 = (S0 + 1d0)*0.622d0*es(T0-273.15d0)/(P0-es(T0-273.15d0))

        y(1) = 0d0
        y(2) = T0
        y(3) = wv0
        y(5) = S0

        !-------------------------------------------
        ! SETUP AEROSOL
        !-------------------------------------------

        !-------------------------------------------
        ! 1) Compute the dry size/# for all the bins
        ! 2) Compute the equil wet size for all bins
        !    - using bisection method
        !-------------------------------------------
        do i = 1, nmodes
            write (6,80001) i, aer_mus(i), aer_sigs(i), aer_Ns(i)
80001 FORMAT ('Mode number - ',I2,' mu: ',F5.3,' sigma: ',F3.1, &
              ' N: ',F7.1)
            call aerosol_bins(aer_mus(i), aer_sigs(i), aer_Ns(i), &
                              r_drys((i-1)*nbins:i*nbins),        &
                              Nis((i-1)*nbins:i*nbins)             )
            kappas((i-1)*nbins:i*nbins) = aer_kaps(i)
        end do

        call equil_bins(T0, S0, r_drys, kappas, r0s)

        wc0 = 0d0
        do j = 1, nbins
            wc0 = wc0 + (4d0*pi/3d0)*rho_w*Nis(j) &
                       *(r0s(j)**3d0 - r_drys(j)**3d0)
        end do
        y(4) = wc0
        y(6:) = r0s

        !do i = 1, nbins
        !    write (*,*) i, r_drys(i), r0s(i), r0s(i)/r_drys(i)
        !end do 
        !write (*,*) "--------------------------"
        write (*,*) "TOTAL COMPARE:", N, sum(Nis)*1d-6
        !write (*,*) wc0

        !-------------------------------------------
        ! SETUP AND RUN DVODE
        !-------------------------------------------
        write (*,*) "Running parcel model"

        open (unit=7, file='parcel.dat')
        ierror = 0
        itask  = 1
        istate = 1    
        neq = 5+nbins


        rtol = 1d-7
        atol(1:5) = (/ 1d-4, 1d-4, 1d-10, 1d-10, 1d-8 /)
        atol(6:neq) = 1d-12

        options = SET_NORMAL_OPTS(dense_j=.true., abserr_vector=atol, &
                                  relerr=rtol)

        tfinal = 1800.0d0
        deltat = 0.1d0
        t = 0.0d0
        tout = t
        ntout  = 1 + int((tfinal - t)/deltat)

        !print *, y(1:5)
        !call parcel_deriv(neq,t,y,ydt)
        !stop

        Smax = S0 - 100d0
        Smax_found = .false.
        do while (tout < tfinal)
            tout = tout + deltat
            call DVODE_F90(parcel_deriv,neq,y,t,tout,itask,istate,options)
            call GET_STATS(rstats,istats)

            if (istate < 0) then ! there was an error
                write (7,90003) istate
                stop
            end if

            ! print the solution and write the plot data file
            !write (7, *) tout, y
            !write (6, 90001) tout, y(5), y(neq)
            !print *, y(5), Smax, Smax_found

            ! Have we found the maximum Smax?
            if (y(5) > Smax) then 
                Smax = y(5)
            else if ((y(5) < Smax) .and. (.not. Smax_found)) then
                Smax_found = .true.
                tfinal = tout+60.d0
                !write (6,90000) tout 
            end if

            if (tout > tfinal) exit

        end do

        write (6,*) "------------------------------------"
        write (6,90004) Smax
        write (6,90002) istats(11), istats(12), istats(13)

90000 FORMAT (' Found Smax at time (t=',F6.1,'), running for 60 more seconds')
90001 FORMAT (' time ',F6.1,' | S = ',F8.6,'  rmax = ', ES7.2E1)
90002 FORMAT (' Steps = ',I10,' f-s = ',I10,' J-s = ',I10)
90003 FORMAT (' An error occurred in VODE_F90. ISTATE = ',I3)   
90004 FORMAT (' Smax = ',ES9.3)

    end subroutine run_model

    subroutine parcel_deriv(neq, t, y, ydot)

        integer, intent (in) :: neq
        real(kind=dp), intent (in) :: t 
        real(kind=dp), intent (in) :: y(neq)
        real(kind=dp), intent (out) :: ydot(neq)

        integer :: nbins
        real (kind=dp) :: dz_dt, dwc_dt, dwv_dt, dT_dt, dS_dt
        real (kind=dp) :: drs_dt(neq-5)
        integer :: i
        real (kind=dp) :: G_a, G_b, G, r_i, r_dry_i, delta_S,    & 
                          kappa_i, dr_dt, Ni, dv_r, ka_r, P_atm, &
                          A, B, Seq_r
        real (kind=dp) :: alpha, gamma

        real(kind=dp) :: z, Te, wv, wc, S
        real(kind=dp) :: rs(neq-5)

        real(kind=dp) :: T_c, pv_sat, wv_sat, Tv, P, rho_air

        nbins = neq-5

        z  = y(1)
        Te = y(2)
        wv = y(3)
        wc = y(4)
        S  = y(5)
        rs(:) = y(6:neq)

        T_c = Te - 273.15d0
        pv_sat = es(T_c)
        wv_sat = wv/(S + 1d0)
        Tv = (1d0 + 0.61d0*wv)*Te
        P = pv_sat*(1d0 + 0.622*(S + 1d0)/wv)
        rho_air = P/Rd/Tv

        ! ---- Compute tendencies

        dwc_dt = 0d0
        do i = 1, nbins
            r_i     = rs(i)
            r_dry_i = r_drys(i)
            kappa_i = kappas(i)
            Ni      = Nis(i)

            dv_r = dv(Te, r_i, P, accom)
            ka_r = ka(Te, r_i, rho_air)

            G_a = (rho_w*R_uni*Te)/pv_sat/dv_r/Mw
            G_b = (L*rho_w*((L*Mw/R_uni/Te)-1d0))/ka_r/Te
            G   = 1d0/(G_a + G_b)

            Seq_r   = Seq(r_i, r_dry_i, Te, kappa_i)
            delta_S = S - Seq_r

            dr_dt = (G/r_i)*delta_S
            dwc_dt = dwc_dt + Ni*r_i*r_i*dr_dt
            drs_dt(i) = dr_dt

            !if (i == nbins) then
            !    print *, 0, Te, r_i, P, accom, rho_air
            !    print *, 1, r_i, r_dry_i, Ni
            !    print *, 2, dv_r, ka_r
            !    print *, 3, G_a, G_b
            !    print *, 4, Seq_r
            !    print *, 5, dr_dt
            !endif 

        end do
        dwc_dt = dwc_dt*4d0*pi*rho_w/rho_air

        dwv_dt = -dwc_dt
        dT_dt  = -grav*V/Cp - L*dwv_dt/Cp
        dz_dt  = V

        alpha = grav*Mw*L/Cp/R_uni/Te/Te - grav*Ma/R_uni/Te
        gamma = P*Ma/Mw/pv_sat + Mw*L*L/Cp/R_uni/Te/Te
        dS_dt = alpha*V - gamma*dwc_dt

        ydot(1) = dz_dt
        ydot(2) = dT_dt
        ydot(3) = dwv_dt
        ydot(4) = dwc_dt
        ydot(5) = dS_dt
        ydot(6:neq) = drs_dt(:)

        !print *, ydot(1:5)

        return

    end subroutine parcel_deriv

    subroutine arg2000( Te, P, accom, Smax )
        real (kind=dp), intent(in) :: Te, P, accom 
        real (kind=dp), intent(out) :: Smax

        real (kind=dp) :: wv_sat, alpha, gamma
        real (kind=dp) :: G, G_a, G_b, G_0, G_ac, G_ac1, dv_r, ka_r

        real (kind=dp) :: am, N, sigma, kappa, fi, gi
        real (kind=dp) :: A, rc_mode, Smi2
        real (kind=dp) :: zeta, etai, Spa, Spb, S_part
        real (kind=dp), dimension(nmodes) :: Smis, Sparts
        integer :: i

        wv_sat = es(Te - 273.15d0)
        alpha = grav*Mw*L/Cp/R_uni/Te/Te - grav*Ma/R_uni/Te
        gamma = R_uni*Te/wv_sat/Mw + Mw*L*L/Cp/Ma/Te/P

        dv_r = 1d-4*(0.211d0/(P*1.01325d-5)*((Te/273d0)**1.94d0))
        ka_r = 1d-3*(4.39d0 + 0.071d0*Te)

        G_a = (rho_w*R_uni*Te)/wv_sat/dv_r/Mw
        G_b = (L*rho_w*((L*Mw/R_uni/Te)-1d0))/ka_r/Te
        G_0 = 1d0/(G_a + G_b)

        do i = 1, nmodes
            am    = aer_mus(i)*1d-6
            N     = aer_Ns(i)*1d6
            sigma = aer_sigs(i)
            kappa = aer_kaps(i)

            fi = 0.5d0*exp(2.5d0*(log(sigma)**2d0))
            gi = 1d0 + 0.25d0*log(sigma)

            A = 2d0*Mw*sigma_w(Te)/R_uni/Te/rho_w
            Smi2    = sqrt( 4d0*(A**3)/27d0/kappa/(am**3) )
            rc_mode = sqrt( 3d0*kappa*(am**3)/A )

            if (accom == 1d0) then
                G = G_0
            else
                dv_r = dv(Te, rc_mode, P, accom)
                G_a  = (rho_w*R_uni*Te)/wv_sat/dv_r/Mw
                G_b  = (L*rho_w*((L*Mw/R_uni/Te)-1d0))/ka_r/Te
                G_ac = 1d0/(G_a + G_b)

                dv_r  = dv(Te, rc_mode, P, 1d0)
                G_a   = (rho_w*R_uni*Te)/wv_sat/dv_r/Mw
                G_ac1 = 1d0/(G_a + G_b)

                G = G_0*G_ac/G_ac1
            end if

            zeta = (2d0/3d0)*A*sqrt(alpha*V/G)
            etai = ((alpha*V/G)**(3d0/2d0))/N/gamma/rho_w/2d0/pi

            Spa = fi*((zeta/etai)**(1.5d0))
            Spb = gi*(((Smi2**2d0)/(etai + 3d0*zeta))**(0.75d0))
            S_part = (1d0/(Smi2**2d0))*(Spa + Spb)           

            Smis(i) = Smi2
            Sparts(i) = S_part
        end do

        Smax = 1d0/sqrt(sum(Sparts))

    end subroutine arg2000

end module model