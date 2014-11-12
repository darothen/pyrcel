module aerosol

    use constants
    use thermo
    use util

    implicit none

contains

    subroutine aerosol_bins( mu, sigma, N, r_drys, Nis )
        real (kind=dp), intent(in) :: mu, sigma, N
        real (kind=dp), dimension(:), intent(out) :: r_drys, Nis

        !------------------

        real (kind=dp) :: lr, rr, a, b
        real (kind=dp), dimension(size(r_drys)+1) :: bins
        integer :: i, nbins

        nbins = size(r_drys)

        ! left-right extent
        !lr = log10(mu/10d0/sigma)
        !rr = log10(mu*10d0*sigma)        
        lr = mu/10d0/sigma
        rr = mu*10d0*sigma

        ! populate the log of the bins
        call logspace ( lr, rr, bins )
        !do i = 1, nbins+1
        !    bins(i) = (rr - lr)*real(i-1, kind=dp) &
        !              / real(nbins-1, kind=dp) + lr
        !end do
        !bins = 10.d0**bins

        ! compute the midpoint size in each bin and the number
        do i = 1, nbins
            a = bins(i)
            b = bins(i+1)
            r_drys(i) = sqrt(a*b)
            Nis(i) = 0.5*(b-a)*(lognorm_pdf(a, mu, sigma, N) + &
                                lognorm_pdf(b, mu, sigma, N))
        end do

        r_drys = r_drys*1d-6
        Nis = Nis*1d6        

    end subroutine

    subroutine equil_bins ( T0, S0, r_drys, kappas, r_wets)
        real (kind=dp), intent(in) :: T0, S0
        real (kind=dp), dimension(:), intent(in) :: r_drys, kappas
        real (kind=dp), dimension(:), intent(out) :: r_wets

        !------------------
        integer :: i, count, nbins
        real (kind=dp) :: a, b, c, fa, fb, fc
        integer, parameter :: max_n = 100
        real (kind=dp), parameter :: tol = 1d-10

        nbins = size(r_drys)

        write (*,"('computing wet equilibrium radii for',I2,' bins')") nbins
        do i = 1, nbins
            !write (*,*) "bin", i

            a = r_drys(i)
            b = kohler_crit(T0, r_drys(i), kappas(i))
            c = (a + b)/2d0

            fa = Seq(a, r_drys(i), T0, kappas(i)) - S0
            fb = Seq(b, r_drys(i), T0, kappas(i)) - S0
            fc = Seq(c, r_drys(i), T0, kappas(i)) - S0

            !write (*,*), fa, fc, fb

            count = 1
            do while (abs(fc) > tol)      
                !write (*,*) "   ", count, a, b, abs(fc)
                if (count > max_n) then
                    write (*,*) "bisection failed on element", count
                    stop
                endif
                
                c = (a + b)/2d0

                fa = Seq(a, r_drys(i), T0, kappas(i)) - S0
                fb = Seq(b, r_drys(i), T0, kappas(i)) - S0
                fc = Seq(c, r_drys(i), T0, kappas(i)) - S0

                if (fc*fa > 0) then
                    a = c
                else
                    b = c
                endif

                count = count + 1
            end do
            r_wets(i) = c

        end do

        return

    end subroutine equil_bins

    real (kind=dp) function lognorm_pdf( x, mu, sigma, N )
        real (kind=dp), intent(in) :: x, mu, sigma, N
        real (kind=dp) :: scaling, exponent

        scaling = N/sqrt(2d0*pi)/log(sigma)
        exponent = (log(x/mu)**2d0)/2d0/(log(sigma)**2d0)
        lognorm_pdf = (scaling/x)*exp(-1d0*exponent)

    end function lognorm_pdf


end module aerosol