module util

    use constants, only: dp

    implicit none

contains

    subroutine linspace(xl, xr, x)
        real (kind=dp), intent(in) :: xl, xr
        real (kind=dp), dimension(:), intent(out) :: x

        integer i, n

        n = size(x)
        if (n == 1) then
            if ( xl /= xr ) then
                print *, "ERROR: n = 1, but xl != xr"
                stop
            else
                x = xl
            end if
        else
            do i = 1, n
                x(i) = (xr - xl)*real(i-1, dp) / real(n-1, dp) + xl
            end do
        end if

        return

    end subroutine linspace

    subroutine logspace(xl, xr, x)
        real (kind=dp), intent(in) :: xl, xr
        real (kind=dp), dimension(:), intent(out) :: x

        if ((size(x) == 1) .and. ( xl /= xr )) then
            print *, "ERROR: n = 1, but xl != xr"
            stop
        end if

        call linspace ( log10(xl), log10(xr), x )
        x = 10.d0**x

        return
        
    end subroutine logspace

end module util