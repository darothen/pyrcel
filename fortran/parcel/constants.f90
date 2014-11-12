module constants

    implicit none
    public

    integer, parameter :: dp = kind(1.d0)

    real (kind=dp), parameter :: grav = 9.81d0
    real (kind=dp), parameter :: Cp = 1004.0d0
    real (kind=dp), parameter :: L = 2.25d6
    real (kind=dp), parameter :: rho_w = 1.0d3
    real (kind=dp), parameter :: Rd = 287.d0
    real (kind=dp), parameter :: Rv = 461.5d0
    real (kind=dp), parameter :: R_uni = 8.314d0
    real (kind=dp), parameter :: Mw = 0.018d0
    real (kind=dp), parameter :: Ma = 0.0289d0
    !real (kind=dp), parameter :: Dv = 3.d-5
    real (kind=dp), parameter :: accom = 1.d0
    !real (kind=dp), parameter :: Ka = 2.d-2
    real (kind=dp), parameter :: at = 0.96d0
    real (kind=dp), parameter :: epsilon = 0.622d0

    real (kind=dp), parameter :: pi = 4 * atan(1.0d0)

end module constants