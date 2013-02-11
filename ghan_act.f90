! Wrapper for Steve Ghan's explicit/parameterized aerosol activation code
! subgrid.f (repackaged as subgrid_mod.f90)
!
! Author - Daniel Rothenberg (darothen@mit.edu)
! Version - 2/10/2013

! Compilation -
! gfortran -c subgrid_mod.f90 -o subgrid_mod.o -ffixed-form
! gfortran -c ghan_act.f90 -o ghan_act.o
! f2py -c subgrid_mod.o -m ghan_act ghan_act.f90

subroutine activation(tair, pres, wbar, explicit, &
                      nmode, na, sig, rhodry_in, hygro_in, rad, &
                      fn, smax)

    use subgrid, only: activate

    implicit none

    ! -- Input variables
    real :: tair ! parcel temperature, K
    real :: pres ! parcel pressure, mb
    real :: wbar ! vertical velocity, m/s
    logical :: explicit ! switch for explicit or Razzak

    integer :: nmode ! number of modes
    real, dimension(nmode) :: na ! mode number concs, #/cm3
    real, dimension(nmode) :: sig ! mode sigmas
    real, dimension(nmode) :: rhodry_in ! mode dry density, g/cm3
    real, dimension(nmode) :: hygro_in ! mode hygroscopicities
    real, dimension(nmode) :: rad ! mode number mode radius, microns

    ! -- Output variables
    real, dimension(nmode), intent(out) :: fn ! number fraction activated
    real, intent(out) :: smax

    ! -- Local variables and parameters
    real :: third = 1./3.
    real :: pi = 3.14159
    real, dimension(nmode) :: alnsign
    real, dimension(nmode) :: dgnum ! mode number mode diameters, microns
    real, dimension(1, nmode) :: mass, vol, rhodry, hygro
    real, dimension(nmode) :: am ! number mode radius (cm)

    real, dimension(nmode) :: fm ! mass fraction activated
    real, dimension(nmode) :: fluxn ! number flux into cloud
    real, dimension(nmode) :: fluxm ! fractional mass flux into cloud

    logical :: top = .true.
    logical :: interior = .false.
    integer :: gaskin = 0 ! disable gas kinetics
    real :: sigw = 0.0
    real :: wdiab = 0.0
    real :: sds = 3
    real :: fmax = 0.99
    real :: eps = 0.3

    integer :: m
    logical :: explicit_activate, razzak
    real :: hydrorad

    ! -- F2PY VARIABLE BINDINGS
        ! f2py intent(in) :: tair, pres, wbar, explicit, na, sig, rhodry_in, hygro_in, rad
        ! f2py intent(out) :: fn, smax
        ! f2py intent(hide) :: nmode

    !-----------------------------------------!
    ! Begin subroutine code
    ! convert shape of input rhodry and hygro
    rhodry(1,:) = rhodry_in(:)
    hygro(1,:) = hygro_in(:)

    ! Toggle explicit or razzak
    explicit_activate = explicit
    razzak = .not.explicit

    ! convert pressure to CGS
    pres = pres*1.e3

    ! Mode mass/volume calculations
    alnsign = alog(sig)
    dgnum = 2*rad
    ! convert radius to ratio of mass of number
    rad = rad*1.e-4
    mass(1,:) = (4./3.)*pi*rhodry(1,:)*na*(rad**3)*exp(4.5*alnsign*alnsign)
    vol(1,:) = mass(1,:)/rhodry(1,:)
    am = exp(-1.5*alnsign*alnsign)*(3.*vol(1,:)/(4.*pi*na))**third
    dgnum = 2*am

    write (*,*) 'tair=', tair, 'pres=', pres

    ! output aerosol definitions
    write (*, 9205), 'Mode    ','N (#/cm3)   ','M (ug/m3)   ', &
       'dgnum (cm)  ','sigma-g     ','rho-dry     ', &
       'hygro          '
9205 format(10a)

    do m = 1,nmode
        write (*, 9210) m,na(m),1.e12*mass(1,m),dgnum(m),sig(m),rhodry(1,m),hygro(1,m)
9210    format(i5,9g12.3)
    enddo
    write (*,*) 'wbar, sigw, wdiab (cm/s) =', wbar,sigw,wdiab


    ! Call the activation routine
    call activate(top,gaskin,eps,wbar,sigw,wdiab,tair,pres,razzak, &
                  explicit_activate,na,1,nmode,mass,sig,hygro, &
                  rhodry,interior,hydrorad,fmax,sds, &
                  fn,fm,fluxn,fluxm,smax)

    !write (*, *) "at the end - smax=", smax, "act_frac=",fn

    return

end subroutine activation