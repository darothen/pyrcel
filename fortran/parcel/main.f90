program main

    use aerosol
    use constants
    use model
    use thermo
    use util

    implicit none

    integer, parameter :: nVs = 1

    !-------------------------------------------    
    real (kind=dp) :: T0, S0, P0, Smax
    integer :: i

    real (kind=dp), dimension(nVs) :: Vs, Smaxes
    
    nbins = 200
    nmodes = 1

!----- Allocate run-time variables ---------
    allocate( aer_mus(nmodes)      )
    allocate( aer_Ns(nmodes)       )
    allocate( aer_sigs(nmodes)     )
    allocate( aer_kaps(nmodes)     )
    allocate( kappas(nbins*nmodes) )
    allocate( Nis(nbins*nmodes)    )
    allocate( r_drys(nbins*nmodes) )
    allocate( r0s(nbins*nmodes)    )

    V = 0.2d0 ! setting module-level variable in model.mod
    T0 = 283d0
    P0 = 85000d0
    S0 = -0.01d0

    ! Aerosol modes
    !aer_mus  (1) = 0.2d0
    !aer_sigs (1) = 1.59d0
    !aer_Ns   (1) = 400d0
    !aer_kaps (1) = 0.507d0
    aer_mus  (1) = 0.2d0
    aer_sigs (1) = 1.59d0
    aer_Ns   (1) = 400d0
    aer_kaps (1) = 0.507d0


    !aer_mus  (2) = 0.055d0
    !aer_sigs (2) = 2.0d0
    !aer_Ns   (2) = 4000d0
    !aer_kaps (2) = 0.25d0

    if (nVs > 1) then
        call logspace ( 1d-1, 10d0, Vs)
        do i = 1, nVs
            V = Vs(i)
            call run_model( T0, S0, P0, Smax)
            Smaxes(i) = Smax
        end do
        
        do i = 1, nVs
            call arg2000(T0, P0, accom, Smax)
            print *, i, Vs(i), Smaxes(i), Smax
        end do
    else
        call run_model( T0, S0, P0, Smax )
        print *, V, Smax
        call arg2000( T0, P0, accom, Smax )
        print *, "ARG", Smax
    end if


!-------------------------------------------
! SIMULATION PARAMETERS
!-------------------------------------------

! ---- Deallocation ------
    deallocate( kappas )
    deallocate( Nis    )
    deallocate( r_drys )
    deallocate( r0s    )
    deallocate( aer_mus  )
    deallocate( aer_Ns   )
    deallocate( aer_sigs )
    deallocate( aer_kaps )

end program main