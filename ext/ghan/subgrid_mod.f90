c gfortran -c subgrid_mod.f90 -o subgrid_mod.mod -ffixed-form

      module subgrid

c--------------------------------------------------------------------------
c this data subprogram loads variables into the internal common
c blocks used by the odepack solvers.  the variables are
c defined as follows..
c   illin  = counter for the number of consecutive times the package
c             was called with illegal input.  the run is stopped when
c             illin reaches 5.
c   ntrep  = counter for the number of consecutive times the package
c             was called with istate = 1 and tout = t.  the run is
c             stopped when ntrep reaches 5.
c   mesflg = flag to control printing of error messages.  1 means print,
c             0 means no printing.
c   lunit  = default value of logical unit number for printing of error
c             messages.
c--------------------------------------------------------------------------
      integer :: illin=0, iduma(9), ntrep=0, idumb(2), iowns(6),
     &           icomm(19),  mesflg=1, lunit=6
      real*8 rowns(209), rcomm(9)

c
c-------------------------- end of block data -----------------------------

      contains

      subroutine activate(top,gaskin,eps,wbar,sigw,wdiab,tair,pres,
     &                    razzak,explicit_activate,na,ntype,nmode,mass,
     &                    sign,hygro,rhodry,interior, hydrorad,
     &                    fmax,sds,fn,fm,fluxn,fluxm,
     &                    smax)


c     cgs units
      integer pmode
      parameter (nx=1000,pmode=40)
      logical top,razzak,explicit_activate, interior
      real wbar    ! grid cell mean vertical velocity (cm/s)
      real sigw    ! subgrid standard deviation of vertical velocity (cm/s)
      real wdiab   ! diabatic vertical velocity (0 if adiabatic)
      real tair    ! air temperature (K)
      real pres    ! air pressure (cgs)
       real hydrorad ! product of number and radius for hydrometeors
      integer ntype      ! number of aerosol types
      integer nmode      ! number of aerosol modes
      real na(pmode)     ! aerosol number concentration (/cc)
      real mass(ntype,pmode) ! dry aerosol mass concentration (g/cc)
      real sign(pmode)    ! geometric standard deviation of size distribution
      real alnsign(pmode) ! ln geometric standard deviation
      integer gaskin ! 1 for gas kinetic effects on vapor diffusivity, 0 for none

      real hygro(ntype,pmode)     ! hygroscopicity
      real mw            ! molecular weight of water (grams/mole)
      real rhow          ! density of liquid water (g/cm3)
      real rhodry(ntype,pmode) ! density of dry aerosol (g/cm3)
      real latent        ! latent heat of evaporation (erg/mole/deg)
      real surften       ! surface tension of water w/respect to air (dyne/cm)
      real np     ! potential number that can be activated, given forcing
c                 ! Ghan et al, Part I, eqn (25), * U
      real amb(pmode), ! (number mode radius (cm))**bexp
     &     am(pmode), ! number mode radius (cm)
     &     e(pmode),  ! Ghan et al, Part II, eqn (27)
     &     u(pmode),  ! Ghan et al, Part II, eqn (18)
     &     h(pmode),  ! Ghan et al, Part II, eqn (15)
     &     h_over_a_twothirds(pmode), ! Ghan et al, Part II, eqn (24)
     &     bexp(pmode),! exponent 4/(sqrt(2 pi) ln sigma)
     &     acb(pmode)  ! (size of small particle activated (cm))**bexp
      real expsign(pmode) ! Ghan et al, Part II, eqn (8)
      real fn(pmode)   ! number fraction activated
      real fm(pmode)   ! mass fraction activated
      real fluxn(pmode)   ! number flux into cloud
      real fluxm(pmode)   ! fractional mass flux into cloud
      real fnold(pmode)   ! number fraction activated
      real fmold(pmode)   ! mass fraction activated
      real integ,pi
      real sumflxn(pmode)
      real sumflxm(pmode)
      real sumfn(pmode)
      real sumfm(pmode)
      real sm(pmode),eta(pmode),lnsm(pmode),b(pmode)
      real sqrtg(pmode),beta(pmode),zeta(pmode)

c     activation coef (cm4/s)
      data c/0.003/,pi/3.14159/
      save pi
      data data pi/3.14159/
      save mw,rhow,ugascon,surften,cp,gasair,gasv
      save latent,grav,p0,t0,epsilon
      data mw/18./,rhow/1./,ugascon/8.3e7/,surften/76./
      data cp/1004.e4/,gasair/287.e4/,gasv/461.e4/,latent/2.5e10/
      data grav/981./,p0/1013.25e3/,t0/273.15/,epsilon/0.622/
      save third,sq2,sqpi
      data third/0.33333333/,sq2/1.414214/,sqpi/1.772454/
      logical new
      real*4 erf,erfc
      real sixth
      data sixth/0.166666667/
      save sixth
c      data accomc/0.042/,accomt/0.96/
      data accomc/1.0/,accomt/0.96/
      !write(*,*) 'enter accomodation coefficient (default is',accomc,')'
      !read(*,*) accomc
      if(nmode.gt.pmode)then
         print *,'nmode=',nmode,' pmode=',pmode
         stop
      endif
      if(explicit_activate)then
c        initialization for explicit activation calculation
         call explactconstant
         presmb=pres*1.e-3
         call explactcoeff(tair,presmb)
      endif

c     effect of organics on surface tension is neglected
      a=2.*mw*surften/(ugascon*tair*rhow)
      print *,'surface tension A=',a
      rhoair=pres/(gasair*tair)
      es=6.11*1.e3*exp(latent/gasv*(1/t0-1/tair))
      qs=epsilon*es/pres
      dqsdt=latent/(gasv*tair*tair)*qs
      dqsdp=-qs/pres
      path=6.6e-6*(p0/pres)*(tair/t0)
      deltat=2.7*path
      deltav=1.3*path
      diff0=0.211*(p0/pres)*(tair/t0)**1.94
      conduct0=(5.69+0.017*(tair-t0))*4.186e2
      g0=1./(rhow/(diff0*rhoair*qs)
     &    +latent*rhow/(conduct0*tair)*(latent/(gasv*tair)-1.))
      if(gaskin.eq.1)then
         diffparam=2.*diff0/accomc*sqrt(2.*pi/(gasv*tair))
	 dropmax=5.e-4
	 dropmin=0.6e-4
	 arg=(dropmin+diffparam)/(dropmax+diffparam)
	 alogarg=alog(arg)
         diff=diff0*(1.+diffparam/(dropmax-dropmin)*alogarg)
	 conduct=conduct0
	 print *,'diffparam,arg,alogarg,diff0,diff=',
     &  	  diffparam,arg,alogarg,diff0,diff
      else
         diff=diff0
	 conduct=conduct0
      endif

c     Razzak parameters
      if(razzak)then
         alpha=grav*(latent/(cp*gasv*tair*tair)-1./(gasair*tair))
         gamma=(1+latent/cp*dqsdt)/(rhoair*qs)
            g=1./(rhow/(diff*rhoair*qs)
     &      +latent*rhow/(conduct0*tair)*(latent/(gasv*tair)-1.))
      endif

      do m=1,nmode
c        internal mixture of aerosols
         vol=0
         sol=0
         tmass=0
         b(m)=0
         do n=1,ntype
            vol=vol+mass(n,m)/rhodry(n,m)
            tmass=tmass+mass(n,m)
            b(m)=b(m)+hygro(n,m)*mass(n,m)/rhodry(n,m)
         enddo
         b(m)=b(m)/(vol)
         print *,'mode=',m,' hygro=',b(m)
         h(m)=sqrt(3.*b(m)/a)
         h_over_a_twothirds(m)=(h(m)/a)**(2./3.)
         alnsign(m)=alog(sign(m))
c        number mode radius (cm)
         am(m)=exp(-1.5*alnsign(m)*alnsign(m))*
     &       (3.*vol/(4.*pi*na(m)))**third
         rm=h(m)*am(m)**1.5
         diff=diff0/(rm/(rm+deltav)
     &      +diff0/(rm*accomc)*sqrt(2*pi/(gasv*tair)))
         diff1=diff0/(rm/(rm+deltav)
     &      +diff0/(rm*1.0)*sqrt(2*pi/(gasv*tair)))
c         diff=diff0
         conduct=conduct0/(rm/(rm+deltat)+
     &        conduct/(rm*accomt*cp*rhoair)*sqrt(2*pi/(gasair*tair)))
c         conduct=conduct0
         g1=1./(rhow/(diff1*rhoair*qs)
     &    +latent*rhow/(conduct*tair)*(latent/(gasv*tair)-1.))
         gac=1./(rhow/(diff*rhoair*qs)
     &    +latent*rhow/(conduct*tair)*(latent/(gasv*tair)-1.))
         g=g0*gac/g1
         sqrtg(m)=sqrt(g)
         beta(m)=4*pi*rhow*g*gamma
         gammastar=4*pi*rhow*gamma
       print *,'G0=',g0
       print *,'G1=',g1
       print *,'G=',g
       print *,'alpha=',alpha
       print *,'gammastar=',gammastar
       print *,'beta=',beta(m)
         u(m)=8.*pi*rhow*g*a/(3.*rhoair)
         earg=-1.5*alnsign(m)/sq2
         e(m)=0.5*exp(1.125*alnsign(m)*alnsign(m))*(1.-erf(earg))-0.5
         expsign(m)=exp(-12.*alnsign(m)/(sq2*sqpi))
         bexp(m)=4/(alnsign(m)*sqrt(2.*pi))
         amb(m)=am(m)**bexp(m)
         sm(m)=2*a/(3.*am(m))*sqrt(a/(3.*b(m)*am(m)))
         lnsm(m)=alog(sm(m))
	 write (6,'(a,e12.4,a,e12.4)')'dry rad=',am(m),' crit S',sm(m)
      enddo

      if(explicit_activate)then
      	call explactstate(nmode,na,alnsign,am,b,tair,presmb)
      endif

      if(top)then
        wmax=0.
        w=min(0.,-wdiab)
      else
        wmax=wbar+sds*sigw
        w=-wdiab
      endif
      w=amax1(w,wbar-sds*sigw)
      wmin=w
      dwmax=eps*sigw
      dw=dwmax
      dfmax=0.2
      dfmin=0.1
      write(*,*) 'wmax=',wmax,' wmin=',wmin
      do m=1,nmode
         sumflxn(m)=0.
         sumflxm(m)=0.
         sumfn(m)=0.
         sumfm(m)=0.
      enddo
c     write(*,9005)
c9005 format(11x,'z',11x,'w',11x,'f',11x,'g',11x,'i')

      if(sigw.gt.1.e-3)then

         if(wmax.le.wmin)go to 20

         z=(w-wbar)/(sigw*sq2)
         sumg=sigw*0.5*sq2*sqpi*(1.+erf(z))
         fold=0
         wold=0
         gold=0


         do n=1,nx

 100        wnuc=w+wdiab


            if(n.eq.1)then
               new=.true.
            else
               new=.false.
            endif

	    if(razzak)then
               alw=alpha*wnuc
               sqrtalw=sqrt(alw)
               do m=1,nmode
                  zeta(m)=2.*sqrtalw*a/(3.*sqrtg(m))
                  print *,'zeta=',zeta(m)
                  eta(m)=2*alw*sqrtalw/(na(m)*beta(m)*
     &                   sqrtg(m))
                  print *,'mode=',m,' eta=',eta
               enddo
               if(interior)then
                   smax=alw/(beta(m)*hydrorad)
               else
                   call activem(zeta,eta,nmode,sm,alnsign,smax)
               endif

               arg=sm(1)/smax
               arg=arg*arg
               arg=amax1(arg,1.e-10)
               arg=amin1(arg,1.e10)
               x=alog(arg)/(3*sq2*alnsign(1))
               fnew=0.5*(1.-erf(x))

	    elseif(explicit_activate)then


	      call explact(wnuc,smax)

              lnsmax=alog(smax)

              x=2*(lnsm(nmode)-lnsmax)/(3*sq2*alnsign(nmode))
              fnew=0.5*(1.-erf(x))

	    else
               print *,'must use either razzak or explicit_activate'
               stop
            endif

            if(fnew-fold.gt.dfmax)then
               dw=0.5*dw
               w=wold+dw
               go to 100
            endif

            if(fnew-fold.lt.dfmin)then
               dw=amin1(2*dw,dwmax)
            endif
            fold=fnew

            z=(w-wbar)/(sigw*sq2)
            g=exp(-z*z)

            fnmin=1.

            do m=1,nmode
               if(razzak.or.explicit_activate)then
                 arg=sm(m)/smax
                 arg=arg*arg
                 arg=amax1(arg,1.e-10)
                 arg=amin1(arg,1.e10)
                 x=alog(arg)/(3*sq2*alnsign(m))
                 fn(m)=0.5*(1.-erf(x))
                 fm(m)=0.5*(1.-erf(x-1.4*sq2*alnsign(m)))
	       endif
               fnmin=amin1(fn(m),fnmin)
c              integration is second order accurate
c              assumes linear variation of f*g with w
               wb=(w+wold)
               fnbar=(fn(m)*g+fnold(m)*gold)
               fmbar=(fm(m)*g+fmold(m)*gold)
               if(w.gt.0.)then
                  sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar
     &                +(fn(m)*g*w+fnold(m)*gold*wold))*dw
                  sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar
     &                +(fm(m)*g*w+fmold(m)*gold*wold))*dw
               endif
               sumfn(m)=sumfn(m)+0.5*fnbar*dw
               sumfm(m)=sumfm(m)+0.5*fmbar*dw
               fnold(m)=fn(m)
               fmold(m)=fm(m)
            enddo

c           print *,'w=',w,' p(w)=',g,' f=',(fn(m),m=1,nmode)

            sumg=sumg+0.5*(g+gold)*dw
            gold=g
            wold=w
            w=w+dw

            if(n.gt.1.and.(w.gt.wmax.or.fnmin.gt.fmax))go to 20

         enddo

         write(*,*) 'nx is too small in activate'
         write(*,*) 'w=',w,' wmax=',wmax
         stop

      else

         n=1
         w=wbar
         wnuc=w+wdiab
         new=.true.

         if(razzak)then
            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            do m=1,nmode
               zeta(m)=2.*sqrtalw*a/(3.*sqrtg(m))
               eta(m)=2*alw*sqrtalw/(na(m)*beta(m)*sqrtg(m))
            enddo

            if(interior)then
                   smax=alw/(beta(m)*hydrorad)
            else
               call activem(zeta,eta,nmode,sm,alnsign,smax)
            endif

	 elseif(explicit_activate)then

        write (100, "(30a)") 'press,temp,supersat,activen,time'
	    call explact(wnuc,smax)

	 endif

         do m=1,nmode
             if(razzak.or.explicit_activate)then
                 arg=sm(m)/smax
                 arg=arg*arg
                 arg=amax1(arg,1.e-10)
                 arg=amin1(arg,1.e10)
                 x=alog(arg)/(3*sq2*alnsign(m))
                 fn(m)=0.5*(1.-erf(x))
                 fm(m)=0.5*(1.-erf(x-1.5*sq2*alnsign(m)))
             endif
             fluxn(m)=fn(m)*w*na(m)
             fluxm(m)=fm(m)*w
	     print *,'smax=',smax
         enddo

         go to 30

      endif

   20 continue

      if(.not.top)then

c        contribution from all updraft stronger than wmax
c        assuming constant f (close to fmax)
         wnuc=w+wdiab
         new=.false.

         if(razzak)then
            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            do m=1,nmode
               zeta(m)=2.*sqrtalw*a/(3.*sqrtg(m))
               eta(m)=2*alw*sqrtalw/(na(m)*beta(m)*sqrtg(m))
            enddo

            if(interior)then
                   smax=alw/(beta(m)*hydrorad)
            else
               call activem(zeta,eta,nmode,sm,alnsign,smax)
            endif

	 elseif(explicit_activate)then

	    call explact(wnuc,smax)

	 endif
         z=(w-wbar)/(sigw*sq2)
         g=exp(-z*z)
         integ=sigw*0.5*sq2*sqpi*(1.-erf(z))
            sumg=sumg+integ

         do m=1,nmode
            if(razzak.or.explicit_activate)then
                 arg=sm(m)/smax
                 arg=arg*arg
                 arg=amax1(arg,1.e-10)
                 arg=amin1(arg,1.e10)
                 x=alog(arg)/(3*sq2*alnsign(m))
                 fn(m)=0.5*(1.-erf(x))
                 fm(m)=0.5*(1.-erf(x-1.5*sq2*alnsign(m)))
            endif
            sumflxn(m)=sumflxn(m)+(wbar*integ+sigw*sigw*g)*fn(m)
            sumflxm(m)=sumflxm(m)+(wbar*integ+sigw*sigw*g)*fm(m)
            sumfn(m)=sumfn(m)+fn(m)*integ
            sumfm(m)=sumfm(m)+fm(m)*integ
         enddo

      endif

c      write(*,9010)z,w,fn,g,integ
c 9010    format(5f12.4)

      do m=1,nmode
         fn(m)=sumfn(m)/(sq2*sqpi*sigw)
         fm(m)=sumfm(m)/(sq2*sqpi*sigw)
         fluxn(m)=sumflxn(m)*na(m)/(sq2*sqpi*sigw)
         fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
      enddo

   30 continue
      write(* ,9030) n, wmax
      write(* ,9035)
c      write(23,9030) n, wmax
c      write(23,9035)
         totfn=0
         totn=0
         totfm=0
         totm=0
      do m=1,nmode
c         write(*,*) 'mode',m,' fn=',fn(m),' fm=',fm(m),
c     &        ' fluxn',fluxn(m),' fluxm',fluxm(m),
c     &        ' from ',n,' bins'
         write(* ,9040) m, fn(m), fm(m), fluxn(m), fluxm(m)
c         write(23,9040) m, fn(m), fm(m), fluxn(m), fluxm(m)
             totfn=totfn+fn(m)*na(m)
             totn=totn+na(m)
             do n=1,ntype
                totfm=totfm+fm(m)*mass(n,m)
                totm=totm+mass(n,m)
             enddo
      enddo
      totfn=totfn/totn
      totfm=totfm/totm
      write (*, *)'totfn=',totfn,' totfm=',totfm
c      write (23,*)'totfn=',totfn,' totfm=',totfm

9030  format( / 'nbins=', i5, 5x, 'wmax=', f13.5 )
9035  format( 'mode', 7x, 'fn', 6x, 7x, 'fm', 6x,5x,
     &  'fluxn', 5x, 4x, 'fluxm/m' )
9040  format( i2, 4(f15.7) )

      print *,'sumg=',sumg
      print *,'sq2*sqpi*sigw=',sq2*sqpi*sigw

      return
      end


c-----------------------------------------------------------------------
c      function erf(x)
c      erf=r_erf(x)
c      return
c      end


c-----------------------------------------------------------------------
      subroutine activem(zeta,eta,nmode,sm,alnsign,smax)

c     calculates size of smallest particle activated for multiple
c     competing aerosol modes.
c     Abdul-Razzak and Ghan, JGR 2000.

      real eta(nmode),! Razzak
     &     sm(nmode), ! Razzak
     &     zeta(nmode), ! Razzak
     &     alnsign(nmode)! ln sigma
      real smax ! maximum supersaturation
      real twothird,sum
      data twothird/0.666666666/
      save twothird
      integer m  ! mode index

      print *,'zeta(1)=',zeta(1)
      print *,'eta(1)=',eta(1)

      do m=1,nmode
c        if(zeta(1).gt.1.e5*eta(1).or.sm(1)*sm(1).gt.1.e5*eta(1))then
         if(zeta(1).gt.1.e7*eta(1).or.sm(1)*sm(1).gt.1.e8*eta(1))then
c           weak forcing. essentially none activated
               smax=10.
         else
c           significant activation of this mode. calc activation all modes.
            go to 1
         endif
      enddo

      return

  1   continue

      sum=0
      do m=1,nmode

            sum=sum+0.5*exp(2.50*alnsign(m)*alnsign(m))/(sm(m)*sm(m))*
     &        (zeta(m)/eta(m))**1.5
     &       +(1.+0.25*alnsign(m))*(sm(m)*sm(m)/(eta(m)+3*zeta(m)))**0.75
     &         /(sm(m)*sm(m))

      enddo

      smax=1./sqrt(sum)

      return


      return
      end


c-----------------------------------------------------------------------
	subroutine explactstate(nmode,totaln,sigma,rmean,
     &                        	hygro,tt,pp)

	parameter (nbin=200)
	parameter (nbinp=nbin+1)

	real*4 totaln(nmode),sigma(nmode),rmean(nmode)
	real*4 hygro(nmode) ! Hygroscopity

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &                denw,delt,delv

        real*8 latentht,g1,g2,d1,d2,cond1,cond2,alpha,beta,cons,ai
	common /explactcoeffs/ latentht,g1,g2,d1,d2,cond1,cond2,
     &        	alpha,beta,cons,ai

        real*8 temp,press,vel,conc(nbin),bn(nbin),rdry(nbin),rwet(nbin)
        real*8 scrit(nbin),ycrit(nbin)
        real*8 wn,an
	common /explactstates/temp,press,vel,conc,bn,rdry,rwet,
     &         	scrit,ycrit,wn,an

	real*8 sum,s,ps,u(nbin),umin(nbin),umax(nbin),ub(nbinp)
	real*4 tt,pp,vol,vn,sqrt2,two
	data two/2./
	save two

c
        temp=tt
        press=pp
	vol=1.0d0
	s=0.0d0

	nbinmode=nbin/nmode
	if(nbinmode*nmode.ne.nbin)then
	   print *,
     &     'error in explactstate: nmode must divide evenly into nbin'
	   print *,'nmode=',nmode,' nbin=',nbin
	   stop
	endif
c
c        calculating dry radius, initial radius, critical radius,
c        no. of moles, and concentrations.
c
c        find boundaries of bins, assuming same number in each bin


c           assume lognormal distribution each mode
           dx=2./nbinmode
           do i=1,nbinmode+1
              xb=1.-(i-1)*dx
              ub(i)=erfinv(xb)
           enddo
           print *,'ub=', ub(1)
c           ub(1)=4.452

c           center of bin

           do i=1,nbinmode
              u(i)=0.5*(ub(i+1)+ub(i))
           enddo

	istart=0
	sumw=0.0d0
        sqrt2=sqrt(two)
        sumn=0.

	do 100 m=1,nmode
          write(*,'(a,f10.2,f10.6,f5.1,f5.3,a,i1)')
     &      'totaln,rmean,sigma,hygro=',
     &         totaln(m),rmean(m),sigma(m),hygro(m),' for mode',m
	   do i=1,nbinmode
	      ii=i+istart
	         conc(ii)=totaln(m)/nbinmode ! same conc in each bin
	      conc(ii) = max( 0.0d0, conc(ii) )
              sumn=sumn+conc(ii)
	   enddo

	   do i=1,nbinmode
	      ii=i+istart
 	         rdry(ii)=rmean(m)*dexp(u(i)*sigma(m)*sqrt2)
	      bn(ii)=hygro(m)*((rdry(ii))**3)
	      ycrit(ii)=dsqrt(3.0d0*bn(ii)/ai)
c             initialize assuming equilibrium size at supersat=0.
	      rwet(ii)=dsqrt(bn(ii)/ai)
	      scrit(ii)=2.0d0*ai/(3.0d0*ycrit(ii))
	      sumw=sumw+conc(ii)*rwet(ii)**3
	   enddo
	   istart=istart+nbinmode
  100   continue

	wn=4.*pi*denw*sumw/(3.0d+06*fmw)
	tempc=temp-273.15d0
	a1=-21.2929759d0
	a2=-7629.01437d0
	a3=-177.013938d0
	a4=10702.3730d0
	a5=1.81207195d0
	a6=2479.85588d0
	a7=-23.1047743d0
	ps=a1+a2/(a3+tempc+a4/(a5+tempc+a6/(a7+tempc)))
	vn=ps*vol*(s+1.0d0)/(ugc*temp*10000.0d0)
	an=press*vol/(ugc*temp*10000.0d0)-vn
c
c	write (60,303)
 303    format ('  class    # part      dry r(mic)',
     1  '   init. r(mic)   crit. r(mic)   crit. s'/)
	do i=1,nbin
	   rdrymic=rdry(i)*10000.0d0
	   ycritmic=ycrit(i)*10000.0d0
	   ymic=rwet(i)*10000.0d0
c           write (60,304) i,conc(i),rdrymic,ymic,ycritmic,scrit(i)
        enddo
 304	format (1x,i3,d14.4,f12.4,3d14.4)

	return
	end
	subroutine explact(vv,smax)

	parameter (nbin=200)
	parameter (nbinp=nbin+1)

c	external fex,jex

        real*8 temp,press,vel,conc(nbin),bn(nbin),rdry(nbin),rwet(nbin)
        real*8 scrit(nbin),ycrit(nbin)
        real*8 wn,an
	common /explactstates/temp,press,vel,conc,bn,rdry,rwet,
     &         	scrit,ycrit,wn,an

	real*8 y(nbinp),smaxdp
	real*4 smax

        integer lrw,liw
        parameter (lrw=2000000,liw=2000)
        real*8 rwork(lrw)
        real*8 atol,rtol,t,tout,dt,smin
	integer iwork(liw),itol
	integer nstepmax,naftermax
	parameter (nstepmax=1000)
	real*4 smaxhist(nstepmax)
	real*8 	activen
	real*4 vv
        real tt,pp

        vel=vv
	iopt=0
	mf=22
	istate=1
	itask=1
	itol=1
	rtol=1.0d-6
	atol=1.0d-10
	t = 0.0d0
	tout = 0.1d0
	dt = 1.0d0
        timef = 500.0d0
	activen=0.0d0
	smaxdp=0.0
	if(vel.le.0.01)then
	   smax=1.e-10
	   return
	endif
	neq=nbin+1

	do i=1,neq-1
	   y(i)=rwet(i)
	enddo
	y(neq)=0. ! initial supersaturation

c	write(60,'(a,f5.1,a,f8.5,a,144f6.3))')'t=',t,' s=',y(neq),
c     &  	' r=',(y(i)*1.e4,i=1,neq-1)

	nstep=0
	naftermax=0

 20     call lsode (fex,(/ neq /),y,t,tout,itol,(/ rtol /),(/ atol /),
     &               itask,istate,iopt,rwork,lrw,iwork,liw,jex,mf)


	if (istate .lt. 0) go to 40

	nstep=nstep+1

	call newstate (neq,y,activen,dt,smaxdp,tout)

        tt=temp
        pp=press
        call explactcoeff(tt,pp)

c	write(60,'(a,f5.1,a,f8.5,a,144f6.3))')'t=',t,' s=',y(neq),
c     &  	' r=',(y(i)*1.e4,i=1,neq-1)

c 	write (6,60) t,y(neq),activen
 60     format (4h t =,f8.2,6h   s =,e12.4,7h  ccn =,f12.4)

	smin=0.01*smaxdp
	smaxhist(nstep)=smaxdp
	if(nstep.gt.nstepmax) go to 51
	if (y(neq) .lt. smin) go to 51
	if (y(neq) .lt. smaxdp) naftermax=naftermax+1
	if (naftermax.gt.10) go to 51
	nstep2=nstep/2
	if(nstep.gt.1)then
 	   if(smaxdp.lt.1.01*smaxhist(nstep2))go to 51 ! smax is increasing slowly, so stop
 	endif
	tout=tout+dt
	if (tout .gt. timef) go to 51

	go to 20

 40     write (6,50) istate
 50     format (///22h error halt.. istate= ,i3)
c
 51     smax=smaxdp
c        write(6,*) 'nstep=',nstep
        write(6,61)vel,nstep,smax,activen
 61     format('v=',f5.1,' nsteps=',i3,' smax=',f8.6,' ccn=',f6.1)

        return
        end


	subroutine fex (neq,y,ydot)

	parameter (nbin=200)
	parameter (nbinp=nbin+1)

        real*8 latentht,g1,g2,d1,d2,cond1,cond2,alpha,beta,cons,ai
	common /explactcoeffs/ latentht,g1,g2,d1,d2,cond1,cond2,
     c 	        alpha,beta,cons,ai

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &                denw,delt,delv


        real*8 temp,press,vel,conc(nbin),bn(nbin),rdry(nbin),rwet(nbin)
        real*8 scrit(nbin),ycrit(nbin)
        real*8 wn,an
	common /explactstates/temp,press,vel,conc,bn,rdry,rwet,
     &         	scrit,ycrit,wn,an

	real*8 y(nbinp),ydot(nbinp)

        real*8 d,g,cond,rsum,ycube,d12,cond12

	rsum=0.0d0
	d12=d1*d2
	cond12=cond1*cond2/c1
	do  i=1,neq-1
	   ycube=y(i)*y(i)*y(i)
	   d=d1/(y(i)/(y(i)+delv)+d12/(y(i)*c2))
	   cond=cond1/(y(i)/(y(i)+delt)+cond12/y(i))
	   g=1.0/(g1/d+g2/cond)
	   ydot(i)=g*(y(neq)-ai/y(i)+bn(i)/ycube)/y(i)
	   rsum=rsum+g*conc(i)*(y(neq)-ai/y(i)+bn(i)/ycube)*y(i)
        enddo

        ydot(neq)=alpha*vel-beta*rsum

	return
	end

	subroutine jex (neq,y,pd,nrpd)

c        not used

	parameter (nbin=200)
	parameter (nbinp=nbin+1)

	real*8 pd(nrpd,nbinp),y(nbinp)

        real*8 latentht,g1,g2,d1,d2,cond1,cond2,alpha,beta,cons,ai
	common /explactcoeffs/ latentht,g1,g2,d1,d2,cond1,cond2,
     &       	        alpha,beta,cons,ai

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &                denw,delt,delv

        real*8 temp,press,vel,conc(nbin),bn(nbin),rdry(nbin),rwet(nbin)
        real*8 scrit(nbin),ycrit(nbin)
        real*8 wn,an
	common /explactstates/temp,press,vel,conc,bn,rdry,rwet,
     &         	scrit,ycrit,wn,an

	real*8 d,cond,g,rsum

	rsum=0.0d0
	do 10 k = 1,neq-1
	   d=d1/(y(k)/(y(k)+delv)+d1*d2/(y(k)*c2))
	   cond=cond1/(y(k)/(y(k)+delt)+cond1*cond2/(y(k)*c1))
	   g=1.0/(g1/d+g2/cond)
	   pd(k,k)=g*(y(k)*(ai/y(K)**2 - 3.0d0*bn(k)/y(k)**4)-(y(neq)
     &             - ai/y(k)+bn(k)/y(k)**3))/(y(k))**2
	   pd(neq,k) = -g*beta*conc(k)*(y(neq)-2.0d0*bn(k)/y(k)**3)
	   rsum = rsum+conc(k)*y(k)
 10     continue
	do 20 kk=1,neq-1
	   pd(kk,neq) = g/y(kk)
 20     continue
	pd(neq,neq) = -g*beta*rsum

	return
	end

	subroutine explactcoeff(temp,press)

	implicit none
	real*4 temp,press

        real*8 latentht,g1,g2,d1,d2,cond1,cond2,alpha,beta,cons,ai
	common /explactcoeffs/ latentht,g1,g2,d1,d2,cond1,cond2,
     &        	alpha,beta,cons,ai

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &             denw,delt,delv

        real*8 tempc,a1,a2,a3,a4,a5,a6,a7,ps,gt,betat,vol

	vol=1.0d0
	ai=3.30d-05/temp
	ai=3.25d-05/temp
	tempc=temp-273.15d0
	latentht = (597.89513d0-0.553821d0*tempc)*4187.0d0
	d1 = (1013.250/press)*(0.20896284d0+0.001585105d0*
     &             tempc+0.0000030892d0*(tempc)**2)
	a1=-21.2929759d0
	a2=-7629.01437d0
	a3=-177.013938d0
	a4=10702.3730d0
	a5=1.81207195d0
	a6=2479.85588d0
	a7=-23.1047743d0
	ps=a1+a2/(a3+tempc+a4/(a5+tempc+a6/(a7+tempc)))
	cond1= 0.024d0+8.0d-5*tempc
	cond2=dsqrt(2.0*pi*ugc*temp/fma)/(cp*press)
	d2=dsqrt(2.0*pi*fmw/(ugc*temp))/100.0d0
	g1=denw*ugc*temp/(100.0d0*ps*fmw)
	g2=latentht*denw*(latentht*fmw/(ugc*temp)- 1.0d0)/(10000.0d0*temp)
	gt=1.0/(g1/d1+g2/cond1)
	alpha = (grav*fmw*latentht/(ugc*cp*temp**2)-
     &               grav*fma/(ugc*temp))/10000.0d0
	beta=4.0*pi*denw*(ugc*temp/(100.0d0*ps*fmw) +
     &             fmw*latentht**2/(100.0d0*cp*press*fma*temp))/(vol)
	betat=beta*gt

	return
	end

	subroutine explactconstant

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &                denw,delt,delv

        pi=3.14159
	cp = 1005.45d0
	denw = 1000.0d0
	grav=980.0d0
	ugc=8.31432d0
	fmw=0.0180153d0
	fma=0.0289644
	c1 = .96d0 ! thermal accom. coeff.
	c2 =1.0d0 ! condensation coeff.
	c3=1.40d0
	delt = 2.16d-05 ! thermal jump length - delta t - (cm)
	delv = 1.096d-05 ! vapor jump length - delta t - (cm)

	return
	end

        function erfinv(x)
	real c,t,e,dt
	integer n
	real erfinv
        real*4 erf
        c=2./sqrt(3.14159)
        t=0.
        n=0
  10    continue
        e=erf(t)
        dt=-(e-x)/(c*exp(-t*t))
        t=t+dt
        if(abs(dt).lt.1.e-4)go to 20
        n=n+1
        go to 10
  20    continue
c        print *,'x,t,e,dt=',x,t,e,dt
        erfinv=t
        return
        end

	subroutine newstate (neq,y,activen,dt,smax,t)

	parameter (nbin=200)
	parameter (nbinp=nbin+1)

	real*8 y(nbinp)
	real*8 dt,activen
	real*8 smax,t

	real*8 pi,cp,grav,fmw,ugc,fma,c1,c2,c3,denw,delt,delv
	common /explactconstants/pi,cp,grav,fmw,ugc,fma,c1,c2,c3,
     &                denw,delt,delv

        real*8 temp,press,vel,conc(nbin),bn(nbin),rdry(nbin),rwet(nbin)
        real*8 scrit(nbin),ycrit(nbin)
        real*8 wn,an
	common /explactstates/temp,press,vel,conc,bn,rdry,rwet,
     &         	scrit,ycrit,wn,an

        real*8 latentht,g1,g2,d1,d2,cond1,cond2,alpha,beta,cons,ai
	common /explactcoeffs/ latentht,g1,g2,d1,d2,cond1,cond2,
     &        	alpha,beta,cons,ai
        real*8 sum,wno,presso,tempo

9850   format(f12.8,',',f12.8,',',e12.6,',',f12.6,','f10.4)

	wno=wn
	wn=0.0d0
	sum=0.d0
	do i=1,neq-1
	   sum=sum+conc(i)*y(i)**3
        enddo
        wn=4.0*pi*denw*sum/(3.0d+06*fmw)
        presso=press
	tempo=temp
        term1=vel*grav*dt/(cp*10000.0d0)
        term2=latentht*fmw*(wn-wno)/(cp*fma*an)
	temp=tempo-term1+term2
	press=presso*dexp(-grav*fma*vel*dt/(ugc*temp*10000.0d0))
	if (y(neq) .gt. smax) then
           activen=0.0
	   do i=1,neq-1
	      if (y(neq) .ge. scrit(i)) activen=activen+conc(i)
           enddo
	   smax=y(neq)
           print *,'press,temp,supersat,activen=',
     *              press,temp,y(neq),activen
           write (100, 9850) press,temp,y(neq),activen,t
	   return
	else
	   return
	endif

	end


c*****************************************************************************


      subroutine lsode (f, neq, y, t, tout, itol, rtol, atol, itask,
     &                  istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      real*8 y, t, tout, rtol, atol, rwork
      dimension neq(1), y(*), rtol(1), atol(1), rwork(lrw), iwork(liw)
c      external prepj, solsy
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     &         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, iflag, imxer, kgo, lf0,
     &         leniw, lenrw, lenwm, ml, mord, mu, mxhnl0, mxstp0
      real*8 rowns,ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,
     &         tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0
c     &         d1mach, vnorm
      dimension mord(2)
      logical ihit
c-----------------------------------------------------------------------
c the following internal common block contains
c (a) variables which are local to any subroutine but whose values must
c      be preserved between calls to the routine (own variables), and
c (b) variables which are communicated between subroutines.
c the structure of the block is as follows..  all real*4 variables are
c listed first, followed by all integers.  within each type, the
c variables are grouped with those local to subroutine lsode first,
c then those local to subroutine stode, and finally those used
c for communication.  the block is declared in subroutines
c lsode, intdy, stode, prepj, and solsy.  groups of variables are
c replaced by dummy arrays in the common declarations in routines
c where those variables are not used.
c-----------------------------------------------------------------------
      common /ls0001/ rowns(209),
     &         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
     &         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     &         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),
     &         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

      character*60 msg
c
      data  mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/
c-----------------------------------------------------------------------
c block a.
c this code block is executed on every call.
c it tests istate and itask for legality and branches appropriately.
c if istate .gt. 1 but the flag init shows that initialization has
c not yet been done, an error return occurs.
c if istate = 1 and tout = t, jump to block g and return immediately.
c-----------------------------------------------------------------------
      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0
c-----------------------------------------------------------------------
c block b.
c the next code block is executed for the initial call (istate = 1),
c or for a continuation call with parameter changes (istate = 3).
c it contains checking of all inputs and various initializations.
c
c first check legality of the non-optional inputs neq, itol, iopt,
c mf, ml, and mu.
c-----------------------------------------------------------------------
      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      meth = mf/10
      miter = mf - 10*meth
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 5) go to 608
      if (miter .le. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue
c next process and check the optional inputs. --------------------------
      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0d0
      hmxi = 0.0d0
      hmin = 0.0d0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0d0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0d0) go to 615
      hmxi = 0.0d0
      if (hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0d0) go to 616
c-----------------------------------------------------------------------
c set work array pointers and check lengths lrw and liw.
c pointers to segments of rwork and iwork are named by prefixing l to
c the name of the segment.  e.g., the segment yh starts at rwork(lyh).
c segments of rwork (in order) are denoted  yh, wm, ewt, savf, acor.
c-----------------------------------------------------------------------
 60   lyh = 21
      if (istate .eq. 1) nyh = n
      lwm = lyh + (maxord + 1)*nyh
      if (miter .eq. 0) lenwm = 0
      if (miter .eq. 1 .or. miter .eq. 2) lenwm = n*n + 2
      if (miter .eq. 3) lenwm = n + 2
      if (miter .ge. 4) lenwm = (2*ml + mu + 1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      if (miter .eq. 0 .or. miter .eq. 3) leniw = 20
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
c check rtol and atol for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
	if (itol .ge. 3) rtoli = rtol(i)
	if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
	if (rtoli .lt. 0.0d0) go to 619
	if (atoli .lt. 0.0d0) go to 620
 70     continue
      if (istate .eq. 1) go to 100
c if istate = 3, set flag to signal parameter changes to stode. --------
      jstart = -1
      if (nq .le. maxord) go to 90
c maxord was reduced below nq.  copy yh(*,maxord+2) into savf. ---------
      do 80 i = 1,n
 80     rwork(i+lsavf-1) = rwork(i+lwm-1)
c reload wm(1) = rwork(lwm), since lwm may have changed. ---------------
 90   if (miter .gt. 0) rwork(lwm) = dsqrt(uround)
      if (n .eq. nyh) go to 200
c neq was reduced.  zero part of yh to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0d0
      go to 200
c-----------------------------------------------------------------------
c block c.
c the next block is for the initial call only (istate = 1).
c it contains all remaining initializations, the initial call to f,
c and the calculation of the initial step size.
c the error weights in ewt are inverted after being loaded.
c-----------------------------------------------------------------------
 100  uround = d1mach(4)
      tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0d0) go to 625
      if (h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0)h0 = tcrit - t
 110  jstart = 0
      if (miter .gt. 0) rwork(lwm) = dsqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0d0
      nqu = 0
      ccmax = 0.3d0
      maxcor = 3
      msbp = 20
      mxncf = 10
c initial call to f.  (lf0 points to yh(*,2).) -------------------------
      lf0 = lyh + nyh
      call f (neq, y, rwork(lf0))
      nfe = 1
c load the initial value vector in yh. ---------------------------------
      do 115 i = 1,n
 115    rwork(i+lyh-1) = y(i)
c load and invert the ewt array.  (h is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0d0
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 120 i = 1,n
	if (rwork(i+lewt-1) .le. 0.0d0) go to 621
 120    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
c-----------------------------------------------------------------------
c the coding below computes the step size, h0, to be attempted on the
c first step, unless the user has supplied a value for this.
c first check that tout - t differs significantly from zero.
c a scalar tolerance quantity tol is computed, as max(rtol(i))
c if this is positive, or max(atol(i)/abs(y(i))) otherwise, adjusted
c so as to be between 100*uround and 1.0e-3.
c then the computed value h0 is given by..
c                                       neq
c   h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
c                                        1
c where   w0     = max ( abs(t), abs(tout) ),
c          f(i)   = i-th component of initial value of f,
c          ywt(i) = ewt(i)/tol  (a weight for y(i)).
c the sign of h0 is inferred from the initial values of tout and t.
c-----------------------------------------------------------------------
      if (h0 .ne. 0.0d0) go to 180
      tdist = dabs(tout - t)
      w0 = dmax1(dabs(t),dabs(tout))
      if (tdist .lt. 2.0d0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = dmax1(tol,rtol(i))
 140  if (tol .gt. 0.0d0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
	if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
	ayi = dabs(y(i))
	if (ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
 150    continue
 160  tol = dmax1(tol,100.0d0*uround)
      tol = dmin1(tol,0.001d0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      sum = 1.0d0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0d0/dsqrt(sum)
      h0 = dmin1(h0,tdist)
      h0 = dsign(h0,tout-t)
c adjust h0 if necessary to meet hmax bound. ---------------------------
 180  rh = dabs(h0)*hmxi
      if (rh .gt. 1.0d0) h0 = h0/rh
c load h with h0 and scale yh(*,2) by h0. ------------------------------
      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270
c-----------------------------------------------------------------------
c block d.
c the next code block is for continuation calls only (istate = 2 or 3)
c and is to check stop conditions before taking a step.
c-----------------------------------------------------------------------
 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0d0 + 100.0d0*uround)
      if ((tp - tout)*h .gt. 0.0d0) go to 623
      if ((tn - tout)*h .lt. 0.0d0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
      if ((tcrit - tout)*h .lt. 0.0d0) go to 625
      if ((tn - tout)*h .lt. 0.0d0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0d0) go to 624
 245  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      if (istate .eq. 2) jstart = -2
c-----------------------------------------------------------------------
c block e.
c the next block is normally executed for all calls and contains
c the call to the one-step core integrator stode.
c
c this is a looping point for the integration steps.
c
c first check for too many steps being taken, update ewt (if not at
c start of problem), check for too much accuracy being requested, and
c check for h below the roundoff level in t.
c-----------------------------------------------------------------------
 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
	if (rwork(i+lewt-1) .le. 0.0d0) go to 510
 260    rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0d0) go to 280
      tolsf = tolsf*2.0d0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      msg='lsode--  warning..internal t (=r1) and h (=r2) are'
      call xerrwv(msg,50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='  such that in the machine, t + h = t on the next step  '
      call xerrwv(msg,60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      (h = step size). solver will continue anyway'
      call xerrwv(msg,50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      msg='lsode--  above warning has been issued i1 times.  '
      call xerrwv(msg,50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      it will not be issued again for this problem'
      call xerrwv(msg,50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
 290  continue
c-----------------------------------------------------------------------
c      call stode(neq,y,yh,nyh,yh,ewt,savf,acor,wm,iwm,f,jac,prepj,solsy)
c-----------------------------------------------------------------------
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),
     1   rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm),
     2   f, jac, prepj, solsy)
      kgo = 1 - kflag
      go to (300, 530, 540), kgo
c--------------------------------------------------------------------------
c block f.
c the following block handles the case of a successful return from the
c core integrator (kflag = 0).  test for stop conditions.
c--------------------------------------------------------------------------
 300  init = 1
      go to (310, 400, 330, 340, 350), itask
c itask = 1.  if tout has been reached, interpolate. -------------------
 310  if ((tn - tout)*h .lt. 0.0d0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
c itask = 3.  jump to exit if tout was reached. ------------------------
 330  if ((tn - tout)*h .ge. 0.0d0) go to 400
      go to 250
c itask = 4.  see if tout or tcrit was reached.  adjust h if necessary.
 340  if ((tn - tout)*h .lt. 0.0d0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0d0 + 4.0d0*uround)
      if ((tnext - tcrit)*h .le. 0.0d0) go to 250
      h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
      jstart = -2
      go to 250
c itask = 5.  see if tcrit was reached and jump to exit. ---------------
 350  hmx = dabs(tn) + dabs(h)
      ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
c--------------------------------------------------------------------------
c block g.
c the following block handles all successful returns from lsode.
c if itask .ne. 1, y is loaded from yh and t is set accordingly.
c istate is set to 2, the illegal input counter is zeroed, and the
c optional outputs are loaded into the work arrays before returning.
c if istate = 1 and tout = t, there is a return with no action taken,
c except that if this has happened repeatedly, the run is terminated.
c--------------------------------------------------------------------------
 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c
 430  ntrep = ntrep + 1
      if (ntrep .lt. 5) return
      msg='lsode--  repeated calls with istate = 1 and tout = t (=r1)'
      call xerrwv(msg, 60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
      go to 800
c--------------------------------------------------------------------------
c block h.
c the following block handles all unsuccessful returns other than
c those for illegal input.  first the error message routine is called.
c if there was an error test or convergence test failure, imxer is set.
c then y is loaded from yh, t is set to tn, and the illegal input
c counter illin is set to 0.  the optional outputs are loaded into
c the work arrays before returning.
c--------------------------------------------------------------------------
c the maximum number of steps was taken before reaching tout. ----------
 500  msg='lsode--  at current t (=r1), mxstep (=i1) steps   '
      call xerrwv(msg,50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      taken on this call before reaching tout     '
      call xerrwv(msg, 50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
      istate = -1
      go to 580
c ewt(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  ewti = rwork(lewt+i-1)
      msg='lsode--  at t (=r1), ewt(i1) has become r2 .le. 0.'
      call xerrwv(msg,50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580
c too much accuracy requested for machine precision. -------------------
 520  msg='lsode--  at t (=r1), too much accuracy requested  '
      call xerrwv(msg, 50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      for precision of machine..  see tolsf (=r2) '
      call xerrwv(msg, 50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580
c kflag = -1.  error test failed repeatedly or with abs(h) = hmin. -----
 530  msg='lsode--  at t(=r1) and step size h(=r2), the error'
      call xerrwv(msg,50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      test failed repeatedly or with abs(h) = hmin'
      call xerrwv(msg, 50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560
c kflag = -2.  convergence failed repeatedly or with abs(h) = hmin. ----
 540  msg='lsode--  at t (=r1) and step size h (=r2), the    '
      call xerrwv(msg,50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      corrector convergence failed repeatedly     '
      call xerrwv(msg,50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      or with abs(h) = hmin   '
      call xerrwv(msg, 30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
c compute imxer if relevant. -------------------------------------------
 560  big = 0.0d0
      imxer = 1
      do 570 i = 1,n
	size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
	if (big .ge. size) go to 570
	big = size
	imxer = i
 570    continue
      iwork(16) = imxer
c set y vector, t, illin, and optional outputs. ------------------------
 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      return
c--------------------------------------------------------------------------
c block i.
c the following block handles all error returns due to illegal input
c (istate = -3), as detected before calling the core integrator.
c first the error message routine is called.  then if there have been
c 5 consecutive such returns just before this call to the solver,
c the run is halted.
c--------------------------------------------------------------------------
 601  msg='lsode--  istate (=i1) illegal '
      call xerrwv(msg,30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
      go to 700
 602  msg='lsode--  itask (=i1) illegal  '
      call xerrwv(msg, 30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
      go to 700
 603  msg='lsode--  istate .gt. 1 but lsode not initialized  '
      call xerrwv(msg, 50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      go to 700
 604  msg='lsode--  neq (=i1) .lt. 1     '
      call xerrwv(msg, 30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
      go to 700
 605  msg='lsode--  istate = 3 and neq increased (i1 to i2)  '
      call xerrwv(msg, 50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 606  msg='lsode--  itol (=i1) illegal   '
      call xerrwv(msg, 30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
      go to 700
 607  msg='lsode--  iopt (=i1) illegal   '
      call xerrwv(msg, 30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
      go to 700
 608  msg='lsode--  mf (=i1) illegal     '
      call xerrwv(msg, 30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
      go to 700
 609  msg='lsode--  ml (=i1) illegal.. .lt.0 or .ge.neq (=i2)'
      call xerrwv(msg, 50, 9, 0, 2, ml, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 610  msg='lsode--  mu (=i1) illegal.. .lt.0 or .ge.neq (=i2)'
      call xerrwv(msg, 50, 10, 0, 2, mu, neq(1), 0, 0.0d0, 0.0d0)
      go to 700
 611  msg='lsode--  maxord (=i1) .lt. 0  '
      call xerrwv(msg, 30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
      go to 700
 612  msg='lsode--  mxstep (=i1) .lt. 0  '
      call xerrwv(msg, 30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
      go to 700
 613  msg='lsode--  mxhnil (=i1) .lt. 0  '
      call xerrwv(msg, 30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
      go to 700
 614  msg='lsode--  tout (=r1) behind t (=r2)      '
      call xerrwv(msg, 40, 14, 0, 0, 0, 0, 2, tout, t)
      msg='      integration direction is given by h0 (=r1)  '
      call xerrwv(msg,50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
      go to 700
 615  msg='lsode--  hmax (=r1) .lt. 0.0  '
      call xerrwv(msg,30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
      go to 700
 616  msg='lsode--  hmin (=r1) .lt. 0.0 '
      call xerrwv(msg ,30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
      go to 700
 617  msg='lsode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)'
      call xerrwv(msg,60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
      go to 700
 618  msg='lsode--  iwork length needed, leniw (=i1), exceeds liw (=i2)'
      call xerrwv(msg,60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
      go to 700
 619  msg='lsode--  rtol(i1) is r1 .lt. 0.0        '
      call xerrwv(msg,40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
      go to 700
 620  msg='lsode--  atol(i1) is r1 .lt. 0.0        '
      call xerrwv(msg,40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      msg='lsode--  ewt(i1) is r1 .le. 0.0         '
      call xerrwv(msg,40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
      go to 700
 622  msg='lsode--  tout (=r1) too close to t(=r2) to start integration'
      call xerrwv(msg,60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  msg='lsode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  '
      call xerrwv(msg,60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  msg='lsode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   '
      call xerrwv(msg,60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  msg='lsode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   '
      call xerrwv(msg,60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  msg='lsode--  at start of problem, too much accuracy   '
      call xerrwv(msg,50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
      msg='      requested for precision of machine..  see tolsf (=r1) '
      call xerrwv(msg,60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
      rwork(14) = tolsf
      go to 700
 627  msg='lsode--  trouble from intdy. itask = i1, tout = r1'
      call xerrwv(msg,50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
c
 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      return
 710  msg='lsode--  repeated occurrences of illegal input    '
      call xerrwv(msg,50, 302, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
c
 800  msg='lsode--  run aborted.. apparent infinite loop     '
      call xerrwv(msg,50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
      return
c-------------------------- end of subroutine lsode -----------------------
      end

      subroutine intdy (t, k, yh, nyh, dky, iflag)
clll. optimize
      integer k, nyh, iflag
      integer iownd, iowns,
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real*8 t, yh, dky
      real*8 rowns,ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 c, r, s, tp
      character*60 msg
      dimension yh(nyh,*), dky(*)
      common /ls0001/ rowns(209), ccmax,
     &    el0, h, hmin, hmxi, hu, rc, tn, uround,iownd(14), iowns(6),
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0d0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0d0) go to 90
c
      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = dble(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
	j = nq - jb
	jp1 = j + 1
	ic = 1
	if (k .eq. 0) go to 35
	jj1 = jp1 - k
	do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = dble(ic)
	do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return
c
 80   msg='intdy--  k (=i1) illegal      '
      call xerrwv(msg,30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
      iflag = -1
      return
 90   msg='intdy--  t (=r1) illegal      '
      call xerrwv(msg, 30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
      msg='  t not in interval tcur - hu (= r1) to tcur (=r2)    '
      call xerrwv(  msg, 60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return
c-------------------------- end of subroutine intdy -----------------------
      end
      subroutine stode (neq, y, yh, nyh, yh1, ewt, savf, acor,
     &    wm, iwm, f, jac, pjac, slvs)
clll. optimize
      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      real*8 y, yh, yh1, ewt, savf, acor, wm
      real*8 conit, crate, el, elco, hold, rmax, tesco,
     &    ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,
c     &    r, rh, rhdn, rhsm, rhup, told, vnorm
     &    r, rh, rhdn, rhsm, rhup, told
       dimension neq(1), y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),
     &    acor(*), wm(*), iwm(*)
      common /ls0001/ conit, crate, el(13), elco(13,12),
     &    hold, rmax, tesco(3,12),
     &    ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),
     &    ialth, ipup, lmax, meo, nqnyh, nslp,
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,       &
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0d0
      rc = 0.0d0
      el0 = 1.0d0
      crate = 0.7d0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140
 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dble(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0d0/dble(l)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
      rh = dmin1(rhdn,1.0d0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = dmin1(rh,dabs(h/hold))
      h = hold
      go to 175
 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dble(nq+2)
      go to (160, 170, 200), iret
 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = dmax1(rh,hmin/dabs(h))
 175  rh = dmin1(rh,rmax)
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)
      r = 1.0d0
      do 180 j = 2,l
	r = r*rh
	do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690
 200  if (dabs(rc-1.0d0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
	i1 = i1 - nyh
cdir$ ivdep
	do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue
 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, y, savf)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250
      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0d0
      nslp = nst
      crate = 0.7d0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0d0
 270  if (miter .ne. 0) go to 350
      do 290 i = 1,n
	savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
	y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400
 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
	acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)
 400  if (m .ne. 0) crate = dmax1(0.2d0*crate,del/delp)
      dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0d0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0d0*delp) go to 410
      delp = del
      call f (neq, y, savf)
      nfe = nfe + 1
      go to 270
 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0d0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
	i1 = i1 - nyh
cdir$ ivdep
	do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (dabs(h) .le. hmin*1.00001d0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25d0
      ipup = miter
      iredo = 1
      go to 170
 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0d0) go to 500
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
	do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700
 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
	i1 = i1 - nyh
cdir$ ivdep
	do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0d0
      if (dabs(h) .le. hmin*1.00001d0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0d0
      go to 540
 520  rhup = 0.0d0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0d0/dble(l+1)
      rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
 540  exsm = 1.0d0/dble(l)
      rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
      rhdn = 0.0d0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0d0/dble(nq)
      rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0d0) rh = 1.0d0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1d0) go to 610
      r = el(l)/dble(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1d0)) go to 610
      if (kflag .le. -2) rh = dmin1(rh,0.2d0)
      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150
 640  if (kflag .eq. -10) go to 660
      rh = 0.1d0
      rh = dmax1(hmin/dabs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, y, savf)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150
c--------------------------------------------------------------------------
c all returns are made through this section.  h is saved in hold
c to allow the caller to change h on the next step.
c--------------------------------------------------------------------------
 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0d0
 700  r = 1.0d0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return
c-------------------------- end of subroutine stode -----------------------
      end
      subroutine cfode (meth, elco, tesco)
clll. optimize
      integer meth
      integer i, ib, nq, nqm1, nqp1
      real*8 elco, tesco
      real*8 agamq, fnq, fnqm1, pc, pint, ragq,
     1   rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
      dimension pc(12)
c
      go to (100, 200), meth
c
 100  elco(1,1) = 1.0d0
      elco(2,1) = 1.0d0
      tesco(1,1) = 0.0d0
      tesco(2,1) = 2.0d0
      tesco(1,2) = 1.0d0
      tesco(3,12) = 0.0d0
      pc(1) = 1.0d0
      rqfac = 1.0d0
      do 140 nq = 2,12
	rq1fac = rqfac
	rqfac = rqfac/dble(nq)
	nqm1 = nq - 1
	fnqm1 = dble(nqm1)
	nqp1 = nq + 1
c form coefficients of p(x)*(x+nq-1). ----------------------------------
	pc(nq) = 0.0d0
	do 110 ib = 1,nqm1
	  i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
	pc(1) = fnqm1*pc(1)
c compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
	pint = pc(1)
	xpin = pc(1)/2.0d0
	tsign = 1.0d0
	do 120 i = 2,nq
	  tsign = -tsign
	  pint = pint + tsign*pc(i)/dble(i)
 120      xpin = xpin + tsign*pc(i)/dble(i+1)
c store coefficients in elco and tesco. --------------------------------
	elco(1,nq) = pint*rq1fac
	elco(2,nq) = 1.0d0
	do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/dble(i)
	agamq = rqfac*xpin
	ragq = 1.0d0/agamq
	tesco(2,nq) = ragq
	if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/dble(nqp1)
	tesco(3,nqm1) = ragq
 140    continue
      return
c
 200  pc(1) = 1.0d0
      rq1fac = 1.0d0
      do 230 nq = 1,5
	fnq = dble(nq)
	nqp1 = nq + 1
c form coefficients of p(x)*(x+nq). ------------------------------------
	pc(nqp1) = 0.0d0
	do 210 ib = 1,nq
	  i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
	pc(1) = fnq*pc(1)
c store coefficients in elco and tesco. --------------------------------
	do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
	elco(2,nq) = 1.0d0
	tesco(1,nq) = rq1fac
	tesco(2,nq) = dble(nqp1)/elco(1,nq)
	tesco(3,nq) = dble(nq+2)/elco(1,nq)
	rq1fac = rq1fac/fnq
 230    continue
      return
c-------------------------- end of subroutine cfode -----------------------
      end
      subroutine prepj (neq, y, yh, nyh, ewt, ftem, savf, wm,iwm,f,jac)
clll. optimize
      external f, jac
      integer neq, nyh, iwm
      integer iownd, iowns,
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, i2, ier, ii, j, j1, jj, lenp,
     &    mba, mband, meb1, meband, ml, ml3, mu, np1
      real*8 y, yh, ewt, ftem, savf, wm
      real*8 rowns,ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
c      real*8 con, di, fac, hl0, r, r0, srur, yi, yj, yjj,vnorm
      real*8 con, di, fac, hl0, r, r0, srur, yi, yj, yjj
       dimension neq(1), y(*), yh(nyh,*), ewt(*), ftem(*), savf(*),
     &    wm(*), iwm(*)
      common /ls0001/ rowns(209), ccmax,
     &    el0, h, hmin, hmxi, hu, rc, tn, uround,iownd(14), iowns(6),
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0d0
      call jac (neq, y, wm(3), n)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnorm (n, savf, ewt)
      r0 = 1000.0d0*dabs(h)*uround*dble(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
	yj = y(j)
	r = dmax1(srur*dabs(yj),r0/ewt(j))
	y(j) = y(j) + r
	fac = -hl0/r
	call f (neq, y, ftem)
	do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
	y(j) = yj
	j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      np1 = n + 1
      do 250 i = 1,n
	wm(j) = wm(j) + 1.0d0
 250    j = j + np1
c do lu decomposition on p. --------------------------------------------
      call dgefa (wm(3), n, n, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(2) = hl0
      r = el0*0.1d0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, y, wm(3))
      nfe = nfe + 1
      do 320 i = 1,n
	r0 = h*savf(i) - yh(i,2)
	di = 0.1d0*r0 - h*(wm(i+2) - savf(i))
	wm(i+2) = 1.0d0
	if (dabs(r0) .lt. uround/ewt(i)) go to 320
	if (dabs(di) .eq. 0.0d0) go to 330
	wm(i+2) = 0.1d0*r0/di
 320    continue
      return
 330  ierpj = 1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0d0
      call jac (neq, y, wm(ml3), meband)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min0(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnorm (n, savf, ewt)
      r0 = 1000.0d0*dabs(h)*uround*dble(n)*fac
      if (r0 .eq. 0.0d0) r0 = 1.0d0
      do 560 j = 1,mba
	do 530 i = j,n,mband
	  yi = y(i)
	  r = dmax1(srur*dabs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
	call f (neq, y, ftem)
	do 550 jj = j,n,mband
	  y(jj) = yh(jj,1)
	  yjj = y(jj)
	  r = dmax1(srur*dabs(yjj),r0/ewt(jj))
	  fac = -hl0/r
	  i1 = max0(jj-mu,1)
	  i2 = min0(jj+ml,n)
	  ii = jj*meb1 - ml + 2
	  do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
	wm(ii) = wm(ii) + 1.0d0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call dgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      if (ier .ne. 0) ierpj = 1
      return
c-------------------------- end of subroutine prepj -----------------------
      end
      subroutine solsy (wm, iwm, x, tem)
clll. optimize
      integer iwm
      integer iownd, iowns,
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, meband, ml, mu
      real*8 wm, x, tem
      real*8 rowns,ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real*8 di, hl0, phl0, r
      dimension wm(*), iwm(*), x(*), tem(1)
      common /ls0001/ rowns(209),ccmax,
     &    el0, h, hmin, hmxi, hu, rc, tn, uround,iownd(14), iowns(6),
     &    icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     &    maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      iersl = 0
      go to (100, 100, 300, 400, 400), miter
 100  call dgesl (wm(3), n, n, iwm(21), x, 0)
      return
c
 300  phl0 = wm(2)
      hl0 = h*el0
      wm(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
	di = 1.0d0 - r*(1.0d0 - 1.0d0/wm(i+2))
	if (dabs(di) .eq. 0.0d0) go to 390
 320    wm(i+2) = 1.0d0/di
 330  do 340 i = 1,n
 340    x(i) = wm(i+2)*x(i)
      return
 390  iersl = 1
      return
c
 400  ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      call dgbsl (wm(3), meband, n, ml, mu, iwm(21), x, 0)
      return
c-------------------------- end of subroutine solsy -----------------------
      end
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
clll. optimize
      integer n, itol
      integer i
      real*8 rtol, atol, ycur, ewt
      dimension rtol(1), atol(1), ycur(n), ewt(n)
c
      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
      return
c-------------------------- end of subroutine ewset -----------------------
      end
      real*8 function vnorm (n, v, w)
clll. optimize
      integer n,   i
      real*8 v, w,   sum
      dimension v(n), w(n)
      sum = 0.0d0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      vnorm = dsqrt(sum/dble(n))
      return
c-------------------------- end of function vnorm -------------------------
      end
      subroutine srcom (rsav, isav, job)
      integer isav, job
      integer ieh, ils
      integer i, lenils, lenrls
      real*8 rsav,   rls
      dimension rsav(*), isav(*)
      common /ls0001/ rls(218), ils(39)
      common /eh0001/ ieh(2)
      data lenrls/218/, lenils/39/
c
      if (job .eq. 2) go to 100
c
      do 10 i = 1,lenrls
 10     rsav(i) = rls(i)
      do 20 i = 1,lenils
 20     isav(i) = ils(i)
      isav(lenils+1) = ieh(1)
      isav(lenils+2) = ieh(2)
      return
c
 100  continue
      do 110 i = 1,lenrls
 110     rls(i) = rsav(i)
      do 120 i = 1,lenils
 120     ils(i) = isav(i)
      ieh(1) = isav(lenils+1)
      ieh(2) = isav(lenils+2)
      return
c-------------------------- end of subroutine srcom -----------------------
      end
      real*8 function d1mach (idum)
      integer idum
      real*8 u, comp
      u = 1.0d0
 10   u = u*0.5d0
      comp = 1.0d0 + u
      if (comp .ne. 1.0d0) go to 10
      d1mach = u*2.0d0
      return
c-------------------------- end of function d1mach ------------------------
      end
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      integer nmes, nerr, level, ni, i1, i2, nr,
     &    i, lun, lunit, mesflg, ncpw, nch, nwds
      real*8 r1, r2
      character*60 msg
      common /eh0001/ mesflg, lunit
      data ncpw/4/
      if (mesflg .eq. 0) go to 100
c get logical unit number. ---------------------------------------------
      lun = lunit
c get number of words in message. --------------------------------------
      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1
c write the message. ---------------------------------------------------
      write (lun, '(a)') (msg)
  10  format(1x,15a4)
c the following is valid for ncpw = 2.
c 10  format(1x,30a2)
c--------------------------------------------------------------------------
      if (ni .eq. 1) write (lun, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)
      if (ni .eq. 2) write (lun, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
      if (nr .eq. 1) write (lun, 40) r1
 40   format(6x,23hin above message,  r1 =,d21.13)
      if (nr .eq. 2) write (lun, 50) r1,r2
 50   format(6x,15hin above,  r1 =,d21.13,3x,4hr2 =,d21.13)
c abort the run if level = 2. ------------------------------------------
 100  if (level .ne. 2) return
      stop
c-------------------------- end of subroutine xerrwv ----------------------
      end
      subroutine xsetf (mflag)
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c-------------------------- end of subroutine xsetf -----------------------
      end
      subroutine xsetun (lun)
c
c this routine resets the logical unit number for messages.
c
      integer lun, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (lun .gt. 0) lunit = lun
      return
c-------------------------- end of subroutine xsetun ----------------------
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      real*8 a(lda,*)
      real*8 t
      integer idamax,j,k,kp1,l,nm1
c
c
c      gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
	 kp1 = k + 1
c
c         find l = pivot index
c
	 l = idamaxfn(n-k+1,a(k,k),1) + k - 1
	 ipvt(k) = l
c
c         zero pivot implies this column already triangularized
c
	 if (a(l,k) .eq. 0.0d0) go to 40
c
c            interchange if necessary
c
	    if (l .eq. k) go to 10
	       t = a(l,k)
	       a(l,k) = a(k,k)
	       a(k,k) = t
   10       continue
c
c            compute multipliers
c
	    t = -1.0d0/a(k,k)
	    call dscal(n-k,t,a(k+1,k),1)
c
c            row elimination with column indexing
c
	    do 30 j = kp1, n
	       t = a(l,j)
	       if (l .eq. k) go to 20
		  a(l,j) = a(k,j)
		  a(k,j) = t
   20          continue
	       call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
	 go to 50
   40    continue
	    info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(*),job
      real*8 a(lda,*),b(*)
c      real*8 ddot,t
      real*8 t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c         job = 0 , solve  a * x = b
c         first solve  l*y = b
c
	 if (nm1 .lt. 1) go to 30
	 do 20 k = 1, nm1
	    l = ipvt(k)
	    t = b(l)
	    if (l .eq. k) go to 10
	       b(l) = b(k)
	       b(k) = t
   10       continue
	    call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c         now solve  u*x = y
c
	 do 40 kb = 1, n
	    k = n + 1 - kb
	    b(k) = b(k)/a(k,k)
	    t = -b(k)
	    call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c         job = nonzero, solve  trans(a) * x = b
c         first solve  trans(u)*y = b
c
	 do 60 k = 1, n
	    t = ddot(k-1,a(1,k),1,b(1),1)
	    b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c         now solve trans(l)*x = y
c
	 if (nm1 .lt. 1) go to 90
	 do 80 kb = 1, nm1
	    k = n - kb
	    b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
	    l = ipvt(k)
	    if (l .eq. k) go to 70
	       t = b(l)
	       b(l) = b(k)
	       b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(*),info
      real*8 abd(lda,*)
      real*8 t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c      zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
	 i0 = m + 1 - jz
	 do 10 i = i0, ml
	    abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c      gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
	 kp1 = k + 1
c
c         zero next fill-in column
c
	 jz = jz + 1
	 if (jz .gt. n) go to 50
	 if (ml .lt. 1) go to 50
	    do 40 i = 1, ml
	       abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c         find l = pivot index
c
	 lm = min0(ml,n-k)
	 l = idamaxfn(lm+1,abd(m,k),1) + m - 1
	 ipvt(k) = l + k - m
c
c         zero pivot implies this column already triangularized
c
	 if (abd(l,k) .eq. 0.0d0) go to 100
c
c            interchange if necessary
c
	    if (l .eq. m) go to 60
	       t = abd(l,k)
	       abd(l,k) = abd(m,k)
	       abd(m,k) = t
   60       continue
c
c            compute multipliers
c
	    t = -1.0d0/abd(m,k)
	    call dscal(lm,t,abd(m+1,k),1)
c
c            row elimination with column indexing
c
	    ju = min0(max0(ju,mu+ipvt(k)),n)
	    mm = m
	    if (ju .lt. kp1) go to 90
	    do 80 j = kp1, ju
	       l = l - 1
	       mm = mm - 1
	       t = abd(l,j)
	       if (l .eq. mm) go to 70
		  abd(l,j) = abd(mm,j)
		  abd(mm,j) = t
   70          continue
	       call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
	 go to 110
  100    continue
	    info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(*),job
      real*8 abd(lda,*),b(1)
c      real*8 ddot,t
      real*8 t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c         job = 0 , solve  a * x = b
c         first solve l*y = b
c
	 if (ml .eq. 0) go to 30
	 if (nm1 .lt. 1) go to 30
	    do 20 k = 1, nm1
	       lm = min0(ml,n-k)
	       l = ipvt(k)
	       t = b(l)
	       if (l .eq. k) go to 10
		  b(l) = b(k)
		  b(k) = t
   10          continue
	       call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c         now solve  u*x = y
c
	 do 40 kb = 1, n
	    k = n + 1 - kb
	    b(k) = b(k)/abd(m,k)
	    lm = min0(k,m) - 1
	    la = m - lm
	    lb = k - lm
	    t = -b(k)
	    call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c         job = nonzero, solve  trans(a) * x = b
c         first solve  trans(u)*y = b
c
	 do 60 k = 1, n
	    lm = min0(k,m) - 1
	    la = m - lm
	    lb = k - lm
	    t = ddot(lm,abd(la,k),1,b(lb),1)
	    b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c         now solve trans(l)*x = y
c
	 if (ml .eq. 0) go to 90
	 if (nm1 .lt. 1) go to 90
	    do 80 kb = 1, nm1
	       k = n - kb
	       lm = min0(ml,n-k)
	       b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
	       l = ipvt(k)
	       if (l .eq. k) go to 70
		  t = b(l)
		  b(l) = b(k)
		  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c      constant times a vector plus a vector.
c      uses unrolled loops for increments equal to one.
c      jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c         code for unequal increments or equal increments
c           not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dy(iy) = dy(iy) + da*dx(ix)
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c         code for both increments equal to 1
c
c
c         clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
	dy(i) = dy(i) + da*dx(i)
	dy(i + 1) = dy(i + 1) + da*dx(i + 1)
	dy(i + 2) = dy(i + 2) + da*dx(i + 2)
	dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine dscal(n,da,dx,incx)
c
c      scales a vector by a constant.
c      uses unrolled loops for increment equal to one.
c      jack dongarra, linpack, 3/11/78.
c
      real*8 da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c         code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
	dx(i) = da*dx(i)
   10 continue
      return
c
c         code for increment equal to 1
c
c
c         clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
	dx(i) = da*dx(i)
	dx(i + 1) = da*dx(i + 1)
	dx(i + 2) = da*dx(i + 2)
	dx(i + 3) = da*dx(i + 3)
	dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      real*8 function ddot(n,dx,incx,dy,incy)
c
c      forms the dot product of two vectors.
c      uses unrolled loops for increments equal to one.
c      jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c         code for unequal increments or equal increments
c           not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	dtemp = dtemp + dx(ix)*dy(iy)
	ix = ix + incx
	iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c         code for both increments equal to 1
c
c
c         clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
	dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      integer function idamaxfn(n,dx,incx)
c
c      finds the index of element having max. absolute value.
c      jack dongarra, linpack, 3/11/78.
c
      real*8 dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c         code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
	 if(dabs(dx(ix)).le.dmax) go to 5
	 idamax = i
	 dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c         code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
	 if(dabs(dx(i)).le.dmax) go to 30
	 idamax = i
	 dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine dcopy(n,sx,incx,sy,incy)
c
c      copies a vector, x, to a vector, y.
c      uses unrolled loops for increments equal to 1.
c      jack dongarra, linpack, 3/11/78.
c
      real*8 sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c         code for unequal increments or equal increments
c           not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
	sy(iy) = sx(ix)
	ix = ix + incx
	iy = iy + incy
   10 continue
      return
c
c         code for both increments equal to 1
c
c
c         clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
	sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
	sy(i) = sx(i)
	sy(i + 1) = sx(i + 1)
	sy(i + 2) = sx(i + 2)
	sy(i + 3) = sx(i + 3)
	sy(i + 4) = sx(i + 4)
	sy(i + 5) = sx(i + 5)
	sy(i + 6) = sx(i + 6)
   50 continue
      return
      end

      end module
