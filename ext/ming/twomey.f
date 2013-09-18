	module INIT
		contains
		subroutine Loading(Mass, DryDp)
			include 'spec.inc'
			double precision, dimension(TY, SZ, SP), intent(inout) :: Mass
			double precision, dimension(TY, SZ), intent(inout) :: DryDp
			double precision, dimension(SP+1) :: Temp
			character no
			character filename*15
			integer i, j, k

			filename='Mass.dat'
			open(11, FILE=filename)

			do i=1,SZ
				read(11,*) Temp
				DryDp(1,i)=Temp(1)
				do j=1,SP
					Mass(1,i,j)=Temp(j+1)
				end do
			end do
			close(11)
		end subroutine Loading

	end module INIT

	module CHEM
		contains
!		subroutine Kohler
!		subroutine Activity

		subroutine Kohler(DryDp, Dens, Mw, Ion, Dpc, Dp0, Sc, AA)
			include 'spec.inc'
			double precision, dimension(TY, SZ), intent(inout) :: DryDp, Dens, Mw, Ion, Dpc, Dp0, Sc
			double precision, intent(inout) :: AA
			double precision A, B
			integer count1, count2

			A = (0.66/T)*(76.1-0.155*(T-273.15))/76.1 !micron
!			A = (0.66/T)*76.1/76.1
			count1=0
			do while (count1<TY)
				count1=count1+1
				count2=0
				do while (count2<SZ)
					count2=count2+1
					B = 3.44e13*(PI/6.)*(DryDp(count1,count2)*1e-6)**3.  &
					*Dens(count1,count2)*Ion(count1, count2)/Mw(count1,count2) !micron^2
					Dpc(count1,count2)=(3.*B/A)**0.5
					Dp0(count1,count2)=(B/A)**0.5
!					Dpc(count1,count2)=((B/A)**0.5+(3.*B/A)**0.5)/2.
					Sc(count1,count2)=EXP((4.*A**3./(27.*B))**0.5)
!					Sc(count1,count2)=(4.*A**3./(27.*B))**0.5
!					print *, Dens(count1, count2)
				end do
			end do
			AA = A*1.e-6 ! m
		end subroutine Kohler

		subroutine CalDenMwIon(Mass, Dens, Mw, Ion)
			include 'spec.inc'
			double precision, dimension(TY, SZ, SP), intent(inout) :: Mass
			double precision, dimension(TY, SZ), intent(inout) :: Dens, Mw, Ion
			integer count1, count2, count3
			double precision sum, sum1, sum2, sum3

			count1=0
			do while (count1<TY)
				count1=count1+1
				count2=0
				do while (count2<SZ)
					count2=count2+1
					count3=0
					sum=0.
					sum2=0.
					sum3=0.
					do while (count3<SP)
						count3=count3+1
						sum=sum+Mass(count1,count2,count3)
						sum2=sum2+Mass(count1,count2,count3)/MolW(count3)
						sum3=sum3+Mass(count1,count2,count3)/MolW(count3)*IonN(count3)
					end do
					count3=0
					sum1=0.
					do while (count3<SP)
						count3=count3+1
						sum1=sum1+Mass(count1,count2,count3)/sum/Density(count3)
					end do
					Dens(count1,count2)=1./sum1
					Mw(count1,count2)=sum/sum2
					Ion(count1,count2)=sum3/sum2
				end do
			end do
		end subroutine CalDenMwIon


		subroutine CalNum(DryDp, Mass, Dens, TMass, Num)
			include 'spec.inc'
			double precision, dimension(TY, SZ, SP), intent(inout) :: Mass
			double precision, dimension(TY, SZ), intent(inout) :: DryDp, Dens, TMass, Num
			integer count1, count2, count3
			double precision sum, sum1, sum2, sum3

			do count1 = 1, TY
1				do count2 = 1, SZ
					sum=0
					do count3 = 1, SP
						sum=sum+Mass(count1, count2, count3)
					end do
					TMass(count1, count2)=sum
					Num(count1, count2)=sum*1e-3/Dens(count1, count2)/(PI/6.*  &
					(DryDp(count1, count2)*1.e-6)**3.)  !cm-3
				end do
			end do
		end subroutine CalNum




		subroutine CalcG(Dp, G)
			include 'spec.inc'
			double precision, intent(inout) :: Dp
			double precision, intent(inout) :: G
			double precision :: rhow = 1.0e3  ! density of water (Kg m-3)
			double precision :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
			double precision :: alpc = 1.0  ! mass accomodation coef.
			double precision :: alpt = 0.97  ! thermal accomodation coef.
!			double precision :: alpc = 0.043  ! thermal accomodation coef.
			double precision :: delt = 2.16e-1 !thermal jump (micron)
			double precision :: delv = 1.096e-1 !vapor jump (micron)
			double precision vpres, Dv, ka, Le, mass, heat


			Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
			Dv = Dv/(Dp/(Dp+delv*2.)+2*Dv/(alpc*(Dp*1e-6))*(2.*PI*Mw/R/T)**0.5)
			ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
			ka = ka/(Dp/(Dp+delt*2.)+2*ka/(alpt*(Dp*1e-6)*1.007e3*P)*(2.*PI*R*T/0.028965)**0.5)
			TC = T-ZERO
			vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
			        +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure(Pa)
			Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)

			mass = rhow*R*T/(vpres*Dv*Mw)
			heat = Le*rhow/(ka*T)*(Le*Mw/T/R-1)
!			print *, Dv, vpres, Mw, rhow, mass, heat
			G = 4./(mass+heat) ! (m2 s-1)
		end subroutine CalcG

		subroutine CalcG2(G)
			include 'spec.inc'
!			double precision, intent(inout) :: Dp
			double precision, intent(inout) :: G
			double precision :: rhow = 1.0e3  ! density of water (Kg m-3)
			double precision :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
			double precision :: alpc = 0.2  ! mass accomodation coef.
			double precision :: alpt = 1.  ! thermal accomodation coef.
			double precision vpres, Dv, ka, Le, mass, heat

			Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!			Dv = Dv/(1+2*Dv/(alpc*(Dp*1e-6))*(2.*PI*Mw/R/T)**0.5)
			ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
!			ka = ka/(1+2*ka/(alpt*(Dp*1e-6)*1.007e3*P)*(2.*PI*R*T/0.028965)**0.5)
			TC = T-ZERO
			vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
			        +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure(Pa)
			Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)

			mass = rhow*R*T/(vpres*Dv*Mw)
			heat = Le*rhow/(ka*T)*(Le*Mw/T/R-1)
!			print *, Dv, vpres, Mw, rhow, mass, heat
			G = 4./(mass+heat) ! (m2 s-1)
		end subroutine CalcG2

		subroutine CalcAlphaGamma(alpha, gamma)
			include 'spec.inc'
			double precision, intent(inout) :: alpha, gamma
			double precision :: rhow = 1.0e3  ! density of water (Kg m-3)
			double precision rhoa ! density of air (Kg m-3)
			double precision :: Cpa = 1.007e3 ! specific heat of air (J Kg-1 K-1)
			double precision :: Mw = 0.018  ! molecular weight of water (Kg mol-1)
			double precision :: Ma = 0.028965  ! molecular weight of air (Kg mol-1)
			double precision :: alpc = 1.  ! mass accomodation coef.
			double precision :: g = 9.815 ! gravitational acceleration (m s-2)
			double precision vpres, Dv, ka, Le

			rhoa = P*Ma/R/T  ! (Kg m-3)
			Dv = 0.211/(P/ATM)*(T/ZERO)**1.94*1e-4  ! diffusivity of water vapor (m2 s-1)
!			Dv = Dv/(1+2*Dv/(alpc*Dp)*(2.*PI*Mw/R/T)**0.5)
			ka = 1e-3*(4.39+0.071*T)  ! thermal conductivity (J m-1 s-1 K-1)
			TC = T-ZERO
			vpres = (6.107799961+TC*(4.436518521e-1+TC*(1.428945805e-2+TC*(2.650648471e-4 &
			        +TC*(3.031240396e-6+TC*(2.034080948e-8+6.136820929e-11*TC))))))*1e2  ! saturated water vapor pressure (Pa)
			Le = 597.3*(ZERO/T)**(0.167+3.67e-4*T)*4.182*1e3  ! latent heat of water (J Kg -1)
			alpha = g*Mw*Le/(Cpa*R*T**2.)-g*Ma/(R*T) ! (m-1)
			gamma = R*T/(vpres*Mw)+Mw*Le**2./(Cpa*P*Ma*T) ! (m3 Kg-1)
		end subroutine CalcAlphaGamma

		subroutine Activity(DryDp, WetDp, Mw, Ion, Dens, Act)
			include 'spec.inc'
			double precision, dimension(TY, SZ), intent(inout) :: DryDp, WetDp, Mw, Ion, Dens, Act
			integer i, j
			double precision temp1, temp2
			do i=1,TY
				do j=1,SZ
					Act(i,j)=1.0
					if(WetDp(i,j)>1e-5) then
						temp1=DryDp(i,j)**3.*Dens(i,j)/Mw(i,j)
						temp2=(WetDp(i,j)**3.-DryDp(i,j)**3.)*1.e3/0.018
						Act(i,j)=temp2/(temp1+temp2)*exp(0.66/T/WetDp(i,j))
					endif
				end do
			end do
		end subroutine Activity

		subroutine f1(sm, x, y)
			double precision, intent(inout) :: sm, x, y
			y = x/2.*(sm**2.-x**2.)**0.5+sm**2./2.*asin(x/sm)
		end subroutine f1

	end module CHEM

	module PARA
		contains
		subroutine Grow(Smax,Updraft,alpha,Dpc,Dp0,Sc,WetDp)
			use CHEM
			include 'spec.inc'
			double precision, dimension(TY, SZ), intent(inout) :: Dpc,Dp0,Sc,WetDp
			double precision, intent(inout) :: Smax,Updraft,alpha
			double precision G, retard

			integer count1, count2
			retard=1.
			do count1=1,TY
				do count2=1,SZ
					WetDp(count1,count2)=0.
					call CalcG(Dpc(count1,count2), G)
					if(Smax>Sc(count1,count2)) then
WetDp(count1,count2)=(Dpc(count1,count2)**2.+1e12*G/(alpha*Updraft)*((Smax-1.)**2.5-(Sc(count1,count2)-1)**2.5))**0.5
!Twomey's formulation
!WetDp(count1,count2)=(1e12*G/(alpha*Updraft)*((Smax-1.)**2.4-(Sc(count1,count2)-1)**2.4))**0.5

					endif


				end do
			end do

		end subroutine Grow
		subroutine Conden(Smax,DryDp,WetDp,Dpc,Num,Act,CondenRate,DropletNum,ActDp)
			use CHEM
			include 'spec.inc'
			double precision, dimension(TY, SZ), intent(inout) :: DryDp,WetDp,Dpc,Num,Act
			double precision, intent(inout) :: Smax, CondenRate, DropletNum,ActDp
			double precision G
			integer count1, count2
			CondenRate=0.
			DropletNum=0.
			do count1=1,TY
				do count2=1,SZ
					call CalcG(WetDp(count1,count2), G)
					if (WetDp(count1,count2)>Dpc(count1,count2)) then
						CondenRate=CondenRate+PI/2.*1.e3*G*(WetDp(count1,count2)*1.e-6)*Num(count1,count2)*1.e6*  &
						(Smax-Act(count1,count2))
!						CondenRate=CondenRate+PI/2.*1.e3*G*(WetDp(count1,count2)*1.e-6)*Num(count1,count2)*1.e6*  &
!						(Smax-1.)
					endif
!						print *, WetDp(count1,count2), Dpc(count1,count2)
					if (WetDp(count1,count2)>Dpc(count1,count2)) then
						DropletNum=DropletNum+Num(count1,count2)
						if (WetDp(count1,count2-1)<Dpc(count1,count2)) then
							ActDp=DryDp(count1,count2)
						endif
					endif
				end do
			end do
		end subroutine Conden
		subroutine ActRatio(Smax,DryDp,Sc,Num,ratio,number,size)
			include 'spec.inc'
			double precision, dimension(TY, SZ), intent(inout) :: DryDp,Sc,Num
			double precision, intent(inout) :: Smax, ratio,number,size
			double precision tmp1, tmp2
			integer i,j

			tmp1=0.
			tmp2=0.
			do i=1, TY
				do j=1,SZ
					if(Sc(i,j)<Smax) then
						tmp1=tmp1+Num(i,j)
					endif
					tmp2=tmp2+Num(i,j)
					if(Sc(i,j-1)>Smax .and. Sc(i,j)<Smax) then
						size=DryDp(i,j)
!						tmp1=tmp1+Num(i,j-1)*(Smax-Sc(i,j))/(Sc(i,j-1)-Sc(i,j))
					endif
				end do
			end do
			ratio=tmp1/tmp2
			number=tmp1
		end subroutine ActRatio
	end module PARA

	program MAIN
		use INIT
		use CHEM
		use PARA
		include 'spec.inc'
		double precision, dimension(TY, SZ, SP) :: Mass
		double precision, dimension(TY, SZ) :: DryDp, WetDp, Dp0
		double precision, dimension(TY, SZ) :: Dens, Mw, Ion, Dpc, Sc, TMass, Num, Act
		double precision :: G, alpha, gamma, Smax, Updraft, CondenRate, DropletNum, ActDp
		double precision Smax1, Smax2, A, factor, Spart, tmp1, tmp2, tmp3
		integer i, count1, count2, m, iter_count
!		double precision, dimension(MD) ::
		double precision B, f, gg, Sm, ita, zeta, Dmi
		double precision up(1) /10./
		double precision ratio,number,size

		Updraft = 0.1

	do m=1,1

		Smax1 = 1.000
		Smax2 = 1.5
		Updraft = up(m)

		print *, "V=",Updraft

		call Loading(Mass, DryDp)
!		print *, DryDp
		call CalDenMwIon(Mass, Dens, Mw, Ion)
		call Kohler(DryDp, Dens, Mw, Ion, Dpc, Dp0, Sc, A)
!		print *, A
		call CalcAlphaGamma(alpha, gamma)
		call CalNum(DryDp, Mass, Dens, TMass, Num)
        Num=Num*1e-6
		call CalcG2(G)
		factor=16.*A**2.*alpha*Updraft/(9.*G)
		do while ((Smax2-Smax1)>1e-7)
			Smax=0.5*(Smax1+Smax2)
			if((Smax-1)**4.>=factor) then
				Spart=1.+(Smax-1.)*(0.5*(1.+(1-factor/(Smax-1)**4.)**0.5))**0.5
			else
				tmp3=2.e7*A/3.*(Smax-1.)**(-0.3824)
				if (tmp3>1.) then
					Spart=Smax
				else
					Spart=1.+(Smax-1.)*tmp3
				endif
			endif

			CondenRate=0.
			do count1=1,TY
				do count2=2,SZ-1
					call CalcG(DryDp(count1,count2), G)
!					call CalcG2(G)
					if(Sc(count1,count2)<Spart) then
						call f1(Smax,Sc(count1,count2),tmp1)
						call f1(Smax,Sc(count1,count2+1),tmp2)
		CondenRate=CondenRate+G**1.5/(alpha*Updraft)**0.5*Num(count1,count2)/(Sc(count1,count2)-Sc(count1,count2+1))  &
		*(tmp1-tmp2)
					else if(Sc(count1,count2)<Smax .and. Sc(count1,count2-1)<Smax) then
		CondenRate=CondenRate+G*(2./3.)*Num(count1,count2)*A/(Sc(count1,count2)-Sc(count1,count2+1))  &
		*log((Sc(count1,count2)-1.)/(Sc(count1,count2+1)-1.))
					else if(Sc(count1,count2)<Smax .and. Sc(count1,count2-1)>Smax) then
		CondenRate=CondenRate+G*(2./3.)*Num(count1,count2)*A/(Sc(count1,count2)-Sc(count1,count2+1))  &
		*log((Smax-1.)/(Sc(count1,count2+1)-1.))
					endif
				end do
			end do
			CondenRate=CondenRate*(PI/2.)*gamma*1.e3*(Smax-1.)/(alpha*Updraft)*1e6
			if(CondenRate>1.) then
				Smax2=Smax
			else
				Smax1=Smax
			endif
		end do

!		print *, Smax, CondenRate, alpha*Updraft/gamma, DropletNum,ActDp
		call ActRatio(Smax,DryDp,Sc,Num,ratio,number,size)
		print *, "nenes",(Smax-1), size,ratio
!		print *, (Smax-1), size,ratio

		Smax1 = 1.000
		Smax2 = 1.10

        iter_count = 1
		do while ((Smax2-Smax1)>1e-7)
!		do while (abs(CondenRate-alpha*Updraft/gamma)>1e-7)
			Smax=0.5*(Smax1+Smax2)
			call Grow(Smax,Updraft,alpha,Dpc,Dp0,Sc,WetDp)
			call Activity(DryDp, WetDp, Mw, Ion, Dens, Act)
			call Conden(Smax,DryDp,WetDp,Dpc,Num,Act,CondenRate,DropletNum,ActDp)
			if(CondenRate<(alpha*Updraft/gamma)) then
				Smax1=Smax
			else
				Smax2=Smax
			endif
			print *, iter_count, Smax, CondenRate, alpha*Updraft/gamma
            iter_count = iter_count+1
		end do
		call ActRatio(Smax,DryDp,Sc,Num,ratio,number,size)
			print *, "ming",(Smax-1),size,ratio
!			print *, Updraft,ratio

!CondenRate/(alpha*Updraft/gamma)

!		print *, asin(0.5)*180/3.14

		A=A/2.
		tmp1=0.
		do i=1,MD
            Dmi = (Dm(i)/2.)*1e-6
			call CalcG(Dm(i), G)
!			call CalcG2(G)
			G=G/4.
			zeta=2.*A/3.*(alpha*Updraft/G)**0.5
			ita=(alpha*Updraft/G)**1.5/(2.*PI*1.e3*gamma*N(i)*1e6)
			f=0.5*exp(2.5*log(Sig(i))**2.)
			gg=1+0.25*log(Sig(i))
			B=IonN1(i)*0.018*Density1(i)/(MolW1(i)*1.e3)
			Sm=(2./(B**0.5))*(A/3./Dmi)**1.5
			tmp1=tmp1+(1./(Sm**2.))*(f*(zeta/ita)**1.5+gg*(Sm**2./(ita+3.*zeta))**0.75)
		end do
		Smax=1.+1/(tmp1**0.5)
		call ActRatio(Smax,DryDp,Sc,Num,ratio,number,size)
		print *, "ghan",(Smax-1.), ratio,size


!		Updraft=Updraft*exp(log(10/0.1)/20)


	end do
		tmp1=0.
		do i=SZ,1,-1
			tmp1=tmp1+Num(1,i)
		end do
		tmp2=0.
		do i=SZ,1,-1
			tmp2=tmp2+Num(1,i)
!			print *, DryDp(i), Dpc(i), tmp2/tmp1
		end do

!		print *, tmp1



	end program MAIN
