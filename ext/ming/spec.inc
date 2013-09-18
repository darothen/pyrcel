	integer, parameter :: TY = 1  !  Number of aerosol types
	integer, parameter :: SZ = 40 !  Number of size bins
	integer, parameter :: SP = 1 !  Number of species
	integer, parameter :: MD = 1 !  Number of log-normal dist.

	double precision :: T = 283.15 !  Temperature (K)
	double precision, parameter :: P = 80000. !  Pressure (Pa)
	double precision, parameter :: R = 8.314  !  Gas constant (J mol-1 K-1)
	double precision, parameter :: ZERO = 273.15 !  Zero degree C (K)
	double precision, parameter :: ATM = 1.01325e5 !  Standard atmosphere pressure (Pa)
	double precision, parameter :: PI = 3.1415926
	double precision  Density (SP) /1.7418e3/ ! Density of species (Kg m-3)
	double precision  MolW (SP) /0.132/ ! Molecular weight of species (Kg mol-1)
	double precision  IonN (SP) /3./ ! Number of ions each molecule
	double precision  N(MD) /200./ ! Total Number Concen (cm-3)
	double precision  Dm(MD) /0.02/ ! Mean diameter (micron)
	double precision  Sig(MD) /2.0/ ! Sigma (St. Dev.)
	double precision  IonN1(MD) /3./ ! Number of ions
	double precision  MolW1(MD) /0.132/ ! Mol. weight (Kg mol-1)
	double precision  Density1(MD) /1.7418e3/ ! Density (Kg m-3)

