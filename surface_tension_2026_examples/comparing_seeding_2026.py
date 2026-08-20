import pyrcel as pm

"""Identical initial conditions to run two parcel simulations"""
V = 1.0
T0 = 283.0    #temperature, K
S0 = -0.02    #initial supersaturation
P0 = 85000.0  #pressure, Pa

"""Identical aerosol properties where only sigma differs to isolate the effect of surface tension on activation"""
SPECIES_NAME = "seeding_agent"
DISTRIBUTION = pm.Lognorm(mu=0.02, sigma=2.0, N=1000.0)
KAPPA = 0.54
BINS = 50

SURFACTANT_SIGMA = 0.060  # J/m^2, vs water's ~0.076 at this temperature


"""Run 1: no sigma override, i.e. pure water surface tension for all bins"""
aerosol_no_surfactant = pm.AerosolSpecies(
    SPECIES_NAME, DISTRIBUTION, kappa=KAPPA, bins=BINS
)

model_no_surfactant = pm.ParcelModel(
    [aerosol_no_surfactant], V=V, T0=T0, S0=S0, P0=P0, console=False
)
out_no_surfactant = model_no_surfactant.run(t_end=300.0, output_dt=10.0, terminate=True)
 
print("=== NO SURFACTANT (pure water surface tension) ===")
print(f"S_max = {out_no_surfactant.summary['S_max'] * 100:.4f} %")
print(f"Nd    = {out_no_surfactant.Nd:.3e} m^-3")
print()

"""Run 2: sigma override set"""
aerosol_with_surfactant = pm.AerosolSpecies(
    SPECIES_NAME, DISTRIBUTION, kappa=KAPPA, sigma=SURFACTANT_SIGMA, bins=BINS
)
 
model_with_surfactant = pm.ParcelModel(
    [aerosol_with_surfactant], V=V, T0=T0, S0=S0, P0=P0, console=False
)
out_with_surfactant = model_with_surfactant.run(t_end=300.0, output_dt=10.0, terminate=True)
 
print("=== WITH SURFACTANT (sigma = %.4f J/m^2) ===" % SURFACTANT_SIGMA)
print(f"S_max = {out_with_surfactant.summary['S_max'] * 100:.4f} %")
print(f"Nd    = {out_with_surfactant.Nd:.3e} m^-3")
print()


"""Comparision of the two runs"""
print("=== COMPARISON (same population, only sigma changed) ===")
delta_smax = out_with_surfactant.summary["S_max"] - out_no_surfactant.summary["S_max"]
delta_nd = out_with_surfactant.Nd - out_no_surfactant.Nd
print(f"Change in S_max: {delta_smax * 100:+.4f} percentage points")
print(f"Change in Nd:    {delta_nd:+.3e} m^-3")
print(f"Surfactant lowered S_max needed:      {delta_smax < 0}")
print(f"Surfactant increased activated number: {delta_nd > 0}")
