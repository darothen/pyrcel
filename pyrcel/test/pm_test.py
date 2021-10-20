#!/usr/env python

import time

import matplotlib.pyplot as plt

import pyrcel as pm
from pyrcel.postprocess import simulation_activation

P0 = 95000.0  # Pressure, Pa
T0 = 290.15  # Temperature, K
S0 = -0.01  # Supersaturation, 1-RH
V = 1.0
accom = 0.1

mu = 0.025
sigma = 1.82
kappa = 0.54
N = 1500.0

aerosol_distribution = pm.Lognorm(mu=mu, sigma=sigma, N=N)
aerosol = pm.AerosolSpecies(
    "test", aerosol_distribution, kappa=kappa, bins=200
)

output_dt = 1.0  # simulation seconds; how frequently should output be saved?
solver_dt = (
    10.0
)  # simulation seconds; how frequently should the solver be re-set?
# Why this change? (1) Adding ice is going to introduce a
# time-splitting operation, so the integration logic will need
# to change accordingly and be more sequential, step-by-step.
# (2) The solver works *much* better in terms of adapting to the
# stiffness of the ODE when it doesn't have to try to predict very
# far in advance where the solution will shoot.
#
# In general, use solver_dt = 10.*output_dt
dZ = 1000.0
t_end = dZ / V

results = {}
initial_aerosols = [aerosol]

## Vanilla parcel model run
model = pm.ParcelModel(
    initial_aerosols, V, T0, S0, P0, console=True, accom=accom
)
# raw_input("Continue? ")
start = time.time()
parcel, aerosols = model.run(
    t_end,
    output_dt,
    solver_dt,  # note new argument order!
    terminate=True,
    terminate_depth=50.0,
    max_steps=2000,
    solver="cvode",
    output="dataframes",
)
end = time.time()

Smax = parcel.S.max()
tmax = parcel.S.argmax()

print("Elapsed time:", end - start)
print("   Smax", Smax)
print("   tmax", tmax)
print("")
print("Computing activation")
acts_total = simulation_activation(model, parcel, aerosols)

fig = plt.figure(2)
ax = fig.add_subplot(111)


def quick_plot(ax):
    plt_meta = parcel.plot("z", "S", zorder=3, ax=ax)
    ylims = ax.get_ylim()

    c = plt_meta.lines[-1].get_color()
    plt.vlines(parcel.z.ix[tmax], ylims[0], Smax, color=c, linestyle="dashed")
    plt.hlines(Smax, 0, parcel.z.iloc[tmax], color=c, linestyle="dashed")
    plt.ylim(S0)


quick_plot(ax)

model.save("test.nc", other_dfs=[acts_total])
