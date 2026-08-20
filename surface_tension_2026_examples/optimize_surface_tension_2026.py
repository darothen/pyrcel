
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

import pyrcel as pm
from pyrcel.equilibrate import kohler_crit_approx

V = 1.0
T0 = 283.0
S0 = -0.02
P0 = 85000.0

SPECIES_NAME = "seeding_agent"
DISTRIBUTION = pm.Lognorm(mu=0.02, sigma=2.0, N=1000.0)
KAPPA = 0.54
BINS = 50

SIGMA_MIN = 0.020
SIGMA_MAX = 0.0744  

ND_WEIGHT = 1.0
SMAX_WEIGHT = 0.0

N_STARTS = 5
XATOL = 1e-7

COLOR_POINTS = "#457B9D"       
COLOR_OPTIMUM = "#E63946"     
COLOR_WATER = "#8D99AE"        
COLOR_ACTIVATED = "#2A9D8F"    
COLOR_HAZE = "#E9C46A"        

_log = []  


def run_model(sigma_value):
    aerosol = pm.AerosolSpecies(
        SPECIES_NAME, DISTRIBUTION, kappa=KAPPA, sigma=float(sigma_value), bins=BINS
    )
    model = pm.ParcelModel([aerosol], V=V, T0=T0, S0=S0, P0=P0, console=False)
    out = model.run(t_end=300.0, output_dt=10.0, terminate=True)
    s_max = float(out.summary["S_max"])
    nd = float(out.Nd)
    _log.append((float(sigma_value), s_max, nd))
    return s_max, nd


def run_model_full(sigma_value):
    """Returns everything needed for both the S(t) plot and the radius-vs-height plot."""
    aerosol = pm.AerosolSpecies(
        SPECIES_NAME, DISTRIBUTION, kappa=KAPPA, sigma=float(sigma_value), bins=BINS
    )
    model = pm.ParcelModel([aerosol], V=V, T0=T0, S0=S0, P0=P0, console=False)
    out = model.run(t_end=300.0, output_dt=10.0, terminate=True)
    return model, out


def classify_bins(model, out, sigma_value):
    aer = model.aerosols[0]
    T_smax = out.summary["T_smax"]
    S_max = out.summary["S_max"]
    _, s_crits = kohler_crit_approx(T_smax, aer.r_drys, aer.kappa, aer.sigma)
    activated_mask = np.asarray(S_max >= s_crits)

    offset = 7 
    radii = out.state[:, offset : offset + aer.nr]
    heights = out.heights
    return radii, heights, activated_mask


def objective(sigma_value):
    s_max, nd = run_model(sigma_value)
    return ND_WEIGHT * (nd / nd_ref) + SMAX_WEIGHT * (s_max / smax_ref)


def neg_objective(sigma_value):
    return -objective(sigma_value)


_, nd_lo_ref = run_model(SIGMA_MIN)
smax_hi_ref, nd_hi_ref = run_model(SIGMA_MAX)
nd_ref = max(abs(nd_lo_ref), abs(nd_hi_ref), 1e-30)
smax_ref = max(abs(smax_hi_ref), 1e-30)


print(f"Running {N_STARTS} independent bounded searches...")
edges = np.linspace(SIGMA_MIN, SIGMA_MAX, N_STARTS + 1)
candidates = []

for i in range(N_STARTS):
    lo, hi = edges[i], edges[i + 1]
    result = minimize_scalar(
        neg_objective, bounds=(lo, hi), method="bounded", options={"xatol": XATOL}
    )
    sigma_i = float(result.x)
    obj_i = -float(result.fun)
    candidates.append((sigma_i, obj_i))
    print(f"  [{lo:.4f}, {hi:.4f}] -> sigma={sigma_i:.6f}, objective={obj_i:.6f}")

optimal_sigma, optimal_obj = max(candidates, key=lambda c: c[1])
optimal_s_max, optimal_nd = run_model(optimal_sigma)

print(f"\nOptimal sigma: {optimal_sigma:.6f} J/m^2")
print(f"  -> S_max = {optimal_s_max * 100:.4f} %")
print(f"  -> Nd    = {optimal_nd:.3e} m^-3")
print(f"Total model evaluations used: {len(_log)}")


print("\nRunning full trajectories for plotting...")
model_water, out_water = run_model_full(SIGMA_MAX)
model_opt, out_opt = run_model_full(optimal_sigma)

radii_water, heights_water, activated_water = classify_bins(model_water, out_water, SIGMA_MAX)
radii_opt, heights_opt, activated_opt = classify_bins(model_opt, out_opt, optimal_sigma)

log_arr = np.array(sorted(_log, key=lambda r: r[0]))

fig, axes = plt.subplots(2, 3, figsize=(16, 9))

ax = axes[0, 0]
ax.scatter(log_arr[:, 0], log_arr[:, 2], color=COLOR_POINTS, s=35, alpha=0.8, edgecolors="none")
ax.axvline(optimal_sigma, color=COLOR_OPTIMUM, linestyle="--", label=f"optimum: {optimal_sigma:.4f}")
ax.set_xlabel("Surface tension, sigma (J/m^2)")
ax.set_ylabel("Activated droplet number, Nd (m^-3)")
ax.set_title("Nd at each evaluated sigma")
ax.legend(frameon=False)
ax.grid(alpha=0.25)

ax = axes[0, 1]
ax.scatter(log_arr[:, 0], log_arr[:, 1] * 100, color=COLOR_POINTS, s=35, alpha=0.8, edgecolors="none")
ax.axvline(optimal_sigma, color=COLOR_OPTIMUM, linestyle="--")
ax.set_xlabel("Surface tension, sigma (J/m^2)")
ax.set_ylabel("Peak supersaturation, S_max (%)")
ax.set_title("S_max at each evaluated sigma")
ax.grid(alpha=0.25)

ax = axes[0, 2]
obj_arr = ND_WEIGHT * (log_arr[:, 2] / nd_ref) + SMAX_WEIGHT * (log_arr[:, 1] / smax_ref)
ax.scatter(log_arr[:, 0], obj_arr, color=COLOR_POINTS, s=35, alpha=0.8, edgecolors="none")
ax.scatter([optimal_sigma], [optimal_obj], color=COLOR_OPTIMUM, s=110, zorder=5, label="optimum")
ax.set_xlabel("Surface tension, sigma (J/m^2)")
ax.set_ylabel("Objective")
ax.set_title("Combined objective")
ax.legend(frameon=False)
ax.grid(alpha=0.25)

ax = axes[1, 0]
ax.plot(out_water.time, out_water.S * 100, color=COLOR_WATER, label=f"pure water", linewidth=2)
ax.plot(out_opt.time, out_opt.S * 100, color=COLOR_OPTIMUM, label=f"optimal sigma", linewidth=2)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Supersaturation, S (%)")
ax.set_title("Supersaturation trajectory")
ax.legend(frameon=False)
ax.grid(alpha=0.25)


def plot_radius_profile(ax, radii, heights, activated_mask, s_max, title):
    """One panel: droplet radius (log-x) vs height, solid = activated, dashed = haze."""
    nr = radii.shape[1]
    for i in range(nr):
        color = COLOR_ACTIVATED if activated_mask[i] else COLOR_HAZE
        style = "-" if activated_mask[i] else "--"
        ax.plot(radii[:, i] * 1e6, heights, color=color, linestyle=style, linewidth=1.3, alpha=0.85)
    ax.set_xscale("log")
    ax.set_xlabel("Droplet radius, um")
    ax.set_ylabel("Height (m)")
    ax.set_title(title)
    ax.grid(alpha=0.25, which="both")
    # Legend proxies (real lines are per-bin, not individually labeled)
    ax.plot([], [], color=COLOR_ACTIVATED, linestyle="-", label="activated")
    ax.plot([], [], color=COLOR_HAZE, linestyle="--", label="haze")
    ax.legend(frameon=False, loc="upper left")


plot_radius_profile(
    axes[1, 1], radii_water, heights_water, activated_water, out_water.summary["S_max"],
    f"Droplet radius profile: pure water\n({int(activated_water.sum())}/{len(activated_water)} bins activated)",
)
plot_radius_profile(
    axes[1, 2], radii_opt, heights_opt, activated_opt, out_opt.summary["S_max"],
    f"Droplet radius profile: optimal sigma\n({int(activated_opt.sum())}/{len(activated_opt)} bins activated)",
)

plt.tight_layout()
plt.savefig("surface_tension_optimization.png", dpi=150)
print("\nSaved plot to surface_tension_optimization.png")
plt.show()
