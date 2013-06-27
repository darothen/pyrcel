import numpy as np
import time

from itertools import product
import pickle

param_ranges = {
    'alphas': (0., 1.),
    'gammas': (0.1, 0.6),
    'ds_dcs': (1.0, 1.6),
    'Ns': (100., 3500.),
    'Ms': (200., 600.),
}
params = ['alphas', 'gammas', 'ds_dcs', 'Ns', 'Ms']

n_samples = 10

## ALL RANDOM SAMPLES
samples = {}
for param in params:
    lo, hi = param_ranges[param]
    samples[param] = np.random.uniform(lo, hi, n_samples)

## All combinations
#params_iter = product(*[samples[key] for key in params]) # if want all combinations

## Just the n_samples concatenated
all_params = np.array([samples[param] for param in params]).T
params_iter = [row for row in all_params]

## WRITE OUT PARAMETERS
#
params_set = []
for p in params_iter:
    params_set.append(p)
    print p

with open("all_params", "wb") as f:
    pickle.dump(params_set, f)
