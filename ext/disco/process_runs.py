import numpy as np
import h5py

import pickle, sys

models = ["explicit", "ARG", "FN"]

## Read in data
f_data = open("run.dat", "rb")
run_data = pickle.load(f_data)
f_data.close()

'''
f_info = open("run.info", "rb")
run_info = pickle.load(f_info)
f_info.close()
'''

## REMOVE ALREADY-SEEN PARAMETERS!
'''
mus = np.linspace(0.01, 0.25, 10)
Ns = np.logspace(2, 4, 10)
Vs = np.logspace(np.log10(0.05), 1, 10)
kappas = np.linspace(0.1, 1.2, 10)
sigmas = np.linspace(1.5, 3.0, 10)
from itertools import product
param_set = product(mus, Ns, Vs, kappas, sigmas)
param_vals = [p for p in param_set]
'''


## Process all the data
parameters = []
data = {
    'Smax': {},
    'N_act': {},
    'act_frac': {},
}
for key in data:
    data_dict = {}
    for model in models: data_dict[model] = []
    data[key] = data_dict

for i, d in enumerate(run_data):
    print i
    ps, results = d
    if results == None:
        print "BAD"
        continue
    parameters.append(ps)


    for i, model in enumerate(models):
        _, Smax, act_frac, N_act, _ = results[i]
        data['Smax'][model].append(Smax)
        data['N_act'][model].append(N_act)
        data['act_frac'][model].append(act_frac)
parameters = np.array(parameters)

## Create HDF5 file
h5 = h5py.File('data.h5', 'w')

## Store parameters
dset = h5.create_dataset("parameters", data=parameters)
dset.attrs["name"] = "Parameter Set"
dset.attrs["parameters"] = ["mu", "N", "V", "kappa", "sigma"]
dset.attrs["n_samples"], dset.attrs["n_params"] = parameters.shape

## Store data
for model in models:
    group = h5.create_group(model)
    for key, data_dict in data.iteritems():
        d = np.array(data_dict[model])
        dset = group.create_dataset(key, data=d)

        dset.attrs["name"] = key

    group.attrs["name"] = model


h5.close()