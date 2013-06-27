import numpy as np
import h5py

from sklearn import manifold, cross_validation

from pylab import *
ion()
from mpl_toolkits.mplot3d import Axes3D

## Read in data from storage
with h5py.File("data.h5", 'r') as f:

    params = f['/parameters'][:]
    mus, Ns, Vs, kappas, sigmas = [params[:,i] for i in xrange(params.shape[1])]

    Smax = f['/explicit/Smax'][:]

n_samples, n_features = params.shape

ss = np.logspace(-4, -1, 1000)

train_inds = sorted(np.random.random_integers(0, n_samples-1, 100))
vals_train = params[train_inds]
data_train = Smax[train_inds]

test_inds = sorted(np.random.random_integers(0, n_samples-1, n_samples-100))
vals_test = params[test_inds]
data_test = Smax[test_inds]

## Gaussian Process Regression
from sklearn.gaussian_process import GaussianProcess
X = vals_train[:, :]
y = data_train

# Instantiate a Gaussian Process model
gp = GaussianProcess(regr='linear', verbose=True)
gp.fit(X, y)

N = 1200.
mu = 0.02
kappa = 0.5
sigma = 2.5

ref = np.array([mu, N, kappa, sigma])

f = lambda a, b, tol=0.25: np.abs(a - b)/b < tol
def good_val(val, ref=ref):
    val_pack = np.array([val[i] for i in [0, 1, 3, 4]])
    return np.all(f(val_pack, ref))

## Select a slice of N/mu as a function of V to test
x_pred = np.array([v for v in vals_test if good_val(v)])
y_pred = np.array([d for v, d in zip(vals_test, data_test) if good_val(v)])

v_lin = np.logspace(-1, np.log10(3.), 1000)
v_x_pred = x_pred[:, 2]


y_pred, MSE = gp.predict(x_pred, eval_MSE=True)
sigma = np.sqrt(MSE)

figure(4); clf()
fig, axes = subplots(1, 1, num=4, figsize=(8, 6))
ax4 = axes

ax4.scatter(v_x, marker='o', color='r', s=15, label="Explicit")
ax4.plot(, y_pred, 'k', label="Predicted")

ax4.set_xlim(0.1, 3.0)
ax4.set_ylim(1e-3, 1e-2)

## Support vector regression
'''
from sklearn import svm
X = vals_train[:, :]
## Take log of mus, 100*Vs
#X[:,0] = np.log(X[:,0])
#X[:,2] = 100.*X[:,0]

y = data_train*1e3

fitter = svm.SVR(C=1e5, verbose=True)
fitter.fit(X, y)

figure(3); clf()
fig, axes = subplots(1, 1, num=3, figsize=(8,6))
ax1 = axes

predicted = fitter.predict(vals_test)

scatter(data_test, predicted*1e-3, marker='x', s=15)
plot(ss, ss, color='k', linestyle='dashed')
loglog()
xlim(1e-4, 1e-1)
ylim(1e-4, 1e-1)
'''

'''
#X_r, err = manifold.locally_linear_embedding(vals_1d, n_neighbors=100, n_components=2)
X_r = manifold.Isomap(n_neighbors=20, n_components=2).fit_transform(vals_train)
#X_r = manifold.SpectralEmbedding(n_components=2, n_neighbors=100).fit_transform(vals_train)

figure(1); clf()
fig, axes = subplots(2, 1, num=1, figsize=(6,12))
ax1, ax2 = axes

ax1 = fig.add_subplot(211, projection='3d')
ax1.scatter(vals_train[:, 0], vals_train[:, 1], vals_train[:, 2], c=data_train, cmap=cm.Spectral)
ax1.set_title("Original_data")

ax2.scatter(X_r[:, 0], X_r[:, 1], c=data_train, cmap=cm.Spectral)
ax2.set_title("Projected data")
'''

'''
from sklearn import linear_model
clf = linear_model.Ridge(alpha=0.5)
clf.fit(vals_train, data_train)


f = lambda vals: np.sum(clf.coef_*vals) + clf.intercept_

xs = data_test
ys = [f(vals) for vals in vals_test]
ss = np.logspace(-4, -2, 1000)

figure(2)
scatter(xs, ys, marker='x', s=15)
plot(ss, ss, linestyle='dashed')
loglog()
xlim(1e-4, 1e-2)
ylim(1e-4, 1e-2)
'''

