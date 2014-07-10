
from numpy.polynomial import polynomial, hermite_e, legendre
import numpy as np
from numpy import array

from parcel_model.micro import kohler_crit
from scipy.special import erfc, erf, erfinv

brackets = {'mu': ('uniform_n', 0.01, 0.25), 
            'sigma': ('uniform_n', 1.2, 3.0), 
            'V': ('weibull_n', 2.0, 1.0), 
            'kappa': ('uniform_n', 0.1, 1.2), 
            'N': ('uniform_n', 100.0, 10000.0)}
def map_to_norm(x, param):
    prior, args = brackets[param][0], brackets[param][1:]

    lo, hi = args

    if prior == "lognormal_n": args = map(np.log, args)
    if prior == "uniform_n": 
        if x <= lo: 
            x = lo*(1. + 1e-4)
        elif x >= hi: 
            x = hi*(1. - 1e-4)
        pass

    inv_func = inv_funcs[prior]
    return inv_func(x, *args)

def inv_uni(x, a0, b0, af=-1, bf=1):
    """
    Map a value from the uniform distribution [a, b]
    to the uniform distribution [af, bf]
    """
    return ((bf-af)/(b0-a0))*(x - a0) + af
def uni_to_norm(x, a, b):
    """
    Map a value from the uniform distribution [a, b] to the normal distribution
    """
    return np.sqrt(2)*erfinv(2*(x-a)/(b-a)  - 1.0)
def lognorm_to_norm(x, mu, sigma):
    """
    Map a value from the lognormal distribution with given mu and sigma to the
    standard normal distribution with mean 0 and std 1
    """
    return (np.log(x)-mu)/sigma
def gamma_to_norm(x, a, b):
    """
    Map a value from the gamma distribution with mean and shape parameters a, b to
    the standard normal distribution with mean 0 and std 1
    """
    return (((x/a/b)**(1./3.)) + 2./9./a - 1.0)/np.sqrt(2./9./a)
def weibull_to_norm(x, k, lam=1):
    exp_var = (x**(1./k))/lam
    return np.sqrt(2.)*erfinv(2.*np.exp(-exp_var*lam) - 1.0)

inv_funcs = {
    "uniform": inv_uni,
    "uniform_n": uni_to_norm,
    "gamma_n": gamma_to_norm,
    "weibull_n": weibull_to_norm,
    "lognormal_n": lognorm_to_norm,
}

def eval_hermite(deg, x, analytic=True):
    if analytic: # Use pre-computed polynomials
        if deg == 0:
            return 1.
        elif deg == 1.:
            return x
        elif deg == 2:
            return x*x - 1.
        elif deg == 3:
            return x*x*x - 3.*x
        elif deg == 4:
            return x*x*x*x - 6.*x*x + 3.

    # If not returned, assume we should make the call to the numerical orthogonalizer
    roots = hermite_e.hermegauss(deg)[0] if deg > 1 else np.array([0, ])
    poly =  np.poly1d(roots, r=True)
    return poly(x)

def eval_legendre(deg, x):
    roots = legendre.leggauss(deg)[0] if deg > 1 else np.array([0, ])
    poly =  np.poly1d(roots, r=True)
    return poly(x)

#@profile
def pce_deg1(V, T, P, aerosol):

    As = array([ 0.00022826, -0.00051625, -0.00020289, -0.00039829, -0.0001025 ,
       -0.0001638 ])

    mu = map_to_norm(aerosol.distribution.mu, 'mu')
    N = map_to_norm(aerosol.distribution.N, 'N')
    V = map_to_norm(V, 'V')
    sigma = map_to_norm(aerosol.distribution.sigma, 'sigma')
    kappa = map_to_norm(aerosol.kappa, 'kappa')

    V_1 = eval_hermite(1, V)
    N_1 = eval_hermite(1, N)
    kappa_1 = eval_hermite(1, kappa)
    mu_1 = eval_hermite(1, mu)
    sigma_1 = eval_hermite(1, sigma)

    Smax = As.dot(np.array([1.0, mu_1, N_1, V_1, kappa_1, sigma_1, ]))

    _, S_mode_crit = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
    ui = 2.*np.log(S_mode_crit/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
    N_act = 0.5*aerosol.distribution.N*erfc(ui)
    act_frac = N_act/aerosol.distribution.N

    return Smax, act_frac

#@profile
def pce_deg2(V, T, P, aerosol):

    As = array([  1.14326819e-03,  -7.17119470e-04,  -2.87207700e-04,
        -6.49925186e-04,  -1.32394919e-04,  -2.03550833e-04,
         4.19256128e-04,   1.95758520e-04,   4.28430253e-04,
         5.43726035e-05,   4.01075063e-05,   1.78061247e-04,
         6.89565065e-04,   1.06794441e-04,   2.34925303e-04,
         5.96974401e-04,   9.09500909e-05,   7.94986009e-05,
         1.49255164e-04,   1.22341600e-04,   5.20229943e-05])

    mu = map_to_norm(aerosol.distribution.mu, 'mu')
    N = map_to_norm(aerosol.distribution.N, 'N')
    V = map_to_norm(V, 'V')
    sigma = map_to_norm(aerosol.distribution.sigma, 'sigma')
    kappa = map_to_norm(aerosol.kappa, 'kappa')

    V_1 = eval_hermite(1, V)
    V_2 = eval_hermite(2, V)
    N_1 = eval_hermite(1, N)
    kappa_2 = eval_hermite(2, kappa)
    N_2 = eval_hermite(2, N)
    sigma_2 = eval_hermite(2, sigma)
    mu_2 = eval_hermite(2, mu)
    kappa_1 = eval_hermite(1, kappa)
    mu_1 = eval_hermite(1, mu)
    sigma_1 = eval_hermite(1, sigma)

    Smax = As.dot(np.array([1.0, mu_1, N_1, V_1, kappa_1, sigma_1, mu_2, N_2, V_2, kappa_2, sigma_2, mu_1*N_1, mu_1*V_1, mu_1*kappa_1, mu_1*sigma_1, N_1*V_1, N_1*kappa_1, N_1*sigma_1, V_1*kappa_1, V_1*sigma_1, kappa_1*sigma_1]))

    _, S_mode_crit = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
    ui = 2.*np.log(S_mode_crit/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
    N_act = 0.5*aerosol.distribution.N*erfc(ui)
    act_frac = N_act/aerosol.distribution.N

    return Smax, act_frac

def pce_deg3(V, T, P, aerosol):


    As = array([  9.14798604e-04,  -9.22633928e-04,  -8.06844177e-04,
        -1.43940090e-03,  -1.77269894e-04,  -2.30463363e-04,
         2.65464696e-04,   1.51795190e-04,   3.10878234e-04,
         2.59723751e-05,   3.18675085e-05,  -8.25970041e-05,
        -6.73965549e-05,  -2.32713884e-04,  -3.38946258e-07,
         9.64888514e-06,   1.70653952e-04,   4.44587400e-04,
         7.76929961e-05,   1.90554293e-04,   4.03803792e-04,
         6.04459159e-05,   9.44462548e-05,   1.07125979e-04,
         1.13169278e-04,   5.70518250e-05,  -9.24557555e-05,
        -3.11843201e-04,  -2.00457356e-05,  -2.11958800e-05,
        -5.27092388e-05,  -4.79579441e-04,  -1.51607324e-05,
        -9.92390215e-06,  -2.32074053e-04,  -2.72153026e-04,
        -1.84122185e-05,   1.02236473e-05,  -2.28193188e-05,
        -2.92921954e-05,  -5.24197442e-05,  -8.52954010e-07,
        -6.09448717e-05,  -4.54888486e-05,   7.46022148e-07,
        -3.49603741e-06,  -2.27098230e-04,  -1.80186021e-05,
         1.12536358e-05,  -3.05864945e-05,   3.41873606e-05,
         2.45392802e-05,  -8.72265689e-05,  -4.30851978e-05,
         1.24256190e-06,   1.07910264e-04])


    mu = map_to_norm(aerosol.distribution.mu, 'mu')
    N = map_to_norm(aerosol.distribution.N, 'N')
    V = map_to_norm(V, 'V')
    sigma = map_to_norm(aerosol.distribution.sigma, 'sigma')
    kappa = map_to_norm(aerosol.kappa, 'kappa')

    V_1 = eval_hermite(1, V)
    kappa_3 = eval_hermite(3, kappa)
    V_2 = eval_hermite(2, V)
    N_1 = eval_hermite(1, N)
    kappa_2 = eval_hermite(2, kappa)
    N_2 = eval_hermite(2, N)
    sigma_2 = eval_hermite(2, sigma)
    mu_2 = eval_hermite(2, mu)
    N_3 = eval_hermite(3, N)
    kappa_1 = eval_hermite(1, kappa)
    mu_3 = eval_hermite(3, mu)
    mu_1 = eval_hermite(1, mu)
    sigma_3 = eval_hermite(3, sigma)
    V_3 = eval_hermite(3, V)
    sigma_1 = eval_hermite(1, sigma)

    Smax = As.dot(np.array([1.0, mu_1, N_1, V_1, kappa_1, sigma_1, mu_2, N_2, V_2, kappa_2, sigma_2, mu_3, N_3, V_3, kappa_3, sigma_3, mu_1*N_1, mu_1*V_1, mu_1*kappa_1, mu_1*sigma_1, N_1*V_1, N_1*kappa_1, N_1*sigma_1, V_1*kappa_1, V_1*sigma_1, kappa_1*sigma_1, mu_1*N_2, mu_1*V_2, mu_1*kappa_2, mu_1*sigma_2, N_1*mu_2, N_1*V_2, N_1*kappa_2, N_1*sigma_2, V_1*mu_2, V_1*N_2, V_1*kappa_2, V_1*sigma_2, kappa_1*mu_2, kappa_1*N_2, kappa_1*V_2, kappa_1*sigma_2, sigma_1*mu_2, sigma_1*N_2, sigma_1*V_2, sigma_1*kappa_2, mu_1*N_1*V_1, mu_1*N_1*kappa_1, mu_1*N_1*sigma_1, mu_1*V_1*kappa_1, mu_1*V_1*sigma_1, mu_1*kappa_1*sigma_1, N_1*V_1*kappa_1, N_1*V_1*sigma_1, N_1*kappa_1*sigma_1, V_1*kappa_1*sigma_1]))

    _, S_mode_crit = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
    ui = 2.*np.log(S_mode_crit/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
    N_act = 0.5*aerosol.distribution.N*erfc(ui)
    act_frac = N_act/aerosol.distribution.N

    return Smax, act_frac

def pce_deg4(V, T, P, aerosol):
    As = array([  7.56052516e-04,  -9.34563322e-04,  -5.60707772e-04,
        -2.81045850e-04,  -8.60005463e-05,  -1.82949921e-04,
         3.64882107e-04,   3.49237541e-04,  -3.86741535e-04,
         9.02462861e-05,   8.08886457e-05,  -1.33596686e-05,
        -3.66592784e-05,   1.84829967e-04,   6.31816281e-06,
         1.43484844e-05,  -2.21552918e-05,   7.80306040e-06,
        -1.63991275e-04,  -3.31969381e-06,  -2.77832570e-06,
         1.31325527e-04,   3.36410592e-04,   1.29343179e-04,
         2.53607557e-04,   4.72821983e-04,   1.10686337e-04,
        -7.76721989e-06,   5.59677016e-05,  -2.03761928e-04,
         1.49402457e-04,  -5.81947331e-05,  -2.87190526e-04,
        -3.32166236e-05,  -3.11984126e-05,  -3.12737856e-05,
        -2.90655230e-04,  -3.87123623e-06,   6.32799902e-06,
        -1.99773540e-04,  -2.14317445e-04,   1.88693178e-06,
         4.01362982e-05,  -2.82253395e-05,  -9.73456375e-06,
         5.36536549e-05,   2.32915895e-06,  -5.70239230e-05,
        -1.01407528e-05,   3.66706996e-05,   1.34141867e-06,
        -2.51711061e-04,  -1.77305074e-05,   1.77902144e-05,
        -1.92713812e-04,  -1.36248503e-04,  -6.39443575e-05,
        -7.64443869e-05,  -7.75791871e-05,   8.91567781e-06,
         7.53736769e-07,   6.30223806e-06,  -1.34681343e-04,
        -4.42559165e-06,  -1.00314540e-05,  -1.41222799e-05,
        -1.02167268e-06,  -1.73354939e-06,  -7.14665553e-06,
        -9.84241564e-06,   6.97345884e-05,   6.02591893e-06,
        -1.40684526e-05,  -6.74866315e-06,   3.57188698e-06,
        -6.56244944e-05,   1.71048935e-06,  -9.80075449e-06,
         4.05502943e-06,  -1.47650326e-04,   4.27525023e-06,
         1.48744812e-05,   8.65829131e-05,   1.58888194e-05,
         1.50391139e-05,   2.04711595e-04,   8.58208991e-06,
         2.56625369e-06,   3.83611977e-05,   4.11932517e-05,
         8.16629817e-06,   7.99776514e-06,  -3.15884420e-06,
        -1.86066573e-05,   3.63951273e-05,   1.75620847e-05,
        -6.68152115e-06,  -1.32085279e-05,   2.22445829e-05,
         1.13519416e-05,  -2.76824182e-05,   6.90459595e-05,
         1.70119847e-05,  -3.14717813e-05,   2.62866496e-05,
        -3.92082853e-05,  -2.89646068e-05,   8.58123514e-05,
        -4.97675282e-06,  -4.48064665e-05,  -1.34016155e-05,
         4.09449969e-06,   3.68871564e-05,   7.32857832e-05,
         3.70926853e-05,  -3.92935660e-06,   2.27675971e-05,
         4.16571690e-05,   5.83095830e-06,   3.52573949e-06,
         1.00259494e-04,  -4.64968600e-05,  -6.63668001e-05,
        -3.33467254e-05,  -4.04622167e-05,   1.43066345e-05])


    mu = map_to_norm(aerosol.distribution.mu, 'mu')
    N = map_to_norm(aerosol.distribution.N, 'N')
    V = map_to_norm(V, 'V')
    sigma = map_to_norm(aerosol.distribution.sigma, 'sigma')
    kappa = map_to_norm(aerosol.kappa, 'kappa')

    V_1 = eval_hermite(1, V)
    mu_4 = eval_hermite(4, mu)
    kappa_3 = eval_hermite(3, kappa)
    V_2 = eval_hermite(2, V)
    N_1 = eval_hermite(1, N)
    kappa_2 = eval_hermite(2, kappa)
    N_2 = eval_hermite(2, N)
    sigma_2 = eval_hermite(2, sigma)
    V_4 = eval_hermite(4, V)
    mu_2 = eval_hermite(2, mu)
    N_3 = eval_hermite(3, N)
    kappa_1 = eval_hermite(1, kappa)
    mu_3 = eval_hermite(3, mu)
    N_4 = eval_hermite(4, N)
    mu_1 = eval_hermite(1, mu)
    sigma_3 = eval_hermite(3, sigma)
    V_3 = eval_hermite(3, V)
    sigma_1 = eval_hermite(1, sigma)
    sigma_4 = eval_hermite(4, sigma)
    kappa_4 = eval_hermite(4, kappa)

    Smax = As.dot(np.array([1.0, mu_1, N_1, V_1, kappa_1, sigma_1, mu_2, N_2, V_2, kappa_2, sigma_2, mu_3, N_3, V_3, kappa_3, sigma_3, mu_4, N_4, V_4, kappa_4, sigma_4, mu_1*N_1, mu_1*V_1, mu_1*kappa_1, mu_1*sigma_1, N_1*V_1, N_1*kappa_1, N_1*sigma_1, V_1*kappa_1, V_1*sigma_1, kappa_1*sigma_1, mu_1*N_2, mu_1*V_2, mu_1*kappa_2, mu_1*sigma_2, N_1*mu_2, N_1*V_2, N_1*kappa_2, N_1*sigma_2, V_1*mu_2, V_1*N_2, V_1*kappa_2, V_1*sigma_2, kappa_1*mu_2, kappa_1*N_2, kappa_1*V_2, kappa_1*sigma_2, sigma_1*mu_2, sigma_1*N_2, sigma_1*V_2, sigma_1*kappa_2, mu_1*N_1*V_1, mu_1*N_1*kappa_1, mu_1*N_1*sigma_1, mu_1*V_1*kappa_1, mu_1*V_1*sigma_1, mu_1*kappa_1*sigma_1, N_1*V_1*kappa_1, N_1*V_1*sigma_1, N_1*kappa_1*sigma_1, V_1*kappa_1*sigma_1, mu_1*N_3, mu_1*V_3, mu_1*kappa_3, mu_1*sigma_3, N_1*mu_3, N_1*V_3, N_1*kappa_3, N_1*sigma_3, V_1*mu_3, V_1*N_3, V_1*kappa_3, V_1*sigma_3, kappa_1*mu_3, kappa_1*N_3, kappa_1*V_3, kappa_1*sigma_3, sigma_1*mu_3, sigma_1*N_3, sigma_1*V_3, sigma_1*kappa_3, mu_2*N_2, mu_2*V_2, mu_2*kappa_2, mu_2*sigma_2, N_2*V_2, N_2*kappa_2, N_2*sigma_2, V_2*kappa_2, V_2*sigma_2, kappa_2*sigma_2, mu_1*N_1*V_2, mu_1*N_1*kappa_2, mu_1*N_1*sigma_2, mu_1*V_1*N_2, mu_1*V_1*kappa_2, mu_1*V_1*sigma_2, mu_1*kappa_1*N_2, mu_1*kappa_1*V_2, mu_1*kappa_1*sigma_2, mu_1*sigma_1*N_2, mu_1*sigma_1*V_2, mu_1*sigma_1*kappa_2, N_1*V_1*mu_2, N_1*V_1*kappa_2, N_1*V_1*sigma_2, N_1*kappa_1*mu_2, N_1*kappa_1*V_2, N_1*kappa_1*sigma_2, N_1*sigma_1*mu_2, N_1*sigma_1*V_2, N_1*sigma_1*kappa_2, V_1*kappa_1*mu_2, V_1*kappa_1*N_2, V_1*kappa_1*sigma_2, V_1*sigma_1*mu_2, V_1*sigma_1*N_2, V_1*sigma_1*kappa_2, kappa_1*sigma_1*mu_2, kappa_1*sigma_1*N_2, kappa_1*sigma_1*V_2, mu_1*N_1*V_1*kappa_1, mu_1*N_1*V_1*sigma_1, mu_1*N_1*kappa_1*sigma_1, mu_1*V_1*kappa_1*sigma_1, N_1*V_1*kappa_1*sigma_1]))

    _, S_mode_crit = kohler_crit(T, aerosol.distribution.mu*1e-6, aerosol.kappa)
    ui = 2.*np.log(S_mode_crit/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.distribution.sigma))
    N_act = 0.5*aerosol.distribution.N*erfc(ui)
    act_frac = N_act/aerosol.distribution.N

    return Smax, act_frac

def MARS_param(V, T, P, aerosol):

    N = aerosol.distribution.N
    mu = aerosol.distribution.mu
    sigma = aerosol.distribution.sigma
    kappa = aerosol.kappa

    h = lambda x: np.max([0, x])

    retval = 0.000800009 + \
            h(mu-0.0579383)*-0.00594626 + \
            h(0.0579383-mu)*0.054127      + \
            h(N-893.673)*-1.34729e-07  + \
            h(893.673-N)*1.26029e-06   + \
            h(sigma-1.57034)*h(0.0579383-mu)*-0.0134238    + \
            h(1.57034-sigma)*h(0.0579383-mu)*0.0580031     + \
            h(V-0.687964)*h(0.0579383-mu)*0.0106313     + \
            h(0.687964-V)*h(0.0579383-mu)*-0.0427569    + \
            h(V-1.16618)*h(893.673-N)*-6.66641e-07  + \
            h(1.16618-V)*h(893.673-N)*-7.5702e-07   + \
            h(kappa-0.396568)*-0.000168919  + \
            h(0.396568-kappa)*0.00188196    + \
            h(sigma-2.24391)*-0.000482217  + \
            h(2.24391-sigma)*0.000439122   + \
            h(N-919.342)*h(mu-0.0579383)*-2.3488e-06   + \
            h(919.342-N)*h(mu-0.0579383)*-4.10455e-05  + \
            h(mu-0.0171867)*h(0.396568-kappa)*-0.00726007   + \
            h(0.0171867-mu)*h(0.396568-kappa)*0.470438      + \
            h(V-0.373096)*0.000195676   + \
            h(0.373096-V)*-0.000635335  + \
            h(mu-0.140341)*h(2.24391-sigma)*-0.00113277+ \
            h(0.140341-mu)*h(2.24391-sigma)*0.0086391+ \
            h(N-3308.2)*h(V-0.373096)*-1.68353e-08  + \
            h(3308.2-N)*h(V-0.373096)*6.78315e-07        + \
            h(1.57551-sigma)*h(0.0171867-mu)*h(0.396568-kappa)*14.5272       + \
            h(mu-0.125655)*h(V-0.373096)*-0.000852333  + \
            h(0.125655-mu)*h(V-0.373096)*0.00425788    + \
            h(N-274.315)*h(3308.2-N)*h(V-0.373096)*-9.03331e-11  + \
            h(274.315-N)*h(3308.2-N)*h(V-0.373096)*4.22941e-09   + \
            h(mu-0.0307511)*h(N-893.673)*2.62773e-06   + \
            h(0.0307511-mu)*h(N-893.673)*7.34263e-06   + \
            h(kappa-0.213236)*h(3308.2-N)*h(V-0.373096)*-2.30155e-07  + \
            h(0.213236-kappa)*h(3308.2-N)*h(V-0.373096)*1.62512e-06   + \
            h(kappa-0.23752)*h(0.140341-mu)*h(2.24391-sigma)*-0.00640169   + \
            h(0.23752-kappa)*h(0.140341-mu)*h(2.24391-sigma)*0.0447306     + \
            h(mu-0.0166662)*h(893.673-N)*3.50406e-05   + \
            h(0.0166662-mu)*h(893.673-N)*0.000776855   + \
            h(sigma-2.0951)*h(mu-0.0307511)*h(N-893.673)*4.12976e-07   + \
            h(2.0951-sigma)*h(mu-0.0307511)*h(N-893.673)*-1.72103e-07  + \
            h(mu-0.0552353)*h(3308.2-N)*h(V-0.373096)*-2.0652e-06   + \
            h(0.0552353-mu)*h(3308.2-N)*h(V-0.373096)*-5.05412e-06  + \
            h(N-4275.75)*h(sigma-1.57034)*h(0.0579383-mu)*-9.59983e-07  + \
            h(4275.75-N)*h(sigma-1.57034)*h(0.0579383-mu)*4.8373e-06    + \
            h(mu-0.137268)*h(mu-0.0579383)*0.0212281     + \
            h(0.137268-mu)*h(mu-0.0579383)*-0.0802494    + \
            h(kappa-0.696631)*h(0.125655-mu)*h(V-0.373096)* -0.00223781   + \
            h(0.696631-kappa)*h(0.125655-mu)*h(V-0.373096)*0.00765882

    return retval, None