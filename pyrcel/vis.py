import numpy as np


def plot_distribution(aer, aer_kwargs={}, ax=None, **kwargs):
    """Generate a comparison plot of a given aerosol or
    droplet distribution

    Parameters
    ----------
    aer : `AerosolSpecies`
        The container class for the aerosol
    aer_kwargs : dict
        A dictionary of arguments to pass to the matplotlib
        function which plots the binned distribution
    ax : Axis
        The axes object to plot on

    """

    if ax is None:
        raise ValueError("Must provide axes instance for plotting.")

    ## Add some basic aer_kwargs if not provided
    if not "color" in aer_kwargs:
        aer_kwargs["color"] = "g"
    if not "alpha" in aer_kwargs:
        aer_kwargs["alpha"] = 0.5

    rl, rr = aer.rs[0], aer.rs[-1]
    r_left = aer.rs[:-1]
    r_width = aer.rs[1:] - r_left
    r_mean = np.sqrt(aer.rs[1:] * r_left)
    bin_height = aer.Nis / 1e6
    bars = ax.bar(r_left, bin_height, width=r_width, **aer_kwargs)

    legend_objects = [(bars[0], "%s bins" % aer.species)]

    handles, labels = list(zip(*legend_objects))
    ax.legend(handles, labels, loc="upper right")

    ax.semilogx()
    ax.set_xlim([rl, rr])

    ax.set_xlabel("$r$ ($\mu$m)")
    ax.set_ylabel("Number Concentration (cm$^{-3}$)")
