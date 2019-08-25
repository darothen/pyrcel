from itertools import product

import numpy as np

## Import unit testing librarys
try:
    import unittest2 as unittest
except ImportError:
    import unittest


def prod_to_array(*iterables):
    return np.array(list(product(*iterables)))
