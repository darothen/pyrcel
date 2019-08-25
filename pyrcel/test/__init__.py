from __future__ import absolute_import

import numpy as np
from numpy.testing import assert_array_equal, assert_allclose

from itertools import product

## Import unit testing librarys
try:
    import unittest2 as unittest
except ImportError:
    import unittest


def prod_to_array(*iterables):
    return np.array(list(product(*iterables)))
