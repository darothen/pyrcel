from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

## Cython modules
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("parcel_aux", ["parcel_aux.pyx"],
                             include_dirs=[numpy.get_include(), ],
                             #extra_compile_args=['-fopenmp', ],
                             #extra_link_args=['-fopenmp', ], 
                             )
                   ]
)
