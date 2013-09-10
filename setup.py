from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

#import cython_gsl as cgsl


## Cython modules
setup(
    #include_dirs = [cgsl.get_include(), ], 
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("parcel_aux", ["parcel_aux.pyx"],
                             include_dirs=[numpy.get_include(), ],
                             #extra_compile_args=['-fopenmp', ],
                             #extra_link_args=['-fopenmp', ], 
                             ),
                   #Extension("parcel_aux_bin", ["parcel_aux_bin.pyx"],
                   #          libraries=cgsl.get_libraries(),
                   #          library_dirs=[cgsl.get_library_dir(), ],
                   #          include_dirs=[numpy.get_include(), cgsl.get_cython_include_dir()],
                   #          ),
                   ]
)
