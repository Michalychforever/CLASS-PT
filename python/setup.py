from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as nm
import os
import subprocess as sbp
import os.path as osp

# Recover the gcc compiler
GCCPATH_STRING = sbp.Popen(
    ['gcc', '-print-libgcc-file-name'],
    stdout=sbp.PIPE).communicate()[0]
GCCPATH = osp.normpath(osp.dirname(GCCPATH_STRING)).decode()

#liblist = ["class","gsl","gslcblas","openblas"]
liblist = ["class"]
MVEC_STRING = sbp.Popen(
    ['gcc', '-lmvec'],
    stderr=sbp.PIPE).communicate()[1]
if b"mvec" not in MVEC_STRING:
    liblist += ["mvec","m"]

# Recover the CLASS version
with open(os.path.join('..', 'include', 'common.h'), 'r') as v_file:
    for line in v_file:
        if line.find("_VERSION_") != -1:
            # get rid of the " and the v
            VERSION = line.split()[-1][2:-1]
            break

setup(
    name='classy',
    version=VERSION,
    description='Python interface to the Cosmological Boltzmann code CLASS',
    url='http://www.class-code.net',
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("classy", ["classy.pyx"],
                           include_dirs=[nm.get_include(), "../include","/Users/gcabass/anaconda3/envs/openblas_test/include"],
                           libraries=liblist,
                           library_dirs=["../", GCCPATH],
                           extra_link_args=['/Users/gcabass/anaconda3/envs/openblas_test/lib/libopenblas.dylib','-lgomp'],
                           )],
    #data_files=[('bbn', ['../bbn/sBBN.dat'])]
)


