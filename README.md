# CLASS-PT: nonlinear perturbation theory extension of the Boltzmann code CLASS

This is a modification of the CLASS code that computes the non-linear power spectra of dark matter and biased tracers in one-loop cosmological perturbation theory, for both Gaussian and non-Gaussian initial conditions.
 
CLASS-PT can be interfaced with the MCMC sampler [MontePython](https://github.com/brinckmann/montepython_public) using the (***new and improved***) custom-built likelihoods found [here](https://github.com/oliverphilcox/full_shape_likelihoods).
 
The code is compatible with both python2 and python3.

# Getting started

Read [these instructions](https://github.com/Michalychforever/CLASS-PT/blob/master/instructions.pdf) for the installation details. See also this [troubleshooting](https://github.com/Michalychforever/CLASS-PT/blob/master/troubleshooting.rst) guide.

The installation instuctions for CLASS can be found on the official code [webpage](https://github.com/lesgourg/class_public).

Once you are all set, check out this [jupyter notebook](https://github.com/Michalychforever/CLASS-PT/blob/master/notebooks/nonlinear_pt.ipynb) for the examples of working sessions. Here's a simple example of computing the galaxy power spectrum multipoles with CLASS-PT:

```python
# Import modules
from classy import Class
import numpy as np

# Set usual CLASS parameters
z_pk = 0.5
cosmo = Class()
cosmo.set({'A_s':2.089e-9,
           'n_s':0.9649,
           'tau_reio':0.052,
           'omega_b':0.02237,
           'omega_cdm':0.12,
           'h':0.6736,
           'YHe':0.2425,
           'N_ur':2.0328,
           'N_ncdm':1,
           'm_ncdm':0.06,
           'z_pk':z_pk
          })  
# Set additional CLASS-PT settings
cosmo.set({'output':'mPk',
           'non linear':'PT',
           'IR resummation':'Yes',
           'Bias tracers':'Yes',
           'cb':'Yes', # use CDM+baryon spectra
           'RSD':'Yes',
           'AP':'Yes', # Alcock-Paczynski effect
           'Omfid':'0.31', # fiducial Omega_m
           'PNG':'No' # single-field inflation PNG
         })
cosmo.compute()

# Define some wavenumbers and compute spectra
khvec = np.logspace(-3,np.log10(1),1000) # array of k in 1/Mpc
cosmo.initialize_output(khvec, z_pk, len(khvec))

# Define nuisance parameters and extract outputs
b1, b2, bG2, bGamma3, cs0, cs2, cs4, Pshot, b4 = 2., -1., 0.1, -0.1, 0., 30., 0., 3000., 10.
pk_g0 = cosmo.pk_gg_l0(b1, b2, bG2, bGamma3, cs0, Pshot, b4)
pk_g2 = cosmo.pk_gg_l2(b1, b2, bG2, bGamma3, cs2, b4)
pk_g4 = cosmo.pk_gg_l4(b1, b2, bG2, bGamma3, cs4, b4)
```

You can also use the Mathematica notebook *'read_tables.nb'* to read the code output. We also provide a technical summary of the fNL implementations [here](https://github.com/Michalychforever/CLASS-PT/blob/master/notebooks/summary_orthogonal_github.ipynb).

# Using the code

You can use CLASS-PT freely, provided that in your publications you cite at least the code paper [arXiv:2004.10607](https://arxiv.org/abs/2004.10607). Feel free to cite the other companion papers devoted to new large-scale structure analysis methodologies! 

# Authors
- [Mikhail (Misha) Ivanov](mailto:ivanov@ias.edu)
- Anton Chudaykin 
- Marko Simonovic
- Oliver Philcox
- Giovanni Cabass
