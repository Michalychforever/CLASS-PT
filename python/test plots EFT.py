import numpy as np
from scipy.interpolate import interp1d

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# construct P(k) at z=0
#z_list = np.array(0,dtype=np.float,ndmin=1)
#kmax = cosmo.pars['P_k_max_h/Mpc'] * cosmo.h()
#k_list = np.logspace(-6,np.log10(kmax),100)
#pk = cosmo.get_pk_array(k_list,z_list,len(k_list),len(z_list),False)

data = np.genfromtxt("class_public-EFT-mark-2/python/planck_pk_nl_bestfit_z086.dat")
data1 = np.genfromtxt("class_public-EDE-PT/python/planck_pk_nl_bestfit_z086.dat")

#f,ax_list = plt.subplots(1,3,figsize=(9,3),constrained_layout=True)
#ax = ax_list.ravel()

#for i, axe in enumerate(ax_list):
    #axe.plot(data[:,0], data[:,i+1])
    #if i > 3:
        #break

plt.plot(data[:,0], data1[:,1]/data[:,1])
plt.show()