from classy import Class
import numpy as np
#import matplotlib.pyplot as plt
import sys
from time import time

t1 = time()

cosmo = Class ()
#cosmo.set({'k_pivot':'0.05','A_s':'2.3e-9','n_s':'1.','alpha_s':'0.','Omega_k':'0.','Omega_fld':'0','YHe':0.25,'z_reio':'10.','T_cmb':'2.726','h':'0.7','Omega_b':'0.05','N_ncdm':'0','N_eff':'3.04','Omega_cdm':'0.25','P_k_max_h/Mpc': '100.','output':'mPk','z_pk':'0.0','non linear':' SPT ','IR resummation':' Yes ','Bias tracers':' Yes '})
#cosmo.set({'P_k_max_h/Mpc': '100.','output':'mPk','z_pk':'0.0','non linear':' halofit '})

cosmo.set({'k_pivot':'0.05','A_s':'2.3e-9','n_s':'1.','alpha_s':'0.','Omega_k':'0.','Omega_fld':'0','YHe':0.25,'z_reio':'10.','T_cmb':'2.726','h':'0.7','Omega_b':'0.05','N_ncdm':'0','N_eff':'3.04','Omega_cdm':'0.25','P_k_max_h/Mpc': '100.','output':'tCl,mPk','z_pk':'0.5','non linear':'SPT','IR resummation':' Yes ','Bias tracers':' Yes '})

cosmo.compute()
t11 = time()

#k=0.01
#print k/h,cosmo.pk(k*h,z)[1]*h**3

h=0.7
N=100
z=0.5
kmin=0.001*h
kmax=1*h

k = np.linspace(np.log(kmin),np.log(kmax),N)
k = np.exp(k)

for i in range(len(k)):
    print k[i]/h,cosmo.pk(k[i],z)[0]*h**3, cosmo.pk(k[i],z)[7]*h**1,cosmo.pk(k[i],z)[1]*h**3,cosmo.pk(k[i],z)[2]*h**3,cosmo.pk(k[i],z)[3]*h**3,cosmo.pk(k[i],z)[4]*h**3,cosmo.pk(k[i],z)[5]*h**3,cosmo.pk(k[i],z)[6]*h**3



#print(cosmo.pk(k1,z))
#print(cosmo.pk_lin(k1,z))

#print(cosmo.pk(k1,z)[0])

t2 = time()
print("elapsed all time=",t2-t1)
print("only nonlinear=",t11-t1)

