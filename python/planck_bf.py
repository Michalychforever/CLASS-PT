from classy import Class 
import numpy as np
import matplotlib.pyplot as plt
from numpy import log, exp
import sys
from time import time 

# z = 0.61
# z = 0.38
z=0.86
cosmo = Class()
cosmo.set({'k_pivot':'0.05',
'ln10^{10}A_s':'3.044',	
# 'A_s':'2.19896453e-9',
'n_s':'0.9649',
'YHe':0.25,
'tau_reio':'0.052',
# 'T_cmb':'2.726',
'h':'0.6737',
'omega_b':'0.02237',
'N_ncdm':'0',
#'m_ncdm':'0.06',
#'N_ur':'2.0328',
'N_eff':'3.046',
'omega_cdm':'0.12',
'P_k_max_h/Mpc': '100.',
'output':'mPk,tCl',
'z_pk':z,
'non linear':'PT',
'IR resummation':'Yes',
'Bias tracers':'Yes',
'RSD':'Yes',
'AP':'Yes',
'Omfid':0.31
# 'SigmaFOG':0.
})
t1 = time()
cosmo.compute() 


h = cosmo.h()
Da =cosmo.angular_distance(z)
print("Da=",Da)
fz = cosmo.scale_independent_growth_factor_f(z)
print("fz=",fz)

# omb = cosmo.omega_b()
# Omm = cosmo.Omega_m()
# #        omm = cosmo.Omega_m()
# #        rat = omb/omm
# #        print("rat",rat)
# print("omega_b =",omb)
# print("Omega_m =",Omm)

#omcdm = Omm * h**2. - omb
# omcdm = cosmo.omegach2()
# print("omega_cdm =",omcdm)

#k1 = 0.7*1.028185622909e-5
#k1 = 1.e-6
#z = 0.0
#k1 = 0.1
#print(cosmo.pk(k1,z))
#print(cosmo.pk_lin(k1,z))

k = np.linspace(log(0.0001),log(50),200)
k = np.exp(k)
testout = [ [0 for x in range(42)] for y in range(len(k))];
for i in range(len(k)):
    testout[i][0] = k[i]
    testout[i][41] = cosmo.pk_lin(k[i]*h,z)*h**3
    for j in range(40):
#        print("j=",j)
        testout[i][j+1] = cosmo.pk(k[i]*h,z)[j]*h**3
#        testout[i][0] = k[i]
#        testout[i][1] = cosmo.pk(k[i],z)[0]
#        testout[i][2] = cosmo.pk(k[i],z)[1]
#        testout[i][3] = cosmo.pk(k[i],z)[2]
#        testout[i][4] = cosmo.pk(k[i],z)[3]
#        testout[i][5] = cosmo.pk(k[i],z)[4]
#        testout[i][6] = cosmo.pk(k[i],z)[5]
#        testout[i][7] = cosmo.pk(k[i],z)[6]
#        testout[i][8] = cosmo.pk(k[i],z)[7]
#        testout[i][9] = cosmo.pk(k[i],z)[8]
#	testout[i][10] = cosmo.pk_lin(k[i],0)
# np.savetxt('planck_pk_nl_bestfit_z038.dat', testout)
np.savetxt('planck_pk_nl_bestfit_z086.dat', testout)
t2 = time()
print("overall elapsed time=",t2-t1)
### everything is in units Mpc! 

#l = np.array(range(2,2501))
#factor = l*(l+1)/(2*np.pi)
#raw_cl = cosmo.raw_cl(2500)
#lensed_cl = cosmo.lensed_cl(2500)
#raw_cl.viewkeys()

#z = 0.2
#raw_pk = np.zeros(len(k))
#for i in range(len(k)):
#    raw_pk[i] = k[i]*cosmo.pk(k[i],z)

#cosmo.set({'P_k_max_h/Mpc': '100.','output':'mPk','z_pk':'0.2','non linear':' SPT ','IR resummation':' Yes ','Bias tracers':' No '})
#cosmo.compute()

#z = 0.2
#pk_lin = np.zeros(len(k))
#for i in range(len(k)):
#    pk_lin[i] = float(cosmo.pk(k[i],z))

#plt.loglog(k,raw_pk)
#plt.xlabel(r"$k$")

#plt.ylabel(r"$P(k)$")
#plt.tight_layout()
#plt.savefig("misha_test.pdf")

