import numpy as np
from numpy.fft import fft, ifft , rfft, irfft , fftfreq
from numpy import exp, log, log10, cos, sin, pi, cosh, sinh , sqrt
from scipy.special import gamma
from scipy.special import hyp2f1
from scipy import interpolate
import sys
from time import time
from scipy.integrate import quad
import scipy.integrate as integrate
from scipy import special
from time import time
from scipy.optimize import fsolve
from scipy.special import factorial


Nmax = 128
b = -0.3
kmax = 1.e2
k0 = 0.00005

Delta = log(kmax/k0) / (Nmax - 1)
jsNm = np.arange(-Nmax/2,Nmax/2+1,1)
etam = b + 2*1j*pi*(jsNm)/Nmax/Delta


def Gamma(x):
	return special.gamma(x)

def J(nu1,nu2):
	return (Gamma(1.5 - nu1) * Gamma(1.5 - nu2) * Gamma(nu1 + nu2 - 1.5) / (Gamma(nu1) * Gamma(nu2) * Gamma(3. - nu2 - nu1))) / (8. * pi**(1.5))

def M13(nu1):
	return  (1 + 9 * nu1) * np.tan(nu1*pi)/(4. * 28. * pi * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.))

def M13v2(nu1):
	return -15. * np.tan(nu1*pi)/(28.* pi * (nu1 + 1.) * nu1 * (nu1 - 1.) * (nu1 - 2.) * (nu1 - 3.))

def M22(nu1,nu2):
	return (1.5 - nu1 - nu2) * (0.5 - nu1 - nu2) * (nu1*nu2*(98.*(nu1+nu2)**2.-14.*(nu1+nu2)+36.) - 91.*(nu1+nu2)**2. + 3.*(nu1 + nu2) + 58.)*J(nu1,nu2)/(196. * nu1 * (1. + nu1) * (0.5 - nu1) * nu2 * (1. + nu2) * (0.5 - nu2))


# def M_Id2d2(nu1,nu2):
# 	return J(nu1,nu2)



m13 = np.zeros((Nmax+1),dtype=complex)
for i in range(Nmax+1):
	m13[i] = M13(-0.5 * etam[i])

m13outfile = np.zeros(2*(Nmax+1));
for i in range(Nmax+1):
        m13outfile[i] = m13.real[i]
for i in range(Nmax+1,2*(Nmax+1)):
	m13outfile[i] = m13.imag[i-Nmax-1];

np.savetxt('M13oneline_N128.dat',m13outfile)


m22mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
for j1 in range(Nmax+1):
	for j2 in range(Nmax+1):
		if j1 - Nmax/2 < 1 :
			m22mat[j1][j2] = M22(-0.5*etam[j1],-0.5*etam[j2])		
		else: 
			m22mat[j1][j2] = np.conjugate(m22mat[Nmax - j1][Nmax - j2])


mout = np.zeros((Nmax+1)*(Nmax+1),dtype=complex)
for i in range(Nmax+1):
    for j in range(Nmax+1):
        mout[i * (Nmax+1)+j] = m22mat[i][j]

mout_red = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1))//2, dtype=complex)

for i in range(Nmax+1):
    for j in range(i+1):
        mout_red[i+(2*(Nmax+1)-1-j)*j//2] = m22mat[i][j]

moutoneline = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1)))
for i in range((Nmax+1)*(Nmax+2)//2):
    moutoneline[i] = mout_red.real[i]
for i in range((Nmax+1)*(Nmax+2)//2,(Nmax+1)*(Nmax+2)):
    moutoneline[i] = mout_red.imag[i-(Nmax+1)*(Nmax+2)//2]

np.savetxt('M22oneline_N128.dat',moutoneline)



# Nmax = 256
b = -1.6000001
# kmax = 1.e2
# k0 = 0.00005

# Delta = log(kmax/k0) / (Nmax - 1)
# jsNm = np.arange(-Nmax/2,Nmax/2+1,1)
etam = b + 2*1j*pi*(jsNm)/Nmax/Delta

m13 = np.zeros((Nmax+1),dtype=complex)
b13 = np.zeros((Nmax+1),dtype=complex)
for i in range(Nmax+1):
	m13[i] = M13v2(-0.5 * etam[i])

m13outfile = np.zeros(2*(Nmax+1));
for i in range(Nmax+1):
        m13outfile[i] = m13.real[i]
for i in range(Nmax+1,2*(Nmax+1)):
	m13outfile[i] = m13.imag[i-Nmax-1];

np.savetxt('IFG2oneline_N128.dat',m13outfile)

mId2d2mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
for j1 in range(Nmax+1):
    for j2 in range(Nmax+1):
            if j1 - Nmax/2 < 1 :
                    mId2d2mat[j1][j2] = J(-0.5*etam[j1],-0.5*etam[j2])
            else:
                    mId2d2mat[j1][j2] = np.conjugate(mId2d2mat[Nmax - j1][Nmax - j2])

mId2d2out_red = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1))//2, dtype=complex)
for i in range(Nmax+1):
    for j in range(i+1):
        mId2d2out_red[i+(2*(Nmax+1)-1-j)*j//2] = mId2d2mat[i][j]

mId2d2outoneline = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1)))
for i in range((Nmax+1)*(Nmax+2)//2):
        mId2d2outoneline[i] = mId2d2out_red.real[i]
for i in range((Nmax+1)*(Nmax+2)//2,(Nmax+1)*(Nmax+2)):
        mId2d2outoneline[i] = mId2d2out_red.imag[i-(Nmax+1)*(Nmax+2)//2]

np.savetxt('M22basiconeline_N128.dat',mId2d2outoneline)
print("PT matrices succesfully recomputed!")
