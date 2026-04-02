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


#Nmax = 128
#Nmax = 256
Nmax = 512
#b = -0.3
b=-0.8 #--> WHY DO WE NEED THIS HERE??? No, I need it, because the matrix is computed at \nu1, \nu2 that depend on b. And then, outside, I multiply by k to the two powers, that also depends on b. This raises an issue -> which issue did I have in mind? Right now I cannot see any problem with this... I will need to change the bias for the biased tracers but apart from that it should be ok... Just make sure that the same kmax and k0 are used in C? No, trust Misha on this...
#--> I think that the issue was just something I thought but then resolved in my mind without writing anything...
#--> the only confusion here was that I don't understand the two names for the matrix of 256 points... I will make another file for the generation of the matrices for biased tracers...
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

def M12(nu1,nu2):
    return J(nu1,nu2)*((1./sin(nu1*pi)/sin(nu2*pi)*(-2.*(60. + 448.*nu1*nu1*nu1*nu1*(-1. + nu2)*nu2 - 3.*nu2*(43. + 4.*nu2*(-15. + 7.*nu2)) + 4.*nu1*nu1*nu1*(-21. + 2.*nu2*(169. + 4.*nu2*(-65. + 28.*nu2))) + 4.*nu1*nu1*(45. + 4.*nu2*(-79. + nu2*(169. + 2.*nu2*(-65. + 14.*nu2)))) + nu1*(-129. - 8.*nu2*(-60. + nu2*(158. + nu2*(-169. + 56.*nu2)))))*cos((nu1 + nu2)*pi) - 1./sin((nu1 + nu2)*pi)*((-60. + 3.*nu1*(43. + 4.*nu1*(-15. + 7.*nu1)) + 209.*nu2 - 8.*nu1*(17. + (-2. + nu1)*nu1*(-43. + 56.*nu1))*nu2 + 4.*(101. + 2.*nu1*(-229. + 2.*nu1*(245. + 4.*nu1*(-43. + 7.*nu1))))*nu2*nu2 + 4.*(-345. + 2.*nu1*(463. + 4.*nu1*(-137. + 42.*nu1)))*nu2*nu2*nu2 + 64.*(16. + 3.*nu1*(-12. + 7.*nu1))*nu2*nu2*nu2*nu2 + 224.*(-1. + 2.*nu1)*nu2*nu2*nu2*nu2*nu2)*sin(2.*nu1*pi) + (-60. + nu1*(209. + 4.*nu1*(101. + nu1*(-345. + 8.*(32. - 7.*nu1)*nu1))) + 129.*nu2 + 8.*nu1*(-17. + nu1*(-229. + nu1*(463. + 8.*nu1*(-36. + 7.*nu1))))*nu2 + 4.*(-45. + 4.*nu1*(-43. + nu1*(245. - 274.*nu1 + 84.*nu1*nu1)))*nu2*nu2 + 4.*(21. + 2.*nu1*(155. + 8.*nu1*(-43. + 21.*nu1)))*nu2*nu2*nu2 + 448.*(-1. + nu1)*nu1*nu2*nu2*nu2*nu2)*sin(2.*nu2*pi))))/(28.*nu1*(-1. + 4.*nu1*nu1)*nu2*(-5. + 2.*nu1 + 2.*nu2)*(-1. + 4.*nu2*nu2)))


# def M_Id2d2(nu1,nu2):
# 	return J(nu1,nu2)

#--> bias is not needed of course...

#m13 = np.zeros((Nmax+1),dtype=complex)
#for i in range(Nmax+1):
#	m13[i] = M13(-0.5 * etam[i])

#m13outfile = np.zeros(2*(Nmax+1));
#for i in range(Nmax+1):
#        m13outfile[i] = m13.real[i]
#for i in range(Nmax+1,2*(Nmax+1)):
#	m13outfile[i] = m13.imag[i-Nmax-1];

#np.savetxt('M13oneline_N128.dat',m13outfile)


#m22mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
#for j1 in range(Nmax+1):
#	for j2 in range(Nmax+1):
#		if j1 - Nmax/2 < 1 :
#			m22mat[j1][j2] = M22(-0.5*etam[j1],-0.5*etam[j2])
#		else:
#			m22mat[j1][j2] = np.conjugate(m22mat[Nmax - j1][Nmax - j2])


#mout = np.zeros((Nmax+1)*(Nmax+1),dtype=complex)
#for i in range(Nmax+1):
#    for j in range(Nmax+1):
#        mout[i * (Nmax+1)+j] = m22mat[i][j]

#mout_red = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1))//2, dtype=complex)

#for i in range(Nmax+1):
#    for j in range(i+1):
#        mout_red[i+(2*(Nmax+1)-1-j)*j//2] = m22mat[i][j]

#moutoneline = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1)))
#for i in range((Nmax+1)*(Nmax+2)//2):
#    moutoneline[i] = mout_red.real[i]
#for i in range((Nmax+1)*(Nmax+2)//2,(Nmax+1)*(Nmax+2)):
#    moutoneline[i] = mout_red.imag[i-(Nmax+1)*(Nmax+2)//2]

#np.savetxt('M22oneline_N128.dat',moutoneline)


m12mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
for j1 in range(Nmax+1):
	for j2 in range(Nmax+1):
		if j1 - Nmax/2 < 1 :
			m12mat[j1][j2] = M12(-0.5*etam[j1],-0.5*etam[j2])
		else:
			m12mat[j1][j2] = np.conjugate(m12mat[Nmax - j1][Nmax - j2])


mout = np.zeros((Nmax+1)*(Nmax+1),dtype=complex)
for i in range(Nmax+1):
    for j in range(Nmax+1):
        mout[i * (Nmax+1)+j] = m12mat[i][j]

mout_red = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1))//2, dtype=complex)

for i in range(Nmax+1):
    for j in range(i+1):
        mout_red[i+(2*(Nmax+1)-1-j)*j//2] = m12mat[i][j]

moutoneline = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1)))
for i in range((Nmax+1)*(Nmax+2)//2):
    moutoneline[i] = mout_red.real[i]
for i in range((Nmax+1)*(Nmax+2)//2,(Nmax+1)*(Nmax+2)):
    moutoneline[i] = mout_red.imag[i-(Nmax+1)*(Nmax+2)//2]

#np.savetxt('M12oneline_N128.dat',moutoneline)
#np.savetxt('M12oneline_N256.dat',moutoneline)
np.savetxt('M12oneline_N512.dat',moutoneline)

#--> BELOW HERE IT IS JUST USELESS STUFF FOR ME...

# Nmax = 256
#b = -1.6000001
# kmax = 1.e2
# k0 = 0.00005

# Delta = log(kmax/k0) / (Nmax - 1)
# jsNm = np.arange(-Nmax/2,Nmax/2+1,1)
#etam = b + 2*1j*pi*(jsNm)/Nmax/Delta

#m13 = np.zeros((Nmax+1),dtype=complex)
#b13 = np.zeros((Nmax+1),dtype=complex)
#for i in range(Nmax+1):
#	m13[i] = M13v2(-0.5 * etam[i])

#m13outfile = np.zeros(2*(Nmax+1));
#for i in range(Nmax+1):
#        m13outfile[i] = m13.real[i]
#for i in range(Nmax+1,2*(Nmax+1)):
#	m13outfile[i] = m13.imag[i-Nmax-1];

#np.savetxt('IFG2oneline_N128.dat',m13outfile)

#mId2d2mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)
#for j1 in range(Nmax+1):
#    for j2 in range(Nmax+1):
#            if j1 - Nmax/2 < 1 :
#                    mId2d2mat[j1][j2] = J(-0.5*etam[j1],-0.5*etam[j2])
#            else:
#                    mId2d2mat[j1][j2] = np.conjugate(mId2d2mat[Nmax - j1][Nmax - j2])

#mId2d2out_red = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1))//2, dtype=complex)
#for i in range(Nmax+1):
#    for j in range(i+1):
#        mId2d2out_red[i+(2*(Nmax+1)-1-j)*j//2] = mId2d2mat[i][j]

#mId2d2outoneline = np.zeros(((Nmax+1)*(Nmax+1)+(Nmax+1)))
#for i in range((Nmax+1)*(Nmax+2)//2):
#        mId2d2outoneline[i] = mId2d2out_red.real[i]
#for i in range((Nmax+1)*(Nmax+2)//2,(Nmax+1)*(Nmax+2)):
#        mId2d2outoneline[i] = mId2d2out_red.imag[i-(Nmax+1)*(Nmax+2)//2]

#np.savetxt('M22basiconeline_N128.dat',mId2d2outoneline)
print("PT matrices succesfully recomputed!")
