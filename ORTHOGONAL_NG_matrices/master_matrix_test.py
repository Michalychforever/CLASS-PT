import numpy as np

from sympy.parsing.mathematica import mathematica
from sympy import *

from mpmath import *

mp.dps = 32
mp.pretty = True

nu1 = var('nu1')
nu2 = var('nu2')

with open("I_equivalent_orthogonal_sympy.txt","r") as file:
    expr = file.read()

mathexpr=expr

M12master=lambdify([nu1,nu2],mathexpr,"mpmath")

def M12(nu1,nu2):
    return M12master(nu1,nu2)

b=-0.8 #--> I am still not sure that this is the correct bias choice... Here I just check that everything works...
kmax = 1.e2
k0 = 0.00005

#Nmaxtab=[128,256,512]

Nmaxtab=[32]

print("----------------------------------------\n")

for Nmax in Nmaxtab:

    Delta = log(kmax/k0) / (Nmax - 1)
    jsNm = np.arange(-Nmax/2,Nmax/2+1,1)
    etam = b + 2*1j*pi*(jsNm)/Nmax/Delta

    print("###################")
    print("Nmax is "+str(Nmax))
    print("###################\n")

    m12mat = np.zeros((Nmax+1,Nmax+1),dtype=complex)

    for j1 in range(Nmax+1):
        print(j1)
        for j2 in range(Nmax+1):
            if j1 - Nmax/2 < 1:
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

    np.savetxt('M12master_oneline_N'+str(Nmax)+'_TEST.dat',moutoneline)
    print("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

print("PT matrices for matter successfully computed!")
