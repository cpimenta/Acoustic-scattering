# -*- coding: utf-8  -*-
"""
===========================================================
Implement ChebyShev basis to solve biharmonic problem.
-----------------------------------------------------------
Cristiano Pimenta
William R. Wolf
===========================================================
"""
import numpy as np
from scipy.linalg import eig, norm
from chebPy import cheb
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import time
from scipy.interpolate import BarycentricInterpolator as bi

# Number of points in domain
N = 320

D,x = cheb(N)

# For chord = 1.0 
fact = 0.5
x = ( x + 1.0 ) * fact
D = D / fact

D2 = np.dot(D,D)

I = np.eye(N+1,N+1)
Z = np.zeros((N+1,N+1))

# Build the matrix -> System of ODE 
L = np.vstack((np.hstack((D2,-I)), np.hstack((Z,D2))))
R = np.vstack((np.hstack((Z,Z)), np.hstack((I,Z))))

Np = len(x)
# Boundary conditions

# Position 
x0 = np.argmin(x)
x1 = np.argmax(x)

# Set Direchlet conditions w1(x0) = u(0)
L[x0,:] = 0.0
R[x0,:] = 0.0
L[x0,x0] = 1.0

# Set Neuman conditions w1'(0) = u'(0)
L[x0+1,:]= 0.0
R[x0+1,:] = 0.0
L[x0+1,:Np] = D[x0,:]

# Set w2(1) = 0
L[x1+Np,:] = 0.0 
R[x1+Np,:] = 0.0 
L[x1+Np,x1+Np] = 1.0 

# set w2'(1) = 0
L[x1+Np-1,:] = 0.0 
R[x1+Np-1,:] = 0.0 

L[x1+Np-1,Np:2*Np] = D[x1,:] 

np.set_printoptions(precision=3)

# Solve the Eigenproblem
lam, mode = eig(L,R)

lam = np.abs(lam) 
mode = np.real(mode)
lam = lam**(0.25)
idx = lam.argsort()
beta = lam[idx]

np.set_printoptions(formatter={'float': '{:.4f}'.format})

m = mode[:,idx]

# Normalise Eigenvalues
nmodes = N+1;
for imode in range(nmodes):
    m[:Np,imode] = m[:Np,imode]/norm(m[:Np,imode])

# Prepare to interpolate modal basis into BEM elements
xBEM,yBEM = np.loadtxt('geometricBEM.txt')
NpsBEM = len(xBEM)
uBEM = np.zeros((NpsBEM, nmodes))

"""
========================================================================
---------------- Barycentric Interpolation -  Burret & Treften ---------
========================================================================
"""
savemodes = 100
for nmodes in range(savemodes): 
    # initialize the object given the weights
    B = bi(x,m[:Np,nmodes])  
    # Modes interpolated
    uBEM[:NpsBEM/2,nmodes] = B(xBEM[:NpsBEM/2])
    uBEM[NpsBEM/2:,nmodes] = B(xBEM[NpsBEM/2:])
    aux = np.sum(uBEM[:NpsBEM/2,nmodes]*uBEM[:NpsBEM/2,nmodes] )
    aux2 = (NpsBEM/2)/aux
    uBEM[:,nmodes] *= np.sqrt(aux2)

np.savetxt('beta.txt', beta[:savemodes],fmt='%0.4f', newline='\n')
np.savetxt('modes.txt', uBEM,fmt='%0.4f', newline='\n')

"""
========================================================================
------------------------ Plot Reseults ---------------------------------
========================================================================
"""
plt.close('all')

#Plot using BEM information
fig = plt.figure(figsize=(15,8))
gs = gridspec.GridSpec(3,4,wspace=0.4,hspace=0.4)
uu = 0.0
for imode in range(12): 
    uu = uBEM[:400,imode]
    # uu = m[Np:2*Np,imode]
    ax = plt.subplot(gs[imode])
    ax.plot(xBEM[:400],uu,'b',linewidth=2.0) 
    plt.title(r'Modo = %i , $\beta$ =%0.4f'%(imode ,beta[imode]))
    plt.axis([0,1,-2.0,2.0])
    plt.grid(True)

plt.show()

