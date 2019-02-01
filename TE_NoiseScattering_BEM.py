# -*- coding: utf-8 -*-

"""
################################################################################
# Numerical solution of acoustic scattering by finite perforated elastic plates
# A. V. G. Cavalieri, W. R. Wolf and J. W. Jaworski - PRSA 2016
# -----------------------------------------------------------------------------
# Develped by Msc. Cristiano Pimenta and Prof. Dr. William Wolf
# Input parameters:
# - Modal basis -> Non-dimensional frequency and vibration modes(displacements) 
# - Acoustic wavenumber
# - alphaH -> Open area fraction
# - Omega -> Vacuum bending wave Mach number
# - epsilon -> Intrinsic fluid-loading parameter
# - field -> Set True to compute the acoustic field
################################################################################
"""

# Import libraries
import sys, os
import numpy as np
from scipy.special import *
from scipy import integrate, optimize
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import argparse
import matrixv1

def main(args):
    k0 = args.k0
    alphaH = args.alphaH
    omega = args.omega
    epsilon = args.epsilon
    field = args.field
    print('Epsilon used is: ' + str(epsilon))
    print('Omega used is: ' + str(omega))
    
    # Constants variables
    R = 1.0e-3          # Radius of porous
    kr = 4.0/np.pi  # Rayleigh conductivity
    gamma = 0.5772156649    # Constant of Euler
    print("alphaH/R = " + str(alphaH/R)) 
    
    # Geometry information
    Npsx = 400  # number of panel created on each surface
    Nps = 2 * Npsx + 2  # Total of panels
    thickness = 0.002 

    c = 1.0     # Chord of the plate
    dx = c / Npsx

    # Source position (x,y)
    z = [1.0, 0.004]

    xp, yp = np.zeros(Nps+1), np.zeros(Nps+1)

    # Coordenates points
    for i in xrange(Nps):
        if i <= Npsx:
            xp[i] = dx * i
            yp[i] = +thickness/2.0
        else:
            xp[i] = xp[Nps - i - 1]
            yp[i] = -thickness/2.0

    xp[Nps] = xp[0]
    yp[Nps] = yp[0]

    class gpanel:
        """This class contains all geometry information about to a panel"""
        def __init__(self,x1,y1,x2,y2):
            ''' Initializes the panel.
            Arguments
            x1, y1 --> Coordinates of the first point of the panel
            x2, y2 --> Coordinates of the first point of the panel
            '''
            self.x1, self.y1 = x1, y1
            self.x2, self.y2 = x2, y2
            self.xc, self.yc = (x1 + x2)/2, (y1 + y2)/2   # Center point of the panel
            self.length = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 )   # length of the panel
            # normal -  Important to check
            self.n1 = (y2 - y1)/self.length
            self.n2 = -(x2 - x1)/self.length


    #defining the panels
    panels = np.empty(Nps,dtype=object)
    for i in xrange(Nps):
        panels[i] = gpanel(xp[i], yp[i], xp[i+1], yp[i+1])

    x1 = [pi.x1 for pi in panels]
    y1 = [pi.y1 for pi in panels]
    x2 = [pi.x2 for pi in panels]
    y2 = [pi.y2 for pi in panels]
    xc = [pi.xc for pi in panels]
    yc = [pi.yc for pi in panels]
    n1 = [pi.n1 for pi in panels]
    n2 = [pi.n2 for pi in panels]
    ds = [pi.length for pi in panels]

    H = np.zeros((Nps, Nps), dtype = complex)
    G = np.zeros((Nps, Nps), dtype = complex)
    A = np.zeros((Nps, Nps), dtype = complex)



    # Fortran Wrapp for the calculus of matrix H, G 
    H,G = matrixv1.mntmat.hgmatrix(k0,x1,y1,x2,y2,xc,yc,n1,n2,ds,Nps)

    """
    ==================================================================================
    ----------------- Import Structural Modal Basis information ----------------------
                            Implement poroelastic materials
    ==================================================================================
    """
    beta = np.loadtxt('modalBasis/beta.txt')   # Read the Non-dimensional frequency
    phi = np.loadtxt('modalBasis/modes.txt')   # Read the modal basis 
    nm = len(beta)   # define number of modes

    nplate = Nps

    # Building the matrix D
    D = matrixv1.mntmat.poroelastic(nplate,k0, alphaH, epsilon, omega, n2, ds, beta, phi,Nps,nm)

    # Compute LHS
    A = H - np.dot(G, D)

    # Compute the RHS of the linear system (Acoustic Source)
    
    def quadrupole(k0, xc, yc, z):
        """ This function computs the model source quadrupole"""
        S = []
        for i in range(len(xc)):
        
            dx = xc[i] - z[0]
            dy = yc[i] - z[1]
            
            arg = k0 * np.sqrt(dx**2 + dy**2)
            
            Ha0 = hankel1(0,arg)  # Hankel of 1 order 0
            Ha1 = hankel1(1,arg)
            Ha2 = hankel1(2,arg)
            
            arg1_y2 = - 0.5 * (Ha0 - Ha2) * ( k0**2 * dx * dy / ( dx**2 + dy**2))
            arg2_y2 =  Ha1 * k0 * dx * dy / (( dx**2 + dy**2)**(3.0/2.0))
            
            ay2_dev =  (1j / 4.0 ) * (arg1_y2 + arg2_y2)
            
            S.append( ay2_dev )
        return S
   
    
    S = quadrupole(k0, xc, yc, z)

    # Solve the linear system and get the source values(pressure fluctuation) for each panel
    pl = np.linalg.solve(A,S)
    
    """
    #########################################################################################
    #  Compute the acoustic field for 360 observers
    #########################################################################################
    """
    
    No = 360   # Number of Observer

    # Position of each observer
    ro = 50.0
    theta_o = np.linspace(0,2.0*np.pi,No)
    center = 1.0

    xo1 = center + ro * np.cos(theta_o)
    xo2 = ro * np.sin(theta_o)
    pobs = np.zeros(No, dtype=complex) # Pressure at observer
    # Far field matrix H, G for the observers
    Hobs, Gobs = matrixv1.mntmat.hgobs(k0,xo1,xo2,x1,y1,x2,y2,n1,n2,ds,No,Nps)

    Sf = np.asarray(quadrupole(k0, xo1, xo2, z))

    pscat  = np.dot((np.dot(Gobs , D) - Hobs ) , pl)  # Scattered pressure
    pinc = -Sf   # Pressure from incident source
    pobs =  pscat + pinc

    outputs = "results"
    if not os.path.isdir(outputs):
        os.makedirs(outputs)
    save = os.path.join(outputs, "k%0.1f_alphaH_R%d_Omega_%0.3f.dat" %(k0,int(alphaH/R),omega))
    f=open(save, "wt")
    for i in range(No):
         f.write(str(np.degrees(theta_o[i])) + " " + str(np.abs(pobs[i])) + "\n")
    f.close()
    """
    ######################################################################################## 
    # Compute the noise scattered for a rectangular mesh to obtain the field around
    # Default is false, if is desired set up True 
    ########################################################################################
    """
    if field:
        # Define the domian
        nx = 200
        ny = 200
        lx = 3.0
        ly = 1.0
        xs = -1.0
        ys = -1.0

        # x = np.zeros((nx,ny))
        # y = np.zeros((nx,ny))

        # x[0,:] = xs
        # y[:,0] = ys
        # dx = lx/(nx-1)
        # dy = ly/(ny-1)
        # for i in range(1,nx):
        #     x[i,:] = x[i-1,:] + dx
        # for j in range(1,ny):
        #     y[:,j] = y[:,j-1] + dy
        

        xx = np.linspace(xs,lx,nx)
        yy = np.linspace(ys,ly,ny)

        x,y = np.meshgrid(xx,yy)

        z1 = z[0]
        z2 = z[1]

        dpdn = np.dot(D, pl)

        pfield= np.zeros((nx,ny),dtype=complex)
        pfield = matrixv1.mntmat.field(k0,z1,z2,x1,y1,x2,y2,n1,n2,ds,dpdn,pl,x,y,nx,ny,Nps)

        field_outs = "field_outputs"
        if not os.path.isdir(field_outs):
            os.makedirs(field_outs)
        savefield = os.path.join(field_outs, "k%0.1f_alphaH_R%d_Omega_%0.3f.dat" %(k0,int(alphaH/R),omega)) 
        f=open(savefield, "wt")  
        # Format tecplot
        f.write('TITLE = "mesh" \n')
        f.write('VARIABLES = "X", "Y", "p"\n')
        f.write('ZONE I = '+ str(nx) + ', J = '+str(ny)+', F=POINT \n')

        for j in range(ny):
            for i in range(nx):
                f.write( "{0:0.8E} \t {1:0.8E} \t {2:0.8E} \n" .format(x[i,j], y[i,j], pfield[i,j].real))

        f.close()
    
    # Plot the Directivity
    f=plt.figure(1)
    ax = f.add_axes([0.1, 0.1, 0.8, 0.8], projection='polar')
    ax.plot(theta_o, np.abs(pobs),linestyle='-', linewidth='2')
    #ax.set_rgrids([20, 40, 60,80], angle=0.)
    #ax.legend(["k = 1.0","k = 5.0","k = 10.0"], loc='lower center', prop={'size':16}, bbox_to_anchor=(1.1,-0.1))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("k0", type=float, help="Acoustic wavenumber")
    parser.add_argument("--alphaH", type=float, help="Open area fraction", default=0.0)
    parser.add_argument("--omega", type=float, help="Vacuum bending wave Mach number", default = 0.1)
    parser.add_argument("--epsilon", type=float, help="Intrinsic fluid-loading parameter", default = 0.0)
    parser.add_argument("--field", help="Flag to compute the field", default = False)
    
    args = parser.parse_args()

    main(args)
