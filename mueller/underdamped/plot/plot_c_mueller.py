#!/usr/bin/env python3
### plot committor values with (x,y,c)
### 2025.4.13


import matplotlib
matplotlib.use('agg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


fig,ax=plt.subplots(1,1)

x_0 = np.loadtxt("../data/C_values.dat",usecols=0)
y_0 = np.loadtxt("../data/C_values.dat",usecols=1)
c_0 = np.loadtxt("../data/C_values.dat",usecols=2)

# Define regions A and B (reactant/product)
reactant_center = [-0.27,1.73]
product_center = [0.84,0.00]
radius = 0.1

biga = [-200.0, -100.0, -170.0, 15.0]
a = [-1.0, -1.0, -6.5, 0.7]
b = [0.0, 0.0, 11.0, 0.6]
c = [-10.0, -10.0,-6.5, 0.7]
x0 = [1.0, 0.0, -0.5, -1.0]
y0 = [0.0, 0.5, 1.5, 1.0]

def muller_pot(x,y):
    V = np.zeros_like(x)
    for Ai, ai, bi, ci, xi, yi in zip(biga, a, b, c, x0, y0):
        V += Ai * np.exp(ai * (x - xi)**2 + bi * (x - xi)*(y - yi) + ci * (y - yi)**2)
    return V

def force(x, u1x):
    ene = np.zeros(4)
    xdif = np.zeros(2)

    u1x[0] = 0.
    u1x[1] = 0.
    for i in range(4):
        xdif[0] = x[0] - x0[i]
        xdif[1] = x[1] - y0[i]
        ene[i] = biga[i]*np.exp(a[i]*xdif[0]*xdif[0] + b[i]*xdif[0]*xdif[1] + c[i]*xdif[1]*xdif[1])
        u1x[0] += ene[i]*(2.*a[i]*xdif[0] + b[i]*xdif[1])
        u1x[1] += ene[i]*(2.*c[i]*xdif[1] + b[i]*xdif[0])

def index(i, j):
    return i * Ny + j

# Grid parameters
Nx, Ny = 1000, 1000
xmin, xmax = -1.7, 1.2
ymin, ymax = -0.4, 2.1
kBT = 10.0
u1x = np.zeros(2)

x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing='ij')
U = muller_pot(X,Y)

A_region = (X - reactant_center[0])**2 + (Y - reactant_center[1])**2 < radius**2
B_region = (X - product_center[0])**2 + (Y - product_center[1])**2 < radius**2

"""
### -------------- ref -----------------------------------
# Initialize matrix and RHS vector
N = Nx * Ny
A_mat = lil_matrix((N, N))
b_vec = np.zeros(N)

# Fill the matrix
for i in range(Nx):
    for j in range(Ny):
        k = index(i, j)

        if A_region[i, j]:
            A_mat[k, k] = 1.0
            b_vec[k] = 0.0
        elif B_region[i, j]:
            A_mat[k, k] = 1.0
            b_vec[k] = 1.0
        else:
            pos = [xmin+dx*i, ymin+dy*j]
            force(pos,u1x)
            ux = u1x[0]
            uy = u1x[1]

            im = max(i-1, 0)
            ip = min(i+1, Nx-1)
            jm = max(j-1, 0)
            jp = min(j+1, Ny-1)

            dx2 = dx**2
            dy2 = dy**2

            if i==0 or i==Nx-1:
                A_mat[k, index(ip,j)] += 2*kBT/dx2 if i==0 else 0
                A_mat[k, index(im,j)] += 2*kBT/dx2 if i==Nx-1 else 0
                A_mat[k,k] -= 2*kBT/dx2
            else:
                A_mat[k,index(ip,j)] = kBT/dx2-ux/(2*dx)
                A_mat[k,index(im,j)] = kBT/dx2+ux/(2*dx)

            if j==0 or j==Ny-1:
                A_mat[k,index(i,jp)] += 2*kBT/dy2 if j==0 else 0
                A_mat[k,index(i,jm)] += 2*kBT/dy2 if j==Ny-1 else 0
                A_mat[k,k] -= 2*kBT/dy2
            else:
                A_mat[k,index(i,jp)] = kBT/dy2-uy/(2*dy)
                A_mat[k,index(i,jm)] = kBT/dy2+uy/(2*dy)

            if i!=0 and i!= Nx-1:
                A_mat[k,k] -= 2*kBT/dx2
            if j!=0 and j!= Ny-1:
                A_mat[k,k] -= 2*kBT/dy2

# Solve the linear system
q_flat = spsolve(A_mat.tocsr(), b_vec)
q = q_flat.reshape((Nx, Ny))
### ------------------------------------------------------
"""

### ----------- C_values ---------------------------------
xi = np.linspace(min(x_0), max(x_0), 200)
yi = np.linspace(min(y_0), max(y_0), 200)
xi, yi = np.meshgrid(xi, yi)  

ci = griddata((x_0, y_0), c_0, (xi, yi), method='linear')
ci_positive = np.where(ci <= 0, np.nan, ci)

vmin = 0.0
vmax = 1.0 
levels = np.linspace(vmin,vmax,100)
#a=plt.contourf(xi,yi,ci.clip(max=200),400,cmap=cm.jet,levels=levels,zorder=1)
plt.contour(X, Y, U, levels=[-100,-50,0,50,100], colors='darkgrey', linewidths=2,alpha=0.7)
a=plt.contourf(xi,yi,ci_positive,400,cmap='coolwarm',levels=levels,zorder=1,vmin=vmin,vmax=vmax)
plt.clim(vmin,vmax)
plt.colorbar(a,ticks=np.linspace(0.0,1.0,11))


#plt.contour(xi,yi,ci,levels=20,linewidths=0.5,colors='black',zorder=1)
iso_c_list = [5e-4,1e-3,1e-2,3e-2,0.1,0.5,0.9,0.95,0.98,0.995]
contour=plt.contour(xi, yi, ci_positive, levels=iso_c_list, colors='purple', linewidths=0.5)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.3e')


# plot R and P
circle = plt.Circle(reactant_center, radius, fill=True, color='blue', 
                   linewidth=3, label=f'Reaction Region (r={radius})')
plt.gca().add_patch(circle)
circle = plt.Circle(product_center, radius, fill=True, color='red', 
                   linewidth=3, label=f'Product Region (r={radius})')
plt.gca().add_patch(circle)

plt.title('Committor function')
###set x and y axis
plt.xlim(-1.7, 1.2)
plt.ylim(-0.4, 2.1)

#plt.xticks(np.arange(-1.5,1.5,0.5),fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('x',fontsize=14)
plt.ylabel('y',fontsize=14)

plt.savefig('./c_mueller',bbox_inches="tight", dpi=600)
plt.close()
