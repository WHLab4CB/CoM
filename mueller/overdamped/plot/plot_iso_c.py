#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.interpolate import griddata

# Grid parameters
Nx, Ny = 1000, 1000
xmin, xmax = -1.7, 1.2
ymin, ymax = -0.4, 2.1
kBT = 10.0

# Define regions A and B (reactant/product)
reactant_center = [-0.27,1.73]
product_center = [0.84,0.00]
radius = 0.1


### mueller potential
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


x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing='ij')
U = muller_pot(X,Y)


A_region = (X - reactant_center[0])**2 + (Y - reactant_center[1])**2 < radius**2
B_region = (X - product_center[0])**2 + (Y - product_center[1])**2 < radius**2

### -------------- ref -----------------------------------
# Initialize matrix and RHS vector
u1x = np.zeros(2)
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

### ------------- nn -------------------------------------
epsilon = 1e-15
class CommittorNet(nn.Module):
    def __init__(self):
        super(CommittorNet, self).__init__()
        self.fc = nn.Sequential(
            nn.Linear(2,8),
            nn.Tanh(),
            nn.Linear(8,8),
            nn.Tanh(),
            nn.Linear(8,8),
            nn.Tanh(),
            nn.Linear(8,1),
            nn.Sigmoid()
        )

        for layer in self.fc:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.constant_(layer.bias, 0.1)

        self.register_buffer('reactant_center',torch.tensor([-0.27,1.73],dtype=torch.float32))
        self.register_buffer('product_center',torch.tensor([0.84,0.00],dtype=torch.float32))
        self.radius = 0.1

    def smooth_transition(self,x,center):
        d2 = torch.sum((x-center)**2, dim=1, keepdim=True)
        r_plus_sq = (self.radius + 0.02)**2
        return 0.5 - 0.5*torch.tanh(1000*(d2-r_plus_sq))

    def forward(self, x):
        rhoA = self.smooth_transition(x, self.reactant_center)
        rhoB = self.smooth_transition(x, self.product_center)

        u = self.fc(x)
        q = (1-rhoA-rhoB)*u + rhoB

        return torch.clamp(q.squeeze(),epsilon,1-epsilon)


model = CommittorNet()
model.load_state_dict(torch.load("../data/best_model.pth"))
XY = torch.tensor(np.stack([X.ravel(), Y.ravel()], axis=1), dtype=torch.float32)

with torch.no_grad():
    q_pred = model(XY).numpy().reshape(X.shape)
### ---------------------------------------------------------------


### ----- plot the committor value of mueller ---------------------
plt.figure(figsize=(8, 6))

### c_values
vmin = 0.0
vmax = 1.0
a=plt.contourf(X, Y, q_pred, levels=100, cmap='coolwarm',vmin=vmin,vmax=vmax)
plt.clim(vmin,vmax)
plt.colorbar(a,ticks=np.linspace(0.0,1.0,11))
plt.contour(X, Y, U, levels=[-100,-50,0,100], colors='darkgrey', linewidths=2,alpha=0.7)
iso_c_list = [3e-4,1e-3,1e-2,3e-2,1e-1,0.5,0.9,0.95,0.98,0.995]
contour=plt.contour(X, Y, q_pred, levels=iso_c_list, colors='green', linewidths=0.5)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.3e')

contour_ref=plt.contour(X, Y, q, levels=iso_c_list, colors='k', linewidths=0.5,zorder=2)
plt.clabel(contour_ref, inline=True, fontsize=9, fmt='%.3e')


# plot R and P
circle = plt.Circle(reactant_center, radius, fill=True, color='blue', 
                   linewidth=1, label=f'Reaction Region (r={radius})')
plt.gca().add_patch(circle)
circle = plt.Circle(product_center, radius, fill=True, color='red', 
                   linewidth=1, label=f'Product Region (r={radius})')
plt.gca().add_patch(circle)

plt.title('Committor function')
plt.xlabel('x')
plt.ylabel('y')
#plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('iso_c.png',dpi=600)


### ----- plot the error of committor value of mueller ---------------------
plt.figure(figsize=(8, 6))

### c_values
vmin = -0.1
vmax = 0.1

diff = q_pred - q
mask = (diff >= vmin) & (diff <= vmax)
diff_masked = np.where(mask, diff, np.nan)

a = plt.contourf(X, Y, diff_masked, levels=200, vmin=vmin, vmax=vmax, cmap='RdBu')
plt.clim(vmin,vmax)
plt.colorbar(a)
plt.contour(X, Y, U, levels=[-100,-50,0,100], colors='darkgrey', linewidths=2,alpha=0.7)

# plot R and P
circle = plt.Circle(reactant_center, radius, fill=True, color='blue', 
                   linewidth=1, label=f'Reaction Region (r={radius})')
plt.gca().add_patch(circle)
circle = plt.Circle(product_center, radius, fill=True, color='red', 
                   linewidth=1, label=f'Product Region (r={radius})')
plt.gca().add_patch(circle)


plt.title('Error committor values')

###set x and y axis
plt.xlim(-1.7, 1.2)
plt.ylim(-0.4, 2.1)

plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig('error_c.png',dpi=600)


