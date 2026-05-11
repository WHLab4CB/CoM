#!/usr/bin/env python3
### h
### 2025.7.19

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import torch
import torch.nn as nn
import torch.optim as optim

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

def index(i, j):
    return i * Ny + j

# Grid parameters
Nx, Ny = 400, 400
xmin, xmax = -1.7, 1.2
ymin, ymax = -0.4, 2.1
u1x = np.zeros(2)

x = np.linspace(xmin, xmax, Nx)
y = np.linspace(ymin, ymax, Ny)
dx = x[1] - x[0]
dy = y[1] - y[0]
X, Y = np.meshgrid(x, y, indexing='ij')
U = muller_pot(X,Y)

# Define regions A and B (reactant/product)
reactant_center = [-0.27,1.73]
product_center = [0.84,0.00]
radius = 0.1

A_region = (X-reactant_center[0])**2 + (Y -reactant_center[1])**2 < radius**2 
B_region = (X-product_center[0])**2 + (Y -product_center[1])**2 < radius**2 

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


### plot the committor value of mueller
plt.figure(figsize=(8, 6))

### c_values
vmin = 0.0
vmax = 1.0
a=plt.contourf(X, Y, q_pred, levels=400, cmap='coolwarm',vmin=vmin,vmax=vmax)
plt.clim(vmin,vmax)
plt.colorbar(a,ticks=np.linspace(0.0,1.0,11))
plt.contour(X, Y, U, levels=[-100,-50,0,50,100], colors='darkgrey', linewidths=2,alpha=0.7)
iso_c_list = [5e-4,1e-3,1e-2,3e-2,0.1,0.5,0.9,0.95,0.98,0.995]
contour=plt.contour(X, Y, q_pred, levels=iso_c_list, colors='green', linewidths=0.5)
plt.clabel(contour, inline=True, fontsize=9, fmt='%.4f')

### plot points
points_list_all = []
points_list_all.append(np.loadtxt("../data/reactant_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_5.00e-04_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_1.00e-03_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_1.00e-02_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_3.00e-02_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_1.00e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_5.00e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_9.00e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_9.50e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_9.80e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/committor_9.95e-01_samples.dat"))
points_list_all.append(np.loadtxt("../data/product_samples.dat"))
n_sample_mile = len(points_list_all)
plt_color_list = ['green','yellow','purple','darkorange','black','chocolate','white','green','yellow','purple','darkorange','black','chocolate','white']
for i in range(n_sample_mile):
    for j in range(len(points_list_all[i])):
        plt.scatter(points_list_all[i][j][0],points_list_all[i][j][1],color=plt_color_list[i],s=1)

#plt.scatter(X[A_region], Y[A_region], label='Reactant region (q=0)')
#plt.scatter(X[B_region], Y[B_region], label='Product region (q=1)')

r_circle=plt.Circle((reactant_center[0], reactant_center[1]), radius, fill=False, edgecolor='orange', linestyle='-', linewidth=1)
p_circle=plt.Circle((product_center[0], product_center[1]), radius, fill=False, edgecolor='orange', linestyle='-', linewidth=1)
plt.gca().add_artist(r_circle)
plt.gca().add_artist(p_circle)

plt.xlabel('x')
plt.ylabel('y')
#plt.legend()
plt.tight_layout()
#plt.show()
plt.savefig('sample_iso_c.png',dpi=600)

