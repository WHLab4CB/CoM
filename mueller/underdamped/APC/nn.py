#!/usr/bin/env python3

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
from sklearn.model_selection import train_test_split

f_best_model = '../data/best_model.pth'
f_network_params = '../data/network_params.npz' 
f_loss = open("../data/loss.txt","w")

def export_model_to_numpy(model,filename=f_network_params):
    params = {}
    for name, param in model.named_parameters():
        params[name] = param.detach().cpu().numpy()
    np.savez(filename, **params)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

data = np.loadtxt('../data/C_values.dat')
x = data[:,:2].astype(np.float32)
c = data[:,2].astype(np.float32)

epsilon = 1e-15
c = np.clip(c,epsilon,1-epsilon)

x_train, x_test, c_train, c_test = train_test_split(x, c, test_size=0.3, random_state=42)

x_train_tensor = torch.from_numpy(x_train).to(device)
c_train_tensor = torch.from_numpy(c_train).to(device)
x_test_tensor = torch.from_numpy(x_test).to(device)
c_test_tensor = torch.from_numpy(c_test).to(device)

class CommittorNN(nn.Module):
    def __init__(self):
        super(CommittorNN, self).__init__()

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

def custom_loss(q_pred, c_true):
    term1 = torch.mean((torch.log10(q_pred+epsilon) - torch.log10(c_true+epsilon))**2)
    term2 = torch.mean((torch.log10(1-q_pred+epsilon) - torch.log10(1-c_true+epsilon))**2)
    return term1 + term2

model = CommittorNN().to(device)
optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
train_dataset = TensorDataset(x_train_tensor, c_train_tensor)
train_loader = DataLoader(train_dataset, batch_size=32,shuffle=True)

best_loss = float('inf')
patience = 50
no_improve = 0

for epoch in range(2000):
    model.train()
    train_loss = 0.0
    for batch_x, batch_c in train_loader:
        optimizer.zero_grad()
        outputs = model(batch_x)
        loss = custom_loss(outputs, batch_c)
        if torch.isnan(loss):
            continue
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
        optimizer.step()
        train_loss += loss.item() * batch_x.size(0)

    model.eval()
    with torch.no_grad():
        test_outputs = model(x_test_tensor)
        test_loss = custom_loss(test_outputs, c_test_tensor).item()

    print(f'Epoch {epoch+1:3d} | Train Loss: {train_loss/len(train_loader.dataset):.6f} | Test Loss: {test_loss:.6f}',flush=True,file=f_loss)

    if test_loss < best_loss:
        best_loss = test_loss
        torch.save(model.state_dict(), f_best_model)
        no_improve = 0
        export_model_to_numpy(model)
    else:
        no_improve += 1
        if no_improve >= patience:
            print(f'Early stopping after {epoch+1} epochs')
            break

print('Training Complete. Best model saved to best_model.pth')
