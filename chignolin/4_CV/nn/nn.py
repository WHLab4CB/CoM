#!/usr/bin/env python3

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
from sklearn.model_selection import train_test_split

def export_model_to_numpy(model,filename='./network_params.npz'):
    params = {}
    for name, param in model.named_parameters():
        params[name] = param.detach().cpu().numpy()
    np.savez(filename, **params)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
torch.set_default_dtype(torch.float64)

f_loss = open("./loss.txt","w")

data = np.loadtxt('./C_values.dat')
x = data[:,:4].astype(np.float64)
c = data[:,5].astype(np.float64)

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
            nn.Linear(4,96),
            nn.Tanh(),
            nn.Linear(96,96),
            nn.Tanh(),
            nn.Linear(96,96),
            nn.Tanh(),
            nn.Linear(96,1),
            nn.Sigmoid()
        )

        for layer in self.fc:
            if isinstance(layer, nn.Linear):
                nn.init.xavier_normal_(layer.weight)
                nn.init.constant_(layer.bias, 0.1)

        self.register_buffer('reactant_center',torch.tensor([15.0],dtype=torch.float64))
        self.register_buffer('product_center',torch.tensor([40.0],dtype=torch.float64))
        #self.radius= 0.2

    def smooth_transition_r(self,x,center):
        sum_x = torch.sum(x,dim=1,keepdim=True)
        d2 = sum_x-center+0.2
        return 0.5 - 0.5*torch.tanh(1000*(d2))

    def smooth_transition_p(self,x,center):
        sum_x = torch.sum(x,dim=1,keepdim=True)
        d2 = sum_x-center+0.2
        return 0.5 + 0.5*torch.tanh(1000*(d2))
    
    def forward(self, x):
        rhoA = self.smooth_transition_r(x, self.reactant_center)
        rhoB = self.smooth_transition_p(x, self.product_center)
        
        x_inv = 1.0/(x+1.0e-15)
        u = self.fc(x_inv)
        q = (1-rhoA-rhoB)*u + rhoB

        return torch.clamp(q.squeeze(),epsilon,1-epsilon)

def custom_loss(q_pred, c_true):
    term1 = torch.mean((torch.log10(q_pred+epsilon) - torch.log10(c_true+epsilon))**2)
    term2 = torch.mean((torch.log10(1-q_pred+epsilon) - torch.log10(1-c_true+epsilon))**2)
#    term1 = torch.mean((q_pred-c_true)**2)
    return term1 + term2

model = CommittorNN().to(device)
optimizer = optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-4)
train_dataset = TensorDataset(x_train_tensor, c_train_tensor)
train_loader = DataLoader(train_dataset, batch_size=128,shuffle=True)

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
        torch.save(model.state_dict(), './best_model.pth')
        no_improve = 0
        export_model_to_numpy(model)
    else:
        no_improve += 1
        if no_improve >= patience:
            print(f'Early stopping after {epoch+1} epochs')
            break

print('Training Complete. Best model saved to best_model.pth')

