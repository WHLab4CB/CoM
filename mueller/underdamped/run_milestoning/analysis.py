#!/usr/bin/env python3


import numpy as np

f_mfpt_out = open("../data/mfpt.txt","w")

method = 1
tlist = np.loadtxt('../data/lifetime.dat')
q_series = [0.0,5e-4,1e-3,1e-2,3e-2,0.1,0.5,0.9,0.95,0.98,0.995,1.0]
nq = len(q_series)


### ------------- calculate A -> B -------------------------
if method == 1:
    K = np.zeros((nq,nq))
    for i in range(nq):
        if i==0:
            K[i,i+1] = 1.
        elif i==nq-1:
            K[i,0] = 1.
        else:
            K[i,i+1] = (q_series[i]-q_series[i-1])/(q_series[i+1]-q_series[i-1])
            K[i,i-1] = 1 - K[i,i+1]
else:
    K = np.loadtxt('../data/Kmat.dat.simulate')

print(K)

KT = K.T
w,v = np.linalg.eig(KT)

counter = 0
for i in range(nq):
    if abs(w[i]-1.0) < 1e-12:
        idx = i
        counter += 1

if counter > 1:
    print("Something is wrong in solving the Eigenequation")
    exit()
print(v[:,idx])
flux = np.real(v[:,idx])
acc_flux = 0.
for i in range(nq-1):
    acc_flux += flux[i]*tlist[i]
mfpt = acc_flux / flux[nq-1]
print(f"The A -> B MFPT is {mfpt:.4e}")

### print mfpt.txt
print("------- k matrix -------",file=f_mfpt_out)
for line in K:
    formatted_line = ' '.join(f"{x:.4f}" for x in line)
    print(formatted_line, file=f_mfpt_out)
    #print(' '.join(map(str, line)),file=f_mfpt_out)
print("------- lifetime -------",file=f_mfpt_out)
print(tlist,file=f_mfpt_out)
print("------- flux -------",file=f_mfpt_out)
print(flux,file=f_mfpt_out)
print(f"The A -> B MFPT is {mfpt:.4e}",file=f_mfpt_out)
print(f"The A -> B log10 MFPT is {np.log10(mfpt):6f}",file=f_mfpt_out)
### --------------------------------------------------------------

### ------- calculate B -> A -----------------------------------
if method == 1:
    q_series_temp = []
    for i in range(nq):
        q_series_temp.append(1.0-q_series[i])
    q_series_ba = q_series_temp[::-1]
    K_ba = np.zeros((nq,nq))
    for i in range(nq):
        if i==0:
            K_ba[i,i+1] = 1.
        elif i==nq-1:
            K_ba[i,0] = 1.
        else:
            K_ba[i,i+1] = (q_series_ba[i]-q_series_ba[i-1])/(q_series_ba[i+1]-q_series_ba[i-1])
            K_ba[i,i-1] = 1 - K_ba[i,i+1]
else:
    K_ba = K[::-1,::-1]
    for i in range(nq):
        K_ba[0][i] = 0
    K_ba[0][1] = 1.0
    for i in range(nq):
        K_ba[int(nq-1)][i] = 0
    K_ba[int(nq-1)][0] = 1.0
tlist_ba = tlist[::-1]

print(K_ba)

KT_ba = K_ba.T
w_ba,v_ba = np.linalg.eig(KT_ba)

counter = 0
for i in range(nq):
    if abs(w_ba[i]-1.0) < 1e-12:
        idx_ba = i
        counter += 1

if counter > 1:
    print("Something is wrong in solving the Eigenequation")
    exit()
print(v_ba[:,idx_ba])
flux_ba = np.real(v_ba[:,idx_ba])
acc_flux_ba = 0.
for i in range(nq-1):
    acc_flux_ba += flux_ba[i]*tlist_ba[i]
mfpt_ba = acc_flux_ba / flux_ba[nq-1]
print(f"The B -> A MFPT is {mfpt_ba:.4e}")

### print mfpt.txt
print("------- k_ba matrix -------",file=f_mfpt_out)
for line in K_ba:
    formatted_line = ' '.join(f"{x:.4f}" for x in line)
    print(formatted_line, file=f_mfpt_out)
    #print(' '.join(map(str, line)),file=f_mfpt_out)
print("------- lifetime_ba -------",file=f_mfpt_out)
print(tlist_ba,file=f_mfpt_out)
print("------- flux_ba -------",file=f_mfpt_out)
print(flux_ba,file=f_mfpt_out)
print(f"The B -> A MFPT is {mfpt_ba:.4e}",file=f_mfpt_out)
print(f"The B -> A log10 MFPT is {np.log10(mfpt_ba):6f}",file=f_mfpt_out)


