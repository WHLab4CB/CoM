#!/usr/bin/env python3


import numpy as np
import heapq
from committor_numpy import NumpyCommittor

muller_params = {
    'A': [-200, -100, -170, 15],
    'a': [-1, -1, -6.5, 0.7],
    'b': [0, 0, 11, 0.6],
    'c': [-10, -10, -6.5, 0.7],
    'x0': [1, 0, -0.5, -1],
    'y0': [0, 0.5, 1.5, 1]
}

def muller_potential(position):
    x,y = position
    V = 0.0
    for i in range(4):
        dx = x - muller_params['x0'][i]
        dy = y - muller_params['y0'][i]
        V += muller_params['A'][i]*np.exp(muller_params['a'][i]*dx**2+ muller_params['b'][i]*dx*dy + muller_params['c'][i]*dy**2)

    return V

def muller_gradient(position):
    x,y = position
    grad = np.zeros(2)
    for i in range(4):
        dx = x - muller_params['x0'][i]
        dy = y - muller_params['y0'][i]
        exp_term = np.exp(muller_params['a'][i]*dx**2 + muller_params['b'][i]*dx*dy + muller_params['c'][i]*dy**2)

        dVdx = muller_params['A'][i]*exp_term * (2*muller_params['a'][i]*dx + muller_params['b'][i]*dy)

        dVdy = muller_params['A'][i]*exp_term * (muller_params['b'][i]*dx + 2*muller_params['c'][i]*dy)

        grad += np.array([dVdx, dVdy])
    return grad 


class ConstrainedSampler:
    def __init__(self, committor_model, k=5e3, k2=2e5, mass=1.0, gamma=10.0, temp=10.0, center=(-0.27,1.73), radius=0.1):
        self.committor = committor_model
        self.k = k
        self.k2 = k2
        self.mass = mass
        self.gamma = gamma
        self.temp = temp
        self.sqrt_2gammaT = np.sqrt(2*gamma*temp)/mass
        self.center = center
        self.radius = radius

    def combined_force(self, x, q0):
        f_muller = -muller_gradient(x)

        q = self.committor.forward(x)
        dq_dx = self.committor.gradient(x)
        if q < 1e-12 or q > (1-1e-12):
            return f_muller

        if q0 < 0.5:
            log_diff = np.log(q) - np.log(q0)
            prefactor = self.k * log_diff / q
            f_constr = -prefactor * dq_dx
        else:
            log_diff = np.log(1-q) - np.log(1-q0)
            prefactor = -self.k * log_diff / (1-q)
            f_constr = -prefactor * dq_dx

        return f_muller + f_constr
    

    def combined_force2(self, x):
        f_muller = -muller_gradient(x)

        cx,cy = self.center
        dx = x[0] - cx
        dy = x[1] - cy
        r = np.sqrt(dx**2 + dy**2)

        dV_dr = self.k2 * (r-self.radius)
        dV_dx = dV_dr * (dx / r)
        dV_dy = dV_dr * (dy / r)

        f_constr = -np.array([dV_dx, dV_dy])

        return f_muller + f_constr

    def minimize_to_surface(self, x0, q0, tol=1e-5, max_iter=500):
        x = x0.copy()
        for i in range(max_iter):
            force = self.combined_force(x, q0)
            step_size = 0.1/(1+i//100)

            x += step_size * force / np.linalg.norm(force+1e-8)

            q = self.committor.forward(x)
#            print(f"log_diff: {i}, {x}, {q:.2e}, {muller_potential(x):.2f}, {0.5*self.k*(np.log(q)-np.log(q0))**2:.2f}")
            if abs(np.log(q) - np.log(q0)) < tol:
                break
        return x

    def minimize_to_surface2(self, x0, tol=1e-5, max_iter=500):
        x = x0.copy()
        for i in range(max_iter):
            force = self.combined_force2(x)
            step_size = 0.1/(1+i//100)

            x += step_size * force / np.linalg.norm(force+1e-8)

            if np.linalg.norm(x-self.center) < tol:
                break
        return x

    def langevin_dynamics(self, x_init, q0, dt=1e-5, steps=5e4, save_interval=5000):
        x = x_init.copy()
        v = np.zeros(2)
        v += np.sqrt(self.temp/self.mass)*np.random.randn(2)
        samples = []
        sqrt_dt = np.sqrt(dt)

        for step in range(int(steps)):
            force = self.combined_force(x, q0)

            x += v*dt
            v += - self.gamma*dt/self.mass*v + force/self.mass*dt + self.sqrt_2gammaT * sqrt_dt *np.random.randn(2)

            if step % save_interval == 0:
                samples.append([x.copy()[0],x.copy()[1]])

        return np.array(samples)

    def langevin_dynamics2(self, x_init, dt=1e-5, steps=5e4, save_interval=5000):
        x = x_init.copy()
        v = np.zeros(2)
        v += np.sqrt(self.temp/self.mass)*np.random.randn(2)
        samples = []
        sqrt_dt = np.sqrt(dt)

        for step in range(int(steps)):
            force = self.combined_force2(x)

            x += v*dt
            v += - self.gamma*dt/self.mass*v + force/self.mass*dt + self.sqrt_2gammaT * sqrt_dt *np.random.randn(2)

            if step % save_interval == 0:
                samples.append([x.copy()[0],x.copy()[1]])
        
        return np.array(samples)

    def run_committor_sampling(self, q_targets, initial_pos):
        all_samples = {}

        for q0 in q_targets:
            print(f"\nTarget q = {q0:.2e}")
            print("Minimizing to surface...", end='')
            x_min = self.minimize_to_surface(initial_pos, q0)
            print(f"Reached {x_min} q={self.committor.forward(x_min):.3e}")

            print("Running dynamics...", end='')
            samples = self.langevin_dynamics(x_min,q0)
            print(f"Generated {len(samples)} samples")

            all_samples[q0] = samples
            filename = f"../data/committor_{q0:.2e}_samples.dat"
            np.savetxt(filename, samples, fmt='%.6f')

        return all_samples
    
    def screen_committor_sampling(self, q_targets, initial_pos, xt_c_list,q_prec=0.1,dist_cut=0.2,target_points=10):
        all_samples = {}

        q_pred = []
        for i_xt in range(len(xt_c_list)):
            q_pred.append(self.committor.forward(xt_c_list[i_xt]))
        q_value_list = np.column_stack([xt_c_list, q_pred])

        for target_idx in range(len(q_targets)):
            q0 = q_targets[target_idx]
            print(f"\nTarget q = {q0:.2e}")
            print("Minimizing to surface...", end='')
            x_min = self.minimize_to_surface(initial_pos, q0)
            print(f"Reached {x_min} q={self.committor.forward(x_min):.3e}")

            print("Screening from C_value.dat...", end='')
            if 0 < target_idx < len(q_targets)-1:
                q_prev = q_targets[target_idx - 1]
                q_next = q_targets[target_idx + 1]
            elif target_idx == 0:
                q_prev = 0.0
                q_next = q_targets[target_idx + 1] if len(q_targets) > 1 else 1.0

            else:
                q_prev = q_targets[target_idx - 1]
                q_next = 1.0
            
            q_min = q0 - q_prec * (q0 - q_prev)
            q_max = q0 + q_prec * (q_next - q0)
            
            dist_list = [] 
            for a in range(len(q_value_list)):
                dist_temp = np.sqrt((q_value_list[a][0]-x_min[0])**2+(q_value_list[a][1]-x_min[1])**2)
                if dist_temp < dist_cut:
                    if q_min < q_value_list[a][2] < q_max:
                        dist_q = np.abs(q_value_list[a][2] - q0)
                        dist_list.append([a,dist_q])
            sample_data = []
            if len(dist_list) < target_points:
                for a_idx in range(len(dist_list)):
                    sample_data.append([q_value_list[dist_list[a_idx][0]][0],q_value_list[dist_list[a_idx][0]][1]])
            else:
                sorted_dist = sorted(dist_list, key=lambda x: x[1])
                for a_idx in range(target_points):
                    sample_data.append([q_value_list[sorted_dist[a_idx][0]][0],q_value_list[sorted_dist[a_idx][0]][1]])

            if len(sample_data) == 0:
                print("Warning: No samples found within the specified q-range!")
                all_samples[q0] = np.array([])
                continue
        
            ### print
            print(f"Generated {len(sample_data)} samples")
            print("\n Saving results...")
            all_samples[q0] = sample_data

            filename_xy = f"../data/committor_{q0:.2e}_samples.dat" 
            np.savetxt(filename_xy, sample_data,fmt=['%.6f', '%.6f'])

        return all_samples

    def run_reac_sampling(self, initial_pos):

        print(f"\nTarget cener={self.center}")
        print("Minimizing to surface...", end='')
        x_min = self.minimize_to_surface2(initial_pos)
        print(f"Reached {x_min} radius={np.linalg.norm(x_min-self.center):.3e}")

        print("Running dynamics...", end='')
        samples = self.langevin_dynamics2(x_min)
        print(f"Generated {len(samples)} samples")

        filename = f"../data/reactant_samples.dat"
        np.savetxt(filename, samples, fmt='%.6f')

    def run_product_sampling(self, initial_pos):

        print(f"\nTarget cener={self.center}")
        print("Minimizing to surface...", end='')
        x_min = self.minimize_to_surface2(initial_pos)
        print(f"Reached {x_min} radius={np.linalg.norm(x_min-self.center):.3e}")

        print("Running dynamics...", end='')
        samples = self.langevin_dynamics2(x_min)
        print(f"Generated {len(samples)} samples")

        filename = f"../data/product_samples.dat"
        np.savetxt(filename, samples, fmt='%.6f')


if __name__ == "__main__":
    committor = NumpyCommittor('../data/network_params.npz')

    ini_pos = np.array([-0.5, 0.5])
    sampler = ConstrainedSampler(committor,k=5e3,k2=2e5,mass=1.0,gamma=10.0,temp=10.0,center=(-0.27,1.73),radius=0.1)
    sampler.run_reac_sampling(ini_pos)

    q_targets = [5e-4,1e-3,1e-2,3e-2,0.1,0.5,0.9,0.95,0.98,0.995]
    #all_samples = sampler.run_committor_sampling(q_targets, ini_pos)
    xt_c_list = np.loadtxt("../data/C_values.dat",usecols=(0,1))
    all_samples = sampler.screen_committor_sampling(q_targets, ini_pos, xt_c_list)
    f_out = open("../data/q0_ini.dat","w")
    for q0 in q_targets:
        x_min = sampler.minimize_to_surface(ini_pos, q0)
        print(f"{q0:.2e} {x_min[0]:.6f} {x_min[1]:.6f}",file=f_out)

    sampler = ConstrainedSampler(committor,k=5e3,k2=2e5,mass=1.0,gamma=10.0,temp=10.0,center=(0.84,0.00),radius=0.1)
    sampler.run_product_sampling(ini_pos)
