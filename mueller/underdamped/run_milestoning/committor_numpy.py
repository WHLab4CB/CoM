#!/usr/bin/env python3


import numpy as np

class NumpyCommittor:
    def __init__(self, param_file):
        params = np.load(param_file)

        self.weights = [
            params['fc.0.weight'],
            params['fc.2.weight'],
            params['fc.4.weight'],
            params['fc.6.weight']
        ]
        self.biases = [
            params['fc.0.bias'],
            params['fc.2.bias'],
            params['fc.4.bias'],
            params['fc.6.bias']
        ]

        self.reactant_center = np.array([-0.27,1.73])
        self.product_center = np.array([0.84,0.00])
        self.radius = 0.1

    def _smooth_transition(self, x, center):
        d2 = np.sum((x - center)**2)
        r_plus_sq = (self.radius + 0.02)**2
        return 0.5 - 0.5*np.tanh(1000 *(d2-r_plus_sq))

    def _d_smooth_transition(self, x, center):
        d = x - center
        d2 = np.sum(d**2)
        r_plus_sq = (self.radius + 0.02)**2
        sech2 = 1 / np.cosh(1000 * (d2-r_plus_sq))**2
        return -1000*d*sech2

    def tanh(self, x):
        return np.tanh(x)

    def sigmoid(self, x):
        return 1/(1+np.exp(-x))

    def forward(self, x):
        x = np.asarray(x, dtype=np.float32)

        rhoA = self._smooth_transition(x, self.reactant_center)
        rhoB = self._smooth_transition(x, self.product_center)

        h = x.copy()
        for i in range(3):
            h = h @ self.weights[i].T + self.biases[i]
            h = self.tanh(h)

        h = h @ self.weights[3].T + self.biases[3]
        u = self.sigmoid(h)[0]

        q = (1-rhoA-rhoB)*u + rhoB
        return q

    def gradient(self, x):
        x = np.asarray(x, dtype=np.float32)
        rhoA = self._smooth_transition(x, self.reactant_center)
        rhoB = self._smooth_transition(x, self.product_center)
        activations = [x.copy()]
        pre_activations = []

        h = x.copy()
        for i in range(3):
            h = h @ self.weights[i].T + self.biases[i]
            pre_activations.append(h)
            h = self.tanh(h)
            activations.append(h)

        h = h @ self.weights[3].T + self.biases[3]
        pre_activations.append(h)
        u = self.sigmoid(h)[0]

        delta = self.sigmoid(h) * (1 - self.sigmoid(h))

        for i in reversed(range(3)):
            tanh_deriv = 1 - activations[i+1]**2
            delta = (delta @ self.weights[i+1]) * tanh_deriv

        dq_dx = delta @ self.weights[0]

        drhoA_dx = self._d_smooth_transition(x, self.reactant_center)
        drhoB_dx = self._d_smooth_transition(x, self.product_center)

        dq_dx = (1 - rhoA - rhoB) * dq_dx + (-drhoA_dx-drhoB_dx) * u + drhoB_dx

        return dq_dx

if __name__ == "__main__":
    model = NumpyCommittor('../data/network_params.npz')
    position = [-1.99160785, 1.54277105]

    q_value = model.forward(position)
    gradient = model.gradient(position)

    print(f"Committor value: {q_value:.2e}")
    print(f"Graident: {gradient}")
