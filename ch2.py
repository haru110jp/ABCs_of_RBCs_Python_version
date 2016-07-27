#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     26/07/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 2,  p32. Python Version.
#-------------------------------------------------------------------------------
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

Kappa = 4
Theta = 0.36
Roe = 0.6 # 0.9 in the original code. But it may be a typo.
N = 200 # The number of repetitions
Runs = 3 # The number of runs

K0 = np.array([0.8*Kappa**(1 / (1 - Theta)), Kappa**(1 / (1 - Theta)), \
    1.2 * Kappa**(1 / (1 - Theta))]) # Expand K0 when RUNS != 3
Lambda = np.ndarray(shape=(N, Runs))
Lambda[0, :] = np.ones(Runs)
K = np.ndarray(shape=(N, Runs))
K[0, :] = K0
Y = np.ndarray(shape=(N, Runs))
Y[0, :] = Lambda[0, :]  * (K[0, :] ** Theta) * (65 ** (1 - Theta)) # The production func.
Itr = [k + 1 for k in range(N-1)]
for k in Itr:
    Lambda_k = (1 - Roe) + Roe * Lambda[k-1, :] + 0.04 * (np.random.rand(Runs) - 0.5) # Not 0,02
    Lambda[k] = Lambda_k
    K_k = Kappa * Lambda_k * K[k-1, :] ** Theta
    K[k] = K_k
    Y_k = Lambda_k * (K[k-1, :] ** Theta) * (65 ** (1 - Theta))
    Y[k] = Y_k

plt.subplot(2, 1, 1)
plt.plot(K)

plt.subplot(2, 1, 2)
plt.plot(Y)

plt.show()
