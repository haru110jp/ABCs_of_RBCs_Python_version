#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     27/07/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 4,  p67. Python Version.
#-------------------------------------------------------------------------------
from __future__ import division
from math import log
import numpy as np
import matplotlib.pyplot as plt

# Set initial conditions.
grids = 200 # The number of grids
vlast = np.zeros(grids)
k0 = np.linspace(0.06, 6, grids)
beta = 0.98
delta = 0.1
theta = 0.36
numits = 240 # The number of iterations
v = np.ndarray(grids)
kt1 = np.ndarray(grids)

# Define the function which calculates the value function.
def valfun(k):
    index_k = np.where(k0 == k)[0][0]
    # Calculate the consumption, see p60.
    ct = kt ** theta - k + (1 - delta) * kt
    if ct <= 0:
        val = -10000000 # To exclude too small ct.
    else:
        val = log(ct) + beta * vlast[index_k]
    return val

# Plot the initial value function.
plt.subplot(2, 1, 1)
plt.plot(k0, vlast)

# Begin the recursive calculations.

for t in range(numits):
    for i in range(grids): # Which value of k we focus on?
        valfun_cands = np.ndarray(grids)
        kt = k0[i]
        for j in range(grids): # Find k which maximizes the utility.
            next_t = k0[j]
            valfun_cands[j] = valfun(next_t)
        val_max_k = np.argmax(valfun_cands)
        ktp1 = k0[val_max_k]
        v[i] = valfun_cands[val_max_k]
        kt1[i] = ktp1
    if t/48 == round(t/48): # Change 48 suitably when numits != 240
        plt.subplot(2, 1, 1)
        plt.plot(k0, v)
    vlast = v

plt.subplot(2, 1, 2)
plt.plot(k0, kt1)
plt.plot(k0, k0) # 45 degrees line

plt.show()



