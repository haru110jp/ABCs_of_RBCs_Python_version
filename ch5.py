#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     30/07/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 5,  p87. Python Version.
#-------------------------------------------------------------------------------
from __future__ import division
from math import log
import numpy as np
import matplotlib.pyplot as plt

# Set initial conditions.
grids = 40 # The number of grids
vlast1 = 20 * np.ones(grids)
vlast2 = 20 * np.ones(grids)
k0 = np.linspace(0.4, 16, grids)
kt11 = np.ndarray(grids)
kt12 = np.ndarray(grids)
beta = 0.98
delta = 0.1
theta = 0.36
A1 = 1.75
p11 = 0.9
p12 = 1 - p11
p21 = 0.4
p22 = 1 - p21
A2 = 0.75
numits = 250 # The number of iterations
v1 = np.ndarray(grids)
v2 = np.ndarray(grids)


# Define the function which calculates the value function.
def valfunsto(k):
    index_k = np.where(k0 == k)[0][0]
    # Calculate the consumption, see p60.
    ct = At * kt ** theta - k + (1 - delta) * kt
    if ct <= 0.001:
        val = log(0.001) + beta * (p1 * vlast1[index_k] + p2 * vlast2[index_k]) \
            + 200 * (ct - 0.001) # To exclude too small ct.
    else:
        val = log(ct) + beta * (p1 * vlast1[index_k] + p2 * vlast2[index_k])
    return val

# Plot the initial value function.
plt.subplot(2, 1, 1)
plt.plot(k0, vlast1)
plt.plot(k0, vlast2)

# Begin the recursive calculations.

for t in range(numits):
    for i in range(grids): # Which value of k we focus on?
        valfun_cands1 = np.ndarray(grids)
        valfun_cands2 = np.ndarray(grids)
        kt = k0[i]
        for j in range(grids): # Find k which maximizes the utility.
            next_k = k0[j]
            At = A1
            p1 = p11
            p2 = p12
            valfun_cands1[j] = valfunsto(next_k)
            At = A2
            p1 = p21
            p2 = p22
            valfun_cands2[j] = valfunsto(next_k)
        val_max1_k = np.argmax(valfun_cands1)
        val_max2_k = np.argmax(valfun_cands2)
        ktp11 = k0[val_max1_k]
        v1[i] = valfun_cands1[val_max1_k]
        kt11[i] = ktp11
        val_max2_k = np.argmax(valfun_cands2)
        ktp12 = k0[val_max2_k]
        v2[i] = valfun_cands2[val_max2_k]
        kt12[i] = ktp12
    if t/50 == round(t/50): # Change 50 suitably when numits != 250
        plt.subplot(2, 1, 1)
        plt.plot(k0, v1)
        plt.plot(k0, v2)
    vlast1 = v1
    vlast2 = v2

plt.subplot(2, 1, 2)
plt.plot(k0, kt11)
plt.plot(k0, kt12)
plt.plot(k0, k0) # 45 degrees line

plt.show()