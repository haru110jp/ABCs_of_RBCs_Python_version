#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     08/08/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 6,  p142. Python Version.
#-------------------------------------------------------------------------------
from __future__ import division
from math import log, sqrt
import numpy as np
import matplotlib.pyplot as plt

# Set initial conditions.
theta = 0.36
delta = 0.025
beta = 0.99
A = 1.72

# Compute hbar.
B = ( 1 -(beta * delta * theta) / (1-beta * (1-delta)))
C = ( 1 / (1 - theta)) * B
D= A * C
E= 1 + D
hbar1=1 / E

# Compute kbar.
F = theta / ( 1 / beta - ( 1 - delta))
kbar1 = hbar1 * F**( 1 / (1-theta))

# Main
rbar= 1 / beta - ( 1 - delta)
ybar = kbar1 ** theta * hbar1** (1-theta)
cbar = ybar - delta * kbar1
A = np.array([0, -kbar1, 0, 0])
B = np.array([0, (1-delta)*kbar1, theta, -1])
C = np.array([[1, -1, -1/(1-hbar1), 0], [ybar, -cbar, 0, 0],[-1, 0, 1-theta, 0], [1, 0, 0, -1]])
D = np.array([0, 0, 1, 0])
F=0
G = 0
H = 0
J = np.array([0, -1, 0, beta*rbar])
K = np.array([0, 1, 0, 0])
L = 0
M = 0
N = 0.95
Cinv= np.linalg.inv(C)
a = F - np.dot (np.dot(J, Cinv), A)
b = - (np.dot(np.dot(J, Cinv), B) - G + np.dot(np.dot(K, Cinv), A))
c = - np.dot(np.dot(K, Cinv), B) + H
P1=(-b+sqrt(b**2-4*a*c))/(2*a)
P2=(-b-sqrt(b**2-4*a*c))/(2*a)
if abs(P1)<1:
    P=P1
else:
    P=P2

R = -np.dot(Cinv, (A*P + B))
Q = (np.dot(np.dot(J, Cinv), D) - L) * N + np.dot(np.dot(K, Cinv), D) - M
QD = np.kron(N, ( F - np.dot(np.dot(J, Cinv), A))) + (np.dot(J, R) + F * P + G - np.dot(np.dot(K, Cinv), A))
Q= Q / QD
S = np.dot(-Cinv, (A*Q+D))
