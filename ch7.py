#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     09/08/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 7,  p178. Python Version.
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
kbar = 12.6695
hbar = 0.3335
ybar = kbar**theta * hbar ** (1-theta)
cbar = ybar - delta * kbar
aa = theta* ybar/ kbar + 1 - delta # "aa" frequently appears in p 152
a = np.zeros((3, 3))
a[0,0] = -1/(2*cbar*cbar)*aa*aa - 1/(2*cbar)*theta*(1-theta)*ybar/(kbar*kbar)
a[0,1] = 1/(2*cbar*cbar)*aa
a[1,0] = 1/(2*cbar*cbar)*aa
a[0,2] = -1/(2*cbar*cbar)*aa*(1-theta)*ybar/hbar
a[0,2] = a[0,2] + 1/(2*cbar)*theta*(1-theta)*ybar/(kbar*hbar)
a[2,0] = a[0,2]
a[1,1] = -1/(2*cbar*cbar)
a[1,2] = 1/(2*cbar*cbar)*(1-theta)*ybar/hbar
a[2,1] = a[1,2]
a[2,2]=-1/(2*cbar*cbar)*(1-theta)*ybar/hbar*(1-theta)*ybar/hbar
a[2,2]=a[2,2]-1/(2*cbar)*theta*(1-theta)*ybar/(hbar*hbar)
a[2,2]=a[2,2]-A/(2*(1-hbar)*(1-hbar))
x = np.array([kbar, kbar, hbar])

m = np.zeros((4, 4))


mm1 = 1/cbar*(theta*ybar/kbar+1-delta)
mm2 = (1-theta)*ybar/(cbar*hbar)-A/(1-hbar)
m[0,0] = log(kbar**theta*hbar**(1-theta) - delta*kbar)+A*log(1-hbar) - mm1*kbar+kbar/cbar-mm2*hbar + \
         np.dot(np.dot(x, a),x)
m[0,1] = mm1/2 - np.dot(a[0:3,0],x)
m[1,0] = m[0,1]
m[0,2] = -1/(2*cbar)-np.dot(a[0:3,1],x)
m[2,0] = m[0,2]
m[0,3] = mm2/2-np.dot(a[0:3,2],x)
m[3,0] = m[0,3]
m[1:4,1:4] = a
AA = np.array([[1, 0],[0, 0]])
B = np.array([[0, 0],[1, 0]])
R=m[0:2, 0:2]
Q=m[2:4, 2:4]
W=m[0:2,2:4].T
P = np.array([[1, 0],[0, 1]])

for i in range(1000):
    z = Q + beta*np.dot(np.dot(B.T, P),B)
    zinv = np.linalg.inv(z)
    z2 = beta * np.dot(np.dot(AA,P),B) + W.T
    P = R + beta*np.dot(np.dot(AA, P), AA) - np.dot(np.dot(z2, zinv), z2.T)


F = -np.dot(zinv, (W +beta* np.dot(np.dot(B.T, P), AA)))
