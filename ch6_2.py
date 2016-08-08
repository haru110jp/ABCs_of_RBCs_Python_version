#-------------------------------------------------------------------------------
# Author:      Atsushi Yamagishi
# Created:     08/08/2016
# Licence:     MIT

# The ABCs of RBSs, Chapter 6,  p143. Python Version.
#-------------------------------------------------------------------------------
#function to approximate the ratio of the variance of a jump
#variables in the Hansen model to the variance of the shock to the
#technology process.  The interation limit is 1000.

#Sample: a = 0.954, b = 0.113, c = 0.205, d = 1.452, roe = 0.95
from __future__ import division

def varratio(a,b,c,d,roe,tol):
    vr = 0
    for i in range(1000):
        shortsum=0
        for j in range(i+1):
            shortsum += (a**j) * (roe**(i-j))

        increment = c * b * shortsum + d * roe**(i+1)
        vr = vr + increment*increment
#        if abs(increment) < tol:
#            print "tol achieved", increment, i
#            return vr + d*d
    return vr+d*d
