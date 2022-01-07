import numpy as np
import sympy as sym

import matplotlib.pyplot as plt
import scipy.integrate   as integrate

def RHS_h(state,time):
    #VOP example from my own mind
    #homogenous equations
    y, y1dot, y2dot, y3dot = state

    return np.array([y1dot,y2dot,y3dot,-y2dot])


def analytic_soln(t,S0):
    y0,d1y0,d2y0,d3y0 = S0
    c1 = -d2y0
    c2 = -d3y0
    c3 = y0+d2y0
    c4 = d1y0+d3y0

    return c1*np.cos(t) + c2*np.sin(t) + c3 + c4*t    

#############################################################
time_span = np.arange(0,10,0.1)
S0        = np.array([1,2,9,0.4])
soln      = integrate.odeint(RHS_h,S0,time_span)

plt.plot(time_span,soln[:,0],'k-',label='numerical')
plt.plot(time_span,analytic_soln(time_span,S0),'r.',label='analytic')
plt.legend()
plt.xlabel('$t$',fontsize=16)
plt.ylabel('$y(t)$',fontsize=16)

plt.show()

print(soln[:,0]/analytic_soln(time_span,S0))