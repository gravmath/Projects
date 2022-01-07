import numpy as np

import matplotlib.pyplot as plt
import scipy.integrate   as integrate

def RHS_h(state,time):
    #VOP example from Problem 10-2 in Battin (page 476)
    #homogenous equations
    y, ydot, yddot = state

    return np.array([ydot,yddot,-4.0*ydot])


def RHS_inh(state,time):
    #VOP example from Problem 10-2 in Battin (page 476)
    #inhomogeneous equation
    y, ydot, yddot = state

    return np.array([ydot,yddot,4.0/np.tan(2*time)-4.0*ydot])

def analytic_soln(t,S0):
    y0,dy0,ddy0 = S0
    c1 = -dy0/2.0
    c2 = -ddy0/4.0
    c3 = (y0 + ddy0/4)
    print(c1,c2,c3)

    return c1*np.cos(2*t) + c2*np.sin(2*t) + c3    

#############################################################
time_span = np.arange(np.pi/4.0,10,0.1)
S0        = np.array([-1.0,2,40.0])
soln      = integrate.odeint(RHS_h,S0,time_span)

plt.plot(time_span,soln[:,0],'k-',label='numerical')
plt.plot(time_span,analytic_soln(time_span,S0),'ro',label='analytic')
plt.legend()
plt.xlabel('$t$',fontsize=16)
plt.ylabel('$y(t)$',fontsize=16)

plt.show()