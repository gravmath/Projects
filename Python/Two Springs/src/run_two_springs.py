#run_two_springs.py script\
"""
Use ODEINT to solve the differential equations defined by the vector field
in two_springs.py.
"""

from   scipy.integrate import odeint
from   matplotlib      import pyplot as plt
import two_springs

# Parameter values
# Masses:
m1 = 1.0
m2 = 1.5
# Spring constants
k1 = 8.0
k2 = 40.0
# Natural lengths
L1 = 0.5
L2 = 1.0
# Friction coefficients
b1 = 0.8
b2 = 0.5

# Initial conditions
# x1 and x2 are the initial displacements; y1 and y2 are the initial velocities
x1 = -1.0
y1 = 0.0
x2 = 1.0
y2 = 0.0

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 20.0
numpoints = 250

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t = [stoptime*float(i)/(numpoints-1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
p = [m1,m2,k1,k2,L1,L2,b1,b2]
w0 = [x1,y1,x2,y2]

# Call the ODE solver.
wsol = odeint(two_springs.vectorfield,w0,t,args=(p,),atol=abserr,rtol=relerr)

plt.plot(t,wsol[:,1],label='m1 pos.')
plt.plot(t,wsol[:,2],label='m2 pos.')
plt.legend
plt.xlabel('Time (sec)')
plt.ylabel('Mass Position (m)')
plt.title('Damped Coupled Masses')
plt.legend(loc=0)
plt.show()
print 'all done'