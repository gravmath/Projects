#
# two_body_gravity.py
#
"""
This module defines the vector field for two point masses attracting 
each other through Newtonian gravity.
"""

from   math  import sqrt
from   math  import pow
import numpy as     np
import scipy as     sp

#
#The running module
#
def run():
    r  = 400000e3
    v  = 1000.0
    G  = 6.6738480e-11
    m1 = 5.9742e24
    m2 = 7.36e22
    
    s  = [0.0,0.0,0.0,0.0,0.0,0.0,r,0.0,0.0,0.0,v/sqrt(2),v/sqrt(2)]
    t  = 0.0
    p  = [m1,m2,G]

    rhs = vectorfield(s,t,p)
    return rhs
    
#
# The vector field.
#
def vectorfield(s,t,p):
    """
    Arguments:
        s :  vector of the state variables:
                  s = [x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2]
        t :  time
        p :  vector of the parameters:
                  p = [m1,m2,G]
    """
    x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, vx2, vy2, vz2 = s
    m1, m2, G = p

    #first construct the position of m2 wrt m1
    rho = [x2 - x1, y2 - y1, z2 - z1]
    
    #find the L2 norm of rho
    rho_mag = sqrt( rho[0]*rho[0] + rho[1]*rho[1] + rho[2]*rho[2] )
    
    #calculate the force coeffcient
    coeff = G*m1*m2/pow(rho_mag,3)
    
    # Create the vector field f = (x1',y1',x2',y2'):
    f = [vx1,
         vy1,
         vz1,
         coeff*rho[0],
         coeff*rho[1],
         coeff*rho[2],
         vx2,
         vy2,
         vz2,
         -coeff*rho[0],
         -coeff*rho[1],
         -coeff*rho[2]]
    return f

print run()
