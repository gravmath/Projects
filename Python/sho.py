#
# sho.py
#
"""
This is the vector field definition
for the simple harmonic oscillator
"""

"""
The vector field for the right-hand side
"""

def sho(state, time):
    """
    Arguments:
      state - current state of the system as a list
      time  - ccurrent time as a double
      parms - list of physical parameters, in this case the natural freq.
    """
    omega     = 2
    dstate    = [0,0]
    dstate[0] = state[1]
    dstate[1] = -omega*omega*state[0]
    return dstate
