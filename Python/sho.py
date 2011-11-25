#
#SHO.py
#

def SHO(state,time):
	dstate = [0,0]
	dstate[0] = state[1]
	dstate[1] = -4*state[0]
	return dstate