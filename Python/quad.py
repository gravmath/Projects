#
# quad.py
#
"""
This module defines a set of functions for performing 
numerical quadrature
"""

from numpy import linspace

#******************************************************************************
# trapezoidal rule
#******************************************************************************
def trap(func, min, max, num):
	fmin = func(min)
	fmax = func(max)
	Itot = fmin + fmax
	h    = (max-min)/1.0/num
	for x in linspace(min+h, max-h, num-1):
		fx = func(x)
		Itot += 2.0 * fx
	Itot = h * Itot / 2.0
	return Itot, num
    
#******************************************************************************
# simspon's rule
#******************************************************************************
def simp(func, min, max, num):
    fmin  = func(min)
    fmax  = func(max)
    Itot  = fmin + fmax
    coeff = [4.0,2.0]

    #trim num as needed
    num    = num - num%2
    #calculate the step size
    h      = (max - min)/1.0/num
    #setup a counter for moving through the coeffcient array
    counter = 0

    for x in linspace(min+h, max-h, num-1):
        fx = func(x)
        Itot += coeff[counter%2]*fx
        counter += 1
    Itot = h * Itot / 3.0
    return Itot, num + 1

#******************************************************************************
# simpson's 3/8 rule
#******************************************************************************
def simp38(func, min, max, num):
    fmin = func(min)
    fmax = func(max)
    Itot = fmin + fmax
    coeff = [3.0, 3.0, 2.0]

    #trim num as needed
    num = (num-1) - (num-1)%3 + 1
    #calculate the step size
    h   = (max-min)/1.0/num
    #setup a counter for moving through the coeffcient array
    counter = 0

    for x in linspace(min+h, max-h, num-1):
        fx       = func(x)
        Itot    += coeff[counter%3] * fx
        counter += 1
    Itot = 3.0 * h * Itot / 8.0
    return Itot, num + 1

#******************************************************************************
# boole's rule
#******************************************************************************     
def boole(func,min,max,num):
	#get the function values at the minimum and maximum range
	fmin    = func(min)
	fmax    = func(max)
	Itot    = 7.0 * ( fmax + fmin )
	#allocate the list of coefficients
	coeff   = [32.0,12.0,32.0,14.0]
	#trim num as needed
	num     = (num-1) - (num-1)%4 + 1
	#detemine the step
	h       = (max - min)/1.0/num
	#setup a counter for moving through the coeffcient array
	counter = 0

	for x in linspace(min+h, max-h, num-1):
		fx       = func(x)
		Itot    += coeff[counter%4] * fx
		counter += 1
	Itot = 2.0 * h * Itot / 45.0
	return Itot, num + 1
    
#******************************************************************************
# generate data for each method
#******************************************************************************      
    def gen_data(func, min, max):
	trap_dat = np.zeros([980,2])
	simp_dat = np.zeros([980,2])
	sp38_dat = np.zeros([980,2])
	bool_dat = np.zeros([980,2])
	for i in range(0,979,1):
		trap_dat[i] = trap(  func, min, max, i+20)
		simp_dat[i] = simp(  func, min, max, i+20)
		sp38_dat[i] = simp38(func, min, max, i+20)
		bool_dat[i] = boole( func, min, max, i+20)
	return trap_dat, simp_dat, sp38_dat, bool_dat