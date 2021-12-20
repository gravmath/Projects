import numpy as np
class Universe:
    #########################################
    # SI Units
    #########################################    
    #kilometers
    km = 1000

    #kilometers per second

    #########################################
    # Physical Constants
    #########################################    
    #speed of light m/s
    c = 299792458
    
    #Newton's constant m^3/kg s^2
    G = 6.67430e-11
    
    #electronic charge
    e_unit = 1.609e-19
    
    #vacuum permittivity epsilon_0 (F/m)
    ep_0 = 8.8541878128e-12
    
    #magnetic constant mu_0 (H/m)
    mu_0 = 1.25663706212e-6
    
    #Boltzman constant k_B (ev/K)
    k_B = 8.617333262145e-5

    #########################################
    # Celestial Object Constants
    #########################################    
    #Earth mu
    mu_earth = 398600.4414
    
    #Sun mu
    mu_sun = 1.32712440018e20
    
    #degrees
    deg = np.pi/180.0

    #radians
    rad = 180.0/np.pi

    #null_vec
    null_vec = np.array([np.nan,np.nan,np.nan])

#########################################
# functions for frame calculations
######################################### 
def calc_spherical_elements(state,space):
    #extract the position
    x,y,z = state[0:3]

    #calculate the distance from the origin
    r     = np.sqrt(x*x + y*y + z*z)

    #calculate the latitude as \lambda = acos(z/rho)
    lamb  = np.arcsin(z/r)*space.rad

    #calculate the longitude as \phi = atan2(y,x)
    phi   = np.arctan2(y,x)*space.rad        

    sph_pos_rep = np.array([r,lamb,phi])
    
    return sph_pos_rep

def calc_local_spherical_frame(state,space):
    
    #allocate a container for the frame
    sph_frame = {}
    
    #convert position into spherical elements: 
    #   r - radius; lamb or l for latitude; phi of p for azimuth
    r, lamb, phi = calc_spherical_elements(state)

    cos_lamb = np.cos(lamb*space.deg)
    sin_lamb = np.sin(lamb*space.deg)
    cos_phi  = np.cos(phi*space.deg)
    sin_phi  = np.sin(phi*space.deg)

    #construct e_r
    sph_frame['e_r'] = np.array([cos_lamb*cos_phi,cos_lamb*sin_phi,sin_lamb])

    #construct e_l
    sph_frame['e_l'] = np.array([-sin_lamb*cos_phi,-sin_lamb*sin_phi,cos_lamb])

    #construct e_p
    sph_frame['e_p'] = np.array([-sin_phi,cos_phi,0.0])
    
    return sph_frame

def calc_local_radial_frame(state,space):
    
    #allocate a container for the frame
    rad_frame = {}

    #extract the position and velocity
    pos   = state[0:3]
    vel   = state[3:6]

    #construct the local frame
    e_R   = pos/np.sqrt(pos.dot(pos))
    h     = np.cross(pos,vel)
    if not (h==space.null_vec).all():  #protect against rectilinear motion
        e_C   = h/np.sqrt(h.dot(h))
        e_I   = np.cross(e_C,e_R)

        rad_frame['e_R'] = e_R
        rad_frame['e_I'] = e_I
        rad_frame['e_C'] = e_C
    else:
        rad_frame['e_R'] = space.null_vec
        rad_frame['e_I'] = space.null_vec
        rad_frame['e_C'] = space.null_vec

    return rad_frame

#########################################
# force model functions
######################################### 
def pointmass_grav_acceleration(state1,state2,mu2):

    #calculate the relative state of object 1 relative to object 2
    rel_state = state1 - state2
    
    #calculates the acceleration due to the force on object1 due to object2 
    pos = rel_state[0:3]
    r   = np.sqrt(pos.dot(pos))
    a   = -mu2*pos/r**3
    return a



#solar functions
def solar_wind_speed(state,space,Sun):
    #assumes that the state has its origin at the sun
    
    theta = calc_spherical_elements(state,space)
    if theta < Sun.current_sheet_tilt:
        V = 400.0*space.kps
    else:
        V = 800.0*space.kps
        
    return V

def b_unit_vector(state,space,Sun):
    
    #find the local frame of the object
    r, lamb, phi = calc_spherical_elements(state)
    e_r, e_l, e_p = calc_local_spherical_frame(state)
    
    #construct the parameter a
    V = solar_wind_speed(Sun)
    a = Sun.Omega

def TwoBody(state,t,particle,central_body):
    import pdb; pdb.set_trace()
    
    particle.set_state(state)
    f = pointmass_grav_acceleration(central_body,particle)
    return np.hstack((particle.get_velocity(),f))

def RHS_dust(state,Sun):
    
    #unpack the state
    vel = state[3:]
    
    #calculate the acceleration due to the gravitational force
    a_grav = pointmass_grav_acceleration(state,Sun)
    
    #calculate the Lorentz force
    #a_L    = Lorentz_force(state,Sun)
    
    #calculate the Poynting-Robertson force
    #a_PR   = Poynting_Robertson_force(state,Sun)
    
    #acc = 
