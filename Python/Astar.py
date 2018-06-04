##### A simple A* package
##### 
##### the basic data structure is a grid_point - a dictionary with 4 keys
#####    1) 'loc'     - a tuple indicating the grid location
#####    2) 'terrain' - a integer with the terrain cost
#####    3) 'cost'    - total cost 'terrain' plus 'dist'
#####    4) 'parent'  - either a tuple like 'loc' indicating the location of the 
#####                   the parent or the string 'self' 

##########################################################################
def grid_point(loc,t,c):
    """function to create the grid_point dictionary
	   
	   Arguments:  loc - tuple with the grid_point's location
	               t   - terrain cost
				   c   - total cost
		
       Returns:    grid_point object

       Use:        my_grid_point = grid_point(loc,t,c)

       Note:       typically used by another funtion	   
	   """
    return {'loc':loc,'terrain':t,'cost':c,'parent':()}

##########################################################################	   
def grid_index(loc,num_y):
    """function to return the index of a grid_point in a 1D list (based 
	   on creation order)
	   
	   Arguments:  loc   - tuple with the grid_point's location
	               num_y - number of 'y' rows 
		
       Returns:    integer index locating the grid_point

       Use:        my_grid_index = grid_index(loc,num_y)

       Note:       primary location within the 1-D list of grid_points   
	   """

    return (loc[0]-1)*num_y + (loc[1]-1)
	
	
##########################################################################	
def grid_dist_traditional(loc,targ):
    """function to return the distance between two locations
	   
	   Arguments:  loc   - tuple with the grid_point's location
	               targ  - tuple with the target's location 
		
       Returns:    distance

       Use:        dist = grid_dist_traditional(loc,targ)

       Note:       traditional video game formulation   
	   """
    k = abs(loc[0]-targ[0])
    l = abs(loc[1]-targ[1])
    
    q = min(k,l)
    w = abs(k-l)
    return q + w

##########################################################################		
def initialize_grid(world,targ):
    """function to return the grid used in A*
	   
	   Arguments:  num_x     - number of 'x' rows
	               num_y     - number of 'y' rows
                   grid_dist - function for computing the distance 
				               between points 				   
		
       Returns:    grid      - a list of grid points

       Use:        my_grid = initialize_grid(num_x,num_y,grid_dist_new)

       Note:       assumes walls, mountains, desert, water, etc. 
	               are defined in world
	   """
    grid = []
    for i in range(1,world['num_x']+1):
        for j in range(1,world['num_y']+1):
            loc = (i,j)
            if loc in world['walls']:
                ter = 1000000
            elif loc in world['mountains']:
                ter = 10
            elif loc in world['desert']:
                ter = 7
            elif loc in world['water']:
                ter = 5
            elif loc in world['forest']:
                ter = 2
            else:
                ter = 1
            dist = world['metric'](loc,targ)
            cost = ter + dist
            grid.append(grid_point(loc,ter,cost))
    grid = sorted(grid,key=lambda k : k['loc'])	
	
    return grid

##########################################################################		
def find_index_of_loc_in_list(loc,lst):
    """function to the index in a list of a grid point given 
	   its loc
	   
	   Arguments:  num_x     - number of 'x' rows
	               num_y     - number of 'y' rows
                   grid_dist - function for computing the distance 
				               between points 				   
		
       Returns:    grid      - a list of grid points

       Use:        my_grid = initialize_grid(num_x,num_y,grid_dist_new)

       Note:       assumes walls, mountains, desert, water, etc. 
	               are defined globally
	   """
    counter = 0
    for item in lst:
        if item['loc'] == loc:
            return counter
        counter += 1
    return False

##########################################################################		
def explore_a_spot(loc,grid,open_list,closed_list,world):
    """primary function - it explores a spot function, opening and 
	   scoring neighbors, closing the current spot.
	   
	   Arguments:  loc         - tuple of the location of the explored spot
	               grid        - list of grid points
                   open_list   - open list
                   closed_list - closed_list 				   
		
       Returns:    nothing

       Use:        explore_a_spot(loc,grid,open_list,closed_list)

       Note:       none
	   """
    #put loc on closed_list and remove loc from the open list
    loc_index = grid_index(loc,world['num_y'])
    closed_list.append(grid[loc_index])
    
    #for each neighbor
    for  i in range(-1,2):
        for j in range(-1,2):
            neighbor_loc   = (loc[0]+i,loc[1]+j)
            if neighbor_loc[0] < 1 or neighbor_loc[1] < 1 or neighbor_loc[0] > world['num_x'] or neighbor_loc[1] > world['num_y'] :
                pass
            else:
                neighbor_index = grid_index(neighbor_loc,world['num_y'])
                neighbor       = grid[neighbor_index]
                if neighbor in open_list or neighbor in closed_list:
                    pass
                else:
                    if grid[neighbor_index]['terrain'] < 1000000:
                        grid[neighbor_index]['parent'] = loc
                        open_list.append(grid[neighbor_index])	

##########################################################################							
def construct_path(closed_list,start,targ):
    """function that constructs the path from start to targ
	   
	   Arguments:  closed_list - closed_list
	               start       - tuple of the starting location
				   targ        - tuple of the end location
		
       Returns:    path

       Use:        path = construct_path(closed_list,start,targ)

       Note:       none
	   """
    path = [targ]
    closed_list_by_loc = sorted(closed_list,key=lambda k : k['loc'])
    curr_index         = find_index_of_loc_in_list(targ,closed_list_by_loc)
    flag = 0 
    while flag != 1:
        prev_loc   = closed_list_by_loc[curr_index]['parent']
        curr_index = find_index_of_loc_in_list(prev_loc,closed_list_by_loc)
        path.append(prev_loc)
        if prev_loc == 'self':
            flag = 1
    path.reverse()
    return path						