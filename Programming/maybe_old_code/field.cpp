//******************************************************************************
//Name:  field.cpp                                                             *
//                                                                             *
//Purpose:  source file for defining the field class                           *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
//                    :  10/24/98 - changed SetTensor to use int indices       *
//                    :             rather than a double position array        *
//                    :  11/12/98 - promoted temporary variables in            *
//                    :             GetInterpTensor to protected member data   *
//                    :             in interest of performance                 *
//                    :  12/10/98 - implemented the user-defined grid smoothing*
//                                  radii
// *****************************************************************************

#include "field.hpp"

//******************************************************************************
//Name:  field constructor                                                     *
//                                                                             *
//Purpose:  construct a default field                                          *
//                                                                             *
//Takes: no args and is the default constructor                                *
//******************************************************************************
field::field(void)
{
  int i;
  Resize(3,4,4,4);

  current_time = 0.0;
  is_time_set = FALSE;

  for(i = 0; i < 6; i++) spatial_extent[i] = -9.99e99;
  is_pos_set  = FALSE;

  Delta_mag = 0.0;
  y = 0.0;
  for(i = 0; i < 3; i++)
  {
	  Deltas[i]        = 0.0;
	  current_pos[i]   = 0.0;
	  reference_pos[i] = 0.0;
     pos_diff[i]      = 0.0;
  }

  particle_smoothing_method = 0;
  grid_smoothing_method = 0;
  out_of_range = FALSE;
  is_periodic = FALSE;
  grid_radius = 0;
  grid_upper_bound = 0;
  grid_lower_bound = 0;
  smoothing_length = 1.0;
  normalization = 0.0;
  series.order = -1;
  is_series_set = FALSE;
  is_deriv_field = FALSE;

}

//******************************************************************************
//Name:  field constructor                                                     *
//                                                                             *
//Purpose:  construct a specified field                                        *
//                                                                             *
//Takes: num indices and a variable argument list                              *
//******************************************************************************
field::field(int num_indices, int range1, ...)
{
    //Construct the range array for Resize
	int *range, i;

	range = new int[num_indices];

	//Initialize variable arg pointer
	va_list arg_pnt;

	va_start(arg_pnt, range1);

	//Pack the range array
	range[0] = range1;
	for(i = 1; i < num_indices; i++)
		range[i] = va_arg(arg_pnt, int);

	//Clean up arg pointer
	va_end(arg_pnt);

	Resize(num_indices, range);

	//Initialize member data
   current_time = 0.0;
   is_time_set = FALSE;

   for(i = 0; i < 6; i++) spatial_extent[i] = -9.99e99;
   is_pos_set  = FALSE;

   Delta_mag = 0.0;
   y = 0.0;
   for(i = 0; i < 3; i++)
   {
	   Deltas[i]        = 0.0;
	   current_pos[i]   = 0.0;
	   reference_pos[i] = 0.0;
      pos_diff[i]      = 0.0;
   }

   particle_smoothing_method = 0;
   grid_smoothing_method = 0;
   out_of_range = FALSE;
   is_periodic = FALSE;
   grid_radius = 0;
   grid_upper_bound = 0;
   grid_lower_bound = 0;
   smoothing_length = 1.0;
   normalization = 0.0;
   series.order = -1;
   is_series_set = FALSE;
   is_deriv_field = FALSE;

   //Clean up
   delete [] range;

}

//******************************************************************************
//Name:  SetSpatialExtent                                                      *
//                                                                             *
//Purpose:  set the spatial extent of the field                                *
//                                                                             *
//Takes: an int arg for the direction id and a double for the size             *
//******************************************************************************
void field::SetSpatialExtent(int id, double value)
{
	spatial_extent[id] = value;

	is_pos_set = TRUE;

	for(int i = 0; i < 6; i++)
	{
		if(spatial_extent[i] == -9.99e99) is_pos_set = FALSE;
	}


	if( is_pos_set == TRUE )
	{
		Deltas[0] = (spatial_extent[X_MAX] - spatial_extent[X_MIN]) /
			         (double)(p->range[0] - 1);

    	Deltas[1] = (spatial_extent[Y_MAX] - spatial_extent[Y_MIN]) /
	               (double)(p->range[1] - 1);

    	Deltas[2] = (spatial_extent[Z_MAX] - spatial_extent[Z_MIN]) /
	               (double)(p->range[2] - 1);

      Delta_mag = sqrt( Deltas[0]*Deltas[0] +
                        Deltas[1]*Deltas[1] +
                        Deltas[2]*Deltas[2]);

      smoothing_length = Deltas[0] * grid_radius;
      if( grid_radius == 0 ) smoothing_length = Delta_mag;
	}
}

//******************************************************************************
//Name:  SetSpatialExtent                                                      *
//                                                                             *
//Purpose:  set the spatial extent of the field                                *
//                                                                             *
//Takes: a pointer to a double array                                           *
//******************************************************************************
void field::SetSpatialExtent(double *value)
{

	for(int i = 0; i < 6; i++) spatial_extent[i] = value[i];

   is_pos_set = TRUE;

	Deltas[0] = (spatial_extent[X_MAX] - spatial_extent[X_MIN]) /
  		         (double)(p->range[0] - 1);

  	Deltas[1] = (spatial_extent[Y_MAX] - spatial_extent[Y_MIN]) /
   				(double)(p->range[1] - 1);

  	Deltas[2] = (spatial_extent[Z_MAX] - spatial_extent[Z_MIN]) /
   				(double)(p->range[2] - 1);

   Delta_mag = sqrt( Deltas[0]*Deltas[0] +
                     Deltas[1]*Deltas[1] +
                     Deltas[2]*Deltas[2]);

  smoothing_length = Deltas[0] * grid_radius;
  if( grid_radius == 0 ) smoothing_length = Delta_mag;
}

//******************************************************************************
//Name:  ResizeField                                                           *
//                                                                             *
//Purpose:  set the number of grid points in the field                         *
//                                                                             *
//Takes: a pointer to an int array                                             *
//******************************************************************************
void field::ResizeField(int *num_grid_pts)
{
   int n, *r, i;

   n = NumIndices();
   r = Ranges();

   for(i = 0; i < 3; i++) r[i] = num_grid_pts[i];

   Resize(n, r);

   delete [] r;
}
//******************************************************************************
//Name:  GetIndices                                                            *
//                                                                             *
//Purpose:  get the indices for the closest grid point                         *
//                                                                             *
//Takes: pointer to a position array and a pointer to an index array           *
//******************************************************************************
void field::GetIndices(double *position, int *indices)
{

//Error check
	if( !is_pos_set )
		error("GetIndices error:","spatial extent not set");

	if( position[0] > spatial_extent[X_MAX] + SPATIAL_TOL ||
		 position[0] < spatial_extent[X_MIN] - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

	if( position[1] > spatial_extent[Y_MAX] + SPATIAL_TOL ||
		 position[1] < spatial_extent[Y_MIN] - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

	if( position[2] > spatial_extent[Z_MAX] + SPATIAL_TOL ||
		 position[2] < spatial_extent[Z_MIN] - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

//set the index array
	indices[0] = (int) ( (position[0] - spatial_extent[X_MIN])           /
		                  (spatial_extent[X_MAX] - spatial_extent[X_MIN]) *
						      (p->range[0] - 1)  );

	indices[1] = (int) ( (position[1] - spatial_extent[Y_MIN])           /
		                  (spatial_extent[Y_MAX] - spatial_extent[Y_MIN]) *
 						      (p->range[1] - 1)  );

	indices[2] = (int) ( (position[2] - spatial_extent[Z_MIN])           /
		                  (spatial_extent[Z_MAX] - spatial_extent[Z_MIN]) *
						      (p->range[2] - 1)  );
}

//******************************************************************************
//Name:  GetIndices                                                            *
//                                                                             *
//Purpose:  get the indices for the closest grid point                         *
//                                                                             *
//Takes: pointer to a position array and a pointer to a double index array     *
//******************************************************************************
void field::GetIndices(double *position, double *indices)
{

//Error check
	if( !is_pos_set )
		error("GetIndices error:","spatial extent not set");

	if( position[0] > spatial_extent[X_MAX] + SPATIAL_TOL ||
		 position[0] < spatial_extent[X_MIN] - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

	if( position[1] > spatial_extent[Y_MAX] + SPATIAL_TOL ||
		position[1] < spatial_extent[Y_MIN]  - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

	if( position[2] > spatial_extent[Z_MAX] + SPATIAL_TOL ||
		position[2] < spatial_extent[Z_MIN]  - SPATIAL_TOL )
	  error("GetIndices error:","position out of field bounds");

//set the index array
	indices[0] = (position[0] - spatial_extent[X_MIN])            /
		          (spatial_extent[X_MAX] - spatial_extent[X_MIN])  *
			       (p->range[0] - 1);

	indices[1] = (position[1] - spatial_extent[Y_MIN])            /
		          (spatial_extent[Y_MAX] - spatial_extent[Y_MIN])  *
                (p->range[1] - 1);

	indices[2] = (position[2] - spatial_extent[Z_MIN])            /
		          (spatial_extent[Z_MAX] - spatial_extent[Z_MIN])  *
				    (p->range[2] - 1);
}

//******************************************************************************
//Name:  GetPosition                                                           *
//                                                                             *
//Purpose:  get the tensor containing the position of a grid based on its      *
//          indices                                                            *
//                                                                             *
//Takes: pointer to an index array                                             *
//******************************************************************************
tensor field::GetPosition(int *indices)
{
   tensor pos;

	int i;

//Error check
	if( !is_pos_set )
		error("GetPosition error:","spatial extent not set");

	for(i = 0; i < 3; i++)
	{
	  if( indices[i] < 0 || indices[i] > p->range[i] )
		error("GetPosition error:","requested index off of grid");
	}

//fill the pos tensor
	for(i =	0; i < 3; i++)
	{
		pos.Set((double)indices[i] * Deltas[i] + spatial_extent[2*i], i);
	}

   return pos;
}

//******************************************************************************
//Name:  GetPosition                                                           *
//                                                                             *
//Purpose:  get the position of a grid point based on it indices               *
//                                                                             *
//Takes: pointer to an index array and a pointer to a position array           *
//******************************************************************************
void field::GetPosition(int *indices, double *position)
{

	int i;

//Error check
	if( !is_pos_set )
		error("GetPosition error:","spatial extent not set");

	for(i = 0; i < 3; i++)
	{
	  if( indices[i] < 0 || indices[i] > p->range[i] )
		error("GetPosition error:","requested index off of grid");
	}

//fill the position array
	for(i =	0; i < 3; i++)
	{
		position[i] = (double)indices[i] * Deltas[i] + spatial_extent[2*i];
	}
}

//******************************************************************************
//Name:  SetSmoother                                                           *
//                                                                             *
//Purpose:  set the smoothing method for the particle or the grid              *
//                                                                             *
//Takes: int id and an int for the smoothing method                            *
//******************************************************************************
void field::SetSmoother(int id, int smoothing_method)
{
   if(id == PARTICLE_SMOOTHING)
	   particle_smoothing_method = smoothing_method;

   if(id == GRID_SMOOTHING)
	   grid_smoothing_method = smoothing_method;
}

//******************************************************************************
//Name:  SetGridRadius                                                         *
//                                                                             *
//Purpose:  set the smoothing radius for grid->particle                        *
//                                                                             *
//Takes: int that is the number of grid spacings                               *
//******************************************************************************
void field::SetGridRadius(int grid_rad)
{
   grid_radius = grid_rad;
   smoothing_length = grid_radius;
   if( grid_radius == 0 )
    {
      grid_lower_bound = 0;
      grid_upper_bound = 2;
    }
   else
    {
      grid_lower_bound = -grid_radius - 1;
      grid_upper_bound = grid_radius + 2;
    }
}


//******************************************************************************
//Name:  SetPeriodic                                                           *
//                                                                             *
//Purpose:  set the flag which determines whether periodic boundary conds.     *
//          are employed                                                       *
//                                                                             *
//Takes: int choice with values either TRUE or FALSE                           *
//******************************************************************************
void field::SetPeriodic(int choice)
{
   is_periodic = choice;
}

//******************************************************************************
//Name:  GetTensor                                                             *
//                                                                             *
//Purpose:  get the tensor portion of the field given the position             *
//                                                                             *
//Takes: position stored as a tensor                                           *
//******************************************************************************
tensor field::GetTensor(const tensor &pos)
{

   int i;
   double position[3];

   tensor temp;

   for( i = 0; i < 3; i++ ) position[i] = pos.Val(i);

   temp <= GetTensor(position);

   return temp;
}

//******************************************************************************
//Name:  GetTensor                                                             *
//                                                                             *
//Purpose:  get the tensor portion of the field given the position             *
//                                                                             *
//Takes: point to a position array                                             *
//******************************************************************************
tensor field::GetTensor(double *position)
{

	int *index;
	int num_indices, f_num_indices;
	int *range;
	int *d_index;
	int i, j, temp;

//set up the index array that indexes into the field and pack the position indices
	index = new int[p->num_indices];
	GetIndices(position, index);

//get the parameters necessary to define the dummy return tensor
	num_indices = p->num_indices - 3;
	range = new int[num_indices];

	for(i = 0; i < num_indices; i++)
		range[i] = p->range[i+3];

//set up the d_index array that indexes into the dummy return tensor
	d_index = new int[num_indices];

//Create the dummy tensor which will hold the return values and resize it
	tensor dummy;

	dummy.Resize(num_indices, range);

//Loop over the number of components
	for(i = 0; i < dummy.p->product; i++)
	{
      //test for range = 1 cases here
      f_num_indices = p->num_indices;
      for(j = 0; j < p->num_indices; j++)
       {
         if(p->range[j] == 1) f_num_indices--;
       }
		//pack the index array (see tensor::Contract)
		temp = 0;
		for(j = 3; j < f_num_indices - 1; j++)
		 {
		  index[j] = (i - i%p->scales[j] - temp)/p->scales[j];
		  temp += index[j] * p->scales[j];
		 }
      //handle range = 1 cases here
      if(f_num_indices < p->num_indices)
       {
        index[f_num_indices - 1] = (i - temp)/p->scales[f_num_indices - 1];
        for(j = f_num_indices; j < p->num_indices; j++) index[j] = 0;
       }
      else
       {
        index[j] = (i - temp)/p->scales[j];
       }

		//strip off only the 'field' indices leaving the indices that specify position
		for(j = 0; j < num_indices; j++)
			d_index[j] = index[j+3];

		//set the value of dummy with the value of the field
        dummy.Set(Val(index),d_index);
	}

//clean up
	delete [] range;
	delete [] index;
	delete [] d_index;

//return
   return dummy;

}

//******************************************************************************
//Name:  GetTensor                                                             *
//                                                                             *
//Purpose:  get the tensor portion of the field given the position             *
//                                                                             *
//Takes: position stored as a tensor                                           *
//******************************************************************************
tensor  field::GetTensor(int *indices)
{
   tensor temp;
   int n, *r, new_n, *new_r;
   int i, j, dummy;

   //Get the total number of indices and ranges of the field
   n = NumIndices();
   r = Ranges();

   //Determine the number of indices on the temp tensor
   new_n = n  - 3;

   //Make a new array for resizing temp
   new_r = new int[new_n];

   //Pack the new_r array with the
   for( i = 0; i < new_n; i++ )
    {
      new_r[i] = r[i + 3];
    }

   //Resize the temp tensor
   temp.Resize(new_n, new_r);

   //Let the r array do double duty and pack the first three elements with the
   //values in the indices array
   for( i = 0; i < 3; i++) r[i] = indices[i];

   //now pack temp with the components of the tensor portion of the field
   for( i = 0; i < temp.p->product; i++)
	 {
      //Let the last indices of the r array hold the
      //index structure for the temp tensor
      dummy = 0;
      for( j = 0; j < new_n; j++)
       {
         r[j + 3] = (i - i%temp.p->scales[j] - dummy)/temp.p->scales[j];
         dummy += r[j + 3] * temp.p->scales[j];
       }

       temp.p->m[i] = Val(r);
    }

  //Clean up
  delete [] r;
  delete [] new_r;

  //Return the tensor
  return temp;
}
//******************************************************************************
//Name:  SetTensor                                                             *
//                                                                             *
//Purpose:  Set the tensor portion of the field                                *
//                                                                             *
//Takes: takes a pointer to a position array and the address of the tensor     *
//******************************************************************************
void field::SetTensor(int *indices, const tensor &rval)
{
   int *index;
	int *r_index;
	int i, j, temp;
	double value;

//Error check
	if( rval.p->num_indices != (p->num_indices - 3) )
		error("SetTensor error:","cannot assign the tensor to a field point");

	for(i = 3; i < p->num_indices; i++)
	{
	  if( rval.p->range[i-3] != p->range[i] )
		  error("SetTensor error:","cannot assign the tensor to a field point");
	}

//set up the index array that indexes into the field and set the specified position
	index = new int[p->num_indices];
   for( i = 0; i < 3; i++) index[i] = indices[i];

//set up the index array that indexes into the rval tensor
	r_index = new int[rval.p->num_indices];

//set the field
	for(i = 0; i < rval.p->product; i++)
	{
		//pack the r_index array (see tensor::Contract)
		temp = 0;
		for(j = 0; j < rval.p->num_indices - 1; j++)
		{
		  r_index[j] = (i - i%rval.p->scales[j] - temp)/rval.p->scales[j];
		  temp += r_index[j] * rval.p->scales[j];
		}

		r_index[j] = (i - temp)/rval.p->scales[j];

		//complete the index array with the tensor indices
		for(j = 3; j < p->num_indices; j++) index[j] = r_index[j - 3];

      value = rval.Val(r_index);

      //put the rval tensors value into the field at specified position
      Set(value,index);
	}

  //Clean up
  delete [] index;
  delete [] r_index;

}

//******************************************************************************
//Name:  SetTensor                                                             *
//                                                                             *
//Purpose:  Set the tensor portion of the field                                *
//                                                                             *
//Takes: takes a pointer to a position array and the address of the tensor     *
//******************************************************************************
void field::SetTensor(double *position, const tensor &rval)
{

	int *index;
	int *r_index;
	int i, j, temp;
	double value;

//Error check
	if( rval.p->num_indices != (p->num_indices - 3) )
		error("SetTensor error:","cannot assign a tensor to a field point");

	for(i = 3; i < p->num_indices; i++)
	{
	  if( rval.p->range[i-3] != p->range[i] )
		  error("SetTensor error:","cannot assign a tensor to a field point");
	}
//set up the index array that indexes into the field and set the specified position
	index = new int[p->num_indices];
	GetIndices(position, index);

//set up the index array that indexes into the rval tensor
	r_index = new int[rval.p->num_indices];

//create a tempory tensor for stupid Microsoft (comment out if not using MS)
	tensor temp_t;
	temp_t = rval;

//set the field
	for(i = 0; i < rval.p->product; i++)
	{
		//pack the r_index array (see tensor::Contract)
		temp = 0;
		for(j = 0; j < rval.p->num_indices - 1; j++)
		{
		  r_index[j] = (i - i%rval.p->scales[j] - temp)/rval.p->scales[j];
		  temp += r_index[j] * rval.p->scales[j];
		}

		r_index[j] = (i - temp)/rval.p->scales[j];

		//complete the index array with the tensor indices
		for(j = 3; j < p->num_indices; j++) index[j] = r_index[j - 3];

        //comment out following line when using Microsoft (in the head)
        //value = rval.Val(r_index);

      //comment out following line when not using Microsoft (in the head)
		value = temp_t.Val(r_index);

		//put the rval tensors value into the field at specified position
      Set(value,index);
	}

  //Clean up
  delete [] index;
  delete [] r_index;
}

//******************************************************************************
//Name:  SetTensor                                                             *
//                                                                             *
//Purpose:  Set the tensor portion of the field                                *
//                                                                             *
//Takes: takes a position in tensor form and the address of the tensor         *
//******************************************************************************
void field::SetTensor(const tensor &pos, const tensor &rval)
{

   int i;
   double position[3];

   for( i = 0; i < 3; i++) position[i] = pos.Val(i);

   SetTensor(position, rval);

}

//******************************************************************************
//Name:  SetCurrentTime                                                        *
//                                                                             *
//Purpose:  Set the current time associated with the field                     *
//                                                                             *
//Takes: double value of the time                                              *
//******************************************************************************
void field::SetCurrentTime(double time)
{
	current_time = time;
	is_time_set = TRUE;
}

//******************************************************************************
//Name:  GetCurrentTime                                                        *
//                                                                             *
//Purpose:  Get the current time associated with the field                     *
//                                                                             *
//Takes: void                                                                  *
//******************************************************************************
double field::GetCurrentTime(void)
{
	if(!is_time_set)
		error("GetCurrentTime error:","field time not set");

	return current_time;
}

//******************************************************************************
//Name:  Dist                                                                  *
//                                                                             *
//Purpose:  get the Euclidean distance associated with a position              *
//                                                                             *
//Takes: pointer to the position array                                         *
//******************************************************************************
double field::Dist(double *position)
{
	if(!is_pos_set)
		error("Dist error:","spatial extent not set");

	double dist_ret;

    dist_ret = sqrt( position[0]*position[0] +
					      position[1]*position[1] +
					      position[2]*position[2] );

    return dist_ret;
}

//******************************************************************************
//Name:  SetSeriesOrder                                                        *
//                                                                             *
//Purpose:  set the order of the taylors series expansion                      *
//                                                                             *
//Takes: int whose value is the order of the series                            *
//******************************************************************************
void field::SetSeriesOrder(int order)
{
    if( order > 2 )
      error("SetSeriesOrder:","order must be less than 3");

    series.order = order;

    //allocate space for the series list (note that order+1 is needed to
    //reflect that the series begins with the zeroth term)
    series.list = new field* [order+1];
    for(int i = 0; i <= order; i++) series.list[i] = NULL;
}

//******************************************************************************
//Name:  SetSeriesList                                                         *
//                                                                             *
//Purpose:  set the elements of the taylors series expansion                   *
//                                                                             *
//Takes: an array of pointers to fields                                        *
//******************************************************************************
void field::SetSeriesList(field *list[])
{
 	int i;

   if( series.order == -1 )
    error("SetSeriesList error:","series order not set");

   for( i = 0; i <= series.order; i++) series.list[i] = list[i];

   is_series_set = TRUE;

   tmp = new tensor[series.order + 1];

   //Now resize the temporary tensors for the interpolation here

   sum_Wm.Resize(1,3);

   if( is_dep_field == TRUE )
    {
      int n, *s, *t;
      n = list[0]->NumIndices() - 3;
      s = list[0]->Ranges();
      t = new int[n];
      for( i = 0; i < n; i++) t[i] = s[i+3];
      sum_a2.Resize(n,t);
      delete [] s;
      delete [] t;

      n = list[1]->NumIndices() - 3;
      s = list[1]->Ranges();
      t = new int[n];
      for( i = 0; i < n; i++) t[i] = s[i+3];
      sum_da1.Resize(n,t);
      d_sum_a2.Resize(n,t);
      delete [] s;
      delete [] t;
    }

}

//******************************************************************************
//Name:  SetDerivField                                                         *
//                                                                             *
//Purpose:  set the derivative field flag                                      *
//                                                                             *
//Takes: an int with the choice                                                *
//******************************************************************************
void field::SetDerivField(int type)
{
 	is_deriv_field = type;

}

//******************************************************************************
//Name:  SetFieldDep                                                           *
//                                                                             *
//Purpose:  set the field dependency flag                                      *
//                                                                             *
//Takes: an int with the choice                                                *
//******************************************************************************
void field::SetFieldDep(int type)
{
 	is_dep_field = type;

}

//******************************************************************************
//Name:  GetInterpTensor                                                       *
//                                                                             *
//Purpose:  get the interpolated tensor at an arbitrary position               *
//                                                                             *
//Takes: pointer to the position as packed in a tensor                         *
//******************************************************************************
tensor field::GetInterpTensor(tensor &pos)
{
   double position[3];
   tensor ret_val;

   for(int i = 0; i < 3; i++)  position[i] = pos.Val(i);

   ret_val <= GetInterpTensor(position);

   return ret_val;
}

//******************************************************************************
//Name:  GetInterpTensor                                                       *
//                                                                             *
//Purpose:  get the interpolated tensor at an arbitrary position               *
//                                                                             *
//Takes: pointer to the position array                                         *
//******************************************************************************
tensor field::GetInterpTensor(double *position)
{
 int i, j, k, m, n;
 double scaled_pos[3];
 int grid_pos[3];
 double W;

 if( is_series_set == FALSE )
  error("GetInterpTensor error:","series not set");

 //initialize member data tensors
 sum_Wm   <= 0.0*sum_Wm;
 sum_a2   <= 0.0*sum_a2;
 sum_da1  <= 0.0*sum_da1;
 d_sum_a2 <= 0.0*d_sum_a2;

 out_of_range = FALSE;

//Set the reference position
  for(i = 0; i < 3; i++) reference_pos[i] = position[i];

//Get the real indices for the scaled position
  GetIndices(reference_pos, scaled_pos);

//Create the return tensor and set it to zero
 ret_dummy <= 0.0*GetTensor(reference_pos);

//Reset the normalization
  normalization = 0.0;

//step through the floor steping
  for( i = grid_lower_bound; i < grid_upper_bound; i++)
   {
	 for( j = grid_lower_bound; j < grid_upper_bound; j++)
	  {
	   for( k = grid_lower_bound; k < grid_upper_bound; k++)
		 {
         //Get the indices of the lowest corner of the cell
			grid_pos[0] = (int)floor( scaled_pos[0] + i );
			grid_pos[1] = (int)floor( scaled_pos[1] + j );
         grid_pos[2] = (int)floor( scaled_pos[2] + k );

         //find the position this grid point
         GetPosition(grid_pos,current_pos);

         //form the difference array and its magnitude (for use by smoother)
         //and the difference position
         y = 0;
         for( m = 0; m < 3; m++ )
			 {
           pos_diff[m] = reference_pos[m] - current_pos[m];
           y += pos_diff[m]*pos_diff[m];
           d.Set(pos_diff[m], m);
			 }
         y = sqrt(y)/smoothing_length;
         if( y > 1.0 + SMOOTHING_TOL ) continue;

         //get all of the tensors requested by the order
         for( m = 0; m <= series.order; m++)
           tmp[m] <= series.list[m]->GetTensor(grid_pos);

         //get the (psuedo) scalar smoothing value
         W = Smoother(CWM_SMOOTHER);

			//if the field is basic (i.e. not a derivative of
			//another field)
			if( is_deriv_field == FALSE )
			 {
           if( series.order == 2 )
            {
              m = tmp[2].NumIndices() - 1;
              tmp[2] <= 0.5 *( tmp[2].Contract(d,m,0) ).Contract(d,m-1,0);
            }
           if( series.order == 1 || series.order == 2 )
            {
              m = tmp[1].NumIndices() - 1;
              tmp[1] <= tmp[1].Contract(d,m,0);
            }
           for( m = 0; m <= series.order; m++)
			      ret_dummy <= ret_dummy + W*tmp[m];
			 }


			//if the field is not basic (i.e. is a derivative of
			//another field
         if( is_deriv_field == TRUE )
			 {
			   //if the field is independent (i.e. its smoothing is not
			   //a derivative of the smoohting of its basic field
			   if( is_dep_field == FALSE )
			    {
				  temp0 <= 0.0*tmp[0];
				  if( series.order == 1 )
				   {
				    m = tmp[1].NumIndices() - 1;
				    temp0 <= tmp[1].Contract(d,m,0);
				   }
				  ret_dummy <= ret_dummy + W*(tmp[0] + temp0);
			    }
			   //if the field is dependent (i.e. it smoothing is
			   //based on the derivative of the smoothing of its basic field
			   else
			    {
				  temp0 <= 0.0*tmp[0];
              if( series.order == 2 )
				   {
                 m = tmp[2].NumIndices() - 1;
                 //form the 1/2 * t,ij*(z-x)_i*(z-x)_j term if neccessary
                 temp0 <= 0.5*( tmp[2].Contract(d,m,0) ).Contract(d,m-1,0);
				   }

				  temp1 <= 0.0*tmp[0];
				  if( series.order ==1 || series.order == 2 )
				   {
                 n = tmp[1].NumIndices() - 1;
				     //form the t,i*(z-x)_i term if neccessary
					  temp1 <= tmp[1].Contract(d,n,0);

					  //form the t,i + t,ij*(z-x)_j
                 temp2 <= 0.0*tmp[1];
                 if( series.order == 2 ) temp2 <=  tmp[2].Contract(d,m,0);
                 sum_da1 <= sum_da1 + W*(tmp[1] + temp2);
 				   }
				  //form the t + t,i*(z-x)_i + 1/2 * t,ij*(z-x)_i*(z-x)_j
              a2 <= tmp[0] + temp1 + temp0;

              sum_a2 <=  sum_a2 + W*a2;

              //get the derivative of the smoothing
              Wm <= Smoother_Deriv(CWM_SMOOTHER);

              //keep running sum
              sum_Wm <= sum_Wm + Wm;

              //form the term labeled dN_D_2 in the MathCAD
              d_sum_a2 <= d_sum_a2 + a2 * Wm;
			    }
          }
	    }
	  }
   }

//Normalize the interpolated tensor
  if( is_dep_field == FALSE ) ret_dummy <= (1.0/normalization)*ret_dummy;
  if( is_dep_field == TRUE  )
  {
    ret_dummy <= normalization * ( sum_da1 + d_sum_a2 ) - sum_a2 * sum_Wm;
    ret_dummy <= 1.0/(normalization*normalization) * ret_dummy;
  }

 return ret_dummy;

}

//******************************************************************************
//Name:  SetInterpTensor                                                       *
//                                                                             *
//Purpose:  set the grid tensors based on interpolation from a tensor at an    *
//          arbitrary position                                                 *
//                                                                             *
//Takes: pointer to the position array and address of the rval tensor          *
//******************************************************************************
void field::SetInterpTensor(double *position, const tensor &rval)
{
 int step;
 int i, j, k;
 int scaled_pos[3];
 int grid_pos[3];
 double pos[3];
 tensor dummy, temp;


 step = 1;
 out_of_range = FALSE;

//Set the reference position
  for(i = 0; i < 3; i++) reference_pos[i] = position[i];

//Get the double indices for the scaled position
  GetIndices(reference_pos, scaled_pos);

//step through the floor steping
 while( !out_of_range )
 {
	 for( i = 0; i < 2; i++)
	 {
		 for( j = 0; j < 2; j++)
		 {
			 for( k = 0; k < 2; k++)
			 {
			   grid_pos[0] = (int)floor( scaled_pos[0] + i * step );
			   grid_pos[1] = (int)floor( scaled_pos[1] + j * step );
            grid_pos[2] = (int)floor( scaled_pos[2] + k * step );
            GetPosition(grid_pos, pos);
            dummy <= GetTensor(pos);
            temp <= rval;
			   temp.ScalarMult(Smoother(PARTICLE_SMOOTHING));
			   dummy <= dummy + temp;
			   if( out_of_range ) break;
			 }
		 }
	 }
//temp hack to keep the interpolation to within one cell
 out_of_range = TRUE;
 }

}

//******************************************************************************
//Name:  Smoother                                                              *
//                                                                             *
//Purpose:  general function that mediates the various smoothings              *
//                                                                             *
//Takes: an int telling what type of smoothing                                 *
//       (grid->particle or particle->grid)                                    *
//******************************************************************************
double field::Smoother(int type)
{
   double W;
 
   switch(type)
   {
     case CWM_SMOOTHER:
                if( y > 1.0 + SMOOTHING_TOL )
                 {
                   W = 0.0;
                 }
                else
                 {
                   W = pow((1.0 - y*y), 3);
                 }
                break;

     case DEFAULT_SMOOTHING:
                W = pow((1.0 - y*y), 3);
                break;
   }
   normalization += fabs(W);
	return W;
}

//******************************************************************************
//Name:  Smoother_Deriv                                                        *
//                                                                             *
//Purpose:  general function that models the derivatives of the various        *
//          smoothings in Smoother                                            *
//                                                                             *
//Takes: an int telling what type of smoothing                                 *
//       (grid->particle or particle->grid)                                    *
//******************************************************************************
tensor field::Smoother_Deriv(int type)
{
   double q;
   tensor Wm;
   int i;

  switch(type)
   {
     case CWM_SMOOTHER:
                if( y > 1.0 + SMOOTHING_TOL )
                 {
                   Wm <= 0.0*Wm;
                 }
                else
                 {
                   for(i = 0; i < 3; i++)
                    {
                      q = -6.0 * pos_diff[i]*(1 - y*y)*(1 - y*y)
                          /(Delta_mag*Delta_mag);
                      Wm.Set(q,i);
                    }
                 }
                break;

     case DEFAULT_SMOOTHING:
                break;
   }
	return Wm;
}

//******************************************************************************
//Name:  GetNumIndices                                                         *
//                                                                             *
//Purpose:  return the number of grid points associated with a given           *
//          grid direction                                                     *
//                                                                             *
//Takes: an int with the id type                                               *
//******************************************************************************
int field::GetNumIndices(int id)
{
	switch(id)
     {
       case X_ID: return p->range[0];

       case Y_ID: return p->range[1];

       case Z_ID: return p->range[2];
     }

  return -9999;
}

