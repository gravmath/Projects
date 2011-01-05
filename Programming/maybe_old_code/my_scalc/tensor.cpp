//******************************************************************************
//Name:  tensor.cpp                                                            *
//                                                                             *
//Purpose:  source file for defining the tensor class                          *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
//								 5/18/99  - added vector magnitude field for K.Watt    *
//                       05/23/99 - updated to use Resize0 for speed           *
//                       06/18/99 - renamed Resize0 to Mimic                   *
//******************************************************************************

#include "tensor.hpp"

//******************************************************************************
//Name:  tensor constructor                                                    *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: default constructor with no args                                      *
//******************************************************************************
tensor::tensor(void)
{
}

//******************************************************************************
//Name:  tensor constructor                                                    *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: variable argument list with the number of indices and range for each  *
//******************************************************************************
tensor::tensor(int num_indices, int range1, ...)
{
  int *range, i;

//Allocate space for the range array
  range = new int[num_indices];

//Initialize the variable argument list
  va_list arg_list;

  va_start(arg_list,range1);

//Pack the range array
  range[0] = range1;
  for(i = 1; i < num_indices; i++)
	  range[i] = va_arg(arg_list,int);

  Resize(num_indices, range);

//Clean up
  va_end(arg_list);
  delete [] range;
}

//******************************************************************************
//Name:  tensor non-referencing assignment                                     *
//                                                                             *
//Purpose:  memcpy the memory of one tensor to another                         *
//                                                                             *
//Takes: address of tensor to be assigned                                      *
//******************************************************************************
tensor& tensor::operator<=(const tensor &rval)
{
  int i;
  int test = 0;

//Check to see if the two arrays have the same size and number of indices
//if not resize this
  if( (p->product != rval.p->product) ||
      (p->num_indices != rval.p->num_indices) )
   {
     Mimic(rval);
   }

//Check to see if the two arrays have the same ranges
//if not resize this
  for( i = 0; i < p->num_indices ; i++ )
    if( p->range[i] != rval.p->range[i] ) test = 1;

  if( test == 1 )
   {
     Mimic(rval);
   }

//Now memcpy the rval arrays values into this
  memcpy( p->m, rval.p->m, p->product*sizeof(double) );

  return *this;
}

//******************************************************************************
//Name:  tensor ScalarMult                                                     *
//                                                                             *
//Purpose:  multiply a tensor by a scalar                                      *
//                                                                             *
//Takes: takes the scalar as an argument                                       *
//******************************************************************************
void tensor::ScalarMult( double scalar )
{
 int i;

 for( i = 0; i < p->product; i++)
   p->m[i] *= scalar;

}
//******************************************************************************
//Name:  tensor operator * for RHS mult. by a scalar                           *
//                                                                             *
//Purpose:  multiply a tensor by a scalar                                      *
//                                                                             *
//Takes: takes the scalar and a tensor as an argument                          *
//******************************************************************************
tensor operator*(double scalar,tensor const &B)
{
 tensor new_tensor;

 new_tensor.Multiply(B,scalar);

 return new_tensor;

}

//******************************************************************************
//Name:  tensor operator * for LHS mult. by a scalar                           *
//                                                                             *
//Purpose:  multiply a tensor by a scalar                                      *
//                                                                             *
//Takes: takes the scalar and a tensor as an argument                          *
//******************************************************************************
tensor operator*(tensor const &B, double scalar)
{
 tensor new_tensor;

 new_tensor.Multiply(B,scalar);

 return new_tensor;

}

//******************************************************************************
//Name:  tensor operator * for LHS mult. by a scalar                           *
//                                                                             *
//Purpose:  multiply a tensor by a scalar                                      *
//                                                                             *
//Takes: takes the scalar and a tensor as an argument                          *
//******************************************************************************
void tensor::Multiply(tensor const &rval, double scalar)
{
 int i;

 Mimic(rval);

 for(i = 0; i < rval.p->product; i++)
    p->m[i] = scalar * rval.p->m[i];
}


//******************************************************************************
//Name:  tensor addition                                                       *
//                                                                             *
//Purpose:  add two tensors together                                           *
//                                                                             *
//Takes: the address of the rval tensor                                        *
//******************************************************************************
tensor tensor::operator+(tensor const &B) const
{
  int i;

  if( p->product != B.p->product )
    error("Addition error:","Tensors must be same size");

  if( p->num_indices != B.p->num_indices)
    error("Addition error:","Tensors must have same index structure");

  for( i = 0; i < p->num_indices ; i++ )
    if( p->range[i] != B.p->range[i] )
      error("Addition error:","Tensor indices must have same ranges");

//Create a dummy tensor then Mimic on the fly
  tensor sum;
  sum.Mimic(B);

//Perform the addition
  for( i = 0; i < p->product ; i++ )
    sum.p->m[i] = p->m[i] + B.p->m[i];

  return sum;
}

//******************************************************************************
//Name:  tensor addition                                                       *
//                                                                             *
//Purpose:  add two tensors together                                           *
//                                                                             *
//Takes: the address of the rval tensor                                        *
//******************************************************************************
/*
tensor operator+(tensor &A, tensor &B)
{
  tensor sum;

  sum <= A.Addition(B);

  return sum;
}
*/
//******************************************************************************
//Name:  tensor subtraction                                                    *
//                                                                             *
//Purpose:  subtract two tensors                                               *
//                                                                             *
//Takes: the address of the rval tensor                                        *
//******************************************************************************
tensor tensor::operator-(tensor const &B) const
{
  int i;

  if( p->product != B.p->product )
    error("Subtraction error:","Tensors must be same size");

  if( p->num_indices != B.p->num_indices)
    error("Subtraction error:","Tensors must have same index structure");

  for( i = 0; i < p->num_indices ; i++ )
    if( p->range[i] != B.p->range[i] )
      error("Subtraction error:","Tensor indices must have same ranges");

//Create a dummy tensor then Mimic on the fly
  tensor sum;
  sum.Mimic(B);

//Perform the subtraction
  for( i = 0; i < p->product ; i++ )
    sum.p->m[i] = p->m[i] - B.p->m[i];

  return sum;
}

//******************************************************************************
//Name:  tensor product                                                        *
//                                                                             *
//Purpose:  form the tensor, outer, or Kronicker product of two tensors        *
//                                                                             *
//Takes: the address of the rval tensor                                        *
//******************************************************************************
tensor tensor::operator*(tensor const &B) const
{
  int i, j , counter;
  int tot_indices;
  int *range;

//Determine the number of indices
  tot_indices = p->num_indices + B.p->num_indices;

//tensor product should result in at least two indices -- if not throw an error
  if( tot_indices < 2)
    error("Tensor product error:","Use ScalarMult");

//Check for scalar and psuedoscalars for both objects here
  if( B.p->product == 1 ) tot_indices = p->num_indices;
  if( p->product   == 1 ) tot_indices = B.p->num_indices;

//allocate temporary range array
  range  =  new int[tot_indices];

//pack the range appropriately to reflect index structure
//of the tensor product object
  for( i = 0; i < p->num_indices; i++)
    range[i]  =  p->range[i];

//test to see if this is a scalar -- if so change the for loop
  if( p->product == 1 )
  {
	  for( i = 0; i < tot_indices; i++)
		  range[i] = B.p->range[i];
  }
  else
  {
      for( i = p->num_indices; i < tot_indices; i++)
          range[i]  = B.p->range[i - p->num_indices];
  }

//Create a dummy tensor then resize on the fly
  tensor op;
  op.Resize0(tot_indices, range);

  delete [] range;

//Perform the product
  counter = 0;
  for( i = 0; i < p->product; i++ )
    for( j = 0; j < B.p->product; j++ )
      {
        op.p->m[counter] = p->m[i]*B.p->m[j];
        counter++;
      }

  return op;
}

//******************************************************************************
//Name:     Contract                                                           *
//                                                                             *
//Purpose:  Form the contraction of the two specified tensors on the selected  *
//          indices                                                            *
//                                                                             *
//Takes:    the of the rval tensor                                             *
//******************************************************************************
tensor tensor::Contract(tensor const &B, int ind1, int ind2) const
{
  int i, counter, r_counter;
  int tot_indices;
  int *range;

//Error check the indices and make sure neither is out of range
  if( ind1 >= p->num_indices )
    error("Contraction error:","LVAL index out of range");

  if( ind2 >= B.p->num_indices )
    error("Contraction error:","RVAL index out of range");

  if( p->range[ind1] != B.p->range[ind2] )
	error("Contraction error:","contracted indices must have the same range");

//Determine the total number of indices on the new tensor
  tot_indices = p->num_indices + B.p->num_indices - 2;

//Allocate space for temporary range
  if( tot_indices == 0 )
   {
     range = new int[1];
   }
  else
   {
     range  =  new int[tot_indices];
   }

  r_counter = 0;

//Pack the range array appropriately
  counter = 0;
  for( i = 0; i < p->num_indices; i++)
   {
      if(counter != ind1)
	  {
        range[r_counter]  =  p->range[i];
		  r_counter++;
	  }
      counter++;
   }

  counter = 0;
  //note that tot_indices + 2 is needed to loop over all the appropriate ranges
  for( i = p->num_indices; i < tot_indices + 2; i++)
    {
      if(counter != ind2)
	  {
        range[r_counter] = B.p->range[i - p->num_indices];
		  r_counter++;
	  }
      counter++;
    }

//create a temp tensor and handle scalar result
  tensor con;
  if( tot_indices == 0 )
   {
     tot_indices = 1;
     range[0] = 1;
   }
  con.Resize0(tot_indices, range);
  delete [] range;

//Now comes the fun part
  int *index, *B_index, *con_index;
  int j;
  double contract;
  int temp;

//Allocate index arrays to hold the current index structure
  index      = new int[p->num_indices];
  B_index    = new int[B.p->num_indices];
  con_index	 = new int[con.p->num_indices];

  for(i = 0; i < con.p->product; i++)
  {

   //pack the con_index array
   temp = 0;
	for(j = 0; j < con.p->num_indices - 1; j++)
    {
      con_index[j] = (i - i%con.p->scales[j] - temp)/con.p->scales[j];
      temp += con_index[j] * con.p->scales[j];
    }
    con_index[j] = (i - temp)/con.p->scales[j];

	//form the index structure for the LVAL object
	counter = 0;
	r_counter = 0;
    for( j = 0; j < p->num_indices; j++)
	{
      if(counter != ind1)
	  {
		index[j]  =  con_index[r_counter];
		r_counter++;
	  }
      counter++;
    }

	//form the index structure for the RVAL object
    counter = 0;
    for( j = 0; j < B.p->num_indices; j++)
	{
      if(counter != ind2)
	  {
      B_index[j] = con_index[r_counter];
		r_counter++;
	  }
      counter++;
    }

	//form the contraction
	contract = 0;
	for( j = 0; j < p->range[ind1]; j++)
	{
	  index[ind1]   = j;
	  B_index[ind2] = j;
	  contract += Val(index)*B.Val(B_index);
	 }
    con.p->m[i] = contract;
  }

  //clean up
  delete [] index;
  delete [] B_index;
  delete [] con_index;

  //return the result
  return con;
}

//******************************************************************************
//Name:     Contract                                                           *
//                                                                             *
//Purpose:  Form the self contraction between a given tensor on the selected   *
//          indices                                                            *
//                                                                             *
//Takes:    the selected indices                                               *
//******************************************************************************
tensor tensor::Contract(int ind1, int ind2)
{
  int i, counter, r_counter;
  int tot_indices;
  int *range;

//Error check the indices and make sure neither is out of range
  if( ind1 >= p->num_indices )
    error("Contraction error:","LVAl index out of range");

  if( ind2 >= p->num_indices )
    error("Contraction error:","RVAL index out of range");

  if( p->range[ind1] != p->range[ind2] )
	error("Contraction error:","contracted indices must have the same range");

//Determine the total number of indices on the new tensor
  tot_indices = p->num_indices - 2;

//Create the range array if needed
  if( tot_indices > 0 )
  {
    //Allocate space for temporary range
    range  =  new int[tot_indices];

    //Pack the range array appropriately
    counter   = 0;
    r_counter = 0;
    for( i = 0; i < p->num_indices; i++)
	  {
       if(counter != ind1 && counter != ind2)
	     {
          range[r_counter]  =  p->range[i];
		    r_counter++;
		  }
       counter++;
     }
  }

//create the result
  tensor con;
  if( tot_indices != 0 )
  {
    con.Resize0(tot_indices, range);
    delete [] range;
  }
  else
  {
    range = new int;
    range[0] = 1;
    con.Resize0(1,range);
    delete [] range;
  }

//Now comes the fun part
  int *index, *con_index;
  int j;
  double contract;
  int temp;

//Allocate index arrays to hold the current index structure
  index      = new int[p->num_indices];
  con_index	 = new int[con.p->num_indices];

  for(i = 0; i < con.p->product; i++)
  {

   //pack the con_index array
   temp = 0;
	for(j = 0; j < con.p->num_indices - 1; j++)
    {
      con_index[j] = (i - i%con.p->scales[j] - temp)/con.p->scales[j];
      temp += con_index[j] * con.p->scales[j];
    }
    con_index[j] = (i - temp)/con.p->scales[j];

	//form the index structure for the tensor
	counter = 0;
	r_counter = 0;
    for( j = 0; j < p->num_indices; j++)
	{
      if(counter != ind1)
	  {
		index[j]  =  con_index[r_counter];
		r_counter++;
	  }
      counter++;
    }

	//form the contraction
	contract = 0;
	for( j = 0; j < p->range[ind1]; j++)
	{
	  index[ind1] = j;
	  index[ind2] = j;
	  contract += Val(index);
	 }
    con.p->m[i] = contract;
  }

  //clean up
  delete [] index;
  delete [] con_index;

  //return the result
  return con;
}
//******************************************************************************
//Name:  tensor operator assignment                                            *
//                                                                             *
//Purpose:  set the memory of one object equal to another                      *
//                                                                             *
//Takes: address of array to be assigned                                       *
//******************************************************************************
/*tensor tensor::operator=(const tensor &rval)
{
  if( --p->n == 0 )
    {
      delete p->m;
      delete p->range;
      delete p->scales;
      delete p;
    }

  rval.p->n++;
  p = rval.p;
  return *this;
}
*/
/*
//******************************************************************************
//Name:  tensor copy constructor                                               *
//                                                                             *
//Purpose:  point copy to the same memory and keep track of duplicates         *
//                                                                             *
//Takes: address of tensor to be copied                                        *
//******************************************************************************
tensor::tensor(tensor &x)
{
  x.p->n++;
  p = x.p;
}
//******************************************************************************
//Name:  tensor destructor                                                     *
//                                                                             *
//Purpose:  delete object but care after the memory                            *
//                                                                             *
//Takes: nothing                                                               *
//******************************************************************************
tensor::~tensor()
{
  if( --p->n == 0)
    {
      delete p->m;
      delete p->range;
      delete p->scales;
      delete p;
    }
}

*/
//******************************************************************************
//Name:     Vmag                                                               *
//                                                                             *
//Purpose:  To calculate the Euclidean norm for a one dimensional array/tensor *
//                                                                             *
//Takes:   const address of the desired tensor                                 *
//******************************************************************************
double tensor::Vmag(void)
{
	int i;
   double ret_val = 0.0;

   if ( p->num_indices != 1 ) error("Vmag error:","only defined for vectors");

   for( i = 0; i < p->range[0]; i++) ret_val = ret_val + p->m[i]*p->m[i];

   ret_val = sqrt(ret_val);

   return ret_val;

}

//******************************************************************************
//Name:     ScalProd                                                           *
//                                                                             *
//Purpose:  To calculate the Euclidean scalar product between two              *
//          one-dimensional tensors/arrays                                     *
//                                                                             *
//Takes:   const address of the desired tensors                                *
//******************************************************************************
double tensor::ScalProd(tensor const &A)
{
   double ret_val;


   if (  ( p->num_indices != 1 ) || ( A.p->num_indices != 1 ) )
      error("Vmag error:","only defined for vectors");

   ret_val =  Contract(A,0,0).Val(0);

   return ret_val;

}


/*
//******************************************************************************
//Name:                                                                        *
//                                                                             *
//Purpose:                                                                     *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************

//******************************************************************************
//Name:                                                                        *
//                                                                             *
//Purpose:                                                                     *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************

//******************************************************************************
//Name:                                                                        *
//                                                                             *
//Purpose:                                                                     *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
*/
