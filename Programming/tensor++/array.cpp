//******************************************************************************
//Name:  array.cpp                                                             *
//                                                                             *
//Purpose:  source file for defining the array class                           *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
//                       05/21/99 - improve speed performance by coding member *
//                                  functions as direct access for 7 or fewer  *
//                                  indices - affected Set, Val,               *
//                       05/21/99 - made the straight Set and Val calls inlines*
//                       05/23/99 - made inline functions to return the number *
//                                  of array elements and the scales array     *
//                                  made the NumIndices inline                 *
//                                  made a Resize0 for speed                   *
//                       06/18/99 - renamed Resize0 to Mimic and moved the data*
//                                  structure to public                        *
//                       06/19/99 - new implementation of Resize0              *
//								 04/05/03 - added path name for array.hpp include
//******************************************************************************          

#include "\school\Programming\tensor++\array.hpp"

//******************************************************************************
//Name:  array constructor                                                     *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: no args and is the default constructor                                *
//******************************************************************************

array::array(void)
{
//Allocate space for the structure that holds the array
  p = new array_rep;

//Default number of indices and product
  p->num_indices = 1;
  p->product = 3;

//Allocate the range array and set its default
  p->range = new int [p->num_indices];
  p->range[0] = 3;

//Allocate the scales array and set its default
  p->scales = new int [p->num_indices];
  p->scales[0] = 1;

//Allocate space for the actual array to accomodate operator = definition
  p->m = new double [p->product];

//Initialize the array to zeros
  memset(p->m,0,p->product*sizeof(double));

//Record number of objects pointing to the same memory
  p->n = 1;
}

//******************************************************************************
//Name:  array constructor                                                     *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: variable argument list with the number of indices and range for each  *
//******************************************************************************
array::array(int num_indices, int range1, ...)
{
  va_list arg_pnt;
  int product;
  int range_i;
  int i;
  int one_already = 0;

//Allocate space for the structure that handles the array and its relations
  p = new array_rep;

//Set the structures number of indices
  p->num_indices = num_indices;

//Allocate space for the relations that allow indexing into the array
  p->range  = new int [p->num_indices];
  p->scales = new int [p->num_indices];

//Set a pointer to the variable list
  va_start(arg_pnt, range1);

//Fill in the range array
  p->range[0] = range1;
  product = range1;
  for( i = 1; i < num_indices ; i++)
   {
     range_i = va_arg(arg_pnt, int);
     p->range[i] = range_i;
     product *= range_i;
   }

//Fill in the scales array
  p->scales[0] = product/p->range[0];

  for( i = 1; i< num_indices; i++)
   {
     p->scales[i] = p->scales[i-1]/p->range[i];
     if( p->scales[i] == 1 && one_already == 1)
       p->scales[i] = 0;
     if( p->scales[i] == 1 && one_already == 0)
       one_already = 1;
   }

//Clean up the variable arg list pointer
  va_end(arg_pnt);

//Set the structure's knowledge of the number of total items in the array
  p->product = product;

//Allocate space for the actual array
  p->m = new double [ product ];

//Initialize the array to zeros
  memset(p->m,0,product*sizeof(double));

//Record number of objects pointing to the same memory
  p->n = 1;
}

//******************************************************************************
//Name:  array copy constructor                                                *
//                                                                             *
//Purpose:  point copy to the same memory and keep track of duplicates         *
//                                                                             *
//Takes: address of array to be copied                                         *
//******************************************************************************
array::array(array &x)
{
  x.p->n++;
  p = x.p;
}

//******************************************************************************
//Name:  array operator assignment                                             *
//                                                                             *
//Purpose:  set the memory of one object equal to another                      *
//                                                                             *
//Takes: address of array to be assigned                                       *
//******************************************************************************
array& array::operator=(const array &rval)
{
  if( --p->n == 0 )
    {
      delete [] p->m;
      delete [] p->range;
      delete [] p->scales;
      delete p;
    }

  rval.p->n++;
  p = rval.p;
  return *this;
}

//******************************************************************************
//Name:  array non-referencing assignment                                      *
//                                                                             *
//Purpose:  memcpy the memory of one array to another                          *
//                                                                             *
//Takes: address of array to be assigned                                       *
//******************************************************************************
array& array::operator<=(const array &rval)
{
  int i;
  int n, *r;
  int test = 0;

//Check to see if the two arrays have the same size and number of indices
//if not resize this
  if( (p->product != rval.p->product) ||
      (p->num_indices != rval.p->num_indices) )
   {
     n = rval.NumIndices();
     r = rval.Ranges();
     Resize(n,r);
	  delete [] r;
   }

//Check to see if the two arrays have the same ranges
//if not resize this
  for( i = 0; i < p->num_indices ; i++ )
    if( p->range[i] != rval.p->range[i] ) test = 1;

  if( test == 1 )
   {
     n = rval.NumIndices();
     r = rval.Ranges();
     Resize(n,r);
	  delete [] r;
   }

//Now memcpy the rval arrays values into this
  memcpy( p->m, rval.p->m, p->product*sizeof(double) );

  return *this;
}

//******************************************************************************
//Name:  array destructor                                                      *
//                                                                             *
//Purpose:  delete object but care after the memory                            *
//                                                                             *
//Takes: nothing                                                               *
//******************************************************************************
array::~array()
{
  if( --(p->n) == 0)
    {
      delete [] p->m;
      delete [] p->range;
      delete [] p->scales;
      delete p;
    }
}

//******************************************************************************
//Name:  array Ranges                                                          *
//                                                                             *
//Purpose:  returns a pointer to an array containing the ranges for each index *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
int* array::Ranges(void) const
{
  int *range;

//Allocate space for the range array
  range = new int[p->num_indices];

//Now memcpy the rval arrays values into this
  memcpy( range, p->range, p->num_indices*sizeof(int) );

  return range;
}

//******************************************************************************
//Name:  array Scales                                                          *
//                                                                             *
//Purpose:  returns a pointer to an array containing the scales for each index *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
int* array::Scales(void) const
{
  int *scales;

//Allocate space for the range array
  scales = new int[p->num_indices];

//Now memcpy the rval arrays values into this
  memcpy( scales, p->scales, p->num_indices*sizeof(int) );

  return scales;
}


//******************************************************************************
//Name:  array SetName                                                         *
//                                                                             *
//Purpose:  set the name in the array with an input string (for debug)         *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
void array::SetName(char *name)
{

//string copy the name into the member data
  strcpy(p->name,name);

}


//******************************************************************************
//Name:  array Resize                                                          *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: number of indices, first range, ...                                   *
//******************************************************************************
int array::Resize(int num_indices, int range1, ...)
{

  int i, *range;

//Allocate space for the range array
  range = new int[num_indices];

  va_list arg_resize;

//Set up the pointer to the variable arg list
  va_start(arg_resize, range1);

//Pack the range array
  range[0] = range1;
  for(i = 1; i < num_indices; i++)
	  range[i] = va_arg(arg_resize,int);

//Call Resize to change size of structure parameters
  if( !Resize(num_indices, range) ) return FALSE;

//Clean up
  va_end(arg_resize);
  delete [] range;

  return TRUE;
}

//******************************************************************************
//Name:  array Resize      -- To Be Called from a Derived Class                *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing       *
//                                                                             *
//Takes: number of indices, range array,                                       *
//******************************************************************************
int array::Resize(int num_indices, int *range)
{
  int product;
  int i;
  int one_already;

  one_already = 0;

//check to make sure nobody else is pointing at the memory
  if( p->n != 1)
     error("Resize error:",
           "Cannot Resize - multiple objects pointing to the same memory");

  delete [] p->m;
  delete [] p->range;
  delete [] p->scales;

//Reallocate the range and scale arrays with the new size
  p->range   = new int [num_indices];
  p->scales  = new int [num_indices];

//Update the structure's knowledge of the number of indices
  p->num_indices = num_indices;

//Setup range array
  p->range[0] = range[0];
  product = range[0];
  for( i = 1; i < num_indices ; i++)
   {
     if( range[i] == 0 )
        error("Resize error:","missing or bad index");

     p->range[i] = range[i];
     product *= range[i];
   }

//Setup scales array
  p->scales[0] = product/p->range[0];
  for( i = 1; i< num_indices; i++)
   {
     p->scales[i] = p->scales[i-1]/p->range[i];
     if( p->scales[i] == 1 && one_already == 1)
       p->scales[i] = 0;
     if( p->scales[i] == 1 && one_already == 0)
       one_already = 1;
   }

//Tell the structure how many array elements there are
  p->product = product;

//Allocate space for the actual array
  p->m = new double [ product ];

//Initialize the array to zeros
  memset(p->m,0,product*sizeof(double));

//Reset the number of things pointing
  p->n = 1;
  return TRUE;
}

//******************************************************************************
//Name:  array Resize0      -- To Be Called from a Derived Class (for speed)   *
//                                                                             *
//Purpose:  allocate flat space for storage and index data for accessing but   *
//          doesn't initialize the data array
//                                                                             *
//Takes: number of indices, range array,                                       *
//******************************************************************************
int array::Resize0(int num_indices, int *range)
{
  int product;
  int i;
  int one_already;

  one_already = 0;

//check to make sure nobody else is pointing at the memory
  if( p->n != 1)
     error("Resize error:",
           "Canot Resize - multiple objects pointing to the same memory");

  delete [] p->m;
  delete [] p->range;
  delete [] p->scales;

//Reallocate the range and scale arrays with the new size
  p->range   = new int [num_indices];
  p->scales  = new int [num_indices];

//Update the structure's knowledge of the number of indices
  p->num_indices = num_indices;

//Setup range array
  p->range[0] = range[0];
  product = range[0];
  for( i = 1; i < num_indices ; i++)
   {
     if( range[i] == 0 )
        error("Resize error:","missing or bad index");

     p->range[i] = range[i];
     product *= range[i];
   }

//Setup scales array
  p->scales[0] = product/p->range[0];
  for( i = 1; i< num_indices; i++)
   {
     p->scales[i] = p->scales[i-1]/p->range[i];
     if( p->scales[i] == 1 && one_already == 1)
       p->scales[i] = 0;
     if( p->scales[i] == 1 && one_already == 0)
       one_already = 1;
   }

//Tell the structure how many array elements there are
  p->product = product;

//Allocate space for the actual array
  p->m = new double [ product ];

//Reset the number of things pointing
  p->n = 1;
  return TRUE;
}

//******************************************************************************
//Name:  array Mimic      -- To Be Called from a Derived Class for speed       *
//                                                                             *
//Purpose:  make the data structure for the lval object mimic the rval without *
//          copying the data values                                            *
//                                                                             *
//Takes: number of indices, range array,                                       *
//******************************************************************************
int array::Mimic(const array &rval)
{

//check to make sure nobody else is pointing at the memory
  if( p->n != 1)
     error("Mimic error:",
           "Cannot Mimic - multiple objects pointing to the same memory");

  delete [] p->m;
  delete [] p->range;
  delete [] p->scales;

//Reallocate the range and scale arrays with the new size
  p->range   = new int [rval.p->num_indices];
  p->scales  = new int [rval.p->num_indices];

//Update the structure's knowledge of the number of indices
  p->num_indices = rval.p->num_indices;
  memcpy( p->range,  rval.p->range,  p->num_indices*sizeof(int) );
  memcpy( p->scales, rval.p->scales, p->num_indices*sizeof(int) );
  p->product = rval.p->product;

//Allocate space for the actual array
  p->m = new double [ p->product ];

//Reset the number of things pointing
  p->n = 1;
  return TRUE;
}


//******************************************************************************
//Name:  array Val                                                             *
//                                                                             *
//Purpose:  returns to value of the array given the indices                    *
//                                                                             *
//Takes: variable list of indices                                              *
//******************************************************************************
double array::Val(int range0,
				      int range1,
   				   int range2,
				      int range3,
    				   int range4,
                  int range5,
                  int range6,
                  int range7, ...) const
{
  double return_val;
  int *range, i;

//Allocate space for the range array
  range = new int[p->num_indices];

//Initialize the variable argument list
  va_list arg_pnt;

  va_start(arg_pnt,range7);

//Pack the range array
  range[0] = range0;
  range[1] = range1;
  range[2] = range2;
  range[3] = range3;
  range[4] = range4;
  range[5] = range5;
  range[6] = range6;
  range[7] = range7;
  for(i = 8; i < p->num_indices; i++)
	  range[i] = va_arg(arg_pnt,int);

//Get the return value
  return_val = Val(range);

//Clean up
  va_end(arg_pnt);
  delete [] range;

  return return_val;
}

//******************************************************************************
//Name:  array Val                                                             *
//                                                                             *
//Purpose:  returns to value of the array given the indices                    *
//                                                                             *
//Takes: array of indices -- for use by derived classes                        *
//******************************************************************************
double array::Val(int *range) const
{
  int i, product;

  product = range[0] * p->scales[0];
  if( p->num_indices == 1 && p->range[0] == 1) product = 0;
  for( i = 1; i < p->num_indices; i++)
    product += range[i] * p->scales[i];

  if( product > p->product )
     error("Value error:","indices out of bounds");

  return p->m[product];
}

//******************************************************************************
//Name:  array Set                                                             *
//                                                                             *
//Purpose:  allows the value of the array to be set given the indices          *
//                                                                             *
//Takes: variable list of indices                                              *
//******************************************************************************
void array::Set(double value,
                int range0,
                int range1,
                int range2,
                int range3,
                int range4,
                int range5,
                int range6,
                int range7, ...)
{
  int *range, i;

//Allocate space for the range array
  range = new int[p->num_indices];

//Initialize the variable argument list
  va_list arg_pnt;

  va_start(arg_pnt,range7);

//Pack the range array
  range[0] = range0;
  range[1] = range1;
  range[2] = range2;
  range[3] = range3;
  range[4] = range4;
  range[5] = range5;
  range[6] = range6;
  range[7] = range7;
  for(i = 8; i < p->num_indices; i++)
   range[i] = va_arg(arg_pnt,int);

//Set the value
  Set(value, range);

//Clean up
  va_end(arg_pnt);
  delete [] range;
}

//******************************************************************************
//Name:  array Set                                                             *
//                                                                             *
//Purpose:  allows the value of the array to be set given the indices          *
//                                                                             *
//Takes: array of indices                                                      *
//******************************************************************************
void array::Set(double value, int *range)
{
  int product, i;

  product = range[0] * p->scales[0];
  if(p->num_indices == 1 && p->range[0] == 1) product = 0;
  for( i = 1; i < p->num_indices; i++)
    product += range[i] * p->scales[i];
  if( product > p->product )
    error("Set error:", "indices out of bound");

  p->m[product] = value;

}

//******************************************************************************
//Name:  array max                                                             *
//                                                                             *
//Purpose:  returns the maximum absolute value of an array                     *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
double array::max(void) const
{
  double ret_val = 0;
  double temp;
  int i;

  for( i = 0; i < p->product; i++)
   {
     temp = fabs(p->m[i]);
     if( temp > ret_val ) ret_val = temp;
   }

  return ret_val;

}
//******************************************************************************
//Name:  array print                                                           *
//                                                                             *
//Purpose:  prints out the value of the array given the indices                *
//                                                                             *
//Takes: a message string and variable list of indices                         *
//******************************************************************************
void array::print(const char *msg, int range1, ...) const
{
  int *range, i;

//Allocate space for the range array
  range = new int[p->num_indices];

//Initialize the variable argument list
  va_list arg_list;

  va_start(arg_list,range1);

//Pack the range array
  range[0] = range1;
  for(i = 1; i < p->num_indices; i++)
	  range[i] = va_arg(arg_list,int);

  double value;

  value = Val(range);

  if(*msg) cout << msg;

  cout << setw(6) << setprecision(3) << value;

//Clean up
  va_end(arg_list);
  delete [] range;
}

//******************************************************************************
//Name:  array print                                                           *
//                                                                             *
//Purpose:  prints out the values in the whole array with indices              *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
void array::print(const char *msg, char **buffer) const
{
  int *index, i, j, temp, counter, num_indices;
  int size_of_name;

  //find the size of the name
  //accounting for the [0,0,..] = ## format
  //using 2*p->num_indices - 1 for the 0,0,... format
  //using 5 for the [ ] space = space format
  //using 30 digits for the output
  size_of_name = strlen(msg);
  size_of_name = size_of_name + 2*p->num_indices - 1 + 5 + 30;

  if( *buffer != NULL ) delete [] *buffer;

  *buffer = new char[p->product*size_of_name + 1];
  (*buffer)[0] = '\0';

  //handle scalars separately
  if(p->product == 1)
   {
    //fill the buffer line
    counter = strlen(*buffer);
    strcpy((*buffer)+counter,msg);
    counter = strlen(*buffer) - 1;
    (*buffer)[counter + 1] = '[';
    (*buffer)[counter + 2] = '0';
    (*buffer)[counter + 3] = ']';
    (*buffer)[counter + 4] = ' ';
    (*buffer)[counter + 5] = '=';
    (*buffer)[counter + 6] = ' ';
    sprintf((*buffer)+ counter + 7,"%g\n",Val(0));
    counter = strlen(*buffer);
    (*buffer)[counter] = '\0';
    return;
   }

  //test for range = 1 indices
  num_indices = p->num_indices;
  for(i = 0; i < p->num_indices; i++)
  {
    if(p->range[i] == 1) num_indices--;
  }

  index  = new int[p->num_indices];

  for(i = 0; i < p->product; i++)
   {
    //pack the index array testing for range = 1 cases
    temp = 0;
	 for(j = 0; j < num_indices - 1; j++)
     {
       index[j] = (i - i%p->scales[j] - temp)/p->scales[j];
       temp += index[j] * p->scales[j];
     }
    //range = 1 testing starts here
    if( num_indices < p->num_indices)
     {
       index[num_indices - 1] = (i - temp)/p->scales[num_indices - 1];
       for(j = num_indices; j < p->num_indices; j++) index[j] = 0;
     }
    else
     {
       index[j] = (i - temp)/p->scales[j];
     }

    //fill the buffer line
    counter = strlen(*buffer);
    strcpy((*buffer)+counter,msg);
    counter = strlen(*buffer) - 1;
    (*buffer)[counter + 1] = '[';
    (*buffer)[counter + 2] = '\0';
    for(j = 0; j < p->num_indices; j++)
      {
        counter = strlen(*buffer);
        sprintf( ((*buffer) + counter),"%d,\0",index[j]);
      }
    counter = strlen(*buffer);
    (*buffer)[counter - 1] = ']';
    (*buffer)[counter]     = ' ';
    (*buffer)[counter + 1] = '=';
    (*buffer)[counter + 2] = ' ';
    sprintf((*buffer)+ counter + 3,"%g\n",Val(index));

   }
   //terminate the string properly
   counter = strlen(*buffer);
   (*buffer)[counter+1] = '\0';

   //clean-up
   delete [] index;

}

//******************************************************************************
//Name:  array print                                                           *
//                                                                             *
//Purpose:  prints out the values in the whole array without indices           *
//                                                                             *
//Takes:                                                                       *
//******************************************************************************
void array::print(char **buffer) const
{
  int *index, i, j, temp, num_indices, counter;

  //test for range = 1 indices
  num_indices = p->num_indices;
  for(i = 0; i < p->num_indices; i++)
  {
    if(p->range[i] == 1) num_indices--;
  }

  index  = new int[p->num_indices];
  if( *buffer != NULL ) delete [] *buffer;
  *buffer = new char[p->product*41 + 1];

  (*buffer)[0] = '\0';

  for(i = 0; i < p->product; i++)
   {
    //pack the index array testing for range = 1 cases
    temp = 0;
	 for(j = 0; j < num_indices - 1; j++)
     {
       index[j] = (i - i%p->scales[j] - temp)/p->scales[j];
       temp += index[j] * p->scales[j];
     }
    //range = 1 testing starts here
    if( num_indices < p->num_indices)
     {
       index[num_indices - 1] = (i - temp)/p->scales[num_indices - 1];
       for(j = num_indices; j < p->num_indices; j++) index[j] = 0;
     }
    else
     {
       index[j] = (i - temp)/p->scales[j];
     }

    //fill the buffer line
    counter = strlen(*buffer);
    sprintf((*buffer)+counter,"%g\n",Val(index));
   }

   //terminate the string properly
   counter = strlen(*buffer);
   (*buffer)[counter+1] = '\0';

   //clean-up
   delete [] index;

}

//******************************************************************************
//Name:  array error                                                           *
//                                                                             *
//Purpose:  prints out an error message                                        *
//                                                                             *
//Takes: two message strings           													 *
//******************************************************************************
void array::error(char *msg1, char *msg2) const
{
  cerr << "\nerror:  "     << msg1
       << " "              << msg2
       << endl;
  exit(1);
}
