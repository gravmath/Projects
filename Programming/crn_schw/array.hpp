//******************************************************************************
//Name:  array.hpp                                                             *
//                                                                             *
//Purpose:  header file for defining the array class                           *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
// *****************************************************************************
#ifndef _ARRAYHPP

#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#define FALSE 0
#define TRUE  1
class array
{
  protected:
      struct array_rep
      {
        double *m;
        int num_indices;
        int *range;
        int n;
        int product;
        int *scales;
        char name[10];
      };

    struct array_rep *p;

    void    error( char *msg1, char *msg2 = "") const;

  public:

   //The big four
   array(void);
   array(int num_indices, int range1, ...);
   array(array &x);
   virtual ~array();

   //The operators
   array& operator =  (const array &rval);
   array& operator <= (const array &rval);

   //Member functions
   int*    Ranges(void) const;
   int*    Scales(void) const;
   int     Resize(int num_indices, int range1, ...);
   int     Resize(int num_indices, int *range1);
   int     Resize0(int num_indices, int *range);
   int     Mimic(const array &rval);
   void    Set(double value,
               int range0,
               int range1,
               int range2,
               int range3,
               int range4,
               int range5,
               int range6,
               int range7, ...);
   void    Set(double value, int *range);
   double  Val(int range0,
			      int range1,
   		      int range2,
				   int range3,
    			   int range4,
               int range5,
               int range6,
               int range7, ...) const;
   double  Val(int *range) const;
   void    SetName(char *name);
   double  max(void) const;
   void    print(const char *msg, int range1, ...) const;
   void    print(const char *msg, char **buffer) const;
   void    print(char **buffer) const;

   friend  class Particle;
   friend  class MSField;


   //Inlines start here
   inline int array::NumIndices(void) const
   {
      return p->num_indices;
   }
   inline int array::NumComponents(void) const
   {
   	return p->product;
   }


   inline double array::Val(int range0) const
   {
      return p->m[ range0 * p->scales[0] ];
   }
   inline double array::Val(int range0,
                            int range1) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1] ];
   }
   inline double array::Val(int range0,
                            int range1,
                            int range2) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1]
                +  range2 * p->scales[2] ];
   }
   inline double array::Val(int range0,
                            int range1,
                            int range2,
                            int range3) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1]
                +  range2 * p->scales[2] + range3 * p->scales[3] ];
   }
   inline double array::Val(int range0,
                            int range1,
                            int range2,
                            int range3,
                            int range4) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1]
                +  range2 * p->scales[2] + range3 * p->scales[3]
                +  range4 * p->scales[4] ];
   }
   inline double array::Val(int range0,
                            int range1,
                            int range2,
                            int range3,
                            int range4,
                            int range5) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1]
                +  range2 * p->scales[2] + range3 * p->scales[3]
                +  range4 * p->scales[4] + range5 * p->scales[5] ];
   }

   inline double array::Val(int range0,
                            int range1,
                            int range2,
                            int range3,
                            int range4,
                            int range5,
                            int range6) const
   {
      return p->m[ range0 * p->scales[0] + range1 * p->scales[1]
                +  range2 * p->scales[2] + range3 * p->scales[3]
                +  range4 * p->scales[4] + range5 * p->scales[5]
                +  range6 * p->scales[6] ];
   }


   inline void array::Set(double value,
                          int range0)
   {
      p->m[ range0 * p->scales[0] ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1] ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1,
                          int range2)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1]
         +  range2 * p->scales[2]                      ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1,
                          int range2,
                          int range3)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1]
         +  range2 * p->scales[2] + range3 * p->scales[3] ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1,
                          int range2,
                          int range3,
                          int range4)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1]
         +  range2 * p->scales[2] + range3 * p->scales[3]
         +  range4 * p->scales[4]                      ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1,
                          int range2,
                          int range3,
                          int range4,
                          int range5)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1]
         +  range2 * p->scales[2] + range3 * p->scales[3]
         +  range4 * p->scales[4] + range5 * p->scales[5] ] = value;
   }
   inline void array::Set(double value,
                          int range0,
                          int range1,
                          int range2,
                          int range3,
                          int range4,
                          int range5,
                          int range6)
   {
      p->m[ range0 * p->scales[0] + range1 * p->scales[1]
         +  range2 * p->scales[2] + range3 * p->scales[3]
         +  range4 * p->scales[4] + range5 * p->scales[5]
         +  range6 * p->scales[6]                      ] = value;
   }

};

#define _ARRAYHPP
#endif