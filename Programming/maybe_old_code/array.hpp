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
    array(void);
    array(int num_indices, int range1, ...);
    array(array &x);
    virtual ~array();
    array& operator =  (const array &rval);
    array& operator <= (const array &rval);

    int     NumIndices(void) const;
    int*    Ranges(void) const;
    int     Resize(int num_indices, int range1, ...);
    int     Resize(int num_indices, int *range1);
    double  Val(int range1, ...) const;
	 double  Val(int *range) const;
    void    Set(double value, int range1, ...);
    void    Set(double value, int *range);
    void    SetName(char *name);
    double  max(void) const;
    void    print(const char *msg, int range1, ...) const;
    void    print(const char *msg, char **buffer) const;
    void    print(char **buffer) const;
};

#define _ARRAYHPP
#endif
