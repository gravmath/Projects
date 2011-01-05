//******************************************************************************
//Name:  tensor.hpp                                                            *
//                                                                             *
//Purpose:  header file for defining the tensor class                          *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
// *****************************************************************************
#ifndef _TENSORHPP

#include "array.hpp"

class tensor : public array
{
  protected:

  public:
    tensor(void);
    tensor(int num_indices, int range1, ...);
    tensor& operator<=(const tensor &rval);
    tensor operator+(tensor const &B) const;
    tensor operator-(tensor const &B) const;
    tensor operator*(tensor const &B) const;
    tensor Contract(tensor &B, int ind1, int ind2);
    tensor Contract(int ind1, int ind2);
    void ScalarMult(double scalar);
    friend class field;
};
    tensor operator*(tensor const &B, double scalar);
    tensor operator*(double scalar, tensor const &B);


#define _TENSORHPP
#endif

