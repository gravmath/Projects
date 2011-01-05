//******************************************************************************
//Name:  tensor.hpp                                                            *
//                                                                             *
//Purpose:  source file for defining the tensor class                          *
//                                                                             *
//Modification History:  10/24/98 - added modification history field           *
//								 5/18/99  - added vector magnitude field for K.Watt    *
//                       05/23/99 - updated to use Resize0 for speed           *
//                       06/18/99 - renamed Resize0 to Mimic                   *
//******************************************************************************
#ifndef _TENSORHPP

#include "array.hpp"

class tensor : public array
{
  protected:

  public:
    tensor(void);
    tensor(int num_indices, int range1, ...);

    tensor& operator<=(const tensor &rval);
    tensor  operator+(tensor const &B) const;
    tensor  operator-(tensor const &B) const;
    tensor  operator*(tensor const &B) const;
    tensor  Contract(tensor const &B, int ind1, int ind2) const;
    tensor  Contract(int ind1, int ind2);
    void    ScalarMult(double scalar);
    void    Multiply(tensor const &rval, double scalar);
    double  Vmag(void);
    double  ScalProd(tensor const &A);

    friend class Particle;
    friend class MSField;
};
    tensor operator*(tensor const &B, double scalar);
    tensor operator*(double scalar, tensor const &B);


#define _TENSORHPP
#endif
