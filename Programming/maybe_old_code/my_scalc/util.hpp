// Util.hpp: interface for the Util class.
//
//////////////////////////////////////////////////////////////////////

#ifndef UTIL_HPP
#define UTIL_HPP
#include "complex.hpp"
#include "array.hpp"
#include "tensor.hpp"
#include "complex.hpp"
#include "constant.hpp"
#include <math.h>
#include <iostream.h>
#define EPS 1.0e-7
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define RTPIO2 1.2533141
#define NUSE1 5
#define NUSE2 5
#define NRANSI

class Util  
{
public:
	Util(void);
	int factorial(int x);
	complex sphericalHarmonicY(int l, int m, double theta, double phi);
	double misnerF(int n, double r, double re, double delta);
	double dMisnerF(int n, double r, double re, double delta);
	double sphBesJ(int n, double x);
	double sphNeuN(int n, double x);
	double legendreP(int l, int m, double x);
	double plgndr(int l, int m, double x);
	void gaussj(tensor a, tensor b);
	void invert(tensor a);
};


#endif // UTIL_HPP
