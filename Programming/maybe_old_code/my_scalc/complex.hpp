// Complex variable class

// From "Introducing C++ for Scientists, Engineers, and Mathematicians"

// D.M. Capper



#ifndef COMPLEX_HPP

#define COMPLEX_HPP



#include <iostream.h>

#include <math.h>



class complex {

		friend complex operator+(double x, const complex &v);

		friend complex operator+(const complex &u, double x);

		friend complex operator+(const complex &u, const complex &v);

		friend complex operator*(double x, const complex &v);

		friend complex operator*(const complex &u, double x);

		friend complex operator*(const complex &u, const complex &v);

		friend complex operator/(const complex &u, const complex &v);

		friend double real(const complex &z);

		friend double imag(const complex &z);

		friend double mod(const complex &z);

		friend complex conj(const complex &z);

		friend complex exp(const complex &z);



public:

		double &real(void);

		double &imag(void);



		complex(void);

		complex(double r, double i);

		complex(const complex &z);



		complex &operator=(double x);

		complex &operator=(const complex &z);



		complex &operator+=(double x);

		complex &operator+=(const complex &z);



		complex operator-() const;



private:

		double re, im;

};



// friend functions:

inline complex operator+(double x, const complex &v)

{

	return complex(x+v.re,v.im);

}



inline complex operator+(const complex &u, double x)

{

	return complex(u.re+x, u.im);

}



inline complex operator+(const complex &u, const complex &v)

{

	return complex(u.re + v.re, u.im + v.im);

}



inline complex operator*(double x, const complex &v)

{

	return complex(x*v.re, x*v.im);

}



inline complex operator*(const complex &u, double x)

{

	return complex(u.re*x, u.im*x);

}



inline complex operator*(const complex &u, const complex &v)

{

	return complex(u.re * v.re - u.im * v.im, u.im * v.re + v.im * u.re);

}


inline complex operator/(const complex &u, const complex &v)

{

/* Derivation of complex division:

	x = r*cos(theta), y = r*sin(theta) where r is modulus and theta is phase of z = x + i*y 
	  => z = r(cos(theta) + i*sin(theta)) = r*exp(i*theta)

    => r = sqrt(x*x + y*y) = mod(z), theta = acos(x/r)

    z3 = z1/z2 = [r1*exp(i*theta1)] / [[r1*exp(i*theta1)]
       = (r1/r2)*exp(i*(theta1-theta2))
	   => x3 = (r1/r2)*cos(theta1-theta2), y3 = (r1/r2)*sin(theta1-theta2)
      
*/

	double r1 = mod(u);
	double theta1 = acos(u.re/r1);
	double r2 = mod(v);
	double theta2 = acos(v.re/r2);

	return complex(r1/r2 * cos(theta1-theta2),r1/r2 * sin(theta1-theta2));

}


inline double real(const complex &z)

{

	return z.re;

}



inline double imag(const complex &z)

{

	return z.im;

}



inline double mod(const complex &z)

{

	return(sqrt(z.re * z.re + z.im * z.im));

}



inline complex conj(const complex &z)

// Complex conjugate of z.

{

	return complex(z.re, -z.im);

}



inline complex exp(const complex &z)

// Exponential function for complex argument.

{

	double temp = exp(z.re);



	return complex(temp*cos(z.im), temp * sin(z.im));

}



// member functions and operators:

inline double &complex::real(void)

{

	return re;

}



inline double &complex::imag(void)

{

	return im;

}



inline complex::complex(void)

{

	re = 0.;

	im = 0.;

}



inline complex::complex(double r, double i)

{

	re = r;

	im = i;

}



inline complex::complex(const complex &z)

{

	re = z.re;

	im = z.im;

}



inline complex &complex::operator=(double x)

{

	re = x;

	im = 0.0;

	return *this;

}



inline complex &complex::operator=(const complex &z)

{

	re = z.re;

	im = z.im;

	return *this;

}



inline complex &complex::operator+=(double x)

{

	re += x;

	return *this;

}



inline complex &complex::operator+=(const complex &z)

{

	re += z.re;

	im += z.im;

	return *this;

}



inline complex complex::operator-() const

{

	return complex(-re, -im);

}





#endif // COMPLEX_HPP

