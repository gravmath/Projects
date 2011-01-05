// Util.cpp: implementation of the Util class.
//
//////////////////////////////////////////////////////////////////////

#include "util.hpp"      //<================= Changed to lower case 

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Util::Util()
{

}

// Member functions
int Util::factorial(int x)
{
  int output=1;

	if (x>0)
	{		
		for (int i=1; i<=x; i++)
		{
			output *= i;
		}
	}
	else if (x==0)
	{
		output = 1;
	}
	else
	{
		cout << "ERROR: Factorial is not defined for x<0.\n";
	}
		
	return output;		
}


complex Util::sphericalHarmonicY(int l, int m, double theta, double phi)
{
	complex output;

	if (l>=m)
	{
	 output = sqrt((2.*(double)l+1.)*
                  (double)factorial(l-m)/
                  (4.*PI*(double)factorial(l+m)))
		 *legendreP(l,m,cos(theta))*exp(complex(0.,m*phi));
	}
	else 
	{
		cout << "ERROR: l must be >= m.\n";
	}

	return output;
}

double Util::legendreP(int l, int m, double x)
{
	if (m>=0)
	{
		return plgndr(l,m,x);
	}
	else
	{
		return pow(-1,-m)*(double)factorial(l+m)/(double)factorial(l-m)*plgndr(l,-m,x);
	}
}
	
double Util::misnerF(int n, double r, double re, double delta)
{
	if (delta==0.)
	{
		cout << "ERROR: delta=0 in routine misnerF.\n";
	}

	if (r==0.)
	{
		cout << "ERROR: r=0 in routine misnerF.\n";
	}

	return sqrt((2.*(double)n+1.)/(2.*delta))*plgndr(n,0,(r-re)/delta)/r;	
}
		
double Util::dMisnerF(int n, double r, double re, double delta)
{
	double output = 0.;

	if (delta==0.)
	{
		cout << "ERROR: delta=0 in routine dMisnerF.\n";
	}

	if (r==0.)
	{
		cout << "ERROR: r=0 in routine dMisnerF.\n";
	}

	if (n==0)
	{
		output = -1./(r*r*sqrt(2.*delta));
	}
	else if (n==1)
	{
		output = re/(r*r) * sqrt(3./(2*delta*delta*delta));
	}
	else if (n==2)
	{
		output = (3.*(r*r-re*re)+delta*delta)/(r*r)*sqrt(5./(8.*delta*delta*delta*delta*delta));
	}
	else 
	{
		cout << "ERROR: dMisnerF has only been calculated to n=2.\n";
	}
	return output;
}
		

double Util::sphBesJ(int n, double x)
{
	double output;

	if (x==0.)
	{
		cout << "ERROR: x=0 in routine sphBesJ.\n";
	}

	if (n==0)
	{
		output = sin(x)/x;
	}
	else if (n==1)
	{
		output = (-cos(x) + sin(x)/x)/x;
	}
	else if (n==2)
	{
		output = (-sin(x) - 3.*cos(x)/x + 3.*sin(x)/(x*x))/x;
	}		
	else if (n==3)
	{
		output = (cos(x) - 6.*sin(x)/x - 15.*cos(x)/(x*x) + 15.*sin(x)/(x*x*x))/x;
	}		
	else
	{
		cout << "ERROR: sphBesJ has not been calculated for l!=(0,3) \n";
	}

	return output;
}

double Util::sphNeuN(int n, double x)
{

	double output;

	if (x==0.)
	{
		cout << "ERROR: x=0 in routine sphNeuJ.\n";
	}

	if (n==0)
	{
		output = -cos(x)/x;
	}
	else if (n==1)
	{
		output = (-sin(x) - cos(x)/x)/x;
	}
	else if (n==2)
	{
		output =(cos(x) - 3.*sin(x)/(x) - 3.*cos(x)/(x*x))/x;
	}
	else if (n==3)
	{
		output = (sin(x) + 6.*cos(x)/x - 15.*sin(x)/(x*x) - 15.*cos(x)/(x*x*x))/x;
	}
	else
	{
		cout << "ERROR: sphNeuN has not been calculated for l!=(0,3) \n";
	}

	return output;
}


// From here down, all functions are modified from Numerical Recipes in C
double Util::plgndr(int l, int m, double x)
{
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;

	if (x*x>1. && m!=0)
	{
		cout << "ERROR: m=" << m << " and x*x=" << x*x << ">1 in routine plgndr. \n";
	}

	if (m < 0 || m > l)
		cout << "Bad arguments in routine plgndr\n";
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt(1.0-x*x);
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=m+2;ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}

void Util::invert(tensor a)
{
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	int *r;
	int num_index = a.NumIndices();

	if (num_index != 2)
	{
		cout << "ERROR (invert): This routine can only invert rank-2 tensors.\n";
	}

	r = a.Ranges();

	int n = r[0];

	array indxc(1,n);
	array indxr(1,n);
	array ipiv(1,n);

	for (j=0;j<n;j++) ipiv.Set(0,j);
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((int)ipiv.Val(j) != 1)
				for (k=0;k<n;k++) {
					if ((int)ipiv.Val(k) == 0) {
						if (fabs(a.Val(j,k)) >= big) {
							big=fabs(a.Val(j,k));
							irow=j;
							icol=k;
						}
					} else if ((int)ipiv.Val(k) > 1) cout << "ERROR (invert): Singular Matrix-1\n";
				}
		ipiv.Set(ipiv.Val(icol)+1,icol);
		if (irow != icol) {
			for (l=0;l<n;l++)
			{
				temp=a.Val(irow,l);
				a.Set(a.Val(icol,l),irow,l);
				a.Set(temp,icol,l);
			}
		}
		indxr.Set(irow,i);
		indxc.Set(icol,i);
		if (a.Val(icol,icol) == 0.0) cout << "ERROR (invert): Singular Matrix-2\n";
		pivinv=1.0/a.Val(icol,icol);
		a.Set(1.0,icol,icol);
		for (l=0;l<n;l++) a.Set(a.Val(icol,l)*pivinv,icol,l);
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a.Val(ll,icol);
				a.Set(0.0,ll,icol);
				for (l=0;l<n;l++) a.Set(a.Val(ll,l)-a.Val(icol,l)*dum,ll,l);
			}
	}
	for (l=n-1;l>=0;l--) {
		if ((int)indxr.Val(l) != (int)indxc.Val(l))
			for (k=0;k<n;k++)
			{
				temp=a.Val(k,(int)indxr.Val(l));
				a.Set(a.Val(k,(int)indxc.Val(l)),k,(int)indxr.Val(l));
				a.Set(temp,k,(int)indxc.Val(l));
			}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software +%X(|y. */

#undef NUSE1
#undef NUSE2
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN
#undef PI
#undef NRANSI
#undef RTPIO2
