#ifndef ASYMPTOTICS_H
#define ASYMPTOTICS_H

#include <complex>
#include <iostream>
#include <cmath>

using namespace std;

namespace asymptotics
{
	const complex<double> one = 1.0;
	const complex<double> i(0, 1);
	const double eps = 1.0e-10;

	inline complex<double> eta(complex<double> nu, complex<double> z)
	{
		complex<double> r = z / nu;
		complex<double> s = sqrt(one+r*r);
//cout << "eta(): " << nu << "   " << z << "   " << r << "   " << s << endl;
		return (s + log(r / (one + s)));
	}

	inline complex<double> I(complex<double> nu, complex<double> z)
	{
		double real_part = 0.0, imag_part = 0.0;

		if (nu.imag() > eps)
		{
			//conjugate nu to utilize asymptotic expansion
			nu = conj(nu);
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);
			//cout << nu << "   " << z << "   " << local_eta << "   " << r << "   " << s << endl;
			//conjugate result to get I for original nu (since z is actually real)

			//construct answer
			real_part = ( conj( exp(nu*local_eta) / ( sqrt(2.0*M_PI*nu*s) ) ) ).real();
			imag_part = ( sin(M_PI*nu)*exp(-nu*local_eta) / sqrt(2.0*M_PI*nu*s) ).imag();
		}
		else if ( abs(nu.imag()) < eps && nu.real() < -eps )	//nu is negative and real
		{
			nu *= -1.0;
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);
			real_part = ( exp(nu*local_eta) / ( sqrt(2.0*M_PI*nu*s) ) ).real()
						+ ( 2.0*sin(M_PI*nu)*exp(-nu*local_eta) / sqrt(2.0*M_PI*nu*s) ).real();
			imag_part = 0.0;
		}
		else
		{
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);
			//cout << nu << "   " << z << "   " << local_eta << "   " << r << "   " << s << endl;

			//construct answer
			real_part = ( exp(nu*local_eta) / ( sqrt(2.0*M_PI*nu*s) ) ).real();
			imag_part = -( sin(M_PI*nu)*exp(-nu*local_eta) / sqrt(2.0*M_PI*nu*s) ).imag();
		}

		return ( real_part + i * imag_part );
	}

	inline complex<double> Iprime(complex<double> nu, complex<double> z)
	{
		double real_part = 0.0, imag_part = 0.0;

		if (nu.imag() > eps)
		{
			//conjugate nu to utilize asymptotic expansion
			nu = conj(nu);
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);

			//construct answer
			real_part = ( conj( sqrt(s)*exp(nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ) ).real();
			imag_part = -( sqrt(s)*sin(M_PI*nu)*exp(-nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ).imag();
		}
		else if ( abs(nu.imag()) < eps && nu.real() < -eps )	//nu is negative and real
		{
			nu *= -1.0;
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);
			real_part = ( sqrt(s)*exp(nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ).real()
						+ ( 2.0*sqrt(s)*sin(M_PI*nu)*exp(-nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ).real();
			imag_part = 0.0;
		}
		else
		{
			complex<double> local_eta = eta(nu, z);
			complex<double> r = z / nu;
			complex<double> s = sqrt(one+r*r);

			//construct answer
			real_part = ( sqrt(s)*exp(nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ).real();
			imag_part = ( sqrt(s)*sin(M_PI*nu)*exp(-nu*local_eta) / ( r*sqrt(2.0*M_PI*nu) ) ).imag();
		}

		return ( real_part + i * imag_part );
	}

	//two versions are the same for k >= k_critical, but differ for k < k_critical:
	// (1) - v1 is just a constant
	// (2) - v2 decays to zero (this one probably makes more sense)
	inline complex<double> csc_h_v1(complex<double> z)
	{
		return ( 2.0*exp(-z.real() ) );
	}

	inline complex<double> csc_h_v2(complex<double> z)
	{
		return ( ( 2.0*exp(-z) ).real() );
	}
}

#endif
