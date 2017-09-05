#ifndef ASYMPTOTICS_H
#define ASYMPTOTICS_H

#include <complex>

namespace asymptotics
{
	inline complex<double> eta(complex<double> nu, complex<double> z)
	{
		complex<double> r = nu / z;
		complex<double> s = sqrt(1.0+r*r);
		return (s + log(r / (1.0 + s)));
	}

	inline complex<double> asymptotic_I(complex<double> nu, complex<double z>)
	{
		complex<double> local_eta = eta(nu, z);
		complex<double> r = nu / z;
		complex<double> s = sqrt(1.0+r*r);
		return ( exp(nu*local_eta) / ( sqrt(2.0*M_PI*nu*s) ) );
	}

	//two versions are the same for k >= k_critical, but differ for k < k_critical:
	// (1) - v1 is just a constant
	// (2) - v2 decays to zero (this one probably makes more sense)
	inline asymptotic_csc_h_v1(complex<double> z)
	{
		return ( 2.0*exp(-z.real() ) );
	}

	inline asymptotic_csc_h_v2(complex<double> z)
	{
		return ( ( 2.0*exp(-z) ).real() );
	}
}

#endif
