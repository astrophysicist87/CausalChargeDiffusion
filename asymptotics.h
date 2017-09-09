#ifndef ASYMPTOTICS_H
#define ASYMPTOTICS_H

#include <complex>
#include <iostream>
#include <cmath>
#include <vector>

#include <gsl/gsl_sf_airy.h>

using namespace std;

#include "./special_function_library/cmlib/complex.h"
#include "./special_function_library/cmlib/protom.h"

extern const int n_integ_besselK_points;
extern vector<double> x_integ_besselK_pts, x_integ_besselK_wts;

namespace asymptotics
{
	//small number to avoid divisions by zero, etc.
	const double eps = 1.0e-10;

	//some useful numerical constants used below
	const complex<double> one = 1.0;
	const complex<double> i(0, 1);
	const complex<double> one_third = 1.0/3.0;
	const complex<double> two_thirds = 2.0/3.0;
	const complex<double> three_to_one_third = 1.4422495703074083823;
	const complex<double> three_to_two_thirds = 2.0800838230519041145;
	const complex<double> three_to_one_sixth = 1.2009369551760027267;
	const complex<double> g_1_3 = 2.6789385347077476337;
	const complex<double> g_2_3 = 1.3541179394264004169;
	const complex<double> m1_to_1_3 = 0.5 + 0.86602540378443864676 * i;
	const complex<double> three_by_two_to_two_by_three = 1.3103706971044483036;

	double nus[4] = { -2.0/3.0, -1.0/3.0, 1.0/3.0, 2.0/3.0 };
	double gamma_nus[4] = { -4.0184078020616214505, -4.0623538182792012508, 2.6789385347077476337, 1.3541179394264004169 };

	///////////////////////////////////////////////////////////////
	inline double I_asym(int nu_index, double z)
	{
		double nu = nus[nu_index];

		double sum1 = 0.0, sum2 = 0.0;
		double oldsum1 = 0.0, oldsum2 = 0.0;
		double next_term = 1.0;

		double sign = 1.0;

		long int k = 0;
		do
		{
			oldsum1 = sum1;
			oldsum2 = sum2;

			sum1 += sign * next_term;
			sum2 += next_term;

			double a = 2.0*(k+1.0) - 1.0;
			next_term *= (4.0*nu*nu - a*a) / (8.0*z*(k+1.0));
			sign *= -1.0;
			k++;
			//cout << "I_asym(): " << k-1 << "   " << oldsum1 << "   " << sum1
			//		<< "   " << oldsum2 << "   " << sum2 << endl;
		} while (
					2.0*abs(sum1-oldsum1)/abs(sum1+oldsum1) >= 1e-20
					|| 2.0*abs(sum2-oldsum2)/abs(sum2+oldsum2) >= 1e-20
					);

		double factor1 = exp(z) / sqrt(2.0*M_PI*z);
		double factor2 = exp(-z) / sqrt(2.0*M_PI*z);
		sign = 1.0;

		return ( factor1*sum1 );
	}

	///////////////////////////////////////////////////////////////
	inline double I_series(int nu_index, double z)
	{
		double nu = nus[nu_index];
		double gamma_nu = gamma_nus[nu_index];

		double Isum = 0.0, oldIsum = 0.0;
		double sign = 1.0;
		double next_term = 1.0/(nu*gamma_nu);
		double arg = 0.25*z*z;

		double k = 0.0;
		do
		{
			oldIsum = Isum;
			Isum += next_term;

			double new_factor = arg / ( (k+1.0) * (k+1.0+nu) );
			next_term *= new_factor;
			sign *= -1.0;
			k++;
		} while ( 2.0*abs(Isum-oldIsum)/abs(Isum+oldIsum) >= 1e-30 );

		return ( pow(0.5*z, nu)*Isum );
	}

	inline double BesselK_1_3(double z)
	{
		double result = 0.0;

		double min = 0.0;
		double max = 3.0 / pow(z, 1.0/3.0);
		double cen = 0.5*(min+max);
		double hw = 0.5*(max-min);

		for (int ix = 0; ix < n_integ_besselK_points; ++ix)
		{
			double x = cen + hw * x_integ_besselK_pts[ix];
			double cf = (1.0 + 4.0*x*x/3.0)*sqrt(1.0+x*x/3.0);
			result += hw * x_integ_besselK_wts[ix] * exp(-z * cf);
		}

		return ( sqrt(3.)*result );
	}

	inline double BesselK_2_3(double z)
	{
		double result = 0.0;

		double min = 0.0;
		double max = 3.0 / pow(z, 1.0/3.0);
		double cen = 0.5*(min+max);
		double hw = 0.5*(max-min);

		for (int ix = 0; ix < n_integ_besselK_points; ++ix)
		{
			double x = cen + hw * x_integ_besselK_pts[ix];
			double cf = (1.0 + 4.0*x*x/3.0)*sqrt(1.0+x*x/3.0);
			result += hw * x_integ_besselK_wts[ix] * exp(-z * cf)
						* ( 3.0 + 2.0*x*x ) / sqrt(1.0+x*x/3.0);
		}

		return ( result/sqrt(3.) );
	}

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	inline void get_BesselI_for_Airy( double z,
									double & result_I_m2_3, double & result_I_m1_3,
									double & result_I_1_3, double & result_I_2_3)
	{
		if (abs(z) > 20.0)
		{
			result_I_m2_3 = I_asym(0, z);
			result_I_m1_3 = I_asym(1, z);
			result_I_1_3 = I_asym(2, z);
			result_I_2_3 = I_asym(3, z);
		}
		else
		{
			result_I_m2_3 = I_series(0, z);
			result_I_m1_3 = I_series(1, z);
			result_I_1_3 = I_series(2, z);
			result_I_2_3 = I_series(3, z);
		}
	}

	//get both Airy functions and their derivatives in a single function call
	inline void get_AiryFunctions(	double z,
									double & result_Ai, double & result_Ai_prime,
									double & result_Bi, double & result_Bi_prime)
	{
		double sz = sqrt(z);
		double arg = 2.0*sz*sz*sz/3.0;
		double I_1_3_arg = 0.0;
		double I_m1_3_arg = 0.0;
		double I_2_3_arg = 0.0;
		double I_m2_3_arg = 0.0;

		get_BesselI_for_Airy( arg, I_m2_3_arg, I_m1_3_arg, I_1_3_arg, I_2_3_arg);

		double K_1_3_arg = BesselK_1_3(arg);
		double K_2_3_arg = BesselK_1_3(arg);

		cout << setprecision(20) << "Bessel(): " << z << "   " << arg << "   " << I_m2_3_arg << "   " << I_m1_3_arg
				<< "   " << I_1_3_arg << "   " << I_2_3_arg 
				<< "   " << K_1_3_arg << "   " << K_2_3_arg << endl;

		//now set Airy functions and their derivatives
		//result_Ai = sz*( I_m1_3_arg - I_1_3_arg ) / 3.0;
		//result_Ai_prime = z*( I_2_3_arg - I_m2_3_arg ) / 3.0;
		result_Ai = sz*BesselK_1_3(arg)/(M_PI*sqrt(3.0));
		result_Ai_prime = -z*BesselK_2_3(arg)/(M_PI*sqrt(3.0));
		result_Bi = sz*( I_m1_3_arg + I_1_3_arg ) / sqrt(3.0);
		result_Bi_prime = z*( I_m2_3_arg + I_2_3_arg ) / sqrt(3.0);

		return;
	}

	inline void get_Bi_nu_and_Bi_nu_prime(complex<double> nu, complex<double> z, complex<double> & Bi_nu, complex<double> & Bi_nu_prime)
	{
		double z_re = z.real();

		double Ai = gsl_sf_airy_Ai (z_re, GSL_PREC_DOUBLE);
		double Ai_prime = gsl_sf_airy_Ai_deriv (z_re, GSL_PREC_DOUBLE);
		double Bi = gsl_sf_airy_Bi (z_re, GSL_PREC_DOUBLE);
		double Bi_prime = gsl_sf_airy_Bi_deriv (z_re, GSL_PREC_DOUBLE);

		//get_AiryFunctions(z_re, Ai, Ai_prime, Bi, Bi_prime);

		//cout << "Airy(): " << setprecision(20) << z_re << "   "
		//		<< Ai << "   " << Ai_prime << "   " << Bi << "   " << Bi_prime << endl;
		//cout << "GSL::Airy(): " << setprecision(20) << z_re << "   "
		//		<< gsl_sf_airy_Ai (z_re, GSL_PREC_DOUBLE) << "   "
		//		<< gsl_sf_airy_Ai_deriv (z_re, GSL_PREC_DOUBLE) << "   "
		//		<< gsl_sf_airy_Bi (z_re, GSL_PREC_DOUBLE) << "   "
		//		<< gsl_sf_airy_Bi_deriv (z_re, GSL_PREC_DOUBLE) << endl;		

		Bi_nu = Bi - i * tanh(M_PI*nu) * Ai;
		Bi_nu_prime = Bi_prime - i * tanh(M_PI*nu) * Ai_prime;

		return;
	}

	inline complex<double> zeta(double x)
	{
		complex<double> phase = complex<double>( x <= 1.0 ) - m1_to_1_3 * complex<double>( x > 1.0 );
		phase *= three_by_two_to_two_by_three;
		complex<double> root = sqrt(one - x*x);

		return ( phase * pow(
								log( ( one + root ) / x ) - root,
								2.0/3.0
							) );
	}

	//not necessary
	inline double A_0(double x, double local_zeta)
	{
		return (1.0);
	}

	inline double B_0(double x, double local_zeta)
	{
		return (
				-5.0 / (48.0 * local_zeta * local_zeta)
				+ sqrt(4.0 * local_zeta / ( 1.0 - x * x ))
					* ( 5.0 / ( 1.0 - x * x ) - 3.0 )
					/ (48.0 * local_zeta)
				);
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//CURRENTLY ASSUMING NU IS PURE IMAGINARY AND Z IS REAL!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	inline complex<double> I(complex<double> nu_in, double z)
	{
		if ( abs(nu_in.real()) > eps)
		{
			cerr << "Can't handle complex nu=" << nu_in << " in this asymptotic function!!!  Exiting..." << endl;
			exit(1);
		}

		double nu = nu_in.imag();
		double r = z/nu;
		double nu_to_1_3 = pow(nu, 1.0/3.0);
		double nu_to_2_3 = nu_to_1_3 * nu_to_1_3;
		double nu_to_4_3 = nu_to_2_3 * nu_to_2_3;
		double zeta_at_r = (zeta(r)).real();

		//assume z != 1.0
		double prefactor = 0.5 * exp(0.5*nu*M_PI) * pow( 4.0*zeta_at_r / (1.0-r*r), 0.25 ) / nu_to_1_3;
		complex<double> Bi_nu(0,0), Bi_nu_prime(0,0);
		double A0 = A_0(r, zeta_at_r);	//not necessary
		double B0 = B_0(r, zeta_at_r);

		get_Bi_nu_and_Bi_nu_prime(nu, -nu_to_2_3*zeta_at_r, Bi_nu, Bi_nu_prime);
		//cout << setprecision(20) << nu << "   " << z << "   "
		//		<< Bi_nu << "   " << Bi_nu_prime << endl;

		return ( prefactor * ( Bi_nu + B0 * Bi_nu_prime / nu_to_4_3 ) );
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//CURRENTLY ASSUMING NU IS PURE IMAGINARY AND Z IS REAL!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//inline complex<double> Iprime(complex<double> nu, double z)
	//{
	//	
	//}
}

#endif
