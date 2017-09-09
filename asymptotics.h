#ifndef ASYMPTOTICS_H
#define ASYMPTOTICS_H

#include <complex>
#include <iostream>
#include <cmath>

using namespace std;

#include "./special_function_library/cmlib/complex.h"
#include "./special_function_library/cmlib/protom.h"

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
	//inline complex<double> I_asym(int nu_index, complex<double> z)
	inline double I_asym(int nu_index, double z)
	{
		double nu = nus[nu_index];

		//complex<double> sum1 = 0.0, sum2 = 0.0;
		//complex<double> oldsum1 = 0.0, oldsum2 = 0.0;
		//complex<double> next_term = 1.0;
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
//cout << "I_asym(): " << k-1 << "   " << oldsum1 << "   " << sum1 << "   " << oldsum2 << "   " << sum2 << endl;
		} while (
					2.0*abs(sum1-oldsum1)/abs(sum1+oldsum1) >= 1e-20
					|| 2.0*abs(sum2-oldsum2)/abs(sum2+oldsum2) >= 1e-20
					);

		//complex<double> factor1 = exp(z) / sqrt(2.0*M_PI*z);
		//complex<double> factor2 = exp(-z) / sqrt(2.0*M_PI*z);
		double factor1 = exp(z) / sqrt(2.0*M_PI*z);
		double factor2 = exp(-z) / sqrt(2.0*M_PI*z);
		double phz = 0.0;
		if (-0.5*M_PI < phz && phz < 1.5*M_PI)
			sign = 1.0;
		else if (-0.5*M_PI < -phz && -phz < 1.5*M_PI)
			sign = -1.0;
		else
			cerr << "Error: phz = " << phz << endl;

		//return ( factor1*sum1 + sign*i*exp(sign*nu*M_PI*i)*factor2*sum2 );
		return ( factor1*sum1 );
	}

	///////////////////////////////////////////////////////////////
	//inline complex<double> I_series(int nu_index, complex<double> z)
	inline double I_series(int nu_index, double z)
	{
		double nu = nus[nu_index];
		double gamma_nu = gamma_nus[nu_index];

		//complex<double> sum = 0.0, oldsum = 0.0;
		double sum = 0.0, oldsum = 0.0;
		double sign = 1.0;
		//complex<double> next_term = 1.0/(nu*gamma_nu);
		double next_term = 1.0/(nu*gamma_nu);
		//complex<double> arg = 0.25*z*z;
		double arg = 0.25*z*z;

		double k = 0.0;
		do
		{
			oldsum = sum;
			sum += next_term;

			//complex<double> new_factor = arg / ( (k+1.0) * (k+1.0+nu) );
			double new_factor = arg / ( (k+1.0) * (k+1.0+nu) );
			next_term *= new_factor;
			k++;
		} while ( 2.0*abs(sum-oldsum)/abs(sum+oldsum) >= 1e-20 );

		return ( pow(0.5*z, nu)*sum );
	}

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	//inline void get_BesselI_for_Airy( complex<double> z, 
	inline void get_BesselI_for_Airy( double z,
									complex<double> & result_I_m2_3, complex<double> & result_I_m1_3,
									complex<double> & result_I_1_3, complex<double> & result_I_2_3)
	{
		if (abs(z) > 10.0)
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
	//inline void get_AiryFunctions(	complex<double> z,
	inline void get_AiryFunctions(	double z,
									complex<double> & result_Ai, complex<double> & result_Ai_prime,
									complex<double> & result_Bi, complex<double> & result_Bi_prime)
	{
		//complex<double> sz = sqrt(z);
		double sz = sqrt(z);
		//complex<double> arg = 2.0*sz*sz*sz/3.0;
		//double arg = (2.0*sz*sz*sz/3.0).real();
		double arg = 2.0*sz*sz*sz/3.0;
		complex<double> I_1_3_arg = 0.0;
		complex<double> I_m1_3_arg = 0.0;
		complex<double> I_2_3_arg = 0.0;
		complex<double> I_m2_3_arg = 0.0;

		get_BesselI_for_Airy( arg, I_m2_3_arg, I_m1_3_arg, I_1_3_arg, I_2_3_arg);

//cout << setprecision(20) << arg << "   " << I_m2_3_arg << "   " << I_m1_3_arg << "   " << I_1_3_arg << "   " << I_2_3_arg << endl;

		//now set Airy functions and their derivatives
		result_Ai = sz*( I_m1_3_arg - I_1_3_arg ) / 3.0;
		result_Ai_prime = z*( I_2_3_arg - I_m2_3_arg ) / 3.0;
		result_Bi = sz*( I_m1_3_arg + I_1_3_arg ) / sqrt(3.0);
		result_Bi_prime = z*( I_m2_3_arg + I_2_3_arg ) / sqrt(3.0);

		return;
	}

	inline void get_Bi_nu_and_Bi_nu_prime(complex<double> nu, complex<double> z, complex<double> & Bi_nu, complex<double> & Bi_nu_prime)
	{
		complex<double> Ai(0,0), Ai_prime(0,0), Bi(0,0), Bi_prime(0,0);

		double z_re = z.real();

		get_AiryFunctions(z_re, Ai, Ai_prime, Bi, Bi_prime);

cout << "Airy(): " << setprecision(20) << z_re << "   " << Ai << "   " << Ai_prime << "   " << Bi << "   " << Bi_prime << endl;

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
	//inline complex<double> A_0(double x, complex<double> local_zeta)
	inline double A_0(double x, double local_zeta)
	{
		return (1.0);
	}

	//inline complex<double> B_0(double x, complex<double> local_zeta)
	inline double B_0(double x, double local_zeta)
	{
		return (
				//-5.0 / (48.0 * local_zeta * local_zeta)
				//+ sqrt(4.0 * local_zeta / ( one - x * x ))
				//	* ( 5.0 / ( one - x * x ) - 3.0 )
				//	/ (48.0 * local_zeta)
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
		//complex<double> zeta_at_r = zeta(r);
		double zeta_at_r = (zeta(r)).real();

		//assume z != 1.0
		//complex<double> prefactor = 0.5 * exp(0.5*nu*M_PI) * pow( 4.0*zeta_at_r / (one-r*r), 0.25 ) / nu_to_1_3;
		double prefactor = 0.5 * exp(0.5*nu*M_PI) * pow( 4.0*zeta_at_r / (1.0-r*r), 0.25 ) / nu_to_1_3;
		complex<double> Bi_nu(0,0), Bi_nu_prime(0,0);
		//complex<double> A0 = A_0(r, zeta_at_r);	//not necessary
		//complex<double> B0 = B_0(r, zeta_at_r);
		double A0 = A_0(r, zeta_at_r);	//not necessary
		double B0 = B_0(r, zeta_at_r);

		get_Bi_nu_and_Bi_nu_prime(nu, -nu_to_2_3*zeta_at_r, Bi_nu, Bi_nu_prime);
cout << setprecision(20) << nu << "   " << z << "   " << Bi_nu << "   " << Bi_nu_prime << endl;

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
