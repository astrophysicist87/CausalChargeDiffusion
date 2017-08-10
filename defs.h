#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "gauss_quadrature.h"
#include "lib.h"

/*USAGE:
debugger(__LINE__, __FILE__);
*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

string truestring = "true";
string falsestring = "false";

bool white_noise = true;

inline string return_boolean_string(bool test){return (test ? truestring : falsestring);}

struct chosen_particle
{
	int index;
	double mass;
};

extern const double hbarC;
extern const double Cem;

extern const double tauQ, vQ2, DQ, tauC;

double fraction_of_evolution;

extern long n_interp;

extern double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, sf;
extern double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

extern double chi_T_T, chi_T_mu, chi_mu_mu;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;
extern double * T_pts;
extern double * running_integral_array;

chosen_particle particle1;
chosen_particle particle2;

int current_ik = -1;

extern double T0, mu0, Tc, Pc, nc, sc, wc, muc;
extern double A0, A2, A4, C0, B, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si;

//general functions
inline double incompleteGamma3(double x)
{
	return (exp(-x)*(2.+x*(2.+x)));
}

inline double incompleteGamma4(double x)
{
	return (exp(-x)*(6. + x*(6. + x*( 3. + x))));
}

inline double incompleteGamma5(double x)
{
	return (exp(-x)*(24.+x*(24.+x*(12.+x*(4.+x)))));
}

//equation of state and other thermodynamic relations
inline double P(double T)
{
	return(
		A4*pow(T, 4.0) - C0*T*T - B
	);
}

inline void set_phase_diagram_and_EOS_parameters()
{
	Nf = 3.0;											//number of massless flavors
	T0 = 170.0 / hbarC;											//phase transition curve, T scale
	mu0 = 1218.48 / hbarC;										//phase transition curve, mu scale
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;				//coeff in P(T) (i.e., EOS)
	A2 = Nf / 18.0;										//same
	A0 = Nf / (324.0 * M_PI * M_PI);					//same
	//C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );		//same
	C0 = 0.0;											//same
	B = 0.8 * pow(T0, 4.0);								//same

	return;
}

// functions to guess seed value for T,
// followed by functions to iteratively solve equations
// for T

inline double guess_T(double tau)
{
	return (Ti * pow(taui / tau, 1.0/3.0));
}

///////////////////////////////////////////////////////////
//two separate definitions of s, for convenience
///////////////////////////////////////////////////////////
inline double s_vs_tau(double tau)
{
	return (si * taui / tau);
}

inline double s_vs_T(double T)
{
	return (
		-2.0 * C0 * T + 4.0 * A4 * T*T*T
	);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

inline double w(double T)
{
	return (T * s_vs_T(T));
}

inline double chi_TT(double T)
{
	return ( 2.0 * ( -C0 + 6.0 * A4 * T * T ) );
}

inline double chi_Tmu(double T)
{
	T;
	return ( 0.0 );
}

inline double chi_mumu(double T)
{
	return ( 4.0 * A2 * T * T );
}

inline double norm_int(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	double mByT = (params->mass) / Tf;
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Fn(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	
	double c1 = 0.0;
	//double c2 = sf * chi_tilde_T_T;
	double c2 = chi_tilde_T_T;

	double mByT = (params->mass) / Tf;

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline complex<double> Ftilde_n(double k, void * p)
{
	return (integrate_1D_FT(Fn, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k, p));
}

inline complex<double> psi_plus(complex<double> k, complex<double> x)
{
	complex<double> lambda = sqrt(0.25 - vQ2*k*k);
	return (
			pow(x,lambda-0.5) * exp(-x)
			* Hypergeometric1F1(lambda+1.5, 2.0*lambda+1.0, x)
			);
}

inline complex<double> psi_minus(complex<double> k, complex<double> x)
{
	complex<double> mlambda = -sqrt(0.25 - vQ2*k*k);
	return (
			pow(x,mlambda-0.5) * exp(-x)
			* Hypergeometric1F1(mlambda+1.5, 2.0*mlambda+1.0, x)
			);
}

inline complex<double> psi_dot_plus(complex<double> k, complex<double> x)
{
	complex<double> lambda = sqrt(0.25 - vQ2*k*k);
	return (
			0.5 * (2.0*lambda-1.0) * pow(x,lambda-1.5) * exp(-0.5*x)
			* Hypergeometric0F1(lambda+1.0, x*x/16.0)
			);
}

inline complex<double> psi_dot_minus(complex<double> k, complex<double> x)
{
	complex<double> mlambda = -sqrt(0.25 - vQ2*k*k);
	return (
			0.5 * (2.0*mlambda-1.0) * pow(x,mlambda-1.5) * exp(-0.5*x)
			* Hypergeometric0F1(mlambda+1.0, x*x/16.0)
			);
}


//define some functions for subtracting off singularities from product of colored Green's functions
inline double g1(double x)
{
	double x2 = x*x;
	return (
				-x2/16.0 - x/2.0 - log(x)/8.0
			); 
}

inline double g2(double x)
{
	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x2*x2;
	double lnx = log(x);
	double lnx2 = lnx*lnx;
	return (
				x4/512.0 + x3/32.0 + x2/16.0 - x/4.0 + x2*lnx/128 + x*lnx/16.0 + lnx2/128.0
			);
}

inline double B1(double x, double xp)
{
	return (
				g1(xp) - g1(x) + (xp + 1.0) / 2.0
			);
}

inline double B2(double x, double xp)
{
	double g1xp = g1(xp);
	double Del_g1 = g1xp - g1(x);
	double Del_g2 = g2(xp) - g2(x);

	return (
			0.5 * Del_g1 * (2.0 * g1xp + xp + 1.0) - Del_g2
			};
}

inline GG_self_correlations(double k, double tau1p, double tau2p)
{
	double xf = tauf / tauQ;
	double x1p = tau1p / tauQ;
	double x2p = tau2p / tauQ;
	double prefactor = sqrt(tauf*tauf / (tau1p*tau2p)
						* exp( 0.5 * (x1p + x2p) - xf )
						/ (2.0*vQ2);
	double vQ = sqrt(vQ2);
	double arg = vQ * k * log(tau2p/tau1p);

	double B1_xf_x1p = B1(xf, x1p);
	double B1_xf_x2p = B1(xf, x2p);
	double B2_xf_x1p = B2(xf, x1p);
	double B2_xf_x2p = B2(xf, x2p);

	double C1 = vQ2*k*k + B1_xf_x1p * B1_xf_x2p - B2_xf_x1p - B2_xf_x2p;
	double C2 = vQ * k * (B1_xf_x1p - B1_xf_x2p);

	return (
				prefactor * ( C1 * cos(arg) + C2 * sin(arg) )
			);
}

inline double mediumN(double tau_loc)
{
	double T_loc = interpolate1D(tau_pts, T_pts, tau_loc, n_tau_pts, 0, false);
	double s_loc = s_vs_T(T_loc);
	double chi_Q = chi_mumu(T_loc);
	double numerator = 2.0 * DQ * chi_Q * T_loc;
	double denominator = tau_loc * s_loc * s_loc;	//assume transverse area A already accounted for...
	return ( numerator / denominator );
}

inline complex<double> Gtilde_n_white(double k, double tau, double taup)
{
	double arg = DQ * k * k * ((1.0/tau) - (1.0/taup));	//since tau==tauf always in this calculation,
														//exp(arg) will always be Gaussian in k
	return ( i * k * exp(arg) ); //dimensionless
}

inline complex<double> Gtilde_n_color(double k, double tau, double taup)
{
	complex<double> psi_plus_at_tau = psi_plus(k, tau / tauQ);
	complex<double> psi_minus_at_tau = psi_minus(k, tau / tauQ);
	complex<double> psi_plus_at_taup = psi_plus(k, taup / tauQ);
	complex<double> psi_minus_at_taup = psi_minus(k, taup / tauQ);
	complex<double> psi_dot_plus_at_taup = psi_dot_plus(k, taup / tauQ);
	complex<double> psi_dot_minus_at_taup = psi_dot_minus(k, taup / tauQ);

	complex<double> numerator = psi_plus_at_tau * psi_dot_minus_at_taup - psi_minus_at_tau * psi_dot_plus_at_taup;
	complex<double> denominator = psi_plus_at_taup * psi_dot_minus_at_taup - psi_minus_at_taup * psi_dot_plus_at_taup;

	return ( i * k * numerator / denominator ); //dimensionless
}

inline complex<double> tau_integration(complex<double> (*Gtilde_X)(double, double, double), complex<double> (*Gtilde_Y)(double, double, double), double k)
{
	complex<double> result(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double tau_loc = tau_pts[it];
		double T_loc = T_pts[it];
		double s_loc = s_vs_T(T_loc);
		double chi_Q = chi_mumu(T_loc);
//cout << tau_loc << "   " << 2.0 * DQ * chi_Q * T_loc * tau_loc << endl;
		complex<double> tmp_result = tau_wts[it] * ( 2.0 * DQ * chi_Q * T_loc / tau_loc )
										* (*Gtilde_X)(k, tauf, tau_loc) * (*Gtilde_Y)(-k, tauf, tau_loc);
		result += tmp_result;
	}

//if (1) exit(0);
//cout << "CHECK: " << k << "   " << result.real() << "   " << result.real() - ( chi_mumu(Tf) * Tf * tauf ) << endl;

	return ( ( chi_mumu(Tf) * Tf * tauf ) - result );
}

inline void set_running_transport_integral(double * run_int_array)
{
	const int n_x_pts = 11;
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	run_int_array[0] = 0.0;
	for (int it = 0; it < n_tau_pts-1; ++it)
	{
		double sum = 0.0;
		double t0 = tau_pts[it], t1 = tau_pts[it+1];
		double cen = 0.5*(t0+t1);
		double hw = 0.5*(t1-t0);
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double t_loc = cen + hw * x_pts[ix];
			sum += x_wts[ix] * hw * exp(t_loc / tauQ) * mediumN(t_loc);
		}
		run_int_array[it+1] = exp(-t1 / tauQ) * ( tauQ * exp(t0 / tauQ) * run_int_array[it] + sum ) / tauQ;	//this array contains eta(x) (defined on my whiteboard)
	}

	delete [] x_pts;
	delete [] x_wts;

	return;
}

inline complex<double> colored_tau_integration(
					complex<double> (*Gtilde_X)(double, double, double),
					complex<double> (*Gtilde_Y)(double, double, double),
					double k)
{
	complex<double> locsum(0,0);
	
	const int n_x_pts = 201;	//try this
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	double delta_tau_lower = -10.0 * tauQ, delta_tau_upper = 10.0 * tauQ;	//this bounds the interval where the integrand is large-ish
	for (int itp = 0; itp < n_tau_pts; ++itp)
	{
		double tX_loc = tau_pts[itp];
		double TX_loc = T_pts[itp];
		double sX_loc = s_vs_T(TX_loc);
		complex<double> factor_X = sX_loc * (*Gtilde_X)(k, tauf, tX_loc);		//extra factor of entropy!!!

		double tau_lower = max(taui, tX_loc + delta_tau_lower);			//if lower limit goes before beginning of lifetime, just start at tau0
		double tau_upper = min(tauf, tX_loc + delta_tau_upper);			//if upper limit goes past end of lifetime, just end at tauf
		double hw_loc = 0.5 * (tau_upper - tau_lower);
		double cen_loc = 0.5 * (tau_upper + tau_lower);

		complex<double> sum_X = 0.0;
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double tY_loc = cen_loc + hw_loc * x_pts[ix];

			//check if we want to subtract self-correlations
			double GG_self_corr = 0.0;
			if (subtract_self_correlations)
				GG_self_corr = GG_self_correlations(k, tX_loc, tY_loc)

			double TY_loc = interpolate1D(tau_pts, T_pts, tY_loc, n_tau_pts, 0, false, 2);
			double sY_loc = s_vs_T(TY_loc);
			complex<double> factor_Y = sY_loc * (*Gtilde_Y)(-k, tauf, tY_loc);	//extra factor of entropy!!!

			double min_tp_tpp = min(tX_loc, tY_loc);
			double eta_at_min_tp_tpp = interpolate1D(tau_pts, running_integral_array, min_tp_tpp, n_tau_pts, 0, false, 2);

			double sum_XY = exp(-abs(tX_loc - tY_loc) / tauQ) * eta_at_min_tp_tpp / (2.0*tauQ);
			sum_X += hw_loc * x_wts[ix] * ( factor_X * factor_Y - GG_self_corr ) * sum_XY;
		}
		locsum += tau_wts[itp] * sum_X;
	}

	delete [] x_pts;
	delete [] x_wts;

	return (locsum);
}

inline complex<double> Ctilde_n_n(double k)
{
	complex<double> sum(0,0);
	if ( white_noise )
		sum = tau_integration(Gtilde_n_white, Gtilde_n_white, k);
	else
		sum = colored_tau_integration(Gtilde_n_color, Gtilde_n_color, k);

	return ( sum / (tauf*tauf) );	//fm^2; note that we must include extra factor of tauf^-2 for consistency with manuscript
}

////////////////////////////////////////////////////////////////////////////////
// Functions/stuff to solve for time-dependence of T and/or Tf
////////////////////////////////////////////////////////////////////////////////

struct rparams
{
	double tau;
};

//compute final time-step Tf and muf

int input_get_Tf_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	const double x0 = gsl_vector_get (x, 0);	//T

	const double y0 = P(x0);						//defines P==0 curve

	gsl_vector_set (f, 0, y0);

	return GSL_SUCCESS;
}

void compute_Tf()
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t iter = 0;

	const size_t n = 1;
	struct rparams p = {taui};
	gsl_multiroot_function f = {&input_get_Tf_f, n, &p};

	double x_init[n] = {Ti};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);

	gsl_T = gsl_multiroot_fsolver_hybrids;
	gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
	gsl_multiroot_fsolver_set (gsl_s, &f, x);

	do
	{
		iter++;

		status = gsl_multiroot_fsolver_iterate (gsl_s);

		if (status)   /* check if solver is stuck */
			break;

		status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	//finally, store results
	Tf = gsl_vector_get (gsl_s->x, 0);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);

	return;
}

//compute T at each time step

int input_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T

	const double y0 = s_vs_T(x0) - s_vs_tau(tau_local);

	gsl_vector_set (f, 0, y0);

	return GSL_SUCCESS;
}

void populate_T_vs_tau()
{
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t iter = 0;

		const size_t n = 1;
		struct rparams p = {tau_pts[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[n] = {guess_T(tau_pts[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//finally, store results
		T_pts[it] = gsl_vector_get (gsl_s->x, 0);

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

double compute_newTf(double new_tauf)
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t iter = 0;

	const size_t n = 1;
	struct rparams p = {new_tauf};
	gsl_multiroot_function f = {&input_f, n, &p};

	double x_init[n] = {guess_T(new_tauf)};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);

	gsl_T = gsl_multiroot_fsolver_hybrids;
	gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
	gsl_multiroot_fsolver_set (gsl_s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (gsl_s);

		if (status)   /* check if solver is stuck */
			break;

		status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	//finally, store result
	double result = gsl_vector_get (gsl_s->x, 0);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);
	
	return (result);
}

// End of file

#endif
