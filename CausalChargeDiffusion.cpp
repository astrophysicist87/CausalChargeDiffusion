#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "defs.h"
#include "asymptotics.h"

//const int particle_to_study
	// 1 - pion
	// 2 - proton
	// 3 - kaon
int particle_to_study;

const double hbarC = 197.33;
const double xi_infinity = 5.0;
const double k_infinity = 50.0;
//const double k_critical = 0.5 / sqrt(vQ2);
const int n_Dy = 501;

//const double tauC = 0.5;	//fm/c
const double DQ = 0.162035;	//fm (rough estimate!)
//const double vQ2 = DQ/tauQ;	//N.B. - must have tauQ > DQ for sub-luminal speed!
//const double vQ2 = 10.0;
//const double tauQ = DQ/vQ2;
double vQ2, tauQ;

long n_interp;

double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, sf;
double alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, chi_mu_mu, chi_T_mu, chi_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

const int n_xi_pts = 5000;
const int n_k_pts = 500;	//# of k points should be even to avoid poles in 1F1, etc.!!!
const int n_tau_pts = 51;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * T_pts;
double * running_integral_array;

//const int n_integ_besselK_points = 101;
//vector<double> x_integ_besselK_pts(n_integ_besselK_points), x_integ_besselK_wts(n_integ_besselK_points);

///////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	//set parameters from command-line input
	// set particles to study...
	particle1.index = atoi(argv[1]);
	particle2.index = atoi(argv[2]);
	//Ti = atoi(argv[2]) / hbarC;		//initial trajectory temperature
	Ti = 250.0 / hbarC;
	fraction_of_evolution = 1.0;

	//set speed of sound and correlation timescale
	vQ2 = atof(argv[3]);
	tauQ = DQ/vQ2;

	set_phase_diagram_and_EOS_parameters();

	//other constants
	double Delta_y_step = 0.01;

	switch (particle1.index)
	{
		case 1:	//pion
			particle1.mass = 139.57 / hbarC;
			break;
		case 2:	//proton
			particle1.mass = 939.0 / hbarC;
			break;
		case 3:	//kaon
			particle1.mass = 493.68 / hbarC;
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}
	switch (particle2.index)
	{
		case 1:	//pion
			particle2.mass = 139.57 / hbarC;
			break;
		case 2:	//proton
			particle2.mass = 939.0 / hbarC;
			break;
		case 3:	//kaon
			particle2.mass = 493.68 / hbarC;
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}

	taui = 0.5;		//fm/c

	si = s_vs_T(Ti);

	compute_Tf();

	sf = s_vs_T(Tf);
	tauf = si * taui / sf;

	// added this to allow "snapshots" throughout evolution
	tauf = taui + fraction_of_evolution * (tauf - taui);
	sf = s_vs_tau(tauf);
	Tf = compute_newTf(tauf);
	// rest of code runs the same as before

    // initialize other parameters
    ds = 2.0;

	//set the susceptibilities
	chi_mu_mu = chi_mumu(Tf);
	chi_T_mu = chi_Tmu(Tf);
	chi_T_T = chi_TT(Tf);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	//output some parameters from calculation
	/*cout 	<< "#########################################################" << endl
			<< "# Using following parameters:" << endl
			<< "# taui = " << taui << " fm/c, tauf = " << tauf << " fm/c" << endl
			<< "# si = " << si << ", sf = " << sf << endl
			<< "# Ti = " << Ti*hbarC << " MeV, Tf = " << Tf*hbarC << " MeV" << endl
			<< "# v_Q^2 = " << vQ2 << ", D_Q = " << DQ << " fm/c, tau_Q = " << tauQ << " fm/c" << endl
			<< "# chi_{T,T} = " << chi_T_T << ", chi_{T,mu} = chi_{mu,T} = " << chi_T_mu << ", chi_{mu,mu} = " << chi_mu_mu << endl
			<< "#########################################################" << endl;*/

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);
    //tmp = gauss_quadrature(n_integ_besselK_points, 1, 0.0, 0.0, -1.0, 1.0, x_integ_besselK_pts, x_integ_besselK_wts);

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_vs_tau();

	//////////////////////////////////////////////////////////
	//Check hypergeometric function implementations
	/*
	for (int ik = 0; ik < n_k_pts; ++ik)
	for (int it = 0; it < n_tau_pts; ++it)
	{
		complex<double> k = k_pts[ik];
		complex<double> x = tau_pts[it] / tauQ;
		complex<double> lambda = sqrt(0.25 - vQ2*k*k);
		complex<double> mlambda = -sqrt(0.25 - vQ2*k*k);

		complex<double> result1 = SFL_Hypergeometric1F1(lambda+0.5, 2.0*lambda+1.0, x);
		complex<double> result2 = SFL_Hypergeometric1F1(mlambda+0.5, 2.0*mlambda+1.0, x);
		complex<double> result3 = SFL_Hypergeometric1F1(lambda+1.5, 2.0*lambda+1.0, x);
		complex<double> result4 = SFL_Hypergeometric1F1(mlambda+1.5, 2.0*mlambda+1.0, x);
	}
	if (1) return (0);
	*/

	/*complex<double> nuBIG = 10.0*i;
	complex<double> nuSMALL = 0.05*i;
	long double z = 30.0;
	complex<double> result1 = asymptotics::I(nuBIG, z);
	cout << setprecision(20) << z << "   " << result1.real() << "   " << result1.imag() << endl;
	complex<double> result2 = asymptotics::Iprime(nuBIG, z);
	cout << setprecision(20) << z << "   " << result2.real() << "   " << result2.imag() << endl;
	cout << setprecision(20) << asymptotics::zeta_prime(0.5) << "   "
			<< asymptotics::zeta_prime(2.0) << endl;
	cout << asymptotics::B_0_prime(0.5, asymptotics::zeta(0.5), asymptotics::zeta_prime(0.5)) << "   "
			<< asymptotics::B_0_prime(2.0, asymptotics::zeta(2.0), asymptotics::zeta_prime(2.0)) << endl;
	for (int iz = 0; iz <= 1000; ++iz)
	{
		long double z = 1.0 + 0.06 * iz;
		complex<double> result1 = asymptotics::I(nuBIG, z);
		complex<double> result2 = asymptotics::I(nuSMALL, z);
		cout << z << "   " << result1.real() << "   " << result1.imag()
				<< "   " << result2.real() << "   " << result2.imag() << endl;
	}
	if (1) return (0);*/

	const double k_critical = 0.5 / sqrt(vQ2);
	/*cerr << "k_c: " << vQ2 << "   " << k_critical << endl;*/
	/*for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		if (k*k < k_critical*k_critical)
			continue;
		double tau1 = 0.5*tauf;
		cout << k << "   " << Gtilde_n_white(k, tauf, tau1).imag() << "   " << Gtilde_n_color(k, tauf, tau1).imag()
				<< "   " << new_asymptotic_Gtilde_n_color(k, tauf, tau1).imag() << endl;
	//if (ik >= 2) return (0);
	}
	if (1) return (0);*/
	
	//////////////////////////////////////////////////////////

	//get the ensemble averaged spectra
	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, &particle1);	//by definition of charge balance function (CBF)

	//vectorize calculations to make them run a little faster
	vector<complex<double> > Ftn_particle1_vec, Ftn_particle2_vec;
	vector<complex<double> > Ctnn_vec, Ctnn_no_SC_vec;

	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		current_ik = ik;
		Ftn_particle1_vec.push_back(Ftilde_n(k, &particle1));
		Ftn_particle2_vec.push_back(Ftilde_n(k, &particle2));
	}

	//compute running integral of transport/medium effects
	running_integral_array = new double [n_tau_pts];
	set_running_transport_integral(running_integral_array);

	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		current_ik = ik;
		vector<complex<double> > results(2);
		Ctilde_n_n(k, &results);
		Ctnn_vec.push_back(results[0]);
		Ctnn_no_SC_vec.push_back(results[1]);
	}

	//start computing actual charge balance functions here
	for (int iDy = 0; iDy < n_Dy; iDy++)
	{
		double Delta_y = (double)iDy * Delta_y_step;

		complex<double> sum(0,0), sum_no_SC(0,0);
		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];

			complex<double> Ftn1 = Ftn_particle1_vec[ik];
			complex<double> Ftn2 = Ftn_particle2_vec[ik];
			complex<double> Ctnn = Ctnn_vec[ik];
			complex<double> Ctnn_no_SC = Ctnn_no_SC_vec[ik];

			sum += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn );
			sum_no_SC += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn_no_SC );
			cerr << k << "   " << Ctnn.real() << "   " << Ctnn_no_SC.real() 
					<< "   " << (Ftn1 * conj(Ftn2)).real()
					<< "   " << (Ftn1 * conj(Ftn2) * Ctnn).real()
					<< "   " << (Ftn1 * conj(Ftn2) * Ctnn_no_SC).real() << endl;
		}
		//if (1) return (0);

		complex<double> result = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum;
		complex<double> result_no_SC = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum_no_SC;
		cout << setprecision(15) << Delta_y << "   " << result.real() << "   " << result.imag() << "   " << result_no_SC.real() << "   " << result_no_SC.imag() << endl;
	}

	return 0;
}
