#include <iostream>
#include <fstream>
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

//const int particle_to_study
	// 1 - pion
	// 2 - proton
	// 3 - kaon
int particle_to_study;

const double hbarC = 197.33;
const double Cem = 2.0 / 3.0;	//my current best guess
const double k_infinity = 10.0;
const double xi_infinity = 5.0;
const int n_Dy = 51;

//const double tauC = 0.5;	//fm/c
const double DQ = 0.162035;	//fm (rough estimate!)
//const double vQ2 = DQ/tauQ;	//N.B. - must have tauQ > DQ for sub-luminal speed!
const double vQ2 = 1.0/3.0;
const double tauC = DQ/vQ2;
const double tauQ = tauC;	//for consistency with manuscript

long n_interp;

double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, sf;
double alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, chi_mu_mu, chi_T_mu, chi_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

const int n_xi_pts = 51;
const int n_k_pts = 50;	//# of k points should be even to avoid poles in 1F1, etc.!!!
const int n_tau_pts = 51;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * T_pts;
double * running_integral_array;

int main(int argc, char *argv[])
{
	//set parameters from command-line input
	// set particles to study...
	particle1.index = atoi(argv[1]);
	particle2.index = atoi(argv[2]);
	//Ti = atoi(argv[2]) / hbarC;		//initial trajectory temperature
	Ti = 250.0 / hbarC;
	fraction_of_evolution = 1.0;

	set_phase_diagram_and_EOS_parameters();

	//other constants
	double Delta_y_step = 0.1;

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
	cout 	<< "#########################################################" << endl
			<< "# Using following parameters:" << endl
			<< "# taui = " << taui << " fm/c, tauf = " << tauf << " fm/c" << endl
			<< "# si = " << si << ", sf = " << sf << endl
			<< "# Ti = " << Ti*hbarC << " MeV, Tf = " << Tf*hbarC << " MeV" << endl
			<< "# chi_{T,T} = " << chi_T_T << ", chi_{T,mu} = chi_{mu,T} = " << chi_T_mu << ", chi_{mu,mu} = " << chi_mu_mu << endl
			<< "#########################################################" << endl;

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

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_vs_tau();

	//get the ensemble averaged spectra
	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, &particle1);	//by definition of charge balance function (CBF)

	//vectorize calculations to make them run a little faster
	vector<complex<double> > Ftn_particle1_vec, Ftn_particle2_vec;
	vector<complex<double> > Ctnn_vec;

//	for (int ixi = 0; ixi < n_xi_pts; ++ixi)
//		cout << xi_pts_minf_inf[ixi] << "   " << Fn(xi_pts_minf_inf[ixi], &particle1) << "   " << Fn(xi_pts_minf_inf[ixi], &particle1) << endl;
//if (1) return (0);

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
		Ctnn_vec.push_back(Ctilde_n_n(k));
	}

	//start computing actual charge balance functions here
	for (int iDy = 0; iDy < n_Dy; iDy++)
	{
		double Delta_y = (double)iDy * Delta_y_step;

		complex<double> sum_lattice(0,0);
		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];

			complex<double> Ftn1 = Ftn_particle1_vec[ik];
			complex<double> Ftn2 = Ftn_particle2_vec[ik];
			complex<double> Ctnn = Ctnn_vec[ik];

			sum_lattice += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn );
			//cout << k << "   " << (Ftn1 * conj(Ftn2)).real() << "   " << Ctnn.real() << "   " << (Ftn1 * conj(Ftn2) * Ctnn).real() << endl;
		}
//if (1) return (0);

		complex<double> result_lattice = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum_lattice;
		cout << setprecision(15) << Delta_y << "   " << result_lattice.real() << "   " << result_lattice.imag() << endl;
	}

	return 0;
}
