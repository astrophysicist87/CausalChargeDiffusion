#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

using namespace std;

#include "defs.h"

//const int particle_to_study
	// 1 - pion
	// 2 - proton
	// 3 - kaon
int particle_to_study;

bool print_dndn_k = true;
bool print_dndn_Dxi = true;
//white noise is default
bool white_noise = true;
bool white_Green = true;

const double hbarC = 197.33;
const double xi_infinity = 5.0;
const double k_infinity = 20.0;

const int n_Dy = 501;
double Delta_y_step = 0.01;

const double DQ = 0.162035;	//fm (rough estimate!)
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
const int n_k_pts = 5000;	//# of k points should be even to avoid poles in 1F1, etc.!!!
const int n_tau_pts = 201;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * T_pts;
double * running_integral_array;

///////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	//set parameters from command-line input
	// set particles to study...
	particle1.index = atoi(argv[1]);
	particle2.index = atoi(argv[2]);
	//Ti = atoi(argv[2]) / hbarC;		//initial trajectory temperature
	Ti = 350.0 / hbarC;
	fraction_of_evolution = 1.0;

	//set speed of sound and correlation timescale
	vQ2 = atof(argv[3]);
	tauQ = DQ/vQ2;

	//update white noise and white Green's functions boolean variables
	//assume both are the same for now...
	white_noise = (vQ2 > 20.0);
	white_Green = (vQ2 > 20.0);

	/*const int n_x_pts = 101;
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 6, 0.0, 0.0, 0.0, 1.0, x_pts, x_wts);
	double test_sum = 0.0;
	for (int ix = 0; ix < n_x_pts; ++ix)
		test_sum += x_wts[ix] * exp(-x_pts[ix]*x_pts[ix]);	//should integrate a Gaussian???
	cout << "integral = " << test_sum << endl;
	if (1) exit (0);*/

	set_phase_diagram_and_EOS_parameters();

	switch (particle1.index)
	{
		case 1:	//pion
			particle1.mass = 139.57 / hbarC;
			particle1.name = "pi";
			break;
		case 2:	//proton
			particle1.mass = 939.0 / hbarC;
			particle1.name = "p";
			break;
		case 3:	//kaon
			particle1.mass = 493.68 / hbarC;
			particle1.name = "K";
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
			particle2.name = "pi";
			break;
		case 2:	//proton
			particle2.mass = 939.0 / hbarC;
			particle2.name = "p";
			break;
		case 3:	//kaon
			particle2.mass = 493.68 / hbarC;
			particle2.name = "K";
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}

	taui = 0.5;		//fm/c

	si = s_vs_T(Ti);

	//fixing Tf explicitly instead of
	//calculating it from P=0 curve
	Tf = 150.0 / hbarC;

	sf = s_vs_T(Tf);
	tauf = si * taui / sf;

	// added this to allow "snapshots" throughout evolution
	tauf = taui + fraction_of_evolution * (tauf - taui);
	sf = s_vs_tau(tauf);
	//Tf = compute_newTf(tauf);
	Tf = guess_T(tauf);
	// rest of code runs the same as before

    // initialize other parameters
    ds = 2.0;

	string filename1, filename2, filename3;
	set_outfilenames(filename1, filename2, filename3);

	//set the susceptibilities
	chi_mu_mu = chi_mumu(Tf);
	chi_T_mu = chi_Tmu(Tf);
	chi_T_T = chi_TT(Tf);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	//output some parameters from calculation
	ofstream output_results;
	output_results.open(filename1.c_str());

	output_results 
			<< "#########################################################" << endl
			<< "# Using following parameters:" << endl
			<< "# taui = " << taui << " fm/c, tauf = " << tauf << " fm/c" << endl
			<< "# si = " << si << ", sf = " << sf << endl
			<< "# Ti = " << Ti*hbarC << " MeV, Tf = " << Tf*hbarC << " MeV" << endl
			<< "# v_Q^2 = " << vQ2 << ", D_Q = " << DQ << " fm/c, tau_Q = " << tauQ << " fm/c" << endl
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
    //tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);
	double inverse_smearing_width = 2.0;
    tmp = gauss_quadrature(n_k_pts, 6, 0.0, 0.0, 0.0, inverse_smearing_width, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);
    //tmp = gauss_quadrature(n_integ_besselK_points, 1, 0.0, 0.0, -1.0, 1.0, x_integ_besselK_pts, x_integ_besselK_wts);

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T for remainder of calculation
	for (int it = 0; it < n_tau_pts; ++it)
		T_pts[it] = guess_T(tau_pts[it]);

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
	//////////////////////////////////////////////////////////

	//get the ensemble averaged spectra
	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, &particle1);	//by definition of charge balance function (CBF)

	//vectorize calculations to make them run a little faster
	vector<complex<double> > Ftn_particle1_vec, Ftn_particle2_vec;
	vector<complex<double> > Ctnn_vec, Ctnn_no_SC_vec, SC_vec;

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


	ofstream output_dndn_k, output_dndn_Dxi;	//files for Fourier output, coordinate-space output
	string dndn_k_filename = filename2;
	string dndn_Dxi_filename = filename3;
	if (print_dndn_k)
		output_dndn_k.open(dndn_k_filename.c_str());
	if (print_dndn_Dxi)
		output_dndn_Dxi.open(dndn_Dxi_filename.c_str());

	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		current_ik = ik;
		vector<complex<double> > results(3);
		Ctilde_n_n(k, &results);
		Ctnn_vec.push_back(results[0]);
		Ctnn_no_SC_vec.push_back(results[1]);
		SC_vec.push_back(results[2]);
		if (print_dndn_k)
			output_dndn_k << k << "   " << setprecision(15)
							<< results[0].real() << "   "
							<< results[1].real() << "   "
							<< results[2].real() << endl;
	}

	//start computing actual charge balance functions here
	for (int iDy = 0; iDy < n_Dy; iDy++)
	{
		double Delta_y = (double)iDy * Delta_y_step;
		double Delta_xi = Delta_y;

		complex<double> sum(0,0), sum_no_SC(0,0);
		complex<double> sum_Dxi(0,0), sum_Dxi_no_SC(0,0), SC_Dxi(0,0);

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];

			complex<double> Ftn1 = Ftn_particle1_vec[ik];
			complex<double> Ftn2 = Ftn_particle2_vec[ik];
			complex<double> Ctnn = Ctnn_vec[ik];
			complex<double> Ctnn_no_SC = Ctnn_no_SC_vec[ik];
			complex<double> SC_loc = SC_vec[ik];

			sum += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn );
			sum_no_SC += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn1 * conj(Ftn2) * Ctnn_no_SC );

			sum_Dxi += k_wts[ik] * exp(-0.01*inverse_smearing_width*inverse_smearing_width*k*k) * exp(i * k * Delta_xi) * Ctnn;
			sum_Dxi_no_SC += k_wts[ik] * exp(-0.01*inverse_smearing_width*inverse_smearing_width*k*k) * exp(i * k * Delta_xi) * Ctnn_no_SC;
			SC_Dxi += k_wts[ik] * exp(-0.01*inverse_smearing_width*inverse_smearing_width*k*k) * exp(i * k * Delta_xi) * SC_loc;

			//cerr << k << "   " << Ctnn.real() << "   " << Ctnn_no_SC.real() 
			//		<< "   " << (Ftn1 * conj(Ftn2)).real()
			//		<< "   " << (Ftn1 * conj(Ftn2) * Ctnn).real()
			//		<< "   " << (Ftn1 * conj(Ftn2) * Ctnn_no_SC).real() << endl;
		}

		//omit smearing function factors to just F.T. back to xi-space
		if (print_dndn_Dxi)
			output_dndn_Dxi << Delta_xi << "   " << setprecision(15)
							<< sum_Dxi.real() << "   "
							<< sum_Dxi_no_SC.real() << "   "
							<< SC_Dxi.real() << endl;

		complex<double> result = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum;
		complex<double> result_no_SC = (ds*tauf*Tf / (4.0*M_PI*M_PI * norm)) * sum_no_SC;

		output_results
				<< setprecision(15) << Delta_y << "   "
				<< result.real() << "   " << result.imag() << "   "
				<< result_no_SC.real() << "   " << result_no_SC.imag() << endl;
	}

	if (print_dndn_k)
		output_dndn_k.close();
	if (print_dndn_Dxi)
		output_dndn_Dxi.close();
	output_results.close();

	return 0;
}
