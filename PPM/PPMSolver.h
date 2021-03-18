#pragma once

#include <vector>

using std::vector;

class PPMSolver
{
public:
	PPMSolver(vector<double> pressures, vector<double> density, double dx, double xmin,
		double xmax, double C, double t, int zoneCount);

	void solve();
	
	void addBoundaryConditions(int nleft, int nright);

	vector<double> flatten();

	vector<vector<double>> paraset(int nmin, int nmax, vector<double> dx, vector<double> xa);

	void parabola(int nmin, int nmax, const vector<vector<double>>& para, vector<double>& a,
		vector<double>& da, vector<double>& a6, vector<double>& al, const vector<double>& flat);

	void states(vector<double>& pl, vector<double>& ul, vector<double>& rl, 
		vector<double>& p6, vector<double>& u6, vector<double>& r6, 
		vector<double>& dp, vector<double>& du, vector<double>& dr, 
		vector<double>& plft, vector<double>& ulft, vector<double>& rlft, 
		vector<double>& prgh, vector<double>& urgh, vector<double>& rrgh,
		double dt);

	void riemann(int lmin, int lmax, double gamma, vector<double>& prgh, vector<double>& urgh, vector<double>& rrgh,
		vector<double>& plft, vector<double>& ulft, vector<double>& rlft,
		vector<double>& pmid, vector<double>& umid);

	void evolve(vector<double>& umid, vector<double>& pmid, double dt);

	void remap(vector<double>& flat);

private:
	double dt;

	vector<double> p, r, u, v, w, f;
	vector<double> e, q;


	double gam, gamm;

	double xmin, xmax, dxfac;

	int nmin, nmax, imax;

	vector<double> xa0; // Eulerian grid
	vector<double> xa;  // Lagrangian grid
	vector<double> dx0; // Eulerian spacing
	vector<double> dx;  // Lagrangian spacing
	vector<double> dvol;

	double svel, xvel, ridt;

	double smallp, smallr;
	double radius;
	double xwig;

	double targetTime;
	double courant;
};

