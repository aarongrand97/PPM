#include "PPMSolver.h"

using std::min;
using std::max;

PPMSolver::PPMSolver(vector<double> pressures, vector<double> density, double dx, double xmin,
	double xmax, double C, double t, int zoneCount)
{

	gam = 1.4;
	gamm = gam - 1.0;

	smallp = 1.0e-15;
	smallr = 1.0e-15;

	radius = 1.0;
	xwig = 0.0;

	nmin = 6;
	nmax = zoneCount + 5;
	imax = zoneCount;


	this->xmin = xmin;
	this->xmax = xmax;
	dxfac = (xmax - xmin) / double(imax);

	xa0 = vector<double>(imax + 12);
	dx0 = vector<double>(imax + 12);
	dvol = vector<double>(imax + 12);

	for (int n = nmin; n <= nmax; n++) {
		xa0[n] = xmin + ((double)n - nmin) * dxfac;
		dx0[n] = dxfac;
	}

	// create the initial conditions
	p = vector<double>(imax + 12);
	r = vector<double>(imax + 12);
	u = vector<double>(imax + 12);
	v = vector<double>(imax + 12);
	w = vector<double>(imax + 12);
	f = vector<double>(imax + 12);
	
	int nmid = imax / 2 + 6;
	for (int n = nmin; n < nmid; n++) {
		p[n] = 1.0;
		r[n] = 1.0;
	}
	for (int n = nmid; n <= nmax; n++) {
		p[n] = 0.1;
		r[n] = 0.125;
	}

	for (int n = nmin; n <= nmax; n++) {
		svel = sqrt(gam * p[n] / r[n]) / dx0[n];
		xvel = abs(u[n]) / dx0[n];
		ridt = max(xvel, max(svel, ridt));
	}

	targetTime = t;
	courant = C;

	dt = courant / ridt;

	e = vector<double>(imax + 12);
	q = vector<double>(imax + 12);

	xa = vector<double>(imax + 12);
	this->dx = vector<double>(imax + 12);

	
}

void PPMSolver::solve() {
	double t0 = 0;
	int iter = 0;

	int ncycle = 0, ncycend = 80;


	while(ncycle < ncycend){
		// NEED TO GET THE TIMESTEP
		if (t0 + dt > targetTime) {
			dt = targetTime - t0;
			ncycend = ncycle - 1;
		}

		// Construct total energy, reset Lagrangian coordinates to Eulerian grid
		for (int n = nmin; n <= nmax; n++) {
			e[n] = p[n] / (r[n] * gamm) + 0.5 * (pow(u[n],2) + pow(v[n], 2) + pow(w[n], 2));
			xa[n] = xa0[n];
			dx[n] = dx0[n];
		}

		// -------------------------------------------
		// THIS IS WHERE PPMLMR WOULD GET CALLED
		// -------------------------------------------

		// ADD BOUNDARY CONDTIONS
		addBoundaryConditions(0, 0);

		// FLATTENING COEFFICIENTS
		vector<double> flat = flatten();

		// GET PARABOLIC COEFFICIENTS
		vector<vector<double>> para = paraset(nmin - 4, nmax + 5, dx, xa);
		
		// INTERPOLATE PARABOLAE FOR FLUID VARIABLES
		vector<double> dp(imax + 12), dr(imax + 12), du(imax + 12);
		vector<double> p6(imax + 12), r6(imax + 12), u6(imax + 12);
		vector<double> pl(imax + 12), rl(imax + 12), ul(imax + 12);

		parabola(nmin - 4, nmax + 4, para, p, dp, p6, pl, flat);
		parabola(nmin - 4, nmax + 4, para, r, dr, r6, rl, flat);
		parabola(nmin - 4, nmax + 4, para, u, du, u6, ul, flat);


		// Integrate parabolae to get input states for Riemann problem
		vector<double> plft(imax + 12), ulft(imax + 12), rlft(imax + 12);
		vector<double> prgh(imax + 12), urgh(imax + 12), rrgh(imax + 12);
		states(pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh, dt);

		// Call the Riemann solver to obtain the zone face averages, umid and pmid
		vector<double> pmid(imax + 12), umid(imax + 12);
		riemann(nmin - 3, nmax + 4, gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid);

		// do lagrangian update using umid and pmid
		evolve(umid, pmid, dt);

		remap(flat);

		for (int i = 0; i < v.size(); i++) {
			if (v[i] != 0.0 || w[i] != 0.0) {
				bool breaker = 1;
			}
		}

		// update time and cycle counters
		t0 += dt;
		ncycle += 1;

		// TIME STEP CONTROL
		ridt = 0;
		for (int n = nmin; n <= nmax; n++) {
			svel = sqrt(gam * p[n] / r[n]) / dx0[n];
			xvel = abs(u[n]) / dx0[n];
			ridt = max(xvel, max(svel, ridt));
		}
		double dtx = courant / ridt;
		double dt3 = 1.1 * dt;
		dt = min(dt3, dtx);
	}

	bool stopper = 1;
}

void PPMSolver::addBoundaryConditions(int nleft, int nright)
{
	if (nleft == 0) {
		for (int n = 1; n <= 6; n++) {
			dx[nmin - n] = dx[nmin + n - 1];
			dx0[nmin - n] = dx0[nmin + n - 1];
			xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
			xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];
			r[nmin - n] = r[nmin + n - 1];
			u[nmin - n] = -u[nmin + n - 1];
			v[nmin - n] = v[nmin + n - 1];
			w[nmin - n] = w[nmin + n - 1];
			p[nmin - n] = p[nmin + n - 1];
			e[nmin - n] = e[nmin + n - 1];
			f[nmin - n] = f[nmin + n - 1];
		}
	}
	else if (nleft == 1) {
		for (int n = 1; n <= 6; n++) {
			dx[nmin - n] = dx[nmin];
			dx0[nmin - n] = dx0[nmin];
			xa[nmin - n] = xa[nmin - n + 1] - dx[nmin - n];
			xa0[nmin - n] = xa0[nmin - n + 1] - dx0[nmin - n];
			r[nmin - n] = r[nmin];
			u[nmin - n] = u[nmin];
			v[nmin - n] = v[nmin];
			w[nmin - n] = w[nmin];
			p[nmin - n] = p[nmin];
			e[nmin - n] = e[nmin];
			f[nmin - n] = f[nmin];
		}
	}

	if (nright == 0) {
		for (int n = 1; n <= 6; n++) {
			dx[nmax + n] = dx[nmax + 1 - n];
			dx0[nmax + n] = dx0[nmax + 1 - n];
			xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
			xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];
			r[nmax + n] = r[nmax + 1 - n];
			u[nmax + n] = -u[nmax + 1 - n];
			v[nmax + n] = v[nmax + 1 - n];
			w[nmax + n] = w[nmax + 1 - n];
			p[nmax + n] = p[nmax + 1 - n];
			e[nmax + n] = e[nmax + 1 - n];
			f[nmax + n] = f[nmax + 1 - n];
		}
	}
	else if (nright == 0) {
		for (int n = 1; n <= 6; n++) {
			dx[nmax + n] = dx[nmax];
			dx0[nmax + n] = dx0[nmax];
			xa[nmax + n] = xa[nmax + n - 1] + dx[nmax + n - 1];
			xa0[nmax + n] = xa0[nmax + n - 1] + dx0[nmax + n - 1];
			r[nmax + n] = r[nmax];
			u[nmax + n] = u[nmax];
			v[nmax + n] = v[nmax];
			w[nmax + n] = w[nmax];
			p[nmax + n] = p[nmax];
			e[nmax + n] = e[nmax];
			f[nmax + n] = f[nmax];
		}
	}
}

vector<double> PPMSolver::flatten()
{
	vector<double> steep(imax + 12), flat(imax+12);

	double omega1 = 0.75;
	double omega2 = 5.0;
	double epsilon = 0.33;

	double delp1, delp2, shock, temp1;
	for (int n = nmin - 4; n <= nmax + 4; n++) {
		delp1 = p[n + 1] - p[n - 1];
		delp2 = p[n + 2] - p[n - 2];
		// TODO if(abs(delp2) < small) delp2 = small
		shock = abs(delp1) / min(p[n + 1], p[n - 1]) - epsilon;
		shock = max(0.0, shock);
		if (shock > 0.0) shock = 1.0;
		if (u[n - 1] < u[n + 1]) shock = 0.0;
		temp1 = (delp1 / delp2 - omega1) * omega2;
		steep[n] = shock * max(0.0, temp1);
	}

	// Set phony boundary conditions for the steepness parameter
	steep[nmin - 5] = steep[nmin - 4];
	steep[nmax + 5] = steep[nmax + 4];

	// Set flattening coefficient based on the steepness in neighboring zones
	double temp2;
	for (int n = nmin - 4; n <= nmax + 4; n++) {
		temp2 = max( steep[n - 1], max(steep[n], steep[n + 1]) );
		flat[n] = max(0.0, min(0.5, temp2));
	}

	// flat(n) should be set to old_flat if no shock in this direction
	double old_flat;
	for (int n = nmin - 3; n <= nmax; n++) {
		old_flat = f[n] - int(f[n]);
		if (flat[n] > 0.0) {
			flat[n] = max(flat[n], old_flat);
			f[n] = max(flat[n] + (2.0 * 1 - 3.0), 0.0); // TODO make sure this is correct
		}
		else {
			f[n] = max(0.0, f[n] - 1.0);
			flat[n] = old_flat;
		}
	}
	return flat;	
}

vector<vector<double>> PPMSolver::paraset(int nmin, int nmax, vector<double> dx, vector<double> xa)
{
	vector<vector<double>> para(imax + 12, vector<double>(5));

	vector<double> a(imax + 12), ai(imax + 12), b(imax + 12), bi(imax + 12), c(imax + 12), ci(imax + 12);
	// main loop
	for (int n = nmin - 2; n <= nmax; n++) {
		a[n] = dx[n] + dx[n + 1];
		ai[n] = 1.0 / a[n];
		b[n] = a[n] + dx[n];
		bi[n] = 1.0 / b[n];
		c[n] = a[n] + dx[n + 1];
		ci[n] = 1.0 / c[n];
	}
	// Deal with last element
	a[nmax+1] = dx[nmax + 1] + 0;
	ai[nmax + 1] = 1.0 / a[nmax + 1];
	b[nmax + 1] = a[nmax + 1] + dx[nmax + 1];
	bi[nmax + 1] = 1.0 / b[nmax + 1];
	c[nmax + 1] = a[nmax + 1] + 0;
	ci[nmax + 1] = 1.0 / c[nmax + 1];

	// constants for equation 1.6
	vector<double> d(imax + 12);
	for (int n = nmin - 1; n <= nmax; n++) {
		d[n] = 1.0 / (a[n - 1] + a[n + 1]);
		para[n][0] = dx[n] * ai[n] + 2.0 * dx[n + 1] * dx[n] * d[n] * ai[n] * (a[n - 1] * bi[n] - a[n + 1] * ci[n]);
		para[n][1] = -d[n] * dx[n]     * a[n - 1] * bi[n];
		para[n][2] =  d[n] * dx[n + 1] * a[n + 1] * ci[n];
	}

	// constants for equation 1.7
	for (int n = nmin - 1; n <= nmax; n++) {
		d[n] = dx[n] / (a[n - 1] + dx[n + 1]);
		para[n][3] = d[n] * b[n - 1] * ai[n];
		para[n][4] = d[n] * c[n] * ai[n - 1];
	}
	// Deal with last element ( nmax+1 )
	d[nmax + 1] = dx[nmax + 1] / (a[nmax + 1 - 1] + 0);
	para[nmax + 1][3] = d[nmax + 1] * b[nmax + 1 - 1] * ai[nmax + 1];
	para[nmax + 1][4] = d[nmax + 1] * c[nmax + 1]     * ai[nmax + 1 - 1];


	return para;
}

void PPMSolver::parabola(int nmin, int nmax, const vector<vector<double>>& para, vector<double>& a,
	vector<double>& da, vector<double>& a6, vector<double>& al, const vector<double>& flat)
{
	vector<double> diffa(imax + 12);
	for (int n = nmin - 2; n <= nmax + 1; n++) {
		diffa[n] = a[n + 1] - a[n];
	}

	// Equation 1.7 of C&W
	for (int n = nmin - 1; n <= nmax + 1; n++) {
		da[n] = para[n][3] * diffa[n] + para[n][4] * diffa[n - 1];
		double temp = min( abs(da[n]), min( 2.0 * abs(diffa[n - 1]), 2.0 * abs(diffa[n])));
		// Eq 1.8
		if (da[n] >= 0) {
			da[n] = abs(temp);
		}
		else {
			da[n] = -abs(temp);
		}
	}

	// zero out da(n) if a(n) is a local max/min
	for (int n = nmin - 1; n <= nmax + 1; n++) {
		if (diffa[n - 1] * diffa[n] < 0.0) {
			da[n] = 0.0;
		}
	}

	// Equation 1.6 of C & W
	vector<double> ar(imax + 12);
	for (int n = nmin - 1; n <= nmax; n++) {
		ar[n] = a[n] + para[n][0] * diffa[n] + para[n][1] * da[n + 1] + para[n][2] * da[n];
		al[n + 1] = ar[n];
	}

	// eqn. 4.1 - flatten interpolation in zones with a shock ( flat(n)->1. )
	double onemfl;
	for (int n = nmin; n <= nmax; n++) {
		onemfl = 1.0 - flat[n];
		ar[n] = flat[n] * a[n] + onemfl * ar[n];
		al[n] = flat[n] * a[n] + onemfl * al[n];
	}

	/*MONOTONICITY constraints :
	compute delta_a, a_6
		MONOT : if a is a local max / min, flaten zone structure ar, al->a.
		MONOT : compute monotonized values using eq. 1.10 of C & W
		if parabola exceeds al / ar, reset ar / al so that slope -> 0.
		Recalculate delta_a and a_6*/
	vector<double> deltaa(imax + 12), scrch1(imax + 12), scrch2(imax + 12), scrch3(imax + 12);
	for (int n = nmin; n <= nmax; n++) {
		deltaa[n] = ar[n] - al[n];
		a6[n] = 6.0 * (a[n] - 0.5 * (al[n] + ar[n]));
		scrch1[n] = (ar[n] - a[n]) * (a[n] - al[n]);
		scrch2[n] = deltaa[n] * deltaa[n];
		scrch3[n] = deltaa[n] * a6[n];
	}

	// Eq 1.10
	for (int n = nmin; n <= nmax; n++) {
		if (scrch1[n] <= 0.0) {
			ar[n] = a[n];
			al[n] = a[n];
		}
		if (scrch2[n] < +scrch3[n]) al[n] = 3.0 * a[n] - 2.0 * ar[n];
		if (scrch2[n] < -scrch3[n]) ar[n] = 3.0 * a[n] - 2.0 * al[n];
	}

	// Eq 1.5
	for (int n = nmin; n <= nmax; n++) {
		deltaa[n] = ar[n] - al[n];
		a6[n] = 6.0 * (a[n] - 0.5 * (al[n] + ar[n]));
	}
}

void PPMSolver::states(vector<double>& pl, vector<double>& ul, vector<double>& rl,
	vector<double>& p6, vector<double>& u6, vector<double>& r6,
	vector<double>& dp, vector<double>& du, vector<double>& dr,
	vector<double>& plft, vector<double>& ulft, vector<double>& rlft,
	vector<double>& prgh, vector<double>& urgh, vector<double>& rrgh,
	double dt)
{
	double fourthd = 4.0 / 3.0;
	double hdt = 0.5 * dt;

	vector<double> Cdtdx(imax + 12), fCdtdx(imax + 12);

	for (int n = nmin - 4; n <= nmax + 4; n++) {
		Cdtdx[n] = sqrt(gam * p[n] / r[n]) / dx[n];
		svel = max(svel, Cdtdx[n]);
		Cdtdx[n] = Cdtdx[n] * hdt;
		fCdtdx[n] = 1.0 - fourthd * Cdtdx[n];
	}

	// WOULD CALL FORCES HERE, NOT GOING TO FOR NOW

	//Obtain averages of rho, u, and P over the domain(+/ -)Cdt
		//lft is the + wave on the left  side of the boundary
		//rgh is the - wave on the right side of the boundary
	int np;
	for (int n = nmin - 4; n <= nmax + 4; n++) {
		np = n + 1;
		plft[np] = pl[n] + dp[n] - Cdtdx[n] * (dp[n] - fCdtdx[n] * p6[n]);
		ulft[np] = ul[n] + du[n] - Cdtdx[n] * (du[n] - fCdtdx[n] * u6[n]);
		rlft[np] = rl[n] + dr[n] - Cdtdx[n] * (dr[n] - fCdtdx[n] * r6[n]);
		plft[np] = max(smallp, plft[np]);
		rlft[np] = max(smallr, rlft[np]);
		ulft[np] = ulft[np] + hdt * 0; // ulft(np) = ulft(np) + hdt*(grav(np)+fict(np))

		prgh[n] = pl[n] + Cdtdx[n] * (dp[n] + fCdtdx[n] * p6[n]);
		urgh[n] = ul[n] + Cdtdx[n] * (du[n] + fCdtdx[n] * u6[n]);
		rrgh[n] = rl[n] + Cdtdx[n] * (dr[n] + fCdtdx[n] * r6[n]);
		prgh[n] = max(smallp, prgh[n]);
		rrgh[n] = max(smallr, rrgh[n]);
		urgh[n] = urgh[n] + hdt * 0; // urgh[n] = urgh[n] + hdt * (grav[n] + fict[n]);
	}
}

void PPMSolver::riemann(int lmin, int lmax, double gamma, vector<double>& prgh, vector<double>& urgh, vector<double>& vrgh,
	vector<double>& plft, vector<double>& ulft, vector<double>& vlft,
	vector<double>& pmid, vector<double>& umid)
{
	double gamfac2 = gamma + 1.0;
	double gamfac1 = 0.5 * (gamfac2) / gamma;

	double smallp = 1.0e-25;
	double tol = 1.0e-5;

	// Obtain first guess for Pmid by assuming Wlft, Wrgh = Clft, Crgh
	vector<double> clft(imax + 12), crgh(imax + 12);
	vector<double> plfti(imax + 12), prghi(imax + 12);
	for (int l = lmin; l <= lmax; l++) {
		clft[l] = sqrt(gamma * plft[l] * vlft[l]);
		crgh[l] = sqrt(gamma * prgh[l] * vrgh[l]);
		vlft[l] = 1.0 / vlft[l];
		vrgh[l] = 1.0 / vrgh[l];
		plfti[l] = 1.0 / plft[l];
		prghi[l] = 1.0 / prgh[l];
		pmid[l] = prgh[l] - plft[l] - crgh[l] * (urgh[l] - ulft[l]);
		pmid[l] = plft[l] + pmid[l] * clft[l] / (clft[l] + crgh[l]);
		pmid[l] = max(smallp, pmid[l]);
	}

	/*Iterate up to 8 times using Newton's method to converge on correct Pmid
!    -use previous guess for pmid to get wavespeeds: wlft, wrgh
     -find the slope in the u-P plane for each state: zlft, zrgh
     -use the wavespeeds and pmid to guess umid on each side: umidl, umidr
     -project tangents from (pmid,umidl) and (pmid,umidr) to get new pmid
     -make sure pmid does not fall below floor value for pressure*/
	vector<double> pmold(imax + 12), wlft(imax + 12), wrgh(imax + 12), zlft(imax + 12), zrgh(imax + 12),
		umidl(imax + 12), umidr(imax + 12);

	for (int l = lmin; l <= lmax; l++) {
		for (int n = 1; n <= 12; n++) {
			pmold[l] = pmid[l];
			wlft[l] = 1.0 + gamfac1 * (pmid[l] - plft[l]) * plfti[l];
			wrgh[l] = 1.0 + gamfac1 * (pmid[l] - prgh[l]) * prghi[l];
			wlft[l] = clft[l] * sqrt(wlft[l]);
			wrgh[l] = crgh[l] * sqrt(wrgh[l]);
			zlft[l] = 4.0 * vlft[l] * wlft[l] * wlft[l];
			zrgh[l] = 4.0 * vrgh[l] * wrgh[l] * wrgh[l];
			zlft[l] = -zlft[l] * wlft[l] / (zlft[l] - gamfac2 * (pmid[l] - plft[l]));
			zrgh[l] = zrgh[l] * wrgh[l] / (zrgh[l] - gamfac2 * (pmid[l] - prgh[l]));
			umidl[l] = ulft[l] - (pmid[l] - plft[l]) / wlft[l];
			umidr[l] = urgh[l] + (pmid[l] - prgh[l]) / wrgh[l];
			pmid[l] = pmid[l] + (umidr[l] - umidl[l]) * (zlft[l] * zrgh[l]) / (zrgh[l] - zlft[l]);
			pmid[l] = max(smallp, pmid[l]);
			if (abs(pmid[l] - pmold[l]) / pmid[l] < tol) break;
		}
	}

	for (int l = lmin; l <= lmax; l++) {
		umidl[l] = ulft[l] - (pmid[l] - plft[l]) / wlft[l];
		umidr[l] = urgh[l] + (pmid[l] - prgh[l]) / wrgh[l];
		umid[l] = 0.5 * (umidl[l] + umidr[l]);
	}
}

void PPMSolver::evolve(vector<double>& umid, vector<double>& pmid, double dt) {
	vector<double> dvol1 = dx;
	vector<double> dm(imax + 12), dtbdm(imax + 12), xa1(imax+12), upmid(imax + 12);
	for (int n = nmin - 3; n <= nmax + 4; n++) {
		dm[n] = r[n] * dvol1[n];
		dtbdm[n] = dt / dm[n];
		xa1[n] = xa[n];
		xa[n] = xa[n] + dt * umid[n] / 1.0;
		upmid[n] = umid[n] * pmid[n];
	}

	xa1[nmin - 4] = xa[nmin - 4];
	xa1[nmax + 5] = xa[nmax + 5];

	vector<double> xa2(imax + 12), xa3(imax + 12);
	for (int n = nmin - 4; n <= nmax + 5; n++) {
		xa2[n] = xa1[n] + 0.5 * dx[n];
		dx[n] = xa[n + 1] - xa[n];
		xa3[n] = xa[n] + 0.5 * dx[n];
	}

	dvol = dx;

	vector<double> amid(imax + 12, 1.0);

	vector<double> uold(imax + 12);
	for (int n = nmin - 3; n <= nmax + 3; n++) {
		// density evolution. lagrangian code, so all we have to do is watch the change in the geometry.
		r[n] = r[n] * (dvol1[n] / dvol[n]);
		r[n] = max(r[n], smallr);

		// velocity evolution due to pressure accelerationand forces.
		uold[n] = u[n];
		u[n] = u[n] - dtbdm[n] * (pmid[n + 1] - pmid[n]) * 0.5 * (amid[n + 1] + amid[n]);

		// total energy evolution
		e[n] = e[n] - dtbdm[n] * (amid[n + 1] * upmid[n + 1] - amid[n] * upmid[n]) + 0.5 * dt * (uold[n] * 0 + u[n] * 0); // 0.5*dt*(uold(n)*grav0(n) + u(n)*grav1(n))
		q[n] = e[n] - 0.5 * (pow(u[n],2) + pow(v[n],2) + pow(w[n],2));
		q[n] = max(q[n], smallp / (gamm * r[n]));
	}
}

void PPMSolver::remap(vector<double>& flat) {
	double fourthd = 4.0 / 3.0;

	vector<vector<double>> para = paraset(nmin - 1, nmax + 1, dx, xa);

	vector<double> dr(imax + 12), du(imax + 12), dv(imax + 12), dw(imax + 12), dq(imax + 12), de(imax + 12);
	vector<double> r6(imax + 12), u6(imax + 12), v6(imax + 12), w6(imax + 12), q6(imax + 12), e6(imax + 12);
	vector<double> rl(imax + 12), ul(imax + 12), vl(imax + 12), wl(imax + 12), ql(imax + 12), el(imax + 12);

	parabola(nmin - 1, nmax + 1, para, r, dr, r6, rl, flat);
	parabola(nmin - 1, nmax + 1, para, u, du, u6, ul, flat);
	parabola(nmin - 1, nmax + 1, para, v, dv, v6, vl, flat);
	parabola(nmin - 1, nmax + 1, para, w, dw, w6, wl, flat);
	parabola(nmin - 1, nmax + 1, para, q, dq, q6, ql, flat);
	parabola(nmin - 1, nmax + 1, para, e, de, e6, el, flat);

	// instead of calling volume
	vector<double> dvol0 = dx0;

	// Calculate the volume of the overlapping subshells (delta)
	vector<double> delta(imax + 12);
	for (int n = nmin; n <= nmax + 1; n++) {
		delta[n] = xa[n] - xa0[n];
	}

	/*Calculate the total mass(fluxr), momentum(fluxu), and energy(fluxe)
		in the subshell created by the overlap of the Lagrangianand Eulerican grids.
		If the zone face has moved to the left(deltx > 0), use the integral from the
		left side of zone n(fluxrr).If the zone face has moved to the right
		(deltx < 0), use the integral from the right side of zone nn = n - 1 (fluxrl).*/
	vector<double> fluxr(imax + 12), fluxu(imax + 12), fluxv(imax + 12), fluxw(imax + 12), fluxe(imax + 12), fluxq(imax + 12);
	
	double deltx, fractn, fractn2;
	int nn;

	for (int n = nmin; n <= nmax + 1; n++) {
		deltx = xa[n] - xa0[n];
		if (deltx >= 0.0) {
			nn = n - 1;
			fractn = 0.5 * deltx / dx[nn];
			fractn2 = 1. - fourthd * fractn;
			fluxr[n] = (rl[nn] + dr[nn] - fractn * (dr[nn] - fractn2 * r6[nn])) * delta[n];
			fluxu[n] = (ul[nn] + du[nn] - fractn * (du[nn] - fractn2 * u6[nn])) * fluxr[n];
			fluxv[n] = (vl[nn] + dv[nn] - fractn * (dv[nn] - fractn2 * v6[nn])) * fluxr[n];
			fluxw[n] = (wl[nn] + dw[nn] - fractn * (dw[nn] - fractn2 * w6[nn])) * fluxr[n];
			fluxe[n] = (el[nn] + de[nn] - fractn * (de[nn] - fractn2 * e6[nn])) * fluxr[n];
			fluxq[n] = (ql[nn] + dq[nn] - fractn * (dq[nn] - fractn2 * q6[nn])) * fluxr[n];
		}
		else {
			fractn = 0.5 * deltx / dx[n];
			fractn2 = 1. + fourthd * fractn;
			fluxr[n] = (rl[n] - fractn * (dr[n] + fractn2 * r6[n])) * delta[n];
			fluxu[n] = (ul[n] - fractn * (du[n] + fractn2 * u6[n])) * fluxr[n];
			fluxv[n] = (vl[n] - fractn * (dv[n] + fractn2 * v6[n])) * fluxr[n];
			fluxw[n] = (wl[n] - fractn * (dw[n] + fractn2 * w6[n])) * fluxr[n];
			fluxe[n] = (el[n] - fractn * (de[n] + fractn2 * e6[n])) * fluxr[n];
			fluxq[n] = (ql[n] - fractn * (dq[n] + fractn2 * q6[n])) * fluxr[n];
		}
	}

	/*Advect mass, momentum, and energy by moving the subshell quantities
	into the appropriate Eulerian zone.*/
	vector<double> dm(imax + 12), dm0(imax + 12);
	for (int n = nmin - 1; n <= nmax + 1; n++) {
		dm[n] = r[n] * dvol[n];
		dm0[n] = (dm[n] + fluxr[n] - fluxr[n + 1]);
		r[n] = dm0[n] / dvol0[n];
		r[n] = max(smallr, r[n]);
		dm0[n] = 1. / (r[n] * dvol0[n]);
		u[n] = (u[n] * dm[n] + fluxu[n] - fluxu[n + 1]) * dm0[n];
		v[n] = (v[n] * dm[n] + fluxv[n] - fluxv[n + 1]) * dm0[n];
		w[n] = (w[n] * dm[n] + fluxw[n] - fluxw[n + 1]) * dm0[n];
		e[n] = (e[n] * dm[n] + fluxe[n] - fluxe[n + 1]) * dm0[n];
		q[n] = (q[n] * dm[n] + fluxq[n] - fluxq[n + 1]) * dm0[n];
	}

	// If flow is highly supersonic remap on internal energy, else on total E
	double ekin;
	for (int n = nmin; n <= nmax; n++) {
		ekin = 0.5 * (pow(u[n], 2) + pow(v[n], 2) + pow(w[n], 2));
		if (ekin / q[n] < 100.0) {
			q[n] = e[n] - ekin;
		}
		p[n] = gamm * r[n] * q[n];
		p[n] = max(smallp, p[n]);
	}

}