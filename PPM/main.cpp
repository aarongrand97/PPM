#include <iostream>

#include <vector>

#include "PPMSolver.h"

using std::vector;

int main() {
	int zoneCount = 100;
	vector<double> pressures(zoneCount);
	vector<double> uVelocity(zoneCount);
	vector<double> density(zoneCount);
	vector<double> vVelocity(zoneCount);
	vector<double> wVelocity(zoneCount);

	for (int i = 0; i < zoneCount / 2; i++) {
		pressures[i] = 500000;
		density[i] = 5.8060411979716511;
	}
	for (int i = zoneCount / 2; i < pressures.size(); i++) {
		pressures[i] = 100000;
		density[i] = 1.1612082395943304;
	}

	double dx = 1 / ((double)zoneCount - 1);
	double xmin = 0.0;
	double xmax = 1.0;
	double C = 0.5;
	double t = 2.5e+01;

	PPMSolver ppmSolver = PPMSolver(pressures, density, dx, xmin, xmax, C, t, zoneCount);
	ppmSolver.solve();

}