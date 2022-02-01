#include <string>
#include <vector>

#include "solver.h"

int main() {

  // physical and mathematical constants

  double pi = 3.14159265359;
  double h = 6.62607004e-34;
  double e = 1.60217662e-19;
  double vF = 8.0e5;
  double g = 0.25;

  // variables
  
  // grid
  int nGrid = 100;
  double xMax = 100.0e-6;

  // physics
  double l = 0.5 * h / (e*e * vF);
  double c = 8.0 * e*e * g*g / (h * vF);
  double r0 = 1.0e9;
  double j0 = 25.0e-6;
  double e0 = 1.0e5;
  double frequency = 10.0e6;
  double omega = 2.0*pi * frequency;

  // time
  double periods = 10.0;
  double tMax = periods / frequency;
  double dt = tMax / 1.0e7;

  // states
  vector<double> jPrevious(nGrid, 0.0);
  vector<double> jCurrent(nGrid, 0.0);

  // io
  string filename = "llm.dat";
  int outputNumber = 100;
  double outputStep = tMax / outputNumber;

  // initialise solver

  Solver solver;
  solver.initGrid(nGrid, xMax);
  solver.initTime(dt, tMax);
  solver.initPhysics(l, c, r0, j0, e0, omega);
  solver.initStates(jPrevious, jCurrent);
  solver.initIO(filename, outputStep);

  // run simulation

  solver.run();

  return 0;
}
