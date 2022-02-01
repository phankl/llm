#include <string>
#include <vector>

#include "solver.h"

int main() {

  // variables
  
  // grid
  int nGrid = 100;
  double xMax = 1.0;

  // time
  double cfl = 0.01;
  double tMax = 10.0;

  // physics
  double l = 1.0;
  double c = 1.0;
  double r0 = 1.0;
  double j0 = 1.0;
  double e0 = 1000.0;
  double omega = 1.0;

  // states
  vector<double> jPrevious(nGrid, 0.0);
  vector<double> jCurrent(nGrid, 0.0);

  // io
  string filename = "llm.dat";
  double outputStep = 0.01;

  // initialise solver

  Solver solver;
  solver.initGrid(nGrid, xMax);
  solver.initTime(cfl, tMax);
  solver.initPhysics(l, c, r0, j0, e0, omega);
  solver.initStates(jPrevious, jCurrent);
  solver.initIO(filename, outputStep);

  // run simulation

  solver.run();

  return 0;
}
