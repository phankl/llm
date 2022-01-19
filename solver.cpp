#include "solver.h"

Solver::Solver():
  gridInitialised(false),
  timeInitialised(false),
  physicsInitialised(false),
  statesInitialised(false),
  ioInitialised(false) {
}

Solver::~Solver() {
  outFile.close();
}

void initGrid(int nGridIn, double xMaxIn) {
  
  nGrid = nGridIn;
  xMax = xMaxIn;

  dx = xMax / (nGrid - 1);

  gridInitialised = true;
}

void initTime(double cflIn, double tMaxIn) {
  
  cfl = cflIn;
  tMax = tMaxIn;

  nStep = 0;
  t = 0.0;

  timeInitialised = true;
}

void initPhysics(double lIn, cIn, r0In, j0In) {
  
  l = lIn;
  c = cIn;
  r0 = r0In;
  j0 = j0In;

  waveSpeed = sqrt(l*c);

  physicsInitialised = true;
}

void initStates(vector<double> jPreviousIn, vector<double> jCurrentIn) {

  jPrevious = jPreviousIn;
  jCurrent = jCurrentIn;

  statesInitialised = true;
}

void initIO(string filenameIn, double outputStepIn) {

  filename = filenameIn;
  outputStep = outputStepIn;

  outFile = ofstream(filename);

  ioInitialised = true;
}

void run() {
  if (!(gridInitialised && timeInitialised && physicsInitialised && statesInitialised && ioInitialised)) {
    cout << "Simulation not initialised properly!" << endl;
    exit(EXIT_FAILURE);
  }
  
  bool end = false;
  while (!end) {
    
    // end condition

    if (t + dt > tMax) {
      dt = tMax - t;
      end = true;
    }
    
    fluxUpdate();
    sourceUpdate();

    t += dt;
    nStep++;
  }
}
