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

void Solver::initGrid(int nGridIn, double xMaxIn) {
  
  nGrid = nGridIn;
  xMax = xMaxIn;

  dx = xMax / (nGrid - 1);

  gridInitialised = true;
}

void Solver::initTime(double dtIn, double tMaxIn) {
  
  dt = dtIn;
  tMax = tMaxIn;

  nStep = 0;
  t = 0.0;

  timeInitialised = true;
}

void Solver::initPhysics(double lIn, double cIn, double r0In, double j0In, double e0In, double omegaIn) {
  
  l = lIn;
  c = cIn;
  r0 = r0In;
  j0 = j0In;
  e0 = e0In;
  omega = omegaIn;

  waveSpeed = sqrt(l*c);

  physicsInitialised = true;
}

void Solver::initStates(vector<double> jPreviousIn, vector<double> jCurrentIn) {

  jPrevious = jPreviousIn;
  jCurrent = jCurrentIn;

  jNext = vector<double>(nGrid, 0.0);
  djCurrent = vector<double>(nGrid, 0.0);

  statesInitialised = true;
}

void Solver::initIO(string filenameIn, double outputStepIn) {

  filename = filenameIn;
  outputStep = outputStepIn;

  outFile = ofstream(filename);

  ioInitialised = true;
}

void Solver::calcGradient() {
  
  // boundary values

  djCurrent[0] = jCurrent[1] / (dx*dx);
  djCurrent[nGrid-1] = jCurrent[nGrid-2] / (dx*dx);

  // cout << "Gradient: ";
  for (int i = 1; i < nGrid-1; i++) {
    djCurrent[i] = (jCurrent[i-1] - 2.0*jCurrent[i] + jCurrent[i+1]) / (dx*dx);
    // cout << djCurrent[i] << " ";
  }
  // cout << endl;
}

void Solver::fluxUpdate() {

  // boundary values

  jNext[0] = 0.0;
  jNext[nGrid-1] = 0.0;

  // cout << "Flux update: ";
  for (int i = 1; i < nGrid-1; i++) {
    jNext[i] = 2.0*jCurrent[i] - jPrevious[i] + cfl*cfl*djCurrent[i];
    // cout << jNext[i] << " ";
  }
  // cout << endl;
}

void Solver::sourceUpdate() {
  
  // boundary values

  jNext[0] = 0.0;
  jNext[nGrid-1] = 0.0;

  e = e0 * sin(omega*t); 
  
  // cout << "Source update: ";
  for (int i = 1; i < nGrid-1; i++) {
    if (e > 0) {
      double temp = l*(jNext[i] + j0) + dt*(e + j0*r0);
      double root = -4.0*j0*l*(dt*e + jNext[i]*l) + temp*temp;
      jNext[i] = 0.5*(temp - sqrt(root)) / l;
    }
    else {
      double temp = l*(jNext[i] - j0) + dt*(e - j0*r0);
      double root = 4.0*j0*l*(dt*e + jNext[i]*l) + temp*temp;
      jNext[i] = 0.5*(temp + sqrt(root)) / l;
    }
    // cout << jNext[i] << " "; 
  }
  // cout << endl;
}

void Solver::printFile() {

  for (int i = 0; i < nGrid; i++) {
    double x = i * dx;
    if (x > xMax) x = xMax;
    outFile << x << " " << jCurrent[i] << endl;
  }

  outFile << endl;
}

void Solver::run() {
  if (!(gridInitialised && timeInitialised && physicsInitialised && statesInitialised && ioInitialised)) {
    cout << "Simulation not initialised properly!" << endl;
    exit(EXIT_FAILURE);
  }

  // check if CFL conditions is satisfied

  cfl = waveSpeed * dt / dx;
  double dtMax = dx / waveSpeed;
  if (cfl > 1.0) {
    cout << "CFL condition not met! CFL number is " << cfl << "." << endl;
    cout << "Maximum timestep is " << dtMax << " seconds." << endl;
    exit(EXIT_FAILURE);
  }

  // inital output

  cout << "Time: " << t << "/" << tMax << endl;
  printFile();

  double outputCounter = 0.0;
  bool end = false;
  while (!end) {
    
    // end condition

    if (t + dt > tMax) {
      dt = tMax - t;
      end = true;
    }
    
    // time update

    calcGradient();
    fluxUpdate();
    t += dt;
    sourceUpdate();

    jPrevious = jCurrent;
    jCurrent = jNext;

    outputCounter += dt;
    nStep++;

    // file output

    if (outputCounter > outputStep)
      outputCounter = 0.0;
    if (outputCounter == 0.0 || end) {
      cout << "Time: " << t << "/" << tMax << endl;
      printFile();
    }
  }

  cout << "Finished simulation after " << nStep << " steps." << endl;
}
