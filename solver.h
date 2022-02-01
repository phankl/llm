#ifndef SOLVER_HEADER
#define SOLVER_HEADER

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class Solver{
  public:
    // constructor and destructor
    
    Solver();
    ~Solver();

    // grid parameters
    
    int nGrid;
    double xMax;
    
    double dx;

    // time parameters

    double cfl;
    double tMax;

    // physical parameters

    double l;
    double c;
    double r0;
    double j0;
    double e0;
    double omega;

    double waveSpeed;

    // states

    vector<double> jCurrent;
    vector<double> jPrevious;

    vector<double> djCurrent;

    // io parameters

    string filename;
    double outputStep;

    // functions

    void initGrid(int, double);
    void initTime(double, double);
    void initPhysics(double, double, double, double, double, double);
    void initStates(vector<double>, vector<double>);
    void initIO(string, double);

    void run();

  private:
    
    // initalisation bools

    bool gridInitialised;
    bool timeInitialised;
    bool physicsInitialised;
    bool statesInitialised;
    bool ioInitialised;

    // internal parameters

    int nStep;
    double t;
    double dt;
    double e;
    
    vector<double> jNext;

    ofstream outFile;

    // run functions

    void calcGradient();
    void fluxUpdate();
    void sourceUpdate();

    void printFile();
};

#endif
