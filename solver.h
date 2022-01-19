#ifndef SOLVER_HEADER
#define SOLVER_HEADER

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class SOLVER{
  public:

    // grid parameters
    
    int nGrid;
    double xMax;
    
    double dx;

    // time parameters

    double cfl;
    double tMax;

    double dt;

    // physical parameters

    double l;
    double c;
    double r0;
    double j0;

    double waveSpeed;

    // states

    vector<double> prevState;
    vector<double> currState;

    // io parameters

    string filename;
    double outputStep;

    // functions

    void initGrid(int, double);
    void initTime(double, double);
    void initPhysics(double, double, double, double);
    void initStates(vector<double>, vector<double>);
    void initIO(string, double);

    void run(double);

  private:

    void fluxUpdate();
    void sourceUpdate();

    void printFile();
}

#endif
