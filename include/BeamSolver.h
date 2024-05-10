#ifndef __GENESIS_BEAMSOLVER__
#define __GENESIS_BEAMSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>



class Field;


#include "Undulator.h"
#include "EFieldSolver.h"
#include "TrackBeam.h"



using namespace std;


class BeamSolver {
public:
    BeamSolver();
    virtual ~BeamSolver();
    void initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr, bool hghgrange, int maxharm);
    void advance(double, Beam *, vector<Field *> *, Undulator *);
    void track(double, Beam *, Undulator *, bool);
    void applyR56(Beam *, Undulator *, double);
    double getSCField(int);
    double getHGHGSCField(int);
    void checkAllocation(unsigned long i);

private:
    complex<double> cpart;
    vector<double> rharm;
    vector<complex<double> > rpart;

    double ez{};
    double xks{}, xku{};

    double theta{}, gamma{}, btpar{};
    double dgamma{};
    double k2gg{}, k2pp{}, k3gg{}, k3pp{};

    bool onlyFundamental;

    void RungeKutta(double);
    void ODE(double, double);

    EFieldSolver efield;
    TrackBeam tracker;

};

inline double BeamSolver::getSCField(int islice){
    return efield.getSCField(islice);
}

inline double BeamSolver::getHGHGSCField(int islice){
    return efield.getHGHGSCField(islice);
}

inline void BeamSolver::initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr, bool hghgrange, int maxharm){
    efield.init(rmax, ngrid, nz, nphi, lambda, longr, hghgrange, maxharm);
}


inline void BeamSolver::track(double dz, Beam *beam, Undulator *und, bool last)
{
    tracker.track(dz,beam,und,last);
}

inline void BeamSolver::applyR56(Beam *beam, Undulator *und, double reflen){
    tracker.applyR56(beam,und,reflen);
}


#endif
