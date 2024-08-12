#ifndef __GENESIS_EFIELDSOLVER__
#define __GENESIS_EFIELDSOLVER__

#include <vector>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <mpi.h>

#include "Particle.h"

class Beam;

using namespace std;

extern const double vacimp;
extern const double eev;



class EFieldSolver {
public:
    EFieldSolver();
    virtual ~EFieldSolver();
    void init(double, int, int, int, double, bool, bool, int, double);
    void shortRange(vector<Particle> *, double, double, int);
    // Arguments: slice, current, slicelength, slicespacing, Ldrift, islice
    void hghgRange(vector<Particle> *, double, double, double, double, int);
    void longRange(Beam *beam, double gamma, double aw);
    double getEField(unsigned long i);
    double getHGHGdgamma(unsigned long i);
    bool hasShortRange() const;
    bool hasHGHGRange() const;
    bool hasLongRange() const;
    void allocateForOutput(unsigned long nslice);
    void allocateForHGHGOutput(unsigned long nslice);
    double getSCField(int);
    double getHGHGSCField(int);

private:
    void analyseBeam(vector<Particle> *beam);
    void analyseBeamHGHG(vector<Particle> *beam);
    void constructLaplaceOperator();
    void tridiag();

    vector<double> work1, work2, fcurrent, fsize;  // used for long range calculation
    vector<complex<double> > cwork;
    vector<int> idxr;
    vector<double> lmid, rlog, vol, ldig;
    vector<complex<double> > csrc, clow, cmid, cupp, celm, gam; // used for tridiag routine
    vector<double> ez,efield;
    // hghgez contains the dgamma for each particle
    // hghgefield contains the sum over all dgamma for each slice
    vector<double> hghgez, hghgefield;

    int nz, nphi, ngrid, rank;
    int maxharm;  // the maximum harmonic to be considered for computing the field (default = 0)
    double scaling;  // the scaling factor for the LSC (default = 0)
    double rmax, ks, xcen, ycen, dr;
    double sigmax, sigmay;
    bool longrange;
    bool hghgrange;

};

inline double EFieldSolver::getSCField(int islice) {
    return efield[islice]*511000;  // convert from Lorentz mass unit to eV /m
}

inline double EFieldSolver::getHGHGSCField(int islice) {
    // Returns the hghgefield for the islice
    return hghgefield[islice];  // in Lorentz for now
}

inline bool EFieldSolver::hasShortRange() const{
    return (nz > 0) & (ngrid > 2);
}

inline bool EFieldSolver::hasLongRange() const{
    return longrange;
}

inline bool EFieldSolver::hasHGHGRange() const{
    return (maxharm > 0) & (scaling != 0.);
}

#endif
