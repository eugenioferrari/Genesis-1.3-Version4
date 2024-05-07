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
    void init(double, int, int, int, double, bool, bool, int);
    void shortRange(vector<Particle> *, double, double, int);
    void hghgRange(vector<Particle> *, double, double, double, int, double);
    void longRange(Beam *beam, double gamma, double aw);
    double getEField(unsigned long i);
    bool hasShortRange() const;
    bool hasHGHGRange() const;
    void allocateForOutput(unsigned long nslice);
    double getSCField(int);
    double getHGHGLSC(unsigned long i);

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
    vector<double> hghgez;

    int nz, nphi, ngrid, rank;
    int maxharm;
    double rmax, ks, xcen, ycen, dr;
    double sigmax, sigmay;
    bool longrange;
    bool hghgrange;

};

inline double EFieldSolver::getSCField(int islice) {
    return efield[islice]*511000;  // convert from Lorentz mass unit to eV /m
}

inline bool EFieldSolver::hasShortRange() const{
    return (nz>0) & (ngrid > 2);
}

inline bool EFieldSolver::hasHGHGRange() const{
    return (maxharm > 0) & hghgrange;
}

#endif
