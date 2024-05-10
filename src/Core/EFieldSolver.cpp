#include "EFieldSolver.h"
#include "Beam.h"



EFieldSolver::EFieldSolver() {
    nz = 0;
    nphi = 0;
    ngrid = 0;
    rmax = 0;
    ks = 1;
    xcen = 0;
    ycen = 0;
    dr = 1;
    longrange = false;
    hghgrange = false;
    maxharm = 0;
    fcurrent.clear();
    fsize.clear();
    efield.clear();
    hghgez.clear();
    hghgefield.clear();
    rank=0;
    if (!MPISingle){
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); // assign rank to node
    }
}

EFieldSolver::~EFieldSolver()= default;

void EFieldSolver::init(double rmax_in, int ngrid_in, int nz_in, int nphi_in, double lambda, bool longr_in, bool hghgrange_in, int maxharm_in) {

    rmax = rmax_in;
    ngrid = ngrid_in;
    nz = nz_in;
    nphi = nphi_in;
    ks = 4 * asin(1) / lambda;
    longrange = longr_in;
    hghgrange = hghgrange_in;
    maxharm = maxharm_in;

    // adjust working arrays
    if (ngrid != csrc.size()) {
        vol.resize(ngrid);
        ldig.resize(ngrid+1);
        rlog.resize(ngrid);
        lmid.resize(ngrid);
        clow.resize(ngrid);
        cmid.resize(ngrid);
        cupp.resize(ngrid);
        celm.resize(ngrid);
        csrc.resize(ngrid);
        gam.resize(ngrid);
    }
}

void EFieldSolver::allocateForOutput(unsigned long nslice){
    if (nslice == efield.size()){ return;}
    efield.resize(nslice, 0);
}

void EFieldSolver::allocateForHGHGOutput(unsigned long nslice){
    if (nslice == hghgefield.size()){ return;}
    hghgefield.resize(nslice, 0);
}

void EFieldSolver::longRange(Beam *beam, double gamma0, double aw) {
    auto nsize = beam->beam.size();
    // check if the beam size has been changed:
    if (nsize != beam->longESC.size()){
        beam->longESC.resize(nsize);
    }

    for (unsigned long i =0; i < nsize; i++){
        beam->longESC[i]=0;
    }
    if (!longrange) {
        return;
    }
    double gamma = gamma0/sqrt(1+aw*aw);

    int MPIsize=1;
    int MPIrank=0;
    if (!MPISingle){
        MPI_Comm_size(MPI_COMM_WORLD, &MPIsize); // assign rank to node
        MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank); // assign rank to node
    }

    // resize if needed also after harmonic conversion.
    if (fcurrent.size()!=(MPIsize*nsize)){
        fcurrent.resize(MPIsize*nsize);
        fsize.resize(MPIsize*nsize);
        work1.resize(nsize);
        work2.resize(nsize);
    }

    // fill local slices
    for (int i=0; i<nsize;i++){
        work1[i]=beam->current[i];
        work2[i]=beam->getSize(i);
    }

    // gathering information on full current and beam size profile
    MPI_Allgather(&work1.front(),static_cast<int>(nsize),MPI_DOUBLE,
                  &fcurrent.front(),static_cast<int>(nsize),MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(&work2.front(),static_cast<int>(nsize),MPI_DOUBLE,
                  &fsize.front(),static_cast<int>(nsize),MPI_DOUBLE, MPI_COMM_WORLD);

    double scl = beam->slicelength/2./asin(1)/299792458.0/2/8.85e-12;  // convert to units of electron rest mass.

    for (int i=0; i < nsize; i++){
        double EFld = 0;
        auto isplit = nsize*MPIrank+i;
        for (int j = 0; j < isplit; j++) {
            double ds = static_cast<double>(j - isplit) * beam->slicelength*gamma;
            double coef = 1 +ds/sqrt(ds*ds+fsize[j]); // ds is negative
            EFld += coef*scl*fcurrent[j]/fsize[j];
        }
        for (auto j = isplit+1; j < MPIsize*nsize; j++) {
            double ds = static_cast<double>(j - isplit) * beam->slicelength*gamma;
            double coef = 1 -ds/sqrt(ds*ds+fsize[j]); // ds is negative
            EFld -= coef*scl*fcurrent[j]/fsize[j];
        }
        beam->longESC[i] = EFld;
    }
}


void EFieldSolver::analyseBeam(vector<Particle> *beam){
/*
 *  calculates the center of the beam slice and its extension (5 times rms radius to construct the space charge grid
 */
    auto npart = beam->size();

    if (npart <1) {return;}

    // allocate new working space
    if (npart> cwork.size()) {
        cwork.resize(npart);
        idxr.resize(npart);
    }

    // calculate center of beam slice
    xcen = 0;
    ycen = 0;

    for (int i = 0; i < npart; i++) {
        xcen += beam->at(i).x;
        ycen += beam->at(i).y;
    }
    xcen /= static_cast<double>(npart);
    ycen /= static_cast<double>(npart);

    // calculate radial extension of particle distribution
    //   d

    double rbound = 0;
    for (int i = 0; i < npart; i++) {
        double tx = beam->at(i).x - xcen;
        double ty = beam->at(i).y - ycen;
        double rad2 = tx * tx + ty * ty;
        if (rad2 > rbound) { rbound = rad2; }
    }
    rbound = sqrt(rbound);

    if (rbound > rmax) {
        rmax = rbound*1.5;  // 20% safety margin
        cout << "*** Info (Rank "<< rank << "): Adjusting radial grid size for space charge calculation" << endl;
        cout << "    New grid size: " << rmax <<  " to hold all electrons" << endl;
    }

    dr = rmax/static_cast<double>(ngrid-1);

    for (int i = 0; i < npart; i++) {
        double tx = beam->at(i).x - xcen;
        double ty = beam->at(i).y - ycen;
        double radi = sqrt(tx * tx + ty * ty);
        cwork[i] = complex<double> (tx/radi,ty/radi);   // save argument e^i phi
        idxr[i] = static_cast<int>(floor(radi/dr));
    }
}


void EFieldSolver::analyseBeamHGHG(vector<Particle> *beam){
/*
 *  calculates the center of the beam slice and its extension (5 times rms radius to construct the space charge grid
 */
    auto npart = beam->size();

    if (npart <1) {return;}

    // calculate center of beam slice
    xcen = 0;
    ycen = 0;

    for (int i = 0; i < npart; i++) {
        xcen += beam->at(i).x;
        ycen += beam->at(i).y;
    }
    xcen /= static_cast<double>(npart);
    ycen /= static_cast<double>(npart);

    // calculate radial extension of particle distribution
    //   d
    double radx2 = 0;
    double rady2 = 0;
    for (int i = 0; i < npart; i++) {
        double tx = beam->at(i).x - xcen;
        double ty = beam->at(i).y - ycen;
        radx2 += tx * tx;
        rady2 += ty * ty;
    }
    radx2 /= static_cast<double>(npart);
    rady2 /= static_cast<double>(npart);

    sigmax = sqrt(radx2);
    sigmay = sqrt(rady2);
}

void EFieldSolver::constructLaplaceOperator(){

    // r_j are the center grid points. The boundaries are given +/-dr/2. Thus r_j=(j+1/2)*dr
    auto pi = 2*asin(1);

    vol[0] = pi *dr*dr;
    rlog[0] = 0.5;
    ldig[0] = 0;
    for (int j = 1; j < ngrid; j++) {
        vol[j] = pi * dr*dr * (2*j+1);
        ldig[j] = 2*pi*j;
        rlog[j] = log(static_cast<double>(j + 1) / static_cast<double>(j));
    }
    ldig[ngrid] = 0;
}

double EFieldSolver::getEField(unsigned long i)
{
    return ez[i];
}

double EFieldSolver::getHGHGdgamma(unsigned long i)
{
    return hghgez[i];
}

void EFieldSolver::shortRange(vector<Particle> *beam, double current, double gz2, int islice) {

    auto npart = beam->size();
    if (npart > ez.size()){
        ez.resize(npart);
    }
    for (int i =0; i < npart; i++){
        ez[i] = 0;
    }
    efield[islice] = 0;

    if (!this->hasShortRange()) { return; }
    // calculate center of beam slice and its extension
    this->analyseBeam(beam);
    if (npart == 0) { return; }

    // calculates the tri-diagonal laplace matrix of space charge field equation
    this->constructLaplaceOperator();

    auto pi = 2*asin(1);
    auto coef = -gz2/ks/ks;      // equivalent to 1/[k^2-(k+ku)^2] from fortran code.
    auto econst = vacimp/eev*current/ks/static_cast<double>(npart);
    for (int m = -nphi; m <=nphi; m++){  // loop over azimutal modes
        for (int i = 0; i < ngrid; i++){
            lmid[i] = -ldig[i] - ldig[i+1] - 2.*pi*m*m*rlog[i];
        }
        lmid[ngrid-1] -=2*pi*ngrid;

        for (int l = 1 ; l <= nz; l++) {   // loop over longitudinal modes

            // clear source term
            for (int i =0; i < ngrid; i++){
                csrc[i] = complex<double> (0,0);
            }

            // build source term

            for (int i=0; i < npart; i++){
                csrc[idxr[i]] += pow(cwork[i],-m)*complex<double>(cos(l*beam->at(i).theta), -sin(l*beam->at(i).theta));
            }

            // finalize Laplace operator and source term
            // note that the entire equation is scaled with gammaz^2/l^2/k^2 -> see Eq. 3.29 - 3.33 of thesis
            for (int i =0 ; i < ngrid; i++){
                csrc[i] *= complex<double>(0,econst/l/vol[i]);
                clow[i] = complex<double>(coef*ldig[i]/l/l/vol[i],0);
                cmid[i] = complex<double>(1+coef*lmid[i]/l/l/vol[i],0);
                cupp[i] = complex<double>(coef*ldig[i+1]/l/l/vol[i],0);
            }
            this->tridiag();

            if (l==1){
                efield[islice] = abs(celm[0]);
            }

            for (int i=0; i < npart; i++){
                complex<double> ctemp = pow(cwork[i],m)*complex<double>(cos(l*beam->at(i).theta), sin(l*beam->at(i).theta));
                ez[i] += 2*real(ctemp*celm[idxr[i]]);
            }
        }
    }
}


void EFieldSolver::hghgRange(vector<Particle> *beam, double current, double slicelength, double slicespacing, double Ldrift, int islice) {

    auto npart = beam->size();
    if (npart > hghgez.size()){
        hghgez.resize(npart);
    }
    for (int ip = 0; ip < npart; ip++){
        hghgez[ip] = 0;
    }
    hghgefield[islice] = 0;

    if (!this->hasHGHGRange()) { return; }
    // calculate center of beam slice and its extension
    // cout << "Analysis" << endl;
    this->analyseBeamHGHG(beam);
    // cout << "B done" << endl;
    // sigmax, sigmay are now updated
    if (npart == 0) { return; }

    // calculate the particle density in the slice
    auto c = 299792458.0; // velocity of light in m/s
    auto Q = 1.602176634e-19; // elementary charge in Coulomb
    auto eps0 = 8.854187817620389e-12; // electric permittivity in vacuum (electric constant)
    auto me_eV = 510998.95; //  # mass of electron in eV
    auto pi = 2 * asin(1);
    auto n0 = (current / c) / (Q * pi * sigmax * sigmay);
    auto lambda_seed = slicespacing;
    auto k_seed = 2 * pi / lambda_seed;

    // Correction for the phase of the current spike
    // compute the bunching at the fundamental
    complex<double> bunching = 0;
    const   complex<double> img(0.0,1.0);
    //
    // bunching = np.sum(np.exp(1j * k_seed * data_space * 1)) / npart
    // bunching_phase = np.angle(bunching)
    // t_phase = bunching_phase * slicespacing / (2 * np.pi)
    // data_space -= t_phase
    double s_now;
    for (int ip = 0; ip < npart; ip++) {
        s_now = beam->at(ip).theta * slicelength / (2 * pi);
        bunching += exp(img * k_seed * s_now * static_cast<double>(1));
    }
    bunching /= static_cast<double>(npart);
    double bunching_phase = arg(bunching);
    double t_phase = bunching_phase * slicespacing / (2 * pi);
    /*
    cout << "t_phase=" << t_phase << endl;
    cout << "bunching " << bunching << endl;
    cout << "n0" << n0 << endl;
    cout << "sigmax" << sigmax << endl;
    cout << "sigmay" << sigmay << endl;
    */
    // double Bh;

    for (int nh = 1; nh <= maxharm; nh++) {
        double Bh;
        bunching = 0;
        for (int ip = 0; ip < npart; ip++) {
            s_now = beam->at(ip).theta * slicelength / (2 * pi) - t_phase;
            bunching += exp(img * k_seed * s_now * static_cast<double>(nh));
        }
        bunching /= static_cast<double>(npart);
        Bh = abs(bunching);
        // cout << "h=" << nh << " b=" << Bh << endl;
        for (int ip = 0; ip < npart; ip++) {
            s_now = beam->at(ip).theta * slicelength / (2 * pi) - t_phase;
            hghgez[ip] += Bh * sin(static_cast<double>(nh) * k_seed * s_now) / static_cast<double>(nh);
        }
    }
    // finally convert to dgamma
    // dgamma = campo * 2 * Q * n0 / (eps0 * k_seed) * Ldrift / me_eV
    // hghgefield[islice] = 0;  // sanity is the accumulation of all the dgamma over the slice
    for (int ip = 0; ip < npart; ip++) {
        hghgez[ip] = hghgez[ip] * 2 * Q * n0 / (eps0 * k_seed) * Ldrift / me_eV;
        hghgefield[islice] += hghgez[ip];
        // cout << "ip=" << ip << " dgamma=" << dgamma << endl;
    }
    // cout << "Sanity check =" << sanity << endl;
}

void EFieldSolver::tridiag(){

    complex<double> bet = cmid[0];
    celm[0] = csrc[0] / bet;
    for (int i = 1; i < ngrid; i++) {
        gam[i] = cupp[i - 1] / bet;
        bet = cmid[i] - clow[i] * gam[i];
        celm[i] = (csrc[i] - clow[i] * celm[i - 1]) / bet;
    }
    for (int i = ngrid - 2; i > -1; i--) {
        celm[i] = celm[i] - gam[i + 1] * celm[i + 1];
    }
}


