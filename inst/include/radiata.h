// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// radiata.h: Rcpp wrapper for SMC library -- TESTING Radiata pine example

#include "smctc.h"

class rad_state 
{
public:
    double alpha, beta, phi;
};

class rad_obs
{
public:
    double data_x, data_y;
};

double logWeight(long lTime, const rad_state & X);

smc::particle<rad_state> fInitialiseRAD(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<rad_state> & p, smc::rng *pRng);
void fMoveRAD(long lTime, smc::particle<rad_state> & pFrom, smc::rng *pRng);

extern double mean_x;

extern std::vector<rad_obs> yRAD;