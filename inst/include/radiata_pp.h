// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// radiata.h: Rcpp wrapper for SMC library -- TESTING Radiata pine example

#include "smctc.h"
#include <RcppArmadillo.h>

class rad_state
{
public:
    double alpha, beta, phi;
};

class rad_obs2
{
public:
    arma::vec data_x, data_y;
};

namespace radiata_pp {
    rad_obs2 y;
	double mean_x;
	unsigned long lNumber;
	long lIterates;
	arma::vec temps;
	
	arma::vec logLikelihood(const std::vector<rad_state> & X);
	double logLikelihood_single(const rad_state & X);
	arma::vec logPrior(const std::vector<rad_state> & X);
	double logPrior_single(const rad_state & X);
	
	smc::particle<rad_state> fInitialise(smc::rng *pRng);
	long fSelect(long lTime, const smc::particle<rad_state> & p, smc::rng *pRng);
	void fMove(long lTime, smc::particle<rad_state> & pFrom, smc::rng *pRng);
	int fMCMC(long lTime, smc::particle<rad_state> & pFrom, smc::rng *pRng);
	
	double integrand_mean_alpha(const rad_state&, void*);
	double integrand_mean_beta(const rad_state&, void*);
	double integrand_mean_phi(const rad_state&, void*);
	double integrand_var_alpha(const rad_state&, void*);
	double integrand_var_beta(const rad_state&, void*);
	double integrand_var_phi(const rad_state&, void*);
	
}