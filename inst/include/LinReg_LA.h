// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA.h: Rcpp wrapper for SMC library -- A simple example for estimating
// the parameters of a linear regression model using likelihood annealing SMC.
//
// Copyright (C) 2017         Dirk Eddelbuettel, Adam Johansen and Leah South
//
// This file is part of RcppSMC.
//
// RcppSMC is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppSMC is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

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

namespace LinReg_LA {
    rad_obs2 y;
	double mean_x;
	unsigned long lNumber;
	long lIterates;
	arma::vec temps;
	
	arma::vec logLikelihood(const std::vector<rad_state> & X);
	double logLikelihood_single(const rad_state & X);
	arma::vec logPrior(const std::vector<rad_state> & X);
	double logPrior_single(const rad_state & X);
	
	smc::population<rad_state> fInitialise(smc::rng *pRng);
	long fSelect(long lTime, const smc::population<rad_state> & p, smc::rng *pRng);
	void fMove(long lTime, smc::population<rad_state> & pFrom, smc::rng *pRng);
	int fMCMC(long lTime, smc::population<rad_state> & pFrom, smc::rng *pRng);
	
	double integrand_mean_alpha(const rad_state&, void*);
	double integrand_mean_beta(const rad_state&, void*);
	double integrand_mean_phi(const rad_state&, void*);
	double integrand_var_alpha(const rad_state&, void*);
	double integrand_var_beta(const rad_state&, void*);
	double integrand_var_phi(const rad_state&, void*);
	
}