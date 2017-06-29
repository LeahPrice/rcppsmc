// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// PMMH.h: Example 3.1 of Andrieu et al. (2010). Implementing particle marginal
// Metropolis-Hastings for an application previous described in Gordon et al. (1993)
// and Kitagawa (1996).
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

class parameters 
{
public:
    double sigv, sigw;
};

namespace PMMH {
    arma::vec y; //data
	arma::vec loglike;
	arma::vec logprior;
	unsigned long lNumber;
	
	double logWeight(long lTime, double & value);
	double logPosterior(long lTime, double & value);
	void fInitialise(smc::rng *pRng, double & value, double & logweight);
	void fMove(long lTime, double & value, double & logweight, smc::rng *pRng);
	
	parameters theta_prop;
	unsigned int accept_count;
}

