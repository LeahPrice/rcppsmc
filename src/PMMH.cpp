// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// PMMH.cpp: Example 3.1 of Andrieu et al. (2010). Implementing particle marginal
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
#include "PMMH.h"
#include "rngR.h"
#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>

namespace PMMH {
    const double a_prior = 0.01;
	const double b_prior = 0.01;
}

using namespace std;
using namespace PMMH;


// PMMH() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List PMMH_cpp(arma::vec data, unsigned long inlNumber, unsigned long lMCMCits) { 	
  
  long lIterates;
  Rcpp::NumericVector sigv(lMCMCits);
  Rcpp::NumericVector sigw(lMCMCits);
  
  try {
    lNumber = inlNumber; // number of particles
    
    lIterates = data.n_rows;
    y = data;
	
	loglike = arma::zeros(lMCMCits);
	logprior = arma::zeros(lMCMCits);
	
	double dRand;
	double loglike_prop;
	double logprior_prop;
	long accept_count = 0;
	
    //Initialise and run the sampler
    smc::sampler<double> Sampler(lNumber, HistoryType::NONE);
	smc::rng* pRng1 = new smc::rng();
	theta_prop.sigv = sqrt(10.0);
	theta_prop.sigw = 1.0;
		
	sigv(0) = theta_prop.sigv;
	sigw(0) = theta_prop.sigw;
	
	
	// Getting a particle filtering estimate of the log likelihood.
	smc::moveset<double> Moveset(fInitialise, fMove, NULL);
    Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
    Sampler.SetMoveSet(Moveset);
    Sampler.Initialise();
	Sampler.IterateUntil(lIterates-1);
	loglike(0) = Sampler.GetLogNCPath();
	
	// Inverse gamma prior
	logprior(0) = 2*a_prior*log(b_prior)-2*lgamma(a_prior)-(a_prior+1)*log(sigv(0))-b_prior/sigv(0)-(a_prior+1)*log(sigw(0))-b_prior/sigw(0);
	
	double MH_ratio;
	for (unsigned int i = 1; i<lMCMCits; i++){
		// RW proposal for parameters
		theta_prop.sigv = pRng1->Normal(sigv(i-1),0.15);
		theta_prop.sigw = pRng1->Normal(sigw(i-1),0.08);
		//theta_prop.sigv = pow(pRng1->Gamma(a_prior,b_prior),-1.0);
		//theta_prop.sigw = pow(pRng1->Gamma(a_prior,b_prior),-1.0);
		
		
		// Getting a particle filtering estimate of the log likelihood.
		Sampler.Initialise();
		Sampler.IterateUntil(lIterates-1);
		loglike_prop = Sampler.GetLogNCPath();
		
		// Inverse gamma prior
		logprior_prop = 2*a_prior*log(b_prior)-2*lgamma(a_prior)-(a_prior+1)*log(theta_prop.sigv)-b_prior/theta_prop.sigv-(a_prior+1)*log(theta_prop.sigw)-b_prior/theta_prop.sigw;
		
		MH_ratio = exp(loglike_prop - loglike(i-1));
		//MH_ratio = exp(loglike_prop - loglike(i-1) + logprior_prop - logprior(i-1));
		dRand = pRng1->Uniform(0,1);
      
		if (MH_ratio>dRand){
			sigv(i) = theta_prop.sigv;
			sigw(i) = theta_prop.sigw;
			loglike(i) = loglike_prop;
			logprior(i) = logprior_prop;
			accept_count++;
			//Rcpp::Rcout << "ACCEPTED value at " << i+1 << " with a value of \t" << sigv(i) << "\t " << sigw(i) << " with loglike of \t" << loglike(i) << std::endl;
		} else {
			sigv(i) = sigv(i-1);
			sigw(i) = sigw(i-1);
			loglike(i) = loglike(i-1);
			logprior(i) = logprior(i-1);
		}
	}
	delete pRng1;
		   
    return Rcpp::DataFrame::create(Rcpp::Named("samples_sigv") = sigv,
								   Rcpp::Named("samples_sigw") = sigw,
                                   Rcpp::Named("loglike") = loglike,
                                   Rcpp::Named("logprior") = logprior);
  }
  catch(smc::exception  e) {
    Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
    //exit(e.lCode);		// we're just called from R so we should not exit
  }
  return R_NilValue;          	// to provide a return 
}

namespace PMMH {
///A function to initialise a particle

/// \param pRng 		A pointer to the random number generator which is to be used
/// \param X			A reference to the empty particle value
/// \param logweight	A reference to the empty particle log weight
void fInitialise(smc::rng *pRng, double & X, double & logweight)
{
	X = pRng->Normal(0.0,sqrt(5.0));
	double mean = pow(X,2)/20.0;
	logweight = -log(theta_prop.sigw) - pow(y(0) - mean,2.0)/(2.0*theta_prop.sigw*theta_prop.sigw) -0.5*log(2.0*M_PI);
}

///The proposal function.

///\param lTim		The sampler iteration.
///\param X			A reference to the current particle value
///\param logweight	A reference to the current particle log weight
///\param pRng		A random number generator.
void fMove(long lTime, double & X, double & logweight, smc::rng *pRng)
{
	X = X/2.0 + 25.0*X/(1+pow(X,2)) + 8*cos(1.2*(lTime+1)) + pRng->Normal(0.0,theta_prop.sigv);
	double mean = pow(X,2)/20.0;
	logweight += -log(theta_prop.sigw) - pow(y(lTime) - mean,2.0)/(2.0*theta_prop.sigw*theta_prop.sigw) -0.5*log(2.0*M_PI);
}

}