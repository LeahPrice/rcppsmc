// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA_auto.cpp: Rcpp wrapper for SMC library -- A simple example for estimating
// the parameters of a linear regression model using likelihood annealing SMC,
// with adaptation of the temperature schedule and the multivariate normal random
// walk covariance matrix.
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
#include "LinReg_LA_auto.h"
#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>

namespace LinReg_LA_auto {
	const double a_prior = 3.0;
	const double b_prior = pow(2.0*300.0*300.0,-1.0);
}
using namespace std;
using namespace LinReg_LA_auto;


// LinRegLA_auto() function callable from R via Rcpp:: 
// [[Rcpp::export]]
Rcpp::List LinRegLA_auto_cpp(arma::mat data, unsigned long inlNumber, double resampTol, double tempTol) { 	


	try {
		lNumber = inlNumber;
		rho = tempTol;
		
		lIterates = data.n_rows;
		y.data_x = data.col(0);
		y.data_y = data.col(1);
		mean_x = arma::sum(y.data_x)/lIterates;
		
		
		myParams = new rad_adapt;
		//Initialise and run the sampler
		smc::sampler<rad_state,smc::staticModelAdapt> Sampler(lNumber, HistoryType::RAM, myParams);
		smc::moveset<rad_state> Moveset(fInitialise, fMove, fMCMC);
		
		Sampler.SetResampParams(ResampleType::SYSTEMATIC, resampTol);
		Sampler.SetMoveSet(Moveset);
		//Sampler.SetAdaptSet(myParams);
		Sampler.Initialise();
		
		std::vector<double> temps;
		temps.push_back(0);
		temps.push_back(myParams->GetParams().GetTemp());
		
		std::vector<double> ESS;
		ESS.push_back(Sampler.GetESS());

		long n = 0;
		Rcpp::Rcout << "Current temperature is : " << myParams->GetParams().GetTemp() << std::endl;
		
		while (myParams->GetParams().GetTemp() != 1){
			n++;
			Sampler.Iterate();
			Rcpp::Rcout << "Current temperature is : " << myParams->GetParams().GetTemp() << std::endl;
			temps.push_back(myParams->GetParams().GetTemp());
			ESS.push_back(Sampler.GetESS());
		}
		
		double logNC = Sampler.GetLogNCPath();
		
		delete myParams;
		
		return Rcpp::List::create(Rcpp::Named("logNC") = logNC,Rcpp::Named("Temps") = temps,Rcpp::Named("ESS") = ESS);
	}
	catch(smc::exception  e) {
		Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o
	}
	return R_NilValue;          	// to provide a return 
}

namespace LinReg_LA_auto {


	///The function corresponding to the log likelihood at specified position
	///  \param value     The state to consider 
	double logLikelihood(const rad_state & value){

		double sigma = pow(expl(value.theta(2)),0.5);
		arma::vec mean_reg = value.theta(0)*arma::ones(lIterates) + value.theta(1)*(y.data_x - mean_x*arma::ones(lIterates));
		return arma::sum(-log(sigma)*arma::ones(lIterates) - pow(y.data_y - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI)*arma::ones(lIterates));	

	}
	///The function corresponding to the log prior at specified position
	///  \param value     The state to consider 
	double logPrior(const rad_state & value){
		return -log(1000.0)- pow(value.theta(0) - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(value.theta(1) - 185.0,2.0)/(2.0*100.0*100.0) + value.theta(2)-1.0/b_prior/expl(value.theta(2)) -value.theta(2)*(a_prior+1.0);
	}

	///A function to initialise a particle

	/// \param value		Reference to the empty particle value
	/// \param logweight	Refernce to the empty particle log weight
	void fInitialise(rad_state & value, double & logweight)
	{
		// drawing from the prior
		value.theta = arma::zeros(3);
		value.theta(0) = R::rnorm(3000.0,1000.0);
		value.theta(1) = R::rnorm(185.0,100.0);
		value.theta(2) = log(pow(R::rgamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
		value.loglike = logLikelihood(value);
		logweight = myParams->GetParams().GetTemp()*value.loglike;
	}

	///The move function.

	///\param lTime			The sampler iteration.
	/// \param value		Reference to the current particle value
	/// \param logweight	Refernce to the current particle log weight
	void fMove(long lTime, rad_state & value, double & logweight)
	{
		logweight += (myParams->GetParams().GetTemp() - myParams->GetParams().GetTempPrevious())*value.loglike;
	}

	///The proposal function.

	///\param lTime The sampler iteration.
	///\param value Reference to the value of the particle being moved
	int fMCMC(long lTime, rad_state & value)
	{
		double MH_ratio;
		double dRand;
		int count = 0;
		rad_state value_prop;
		double logprior_curr = logPrior(value);
		double logprior_prop;
		arma::mat randPrep(3,3);
		
		for (unsigned int j=0; j<10; j++){
			value_prop.theta = value.theta + myParams->GetParams().GetCholCov()*Rcpp::as<arma::vec>(Rcpp::rnorm(3));
			value_prop.loglike = logLikelihood(value_prop);
			logprior_prop = logPrior(value_prop);
			
			MH_ratio = exp(myParams->GetParams().GetTemp()*(value_prop.loglike - value.loglike) + logprior_prop - logprior_curr);
			dRand = R::runif(0,1);
			
			if (MH_ratio>dRand){
				value = value_prop;
				logprior_curr = logprior_prop;
				count++;
			}
		}
		return count; 
	}
}