// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg.h: Rcpp wrapper for SMC library -- A simple example for estimating
// the parameters of a linear regression model using data annealing SMC.
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



namespace LinReg {
	unsigned long lNumber;
	
	class rad_state 
	{
	public:
		double alpha, beta, phi;
	};

	class rad_obs
	{
	public:
		arma::vec data_x, data_y;
	};

	class rad_params:
	public smc::algParam<rad_state>
	{
	private:
		double std_alpha,std_beta,std_phi;
	public:
		
		///Returns a random number generated from a Beta distribution with the specified parameters.
		void updateForMCMCExtra(const smc::population<rad_state> & pop) {
			arma::vec vAlpha(lNumber), vBeta(lNumber), vPhi(lNumber);
			arma::vec weights = pop.GetWeight();
			for (unsigned int i=0; i<lNumber; i++){
				vAlpha(i) = pop.GetValueN(i).alpha;
				vBeta(i) = pop.GetValueN(i).beta;
				vPhi(i) = pop.GetValueN(i).phi;
			}
			std_alpha = sqrt(arma::sum(weights%pow(vAlpha-arma::sum(vAlpha%weights),2)));
			std_beta = sqrt(arma::sum(weights%pow(vBeta-arma::sum(vBeta%weights),2)));
			std_phi = sqrt(arma::sum(weights%pow(vPhi-arma::sum(vPhi%weights),2)));
		}
		
		double GetStdAlpha(void){return std_alpha;}
		double GetStdBeta(void){return std_beta;}
		double GetStdPhi(void){return std_phi;}
	};
	
	
	rad_params * myParams;
	rad_obs y;
	double mean_x;
	
	double logWeight(long lTime, const rad_state & value);
	double logPosterior(long lTime, const rad_state & value);
	void fInitialise(rad_state & value, double & logweight);
	void fMove(long lTime, rad_state & value, double & logweight);
	int fMCMC(long lTime, rad_state & value);

	double integrand_mean_alpha(const rad_state&, void*);
	double integrand_mean_beta(const rad_state&, void*);
	double integrand_mean_phi(const rad_state&, void*);
	double integrand_var_alpha(const rad_state&, void*);
	double integrand_var_beta(const rad_state&, void*);
	double integrand_var_phi(const rad_state&, void*);
}
