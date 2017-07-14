// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA_auto.h: Rcpp wrapper for SMC library -- A simple example for estimating
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


namespace LinReg_LA_auto {
	
	class rad_state
	{
	public:
		arma::vec theta; //alpha, beta, phi
		double loglike;
	};

	class rad_obs
	{
	public:
		arma::vec data_x, data_y;
	};
	
	rad_obs y;
	double mean_x;
	unsigned long lNumber;
	long lIterates;
	double rho;
	
	double logLikelihood(const rad_state & value);
	double logPrior(const rad_state & value);
	
	void fInitialise(rad_state & value, double & logweight);
	void fMove(long lTime, rad_state & value, double & logweight);
	int fMCMC(long lTime, rad_state & value);
	
	class rad_adapt:
	public smc::algParam<rad_state,smc::staticModelAdapt>
	{
	public:
		
		void updateForMove(const smc::population<rad_state> & pop) {
			arma::vec loglike(pop.GetN());
			for (int i=0; i<pop.GetN(); i++){
				loglike(i) = pop.GetValueN(i).loglike;
			}
			param.ChooseTemp(pop.GetLogWeight(),loglike,rho*lNumber);
		}
		
		void updateForMCMC(const smc::population<rad_state> & pop) {
			arma::mat thetaMat(pop.GetN(),3);
			for (long i=0; i<pop.GetN(); i++){
				thetaMat.row(i) = pop.GetValueN(i).theta.t();
			}
			param.calcCholCov(thetaMat,pop.GetLogWeight());
		}
		
		void updateEnd(const smc::population<rad_state> & pop) {
			param.SetTempPrevious(param.GetTemp());
		}
		
		~rad_adapt() {};
		
	};
	
	smc::algParam<rad_state,smc::staticModelAdapt> * myParams;
	
	
	
	std::vector<int> mcmcRepeats;
	
}