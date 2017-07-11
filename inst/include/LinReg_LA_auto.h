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
#include <RcppArmadillo.h>


namespace LinReg_LA_auto {
	
	class rad_state
	{
	public:
		arma::vec theta; //alpha, beta, phi
		//arma::vec::fixed<3> theta; //alpha, beta, phi
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
	
	class rad_params:
	public smc::algParam<rad_state>
	{
	private:
		double temp_curr;
		double temp_previous;
		bool done;
		arma::mat cholCov;
	public:
		
		rad_params(long n, ResampleType::Enum restype, double resthresh) : algParam(n,restype,resthresh) {temp_curr = 0; temp_previous = 0; done = 0;}
		
		///Returns a random number generated from a Beta distribution with the specified parameters.
		void updateForMoveExtra(const smc::population<rad_state> & pop) {
			if (CESSdiff(pop,1-temp_curr,rho*this->GetN())>=0){
				temp_curr = 1;
				done = 1;
			} else {
				temp_curr = bisection(temp_curr, temp_curr, 1, pop, rho*this->GetN());
			}
		}
		
		///Returns a random number generated from a Beta distribution with the specified parameters.
		void updateForMCMCExtra(const smc::population<rad_state> & pop) {
			calcCholCov(pop);
		}
		
		///Returns a random number generated from a Beta distribution with the specified parameters.
		void updateEndExtra(const smc::population<rad_state> & pop) {
			temp_previous = temp_curr;
		}
		
		double GetTemp(void){return temp_curr;}
		double GetTempPrevious(void){return temp_previous;}
		double GetDone(void){return done;}
		arma::mat GetCholCov(void){return cholCov;}
		
		///Free the workspace allocated for the algorithm parameters
		~rad_params() {};
		
		double bisection(double curr, double a, double b, const smc::population<rad_state> & pop, const double & desiredESS){
			double f_a = CESSdiff(pop,a-curr,desiredESS);
			double f_b = CESSdiff(pop,b-curr,desiredESS);
			if (f_a*f_b>0){
				Rcpp::stop("Bisection method to choose next temperature failed");
			} else{
				double m, f_m, err;
				m = (a+b)/2.0; 
				f_m = CESSdiff(pop,m-curr,desiredESS);
				//Rcpp::Rcout << "a=" << a << "\t b=" << b  << "\t m=" << m << "\t f(a)=" << f_a << "\t f(b)=" << f_b << "\t f(m)=" << f_m << std::endl;
				err = 10;
				while (err > 0.01){
					if (f_m<0){
						//if (f_a*f_b<0){
						b = m;
						f_b = f_m;
					} else{
						a = m;
						f_a = f_m;
					}
					m = (a+b)/2.0;
					f_m = CESSdiff(pop,m-curr,desiredESS);
					err = abs(f_m);
					//Rcpp::Rcout << "a=" << a << "\t b=" << b  << "\t m=" << m << "\t f(a)=" << f_a << "\t f(b)=" << f_b << "\t f(m)=" << f_m << std::endl;
				}
				return m;
			}
		}

		/// Computes the conditional ESS minus the desired CESS
		double CESSdiff(const smc::population<rad_state> & pop, const double & multiplier, const double & desiredCESS){
			long double sum1 = 0;
			long double sum2 = 0;

			for(int i = 0; i < this->GetN(); i++) {
				sum1 += expl(pop.GetLogWeightN(i) + multiplier*pop.GetValueN(i).loglike);
				sum2 += expl(pop.GetLogWeightN(i) + 2*multiplier*pop.GetValueN(i).loglike);
			}

			return expl(log(this->GetN()) + 2*log(sum1) - log(sum2)) - desiredCESS;
		}
		
		/// Gets the cholesky decomposition of the current covariance matrix
		void calcCholCov(const smc::population<rad_state> & pop){
			long N = this->GetN();
			arma::vec normWeights = exp(pop.GetLogWeight() - log(sum(exp(pop.GetLogWeight()))));
			arma::mat thetaMat(N,3);
			for (long i=0; i<N; i++){
				thetaMat.row(i) = pop.GetValueN(i).theta.t();
			}
			
			arma::mat diff = thetaMat - arma::ones(N,1)*arma::mean(thetaMat,0);
			arma::mat emp_cov = diff.t()*diagmat(normWeights)*diff;
			
			//emp_cov.print("The covariance is: ");
			
			cholCov = arma::chol(emp_cov);
		}
		
	};
	
	rad_params * myParams;
	
	
	
	
}