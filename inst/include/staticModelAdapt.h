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

#include <RcppArmadillo.h>


namespace smc {
	
class staticModelAdapt{
	private:
		double temp_curr;
		double temp_previous;
		arma::mat cholCov;
		
		/// Computes the conditional ESS given the specified temperature differenceminus the desired conditional ESS
		double CESSdiff(const arma::vec & logweight, const arma::vec & loglike, double tempDiff, double desiredCESS){
			double sum1 = arma::sum(exp(logweight + tempDiff*loglike));
			double sum2 = arma::sum(exp(logweight + 2*tempDiff*loglike));

			return expl(log(logweight.n_rows) + 2*log(sum1) - log(sum2)) - desiredCESS;
		}
		
		///Performs the bisection method 
		double bisection(double curr, const arma::vec & logweight, const arma::vec & loglike, const double & desiredCESS){
			double a = curr;
			double b = 1;
			double f_a = CESSdiff(logweight,loglike,a-curr,desiredCESS);
			double f_b = CESSdiff(logweight,loglike,b-curr,desiredCESS);
			if (f_a*f_b>0){
				Rcpp::stop("Bisection method to choose next temperature failed");
			} else{
				double m, f_m, err;
				m = (a+b)/2.0; 
				f_m = CESSdiff(logweight,loglike,m-curr,desiredCESS);
				//Rcpp::Rcout << "a=" << a << "\t b=" << b  << "\t m=" << m << "\t f(a)=" << f_a << "\t f(b)=" << f_b << "\t f(m)=" << f_m << std::endl;
				err = 10;
				while (err > 0.01){
					if (f_m<0){
						b = m;
						f_b = f_m;
					} else{
						a = m;
						f_a = f_m;
					}
					m = (a+b)/2.0;
					f_m = CESSdiff(logweight,loglike,m-curr,desiredCESS);
					err = abs(f_m);
					//Rcpp::Rcout << "a=" << a << "\t b=" << b  << "\t m=" << m << "\t f(a)=" << f_a << "\t f(b)=" << f_b << "\t f(m)=" << f_m << std::endl;
				}
				return m;
			}
		}
		
	public:
		///Free the workspace allocated for the algorithm parameters
		~staticModelAdapt() {};
		
		staticModelAdapt() {temp_curr = 0; temp_previous = 0;}
		
		/// Determines the next temperature using the bisection method
		void ChooseTemp(const arma::vec & logweight, const arma::vec & loglike, double desiredCESS) {
			if (CESSdiff(logweight,loglike,1-temp_curr,desiredCESS)>=0){
				temp_curr = 1;
			} else {
				temp_curr = bisection(temp_curr, logweight, loglike, desiredCESS);
			}
		}
		
		
		/// Gets the cholesky decomposition of the current covariance matrix. Theta should be stored in a [Nxd] matrix where d is the dimension of the parameter
		void calcCholCov(const arma::mat & theta, const arma::vec logweight){
			long N = logweight.n_rows;
			arma::vec normWeights = exp(logweight - log(sum(exp(logweight))));
						
			arma::mat diff = theta - arma::ones(N,1)*arma::mean(theta,0);
			arma::mat emp_cov = diff.t()*diagmat(normWeights)*diff;
			
			cholCov = arma::chol(emp_cov);
		}
		
		double GetTemp(void){return temp_curr;}
		void SetTemp(double tempin){temp_curr = tempin;}
		
		double GetTempPrevious(void){return temp_previous;}
		void SetTempPrevious(double tempin){temp_previous = tempin;}
		
		arma::mat GetCholCov(void){return cholCov;}
	};

	
}