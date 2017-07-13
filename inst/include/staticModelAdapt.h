// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// staticModelAdapt.h: A class containing parameters and functions to update
// these parameters in the context of static Bayesian models. The methods to
// estimate the empirical covariance and its Cholesky decomposition are applicable
// for all applications where the particle values are a vector of doubles. The methods
// to compute the next 'temperature' are relevant to likelihood annealing SMC where the
// power posteriors are defined by P_t(theta|y) \propto P(y|theta)^\gamma_t P(theta).
// Here y is observed data, theta denotes the parameters of the model and \gamma_t
// denotes the 'temperatures' which start at 0 and finish at 1 (the posterior).
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

/// A class containing parameters and functions to update these parameters in the context of static Bayesian models
class staticModelAdapt{
	private:
		/// The current temperature
		double temp_curr;
		/// The previous temperature
		double temp_previous;
		/// The Cholesky decomposition of the empirical covariance matrix estimated from the population of particles
		arma::mat cholCov;
		
		/// Computes the difference between the conditional ESS given the specified temperature difference and the desired conditional ESS
		double CESSdiff(const arma::vec & logweight, const arma::vec & loglike, double tempDiff, double desiredCESS){
			double sum1 = arma::sum(exp(logweight + tempDiff*loglike));
			double sum2 = arma::sum(exp(logweight + 2*tempDiff*loglike));

			return expl(log(logweight.n_rows) + 2*log(sum1) - log(sum2)) - desiredCESS;
		}
		
		///Performs the bisection method to find the temperature which gives the desired conditional ESS with initial limits [temp_curr,1]
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
		
		///The class constructor which sets the current and previous temperatures to zero.
		staticModelAdapt() {temp_curr = 0; temp_previous = 0;}
		
		/// Determines the next temperature by using the bisection method to find the temperature which gives the desired conditional ESS
		///
		/// \param logweight An armadillo vector with the logarithm of the current particle weights
		/// \param loglike An armadillo vector containing the log likelihood of the current particle values
		/// \param desiredCESS The target conditional ESS for the next temperature (generally fixed)
		void ChooseTemp(const arma::vec & logweight, const arma::vec & loglike, double desiredCESS) {
			if (CESSdiff(logweight,loglike,1-temp_curr,desiredCESS)>=0){
				temp_curr = 1;
			} else {
				temp_curr = bisection(temp_curr, logweight, loglike, desiredCESS);
			}
		}
		
		
		/// Calculates the Cholesky decomposition of the empirical covariance matrix based on the current weighted particle set
		///
		/// \param theta An [Nxd] armadillo matrix of doubles for the current particle values, where N is
		/// the number of particles and d is the dimension of the parameter
		/// \param logweight An armadillo vector of the logarithm of the current particle weights
		void calcCholCov(const arma::mat & theta, const arma::vec logweight){
			long N = logweight.n_rows;
			arma::vec normWeights = exp(logweight - log(sum(exp(logweight))));
						
			arma::mat diff = theta - arma::ones(N,1)*arma::mean(theta,0);
			arma::mat emp_cov = diff.t()*diagmat(normWeights)*diff;
			
			cholCov = arma::chol(emp_cov);
		}
		
		/// Returns the current temperature
		double GetTemp(void){return temp_curr;}
		/// Sets the current temperature
		void SetTemp(double tempin){temp_curr = tempin;}
		
		/// Returns the previous temperature
		double GetTempPrevious(void){return temp_previous;}
		/// Sets the previous temperature
		void SetTempPrevious(double tempin){temp_previous = tempin;}
		
		/// Returns the Cholesky decomposition of the empirical covariance matrix based on the current weighted particle set
		arma::mat GetCholCov(void){return cholCov;}
	};

	
}