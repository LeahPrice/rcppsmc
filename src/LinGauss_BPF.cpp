// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_BPF.cpp: Example 5.2 of Guarniero, Johansen and Lee (2016):
// Estimating the transition matrix for a multidimensional linear Gaussian
// model. A bootstrap particle filter is used to estimate the likelihood
// for use in particle marginal Metropolis Hastings.
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

#include "LinGauss_BPF.h"

using namespace std;
using namespace linGauss_BPF;


// linGauss_BPF_impl() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List linGauss_BPF_impl(arma::mat data, arma::rowvec initial, unsigned long lNumber, unsigned long lMCMCits) {    

    try {
        long lIterates = data.n_rows;
        lDim = data.n_cols;
        A.zeros(lDim,lDim);
        lLength = lDim*(lDim+1)/2;
        y = data;
        arma::mat theta(lMCMCits,lLength);
        
        I = arma::eye(lDim,lDim);
        
        arma::vec loglike = arma::zeros(lMCMCits);
        
        double loglike_prop;
        long accept_count = 0;
        
        //Initialise and run the sampler
        smc::sampler<arma::vec,smc::nullParams> Sampler(lNumber, HistoryType::NONE);
        smc::moveset<arma::vec,smc::nullParams> Moveset(fInitialise, fMove, NULL);
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 1.01*lNumber);
        Sampler.SetMoveSet(Moveset);
        
        theta_prop = initial; 
        theta.row(0) = theta_prop;
        GenA(A,theta_prop);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike(0) = Sampler.GetLogNCPath();
        
        Rcpp::NumericVector unifRands = Rcpp::runif(lMCMCits);
        
        double MH_ratio;
        for (unsigned int i = 1; i<lMCMCits; i++){
            // RW proposal for parameters
            theta_prop = theta.row(i-1);
            theta_prop((i-1)%15) += R::rnorm(0.0,sqrt(0.1));
            GenA(A,theta_prop);
            
            // Evaluating prior
            if (any(abs(theta_prop)>5)){
                theta.row(i) = theta.row(i-1);
                loglike(i) = loglike(i-1);
            } else{
                // Getting a particle filtering estimate of the log likelihood.
                Sampler.Initialise();
                Sampler.IterateUntil(lIterates-1);
                loglike_prop = Sampler.GetLogNCPath();
                
                MH_ratio = exp(loglike_prop - loglike(i-1));
                
                if (MH_ratio>unifRands(i-1)){
                    theta.row(i) = theta_prop;
                    loglike(i) = loglike_prop;
                    accept_count++;
                } else {
                    theta.row(i) = theta.row(i-1);
                    loglike(i) = loglike(i-1);
                }
            }
            
        }
        
        return Rcpp::List::create(Rcpp::Named("samples_theta") = theta,
        Rcpp::Named("loglike") = loglike);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;           // not cerr, as R doesn't like to mix with i/o
    }
    return R_NilValue;              // to provide a return 
}

namespace linGauss_BPF {
    
    /// Calculates the log of the multivariate normal density

    /// \param obs      The point at which the density is calculated
    /// \param mu       The mean vector
    /// \param cov      The covariance
    double logMvnPdf(arma::vec obs, arma::vec mu, arma::mat cov){
        double det, sign;
        arma::log_det(det, sign, cov);
        return ( -lDim*0.5*log(2.0*M_PI) - 0.5*det - 0.5*arma::as_scalar((obs-mu).t()*arma::inv(cov)*(obs-mu)));
    }
    
    /// Updates the matrix A based on the proposals

    /// \param A        A reference to the current A matrix
    /// \param a        A vector of the estimated lower triangular elements of A (by row)
    void GenA(arma::mat & A, const arma::rowvec & a)
    {
        // Starting with an upper triangular so that the elements are filled in the correct order.
        A.ones();
        A = trimatu(A);
        A.elem( find(A == 1) ) = a;
        A = A.t();
    }

    ///A function to initialise a particle

    /// \param X            A reference to the empty particle value
    /// \param logweight    A reference to the empty particle log weight
    /// \param params       A reference to the (null) parameter values
    void fInitialise(arma::vec & X, double & logweight, smc::nullParams & params)
    {
        X = Rcpp::as<arma::vec>(Rcpp::rnorm(lDim)); //mu_1
        logweight = logMvnPdf(y.row(0).t(),X,0.25*arma::eye(lDim,lDim)); //g_1
    }

    ///The proposal function.

    ///\param lTime     The sampler iteration.
    ///\param X         A reference to the current particle value
    ///\param params    A reference to the (null) parameter values
    void fMove(long lTime, arma::vec & X, double & logweight, smc::nullParams & params)
    {
        X = A*X + Rcpp::as<arma::vec>(Rcpp::rnorm(lDim)); //f_t
        logweight += logMvnPdf(y.row(lTime).t(),X,0.25*arma::eye(lDim,lDim)); //g_t
    }
}
