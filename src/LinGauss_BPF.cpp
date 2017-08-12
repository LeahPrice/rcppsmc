// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// ADD DESCRIPTION HERE
//
//
//
//
//
//
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


// linGauss_BPF() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::DataFrame linGauss_BPF_impl(arma::mat data, arma::rowvec initial, unsigned long lNumber, unsigned long lMCMCits) {    

    try {
        long lIterates = data.n_rows;
        lDim = data.n_cols;
        A.zeros(lDim,lDim);
        lLength = lDim*(lDim+1)/2;
        y = data;
        arma::mat theta(lMCMCits,lLength);
        
        arma::vec loglike = arma::zeros(lMCMCits);
        
        double loglike_prop;
        long accept_count = 0;
        
        //Initialise and run the sampler
        smc::sampler<arma::vec,smc::nullParams> Sampler(lNumber, HistoryType::NONE);
        smc::moveset<arma::vec,smc::nullParams> Moveset(fInitialise, fMove, NULL);
        Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 1.01*lNumber);
        Sampler.SetMoveSet(Moveset);
        
        // Getting a particle filtering estimate of the log likelihood.
        theta_prop = initial; //theta_prop = {0.9,0.3,0.7,0.1,0.2,0.6,0.4,0.1,0.1,0.3,0.1,0.2,0.5,0.2,0};
        theta.row(0) = theta_prop;
        GenA(A,theta_prop);
        Sampler.Initialise();
        Sampler.IterateUntil(lIterates-1);
        loglike(0) = Sampler.GetLogNCPath();
        
        Rcpp::NumericVector unifRands = Rcpp::runif(lMCMCits);
        
        double MH_ratio;
        for (unsigned int i = 1; i<lMCMCits; i++){
            // RW proposal for parameters
            theta_prop = theta.row(i-1) + Rcpp::as<arma::rowvec>(Rcpp::rnorm(lLength,0,0.1));
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
        
        return Rcpp::DataFrame::create(Rcpp::Named("samples_theta") = theta,
        Rcpp::Named("loglike") = loglike);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;           // not cerr, as R doesn't like to mix with i/o
    }
    return R_NilValue;              // to provide a return 
}

namespace linGauss_BPF {
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
        X = Rcpp::as<arma::vec>(Rcpp::rnorm(lDim,0.0,1.0)); //mu_1
        logweight = arma::sum(-log(sqrt(0.25)) - pow(y.row(0).t() - X,2.0)/(2.0*0.25) -0.5*log(2.0*M_PI)); //g_1
    }

    ///The proposal function.

    ///\param lTime     The sampler iteration.
    ///\param X         A reference to the current particle value
    ///\param params    A reference to the (null) parameter values
    void fMove(long lTime, arma::vec & X, double & logweight, smc::nullParams & params)
    {
        X = A*X + Rcpp::as<arma::vec>(Rcpp::rnorm(lDim,0.0,1.0)); //f_t
        logweight += arma::sum(-log(sqrt(0.25)) - pow(y.row(lTime).t() - X,2.0)/(2.0*0.25) -0.5*log(2.0*M_PI)); //g_t
    }
}
