// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_fAPF.h: Example 5.2 of Guarniero, Johansen and Lee (2016):
// Estimating the transition matrix for a multidimensional linear Gaussian
// model. The fully adapted (auxiliary) particle filter of Pitt and Shepherd
// (1999) is used to estimate the likelihood for use in particle marginal
// Metropolis Hastings.
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

namespace linGauss_fAPF {
        
    arma::mat y; //data
    long lDim; //the dimension of the matrix
    long lLength; //the number of parameters (number of elements on lower diagonal of lDim*lDim matrix)
    arma::mat A;
    arma::rowvec theta_prop;
    long lIterates;
    
    double logMvnPdf(arma::vec obs, arma::vec mu, arma::mat cov);
    
    void GenA(arma::mat & A, const arma::rowvec & a);
    
    void fInitialise(arma::vec & X, double & logweight, smc::nullParams & param);
    void fMove(long lTime, arma::vec & X, double & logweight, smc::nullParams & param);
    
    arma::mat I, prep, SIG, cholSIG;
}
