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


#include "smctc.h"
#include <RcppArmadillo.h>

namespace linGauss_BPF {
    arma::mat y; //data
    long lDim; //the dimension of the matrix
    long lLength; //the number of parameters (number of elements on lower diagonal of lDim*lDim matrix)
    arma::mat A;
    arma::rowvec theta_prop;
    
    void GenA(arma::mat & A, const arma::rowvec & a);
    
    void fInitialise(arma::vec & value, double & logweight, smc::nullParams & param);
    void fMove(long lTime, arma::vec & value, double & logweight, smc::nullParams & param);
}
