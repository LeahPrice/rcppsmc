// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// pfnonlinbs.cpp: Rcpp integration of SMC library -- PF Nonlinear Bootstrp
//
//    The principle source file for an implementation of the bootstrap
//    particle filter of "Novel approaches to nonlinear non-Gaussian
//    Bayesian state estimation", Gordon Salmond and Smith, 
//    IEE PROCEEDINGS-F 140(2):107-113, 1993
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012 - 2013  Dirk Eddelbuettel and Adam Johansen
// Copyright (C) 2017		  Dirk Eddelbuettel, Adam Johansen and Leah South
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
#include "pfnonlinbs.h"
#include <RcppArmadillo.h>

#include <cstdlib>
#include <cmath>

namespace nonlinbs {
const double std_x0 = 2;
const double var_x0 = std_x0 * std_x0;
const double std_x  = sqrt(10.0);
const double var_x  = std_x * std_x;
const double var_y  = 1.0;

const double scale_y = 1.0 / 20.0;

///The observations
arma::vec y;
}

using namespace std;
using namespace nonlinbs;

// [[Rcpp::export]]
Rcpp::List pfNonlinBS_cpp(arma::vec data, long inlNumber) {
  lNumber = inlNumber;
  
  y = data;
  long lIterates = y.n_rows;
  
  //Initialise and run the sampler 
  smc::sampler<double> Sampler(lNumber, HistoryType::NONE);  
  smc::moveset<double> Moveset(fInitialise, fMove, NULL);
  
  Sampler.SetResampleParams(ResampleType::MULTINOMIAL, 1.01 * lNumber);
  Sampler.SetMoveSet(Moveset);
  Sampler.Initialise();
  
  Rcpp::NumericVector resMean = Rcpp::NumericVector(lIterates);
  Rcpp::NumericVector resSD   = Rcpp::NumericVector(lIterates);
  for(int n=0 ; n < lIterates ; ++n) {
    if(n > 0)
      Sampler.Iterate();
    
    resMean(n) = Sampler.Integrate(integrand_mean_x,NULL);      
    resSD(n)  = sqrt(Sampler.Integrate(integrand_var_x, (void*)&resSD(n)));      
  }

  double logNC = Sampler.GetLogNCPath();
  
  return Rcpp::List::create(Rcpp::_["mean"] = resMean,
                            Rcpp::_["sd"] = resSD,
                            Rcpp::_["logNC"] = logNC);
}

namespace nonlinbs {
///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param value     The state to consider 
double logLikelihood(long lTime, const double & value) {
  return -0.5 * pow(y(int(lTime)) - value*value*scale_y,2) / var_y;;
}

///A function to initialise a particle

/// \param pRng A pointer to the random number generator which is to be used
void fInitialise(smc::rng *pRng, double & value, double & logweight) {
	value = pRng->Normal(0,std_x0);
	logweight = -0.5 * pow(y(int(0)) - value*value*scale_y,2) / var_y;
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The population to move.
///\param pRng  A random number generator.
void fMove(long lTime, double & value, double & logweight, smc::rng *pRng) {
  value = 0.5 * value + 25.0*value / (1.0 + value * value) + 8.0 * cos(1.2  * ( lTime)) + pRng->Normal(0.0,std_x);
  logweight += logLikelihood(lTime, value);
}


double integrand_mean_x(const double& x, void *) {
  return x;
}

double integrand_var_x(const double& x, void* vmx) {
  double* dmx = (double*)vmx;
  double d = (x - (*dmx));
  return d*d;
}
}
