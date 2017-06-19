// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg.cpp: A simple example for estimating the parameters of a
// linear regression model using data annealing SMC.
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
#include "LinReg.h"
#include "rngR.h"
#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>

namespace LinReg {
const double std_alpha = 50.0;
const double std_beta = 11.0;
const double std_phi = 0.2;

const double a_prior = 3.0;
const double b_prior = pow(2.0*300.0*300.0,-1.0);
}

using namespace std;
using namespace LinReg;


// LinRegBS() function callable from R via Rcpp::
// [[Rcpp::export]]
Rcpp::List LinRegBS_cpp(arma::mat data, unsigned long inlNumber) { 	
  
  long lIterates;
  
  try {
    lNumber = inlNumber; // number of particles
    
    // Load observations -- or rather copy them in from R
    lIterates = data.n_rows;
    y.data_x = data.col(0);
    y.data_y = data.col(1);
    mean_x = arma::sum(y.data_x)/lIterates;
    
    //Initialise and run the sampler
    smc::sampler<rad_state> Sampler(lNumber, SMC_HISTORY_NONE);  
    smc::moveset<rad_state> Moveset(fInitialise, fMove, fMCMC);
    
    Sampler.SetResampleParams(SMC_RESAMPLE_SYSTEMATIC, 0.5);
    Sampler.SetMoveSet(Moveset);
    Sampler.Initialise();
    
    Rcpp::NumericVector alpham(lIterates), alphav(lIterates), betam(lIterates), betav(lIterates), phim(lIterates), phiv(lIterates), ESS(lIterates);
    
    alpham(0) = Sampler.Integrate(integrand_mean_alpha, NULL);
    alphav(0) = Sampler.Integrate(integrand_var_alpha, (void*)&alpham(0));
    betam(0) = Sampler.Integrate(integrand_mean_beta, NULL);
    betav(0) = Sampler.Integrate(integrand_var_beta, (void*)&betam(0));
    phim(0) = Sampler.Integrate(integrand_mean_phi, NULL);
    phiv(0) = Sampler.Integrate(integrand_var_phi, (void*)&phim(0));
    ESS(0) = Sampler.GetESS();
    
    
    for(int n=1; n < lIterates; ++n) {
      Sampler.Iterate();
      
      alpham(n) = Sampler.Integrate(integrand_mean_alpha, NULL);
      alphav(n) = Sampler.Integrate(integrand_var_alpha, (void*)&alpham(n));
      betam(n) = Sampler.Integrate(integrand_mean_beta, NULL);
      betav(n) = Sampler.Integrate(integrand_var_beta, (void*)&betam(n));
      phim(n) = Sampler.Integrate(integrand_mean_phi, NULL);
      phiv(n) = Sampler.Integrate(integrand_var_phi, (void*)&phim(n));
      ESS(n) = Sampler.GetESS();
    }
    
	double logNC = Sampler.GetLogNCPath();
	
    //return Rcpp::DataFrame::create		   
    return Rcpp::List::create(Rcpp::Named("alpham") = alpham,
                                   Rcpp::Named("alphav") = alphav,
                                   Rcpp::Named("betam") = betam,
                                   Rcpp::Named("betav") = betav,
                                   Rcpp::Named("phim") = phim,
                                   Rcpp::Named("phiv") = phiv,
                                   Rcpp::Named("ESS") = ESS,
								   Rcpp::Named("logNC") = logNC);
  }
  catch(smc::exception  e) {
    Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
    //exit(e.lCode);		// we're just called from R so we should not exit
  }
  return R_NilValue;          	// to provide a return 
}

namespace LinReg {
double integrand_mean_alpha(const rad_state& s, void *){ return s.alpha;}
double integrand_mean_beta(const rad_state& s, void *){ return s.beta;}
double integrand_mean_phi(const rad_state& s, void *){ return s.phi;}

double integrand_var_alpha(const rad_state& s, void* vmx){
  double* dmx = (double*)vmx;
  double d = (s.alpha - (*dmx));
  return d*d;
}

double integrand_var_beta(const rad_state& s, void* vmy){
  double* dmy = (double*)vmy;
  double d = (s.beta - (*dmy));
  return d*d;
}

double integrand_var_phi(const rad_state& s, void* vmx){
  double* dmx = (double*)vmx;
  double d = (s.phi - (*dmx));
  return d*d;
}

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
arma::vec logWeight(long lTime, const std::vector<rad_state> & X){
  arma::vec log_normpdf(lNumber);
  double mean_reg;
  double sigma;
  
  for (unsigned int i=0; i<lNumber; i++){
    mean_reg = X[i].alpha + X[i].beta*(y.data_x(lTime) - mean_x);
    sigma = pow(expl(X[i].phi),0.5);
    log_normpdf(i) = -log(sigma) - pow(y.data_y(lTime) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
  }
  
  return log_normpdf;
}

///The function corresponding to the log posterior at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param X     The state to consider 
double logPosterior(long lTime, const rad_state & X){
  
  double log_prior = -log(1000.0)- pow(X.alpha - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(X.beta - 185.0,2.0)/(2.0*100.0*100.0) + X.phi-1.0/b_prior/expl(X.phi) -X.phi*(a_prior+1.0);
  
  double sigma = pow(expl(X.phi),0.5);
  
  double log_normpdf;
  
  if (lTime==0){
    double mean_reg = X.alpha + X.beta*(y.data_x(0) - mean_x);
    log_normpdf = -log(sigma) - pow(y.data_y(0) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
  } else{
    arma::vec mean_reg = X.alpha*arma::ones(lTime) + X.beta*(y.data_x.rows(0,lTime-1) - mean_x*arma::ones(lTime));
    log_normpdf = arma::sum(-log(sigma)*arma::ones(lTime) - pow(y.data_y.rows(0,lTime-1) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI)*arma::ones(lTime));
  }
  
  return (log_normpdf + log_prior);
}

///A function to initialise a particle

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<rad_state> fInitialise(smc::rng *pRng)
{
  rad_state value;
  // drawing from the prior
    value.alpha = pRng->Normal(3000.0,1000.0);
    value.beta = pRng->Normal(185.0,100.0);
    value.phi = log(pow(pRng->Gamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
	
      double mean_reg = value.alpha + value.beta*(y.data_x(0) - mean_x);
    double sigma = pow(expl(value.phi),0.5);
    double log_normpdf = -log(sigma) - pow(y.data_y(0) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);

  return smc::particle<rad_state>(value,log_normpdf);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The population to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::population<rad_state > & pFrom, smc::rng *pRng)
{
  std::vector<rad_state> * cv_to = pFrom.GetValuePointer();  
  
  pFrom.AddToLogWeight(logWeight(lTime, *cv_to));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The population to move.
///\param pRng  A random number generator.
int fMCMC(long lTime, smc::population<rad_state > & pFrom, smc::rng *pRng)
{
  //Rcpp::Rcout << "lTime is " << lTime << std::endl;
  double MH_ratio;
  double dRand;
  int count = 0;
  
  rad_state * cv_to;
  rad_state cv_to_new;
  
  for (unsigned int i=0; i<lNumber; i++){
    
    for (unsigned int j=0; j<10; j++){
      cv_to = pFrom.GetValuePointerN(i);
      // Some code which could be useful if adding an MCMC step later
      cv_to_new.alpha = pRng->Normal(cv_to->alpha,std_alpha);
      cv_to_new.beta = pRng->Normal(cv_to->beta,std_beta);
      cv_to_new.phi = pRng->Normal(cv_to->phi,std_phi);
      
      MH_ratio = exp(logPosterior(lTime, cv_to_new) - logPosterior(lTime, *cv_to));
      dRand = pRng->Uniform(0,1);
      
      if (MH_ratio>dRand){
        pFrom.SetValueN(cv_to_new,i);
        count++;
      }
    }
  }
  return count;
}
}