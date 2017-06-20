// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// LinReg_LA.cpp: A simple example for estimating the parameters of a
// linear regression model using likelihood annealing SMC.
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
#include "LinReg_LA.h"
#include "rngR.h"
#include <RcppArmadillo.h>

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>

namespace LinReg_LA {	
const double std_alpha = 50.0;
const double std_beta = 11.0;
const double std_phi = 0.2;
const double a_prior = 3.0;
const double b_prior = pow(2.0*300.0*300.0,-1.0);
}
using namespace std;
using namespace LinReg_LA;


// LinRegPPBS() function callable from R via Rcpp:: 
// [[Rcpp::export]]
Rcpp::List LinRegLABS_cpp(arma::mat data, arma::vec intemps, unsigned long inlNumber) { 	
  
  
  try {
	temps = intemps;
    lNumber = inlNumber;
	
    lIterates = data.n_rows;
    y.data_x = data.col(0);
    y.data_y = data.col(1);
    mean_x = arma::sum(y.data_x)/lIterates;
    
    long lTemps = temps.n_rows;
    
    
    //Initialise and run the sampler
    smc::sampler<rad_state> Sampler(lNumber, SMC_HISTORY_RAM);
    smc::moveset<rad_state> Moveset(fInitialise, fMove, fMCMC);
    
    Sampler.SetResampleParams(SMC_RESAMPLE_SYSTEMATIC, 0.5);
    Sampler.SetMoveSet(Moveset);
	Sampler.Initialise();
    
    Rcpp::NumericVector alpham(lTemps), alphav(lTemps), betam(lTemps), betav(lTemps), phim(lTemps), phiv(lTemps), ESS(lTemps);
    
    alpham(0) = Sampler.Integrate(integrand_mean_alpha, NULL);
    alphav(0) = Sampler.Integrate(integrand_var_alpha, (void*)&alpham(0));
    betam(0) = Sampler.Integrate(integrand_mean_beta, NULL);
    betav(0) = Sampler.Integrate(integrand_var_beta, (void*)&betam(0));
    phim(0) = Sampler.Integrate(integrand_mean_phi, NULL);
    phiv(0) = Sampler.Integrate(integrand_var_phi, (void*)&phim(0));
    ESS(0) = Sampler.GetESS();
    
    
    for(int n=1; n < lTemps; ++n) {
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

namespace LinReg_LA {
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


///The function corresponding to the log likelihood at specified position
///  \param value     The state to consider 
double logLikelihood(const rad_state & value){
  
  double sigma = pow(expl(value.phi),0.5);
  arma::vec mean_reg = value.alpha*arma::ones(lIterates) + value.beta*(y.data_x - mean_x*arma::ones(lIterates));
  return arma::sum(-log(sigma)*arma::ones(lIterates) - pow(y.data_y - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI)*arma::ones(lIterates));	
  
}
///The function corresponding to the log prior at specified position
///  \param value     The state to consider 
double logPrior(const rad_state & value){
  return -log(1000.0)- pow(value.alpha - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(value.beta - 185.0,2.0)/(2.0*100.0*100.0) + value.phi-1.0/b_prior/expl(value.phi) -value.phi*(a_prior+1.0);
}

///A function to initialise a particle

/// \param pRng A pointer to the random number generator which is to be used
void fInitialise(smc::rng *pRng, rad_state & value, double & logweight)
{
	// drawing from the prior
    value.alpha = pRng->Normal(3000.0,1000.0);
    value.beta = pRng->Normal(185.0,100.0);
    value.phi = log(pow(pRng->Gamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
    logweight = temps(0)*logLikelihood(value);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param value The value of the particle being moved
///\param logweight The log weight of the particle being moved
///\param pRng  A random number generator.
void fMove(long lTime, rad_state & value, double & logweight, smc::rng *pRng)
{
	logweight += (temps(lTime) - temps(lTime-1))*logLikelihood(value);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param value The value of the particle being moved
///\param pRng  A random number generator.
int fMCMC(long lTime, rad_state & value, smc::rng *pRng)
{
  double MH_ratio;
  double dRand;
  int count = 0;
  rad_state value_prop;
  double loglike_curr = logLikelihood(value);
  double logprior_curr = logPrior(value);
  double loglike_prop;
  double logprior_prop;
    
    for (unsigned int j=0; j<10; j++){
      value_prop.alpha = pRng->Normal(value.alpha,std_alpha);
      value_prop.beta = pRng->Normal(value.beta,std_beta);
      value_prop.phi = pRng->Normal(value.phi,std_phi);
	  
	  loglike_prop = logLikelihood(value_prop);
	  logprior_prop = logPrior(value_prop);
      
      MH_ratio = exp(temps(lTime)*(loglike_prop - loglike_curr) + logprior_prop - logprior_curr);
      dRand = pRng->Uniform(0,1);
      
      if (MH_ratio>dRand){
        value = value_prop;
		loglike_curr = loglike_prop;
		logprior_curr = logprior_prop;
        count++;
      }
    }
  return count;
}
}