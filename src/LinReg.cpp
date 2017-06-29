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
    smc::sampler<rad_state> Sampler(lNumber, HistoryType::RAM);  
    smc::moveset<rad_state> Moveset(fInitialise, fMove, fMCMC);
    
    Sampler.SetResampleParams(ResampleType::SYSTEMATIC, 0.5);
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
	
	//Sampler.OstreamMCMCRecordToStream(Rcpp::Rcout);
	//Sampler.OstreamResamplingRecordToStream(Rcpp::Rcout);
	//Rcpp::Rcout << Sampler << std::endl;
	
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
///  \param value     The state to consider 
double logWeight(long lTime, const rad_state & value){
 
    double mean_reg = value.alpha + value.beta*(y.data_x(lTime) - mean_x);
    double sigma = pow(expl(value.phi),0.5);
    return -log(sigma) - pow(y.data_y(lTime) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);

}

///The function corresponding to the log posterior at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the index of the current distribution)
///  \param value     The state to consider 
double logPosterior(long lTime, const rad_state & value){
  
  double log_prior = -log(1000.0)- pow(value.alpha - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(value.beta - 185.0,2.0)/(2.0*100.0*100.0) + value.phi-1.0/b_prior/expl(value.phi) -value.phi*(a_prior+1.0);
  
  double sigma = pow(expl(value.phi),0.5);
  
  double log_normpdf;
  
  if (lTime==0){
    double mean_reg = value.alpha + value.beta*(y.data_x(0) - mean_x);
    log_normpdf = -log(sigma) - pow(y.data_y(0) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
  } else{
    arma::vec mean_reg = value.alpha*arma::ones(lTime) + value.beta*(y.data_x.rows(0,lTime-1) - mean_x*arma::ones(lTime));
    log_normpdf = arma::sum(-log(sigma)*arma::ones(lTime) - pow(y.data_y.rows(0,lTime-1) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI)*arma::ones(lTime));
  }
  
  return (log_normpdf + log_prior);
}

///A function to initialise a particle

/// \param pRng A pointer to the random number generator which is to be used
void fInitialise(smc::rng *pRng, rad_state & value, double & logweight)
{
  // drawing from the prior
    value.alpha = pRng->Normal(3000.0,1000.0);
    value.beta = pRng->Normal(185.0,100.0);
    value.phi = log(pow(pRng->Gamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
	
    double mean_reg = value.alpha + value.beta*(y.data_x(0) - mean_x);
    double sigma = pow(expl(value.phi),0.5);
    logweight = -log(sigma) - pow(y.data_y(0) - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The population to move.
///\param pRng  A random number generator.
void fMove(long lTime, rad_state & value, double & logweight, smc::rng *pRng)
{
  logweight += logWeight(lTime, value);
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The population to move.
///\param pRng  A random number generator.
int fMCMC(long lTime, rad_state & value, smc::rng *pRng)
{
  //Rcpp::Rcout << "lTime is " << lTime << std::endl;
  double MH_ratio;
  double dRand;
  int count = 0;
  
  rad_state value_prop;
  double logposterior_curr = logPosterior(lTime, value);
  double logposterior_prop;
  
    for (unsigned int j=0; j<10; j++){
		
      value_prop.alpha = pRng->Normal(value.alpha,std_alpha);
      value_prop.beta = pRng->Normal(value.beta,std_beta);
      value_prop.phi = pRng->Normal(value.phi,std_phi);
	  
      logposterior_prop = logPosterior(lTime, value_prop);
	  
      MH_ratio = exp(logposterior_prop - logposterior_curr);
      dRand = pRng->Uniform(0,1);
      
      if (MH_ratio>dRand){
        value = value_prop;
		logposterior_curr = logposterior_prop;
        count++;
      }
    }
  return count;
}
}

namespace std {
  /// Produce a human readable display of an the particle values using the standard stream operators

  /// \param os The output stream to which the display should be made.
  /// \param p  The particle which is to be displayed.
  //template <class rad_state>
  std::ostream & operator << (std::ostream & os, rad_state & value)
  {
	double alpha = value.alpha;
	double beta = value.beta;
	double phi = value.phi;
    os << "(" << alpha << ", " << beta << ", " << phi << ")";
    return os;
  }
}