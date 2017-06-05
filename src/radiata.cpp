// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// radiata.cpp: Rcpp wrapper for SMC library -- TESTING Radiata pine example

#include <Rcpp.h>

#include "smctc.h"
#include "radiata.h"
#include "rngR.h"

#include <cstdio> 
#include <cstdlib>
#include <cstring>
#include <math.h>

#include <iostream>
#include <cmath>
//#include <gsl/gsl_randist.h>

using namespace std;

///The observations
//rad_obs * y;
//Rcpp::NumericMatrix y;
std::vector<rad_obs> yRAD;
double mean_x;

double var_alpha = 2000.0;
double var_beta = 120.0;
double var_phi = 0.05;

double a_prior = 3.0;
double b_prior = pow(2.0*300.0*300.0,-1.0);

double integrand_mean_alpha(const rad_state&, void*);
double integrand_mean_beta(const rad_state&, void*);
double integrand_mean_phi(const rad_state&, void*);
double integrand_var_alpha(const rad_state&, void*);
double integrand_var_beta(const rad_state&, void*);
double integrand_var_phi(const rad_state&, void*);

// pf() function callable from R via Rcpp:: essentially the same as main() from pf.cc 
// minor interface change to pass data down as matrix, rather than a filename
extern "C" SEXP radiataBS(SEXP dataS, SEXP partS) { 	

    long lIterates;

    try {

        //std::string filename = Rcpp::as<std::string>(fileS);
        unsigned long lNumber = Rcpp::as<unsigned long>(partS); // number of particles

        // Load observations -- or rather copy them in from R
        Rcpp::NumericMatrix dat = Rcpp::NumericMatrix(dataS); // so we expect a matrix
        lIterates = dat.nrow(); // number of temperatures after 0
        yRAD.reserve(lIterates);
		double sum_x = 0.0;
        for (long i = 0; i < lIterates; ++i) {
            yRAD[i].data_x = dat(i,0);
            yRAD[i].data_y = dat(i,1);
			sum_x += yRAD[i].data_x;
        }
		mean_x = 1.0 * sum_x/lIterates; 

        //Initialise and run the sampler
        smc::sampler<rad_state> Sampler(lNumber, SMC_HISTORY_NONE);  
        smc::moveset<rad_state> Moveset(fInitialiseRAD, fMoveRAD, NULL);

        Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();
		
        Rcpp::NumericVector alpham(lIterates), alphav(lIterates), betam(lIterates), betav(lIterates), phim(lIterates), phiv(lIterates), IntermNC(lIterates), myRatioNC(lIterates), myESS(lIterates);

        alpham(0) = Sampler.Integrate(integrand_mean_alpha, NULL);
        alphav(0) = Sampler.Integrate(integrand_var_alpha, (void*)&alpham(0));
        betam(0) = Sampler.Integrate(integrand_mean_beta, NULL);
        betav(0) = Sampler.Integrate(integrand_var_beta, (void*)&betam(0));
        phim(0) = Sampler.Integrate(integrand_mean_phi, NULL);
        phiv(0) = Sampler.Integrate(integrand_var_phi, (void*)&phim(0));
        IntermNC(0) = Sampler.GetLogNCStep();
        myRatioNC(0) = Sampler.GetLogNCPath();
        myESS(0) = Sampler.GetESS();
		

        for(int n=1; n < lIterates; ++n) {
            Sampler.Iterate();
      
            alpham(n) = Sampler.Integrate(integrand_mean_alpha, NULL);
            alphav(n) = Sampler.Integrate(integrand_var_alpha, (void*)&alpham(n));
            betam(n) = Sampler.Integrate(integrand_mean_beta, NULL);
            betav(n) = Sampler.Integrate(integrand_var_beta, (void*)&betam(n));
            phim(n) = Sampler.Integrate(integrand_mean_phi, NULL);
            phiv(n) = Sampler.Integrate(integrand_var_phi, (void*)&phim(n));
            IntermNC(n) = Sampler.GetLogNCStep();
            myRatioNC(n) = Sampler.GetLogNCPath();
			      myESS(n) = Sampler.GetESS();
        }

		//return Rcpp::List::create		   
        return Rcpp::DataFrame::create(Rcpp::Named("alpham") = alpham,
                                       Rcpp::Named("alphav") = alphav,
                                       Rcpp::Named("betam") = betam,
                                       Rcpp::Named("betav") = betav,
                                       Rcpp::Named("phim") = phim,
                                       Rcpp::Named("phiv") = phiv,
                                       Rcpp::Named("IntermNC") = IntermNC,
                                       Rcpp::Named("MyRatioNC") = myRatioNC,
                                       Rcpp::Named("MyESS") = myESS);
    }
    catch(smc::exception  e) {
        Rcpp::Rcout << e;       	// not cerr, as R doesn't like to mix with i/o 
        //exit(e.lCode);		// we're just called from R so we should not exit
    }
    return R_NilValue;          	// to provide a return 
}

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
double logWeight(long lTime, const rad_state & X)
{
		//double log_prior = -log(1000.0)- pow(X.alpha - 3000.0,2.0)/(2.0*1000.0*1000.0) -log(100.0)- pow(X.alpha - 185.0,2.0)/(2.0*100.0*100.0) + X.phi-1.0/b_prior/expl(X.phi)-X.phi*(a_prior+1.0);
		
		double mean_reg = X.alpha + X.beta*(yRAD[lTime].data_x - mean_x);
		double sigma = pow(expl(X.phi),0.5);
		double log_normpdf = -log(sigma) - pow(yRAD[lTime].data_y - mean_reg,2.0)/(2.0*sigma*sigma) -0.5*log(2.0*M_PI);
				
		return log_normpdf;
		
	
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<rad_state> fInitialiseRAD(smc::rng *pRng)
{
  rad_state value;
  
  // drawing from the prior
  value.alpha = pRng->Normal(3000.0,1000.0);
  value.beta = pRng->Normal(185.0,100.0);
  value.phi = log(pow(pRng->Gamma(3,pow(2.0*300.0*300.0,-1.0)),-1.0));
  
  return smc::particle<rad_state>(value,logWeight(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMoveRAD(long lTime, smc::particle<rad_state > & pFrom, smc::rng *pRng)
{
  rad_state * cv_to = pFrom.GetValuePointer();  
  
  // Some code which could be useful if adding an MCMC step later
  //cv_to->alpha += pRng->Normal(0,sqrt(var_alpha));
  //cv_to->beta += pRng->Normal(0,sqrt(var_beta));
  //cv_to->phi += pRng->Normal(0,sqrt(var_phi));
  pFrom.AddToLogWeight(logWeight(lTime, *cv_to));
}
