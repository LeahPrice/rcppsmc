#include "SMC2_Ricker_Chopin.h"
#include <cmath>

using namespace std;
using namespace SMC_Ricker;

// [[Rcpp::export]]
Rcpp::List SMC2_Ricker(arma::vec data, unsigned long num_outer, unsigned long num_inner) {
  
  try {
    inN = num_inner;
    y = data;
    
    lIterates_outer = data.n_rows;
    
    boundaries(0,0) = 2; boundaries(0,1) = 5;
    boundaries(1,0) = 1.61; boundaries(1,1) =3;
    boundaries(2,0) = -3; boundaries(2,1) = -0.22;
    
    // The outer sampler
    smc::sampler<Ricker,smc::nullParams> Sampler_outer(num_outer, HistoryType::NONE);
    smc::moveset<Ricker,smc::nullParams> Moveset_outer(fInitialise_outer, fMove_outer, fMCMC_outer);
    Sampler_outer.SetResampleParams(ResampleType::MULTINOMIAL, 0.5);
    Sampler_outer.SetMcmcRepeats(1);
    //Sampler_outer.SetMcmcRepeats(10);
    Sampler_outer.SetMoveSet(Moveset_outer);
    
    // Running the main sampler
    Sampler_outer.Initialise();
    //Sampler_outer.IterateUntil(lIterates_outer-1);
    for(int n=0; n < lIterates_outer; ++n){
      Sampler_outer.Iterate();
      
      Rcpp::Rcout << "at observation " << n+1 << " the log evidence is " << Sampler_outer.GetLogNCPath() << std::endl;
    }
    
    arma::mat theta(num_outer,3);
    arma::vec weights = Sampler_outer.GetParticleWeight();
    
    for (unsigned int i = 0; i<num_outer; i++){
      theta.row(i) = Sampler_outer.GetParticleValueN(i).theta.t();
    }
    
    double logNC = Sampler_outer.GetLogNCPath();
    
    return Rcpp::List::create(Rcpp::Named("samples") = theta,
                              Rcpp::Named("weights") = weights,
                              Rcpp::Named("logZ") = logNC);
  }
  catch(smc::exception  e) {
    Rcpp::Rcout << e;
  }
  return R_NilValue;            // to provide a return
}

namespace SMC_Ricker {

double logPrior(const arma::vec & current)
{
  if ( (sum(current>=boundaries.col(0))==3) && (sum(current<=boundaries.col(1))==3) )
    return 0.0;
  else
    return -std::numeric_limits<double>::infinity();
}

void fInitialise_outer(Ricker & current, double & logweight, smc::nullParams & param)
{
  current.theta = boundaries.col(0) + (boundaries.col(1) - boundaries.col(0)) % Rcpp::as<arma::vec>(Rcpp::runif(3));
  current.logprior = logPrior(current.theta);
  current.loglike = 0;
  
  // The inner sampler
  current.pf = new smc::sampler<double,arma::vec>(inN, HistoryType::NONE);
  current.pf->SetResampleParams(ResampleType::MULTINOMIAL, inN);
  current.pf->SetMoveSet(Moveset_inner);
  current.pf->SetAlgParam(current.theta);
  current.pf->Initialise();
  
  logweight = current.logprior;
}


void fMove_outer(long lTime, Ricker & current, double & logweight, smc::nullParams & param)
{
  double before = current.pf->GetLogNCPath();
  current.pf->Iterate();
  current.loglike = current.pf->GetLogNCPath();
  logweight += current.loglike - before;
}


bool fMCMC_outer(long lTime, Ricker & current, double & logweight, smc::nullParams & param)
{
  
  arma::vec theta_prop = current.theta + 0.25*Rcpp::as<arma::vec>(Rcpp::rnorm(3)); // + cholCovRW*Rcpp::as<arma::vec>(Rcpp::rnorm(3));
  double logprior_prop = logPrior(theta_prop);
  if (isinf(logprior_prop)==false){
    smc::sampler<double,arma::vec> pf_prop(inN, HistoryType::NONE);
    pf_prop.SetResampleParams(ResampleType::MULTINOMIAL, inN);
    pf_prop.SetMoveSet(Moveset_inner);
    pf_prop.SetAlgParam(theta_prop);
    pf_prop.Initialise();
    if (lTime>0){
      pf_prop.IterateUntil(lTime);
    }
    
    double loglike_prop = pf_prop.GetLogNCPath();
    double MH_ratio = exp(loglike_prop - current.loglike + logprior_prop - current.logprior);
    
    if (MH_ratio>R::runif(0,1)){
      current.theta = theta_prop;
      current.loglike = loglike_prop;
      current.logprior = logprior_prop;
      *current.pf = pf_prop;
      return TRUE;
    }
  }
  return FALSE;
}



void fInitialise_inner(double & N, double & logweight, arma::vec & param)
{
  N = 1;
  logweight = 0;
}

void fMove_inner(long lTime, double & N, double & logweight, arma::vec & param)
{
  double r = exp(param(0));
  double phi = exp(param(1));
  double sigmae = exp(param(2));
  N = r*N*exp(-N+sigmae*R::rnorm(0,1));
  double mu = phi*N; //lambda for poisson
  logweight += -lgamma(y(lTime-1)+1) - mu + y(lTime-1)*log(mu);
}

}
