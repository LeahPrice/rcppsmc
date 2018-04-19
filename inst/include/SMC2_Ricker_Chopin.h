
#include "smctc.h"
#include <RcppArmadillo.h>

namespace SMC_Ricker {


long lIterates_outer;

class Ricker{
public:
  arma::vec theta;
  smc::sampler<double,arma::vec> * pf;
  double loglike;
  double logprior;
  
  // destructor
  ~Ricker(){
    if (pf)
      delete pf;
  }
  
  // basic constructor
  Ricker(){
    pf = NULL;
  }
  
  // // copy
  // Ricker(const Ricker & rFrom){
  // 	this->theta = rFrom.theta;
  // 	this->loglike = rFrom.loglike;
  // 	this->logprior = rFrom.logprior;
  // 	//*(this->pf) = *(rFrom.pf);
  // 	this->pf = new smc::sampler<double, arma::vec>(*(rFrom.pf));
  // }
  
  // assignment
  Ricker & operator=(const Ricker & rFrom){
    this->theta = rFrom.theta;
    this->loglike = rFrom.loglike;
    this->logprior = rFrom.logprior;
    
    if(! this->pf) {
      this->pf = new smc::sampler<double,arma::vec>(*(rFrom.pf));
    } else {
      *(this->pf) = *(rFrom.pf);
    }
    
    //*(this->pf) = *(rFrom.pf); // this also works
    
    // if(this->pf) delete this->pf;
    // this->pf = new smc::sampler<double,arma::vec>(*(rFrom.pf));
    
    return *this;
  }
  
};

arma::vec y;
unsigned long inN;

arma::mat boundaries(3,2);

double logPrior(const arma::vec & current);

void fInitialise_outer(Ricker & current, double & logweight, smc::nullParams & param);
void fMove_outer(long lTime, Ricker & current, double & logweight, smc::nullParams & param);
bool fMCMC_outer(long lTime, Ricker & current, double & logweight, smc::nullParams & param);

void fInitialise_inner(double & N, double & logweight, arma::vec & param);
void fMove_inner(long lTime, double & N, double & logweight, arma::vec & param);
smc::moveset<double,arma::vec> Moveset_inner(fInitialise_inner, fMove_inner, NULL);


}
