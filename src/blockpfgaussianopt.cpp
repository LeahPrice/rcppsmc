// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// blockpfgaussianopt.cpp: Rcpp integration of SMC library -- Block PF Gaussian
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2012         Dirk Eddelbuettel and Adam Johansen
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
#include "blockpfgaussianopt.h"
#include "rngR.h"

#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;
using namespace BSPFG;

extern "C" SEXP blockpfGaussianOpt(SEXP dataS, SEXP partS, SEXP lagS)
{
  long lIterates;
  lNumber = Rcpp::as<long>(partS);
  lLag = Rcpp::as<long>(lagS);
  
  y = Rcpp::NumericVector(dataS);
  lIterates = y.size();
  
  //Initialise and run the sampler
  smc::sampler<vector<double> > Sampler(lNumber, SMC_HISTORY_NONE);  
  smc::moveset<vector<double> > Moveset(fInitialise, fMove, NULL);
  
  Sampler.SetResampleParams(SMC_RESAMPLE_SYSTEMATIC, 0.5);
  Sampler.SetMoveSet(Moveset);
  
  Sampler.Initialise();
  Sampler.IterateUntil(lIterates - 1);
  
  //Generate results
  Rcpp::NumericMatrix resValues = Rcpp::NumericMatrix(lNumber,lIterates);
  arma::vec resWeights(lNumber);
  double logNC = Sampler.GetLogNCPath();
  
  std::vector<vector<double> > pValue = Sampler.GetParticleValue();
  for(int i = 0; i < lNumber; ++i) 
  {
    for(int j = 0; j < lIterates; ++j) {
      resValues(i,j) = pValue.at(i).at(j);
    }
  }
  resWeights = Sampler.GetParticleWeight();
  
  return Rcpp::List::create(Rcpp::_["weight"] = resWeights, Rcpp::_["values"] = resValues, Rcpp::_["logNC"] = logNC);
}

using namespace std;

namespace BSPFG {

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<vector<double> > fInitialise(smc::rng *pRng)
{
  std::vector<vector<double> > value(lNumber);
  for (int k=0; k<lNumber; k++){
    value[k].push_back(pRng->Normal(0.5 * y[0],1.0/sqrt(2.0)));
  }
  
  return smc::particle<vector<double> >(value,arma::ones(lNumber));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<vector<double> > & pFrom, smc::rng *pRng)
{
  std::vector<vector<double> > * cv_to = pFrom.GetValuePointer();
  arma::vec logtarget(lNumber);
  if(lTime == 1) {
    for (int k=0; k<lNumber; k++){
      cv_to->at(k).push_back((cv_to->at(k).at(lTime-1) + y[int(lTime)])/2.0 + pRng->Normal(0.0,1.0/sqrt(2.0)));
      logtarget(k) = -0.25*(y[int(lTime)] - cv_to->at(k).at(lTime-1))*(y[int(lTime)]-cv_to->at(k).at(lTime-1));
    }
    pFrom.AddToLogWeight(logtarget);
    
    return;
  }
  
  long lag = min(lTime,lLag);
  
  //These structures should really be made static 
  std::vector<double> mu(lag+1);
  std::vector<double> sigma(lag+1);
  std::vector<double> sigmah(lag+1);
  std::vector<double> mub(lag+1);
  
  for (int k=0; k<lNumber; k++){
    // Forward filtering
    mu[0] = cv_to->at(k).at(lTime-lag);
    sigma[0] = 0;
    for(int i = 1; i <= lag; ++i)
    {
      sigmah[i] = sigma[i-1] + 1;
      
      mu[i] = (sigmah[i] * y[int(lTime-lag+i)] +  mu[i-1]) / (sigmah[i] + 1);
      sigma[i] = sigmah[i] / (sigmah[i] + 1);
    }
    // Backward smoothing
    mub[lag] = mu[lag];
    cv_to->at(k).push_back(pRng->Normal(mub[lag],sqrt(sigma[lag])));
    for(int i = lag-1; i; --i)
    {
      mub[i] = (sigma[i]*cv_to->at(k).at(lTime-lag+i+1) + mu[i]) / (sigma[i]+1);
      cv_to->at(k).at(lTime-lag+i) = pRng->Normal(mub[i],sqrt(sigma[lag]/(sigma[lag] + 1)));
    }
    logtarget(k) = -0.5 * pow(y[int(lTime)] - mu[lag-1],2.0) / (sigmah[lag]+1) ;
  }
  
  // Importance weighting
  pFrom.AddToLogWeight(logtarget);
  
}
}
