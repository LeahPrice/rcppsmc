//   SMCTC: sampler.hh  
//
//   Copyright Adam Johansen, 2008-2009.
// 
//   This file is part of SMCTC.
//
//   SMCTC is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   SMCTC is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with SMCTC.  If not, see <http://www.gnu.org/licenses/>.

//! \file
//! \brief Defines the overall sampler object.
//!
//! This file defines the smc::sampler class which is used to implement entire particle systems.

#ifndef __SMC_SAMPLER_HH

#define __SMC_SAMPLER_HH 1.0

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <RcppArmadillo.h>
#include <typeinfo>

#include "rngR.h"
#include "history.h"
#include "moveset.h"
#include "particle.h"
#include "smc-exception.h"

///Specifiers for various resampling algorithms:
enum ResampleType { SMC_RESAMPLE_MULTINOMIAL = 0, 
                    SMC_RESAMPLE_RESIDUAL, 
                    SMC_RESAMPLE_STRATIFIED, 
                    SMC_RESAMPLE_SYSTEMATIC };

///Storage types for the history of the particle system. 
enum HistoryType { SMC_HISTORY_NONE = 0, 
                   SMC_HISTORY_RAM };	
				   
namespace smc {

/// A template class for an interacting particle system suitable for SMC sampling
template <class Space> 
class sampler
{
protected:
  ///A random number generator.
  rng* pRng;
  
  ///Number of particles in the system.
  long N;
  ///The current evolution time of the system.
  long T;
  
  ///The resampling mode which is to be employed.
  ResampleType rtResampleMode;
  ///The effective sample size at which resampling should be used.
  double dResampleThreshold;
  ///Structure used internally for resampling.
  arma::vec dRSWeights;
  ///Structure used internally for resampling.
  arma::Col<unsigned int> uRSCount;
  ///Structure used internally for resampling.
  arma::Col<unsigned int> uRSIndices;
  
  ///The particles within the system.
  particle<Space> pParticles;
  ///The set of moves available.
  moveset<Space> Moves;
  
  ///The number of MCMC moves which have been accepted during this iteration
  int nAccepted;
  ///A flag which tracks whether the ensemble was resampled during this iteration
  int nResampled;
  
  ///A mode flag which indicates whether historical information is stored
  HistoryType htHistoryMode;
  ///The historical process associated with the particle system.
  history<particle<Space> > History;
  ///An estimate of the overall ratio of normalising constants
  double dlogNCPath;
  ///An estimate of the latest iteration's ratio of normalising constants
  double dlogNCIt;
  
public:
  ///Create an particle system containing lSize uninitialised particles with the specified mode.
  sampler(long lSize, HistoryType htHistoryMode);
  ///Create an particle system constaining lSize uninitialised particles with the specified mode and random number generator.
  // -- no GSL  sampler(long lSize, HistoryType htHistoryMode, const gsl_rng_type* rngType, unsigned long nSeed);
  ///Dispose of a sampler.
  ~sampler();
  ///Calculates and Returns the Effective Sample Size.
  double GetESS(void) const;
  ///Returns a pointer to the History of the particle system
  const history<particle<Space> > * GetHistory(void) const { return &History; }																												   
  ///Returns the current estimate of the log normalising constant ratio over the entire path
  double GetLogNCPath(void) const { return dlogNCPath; }
  ///Returns the current estimate of the log normalising constant ratio over the last step
  double GetLogNCStep(void) const { return dlogNCIt; }
  ///Returns the current estimate of the normalising constant ratio over the entire path
  double GetNCPath(void) const { return exp(dlogNCPath); }
  ///Returns the current estimate of the normalising constant ratio over the last step
  double GetNCStep(void) const { return exp(dlogNCIt); }
  ///Returns the number of particles within the system.
  long GetNumber(void) const {return N;}
  ///Return the value of particles
  const std::vector<Space> &  GetParticleValue(void) { return pParticles.GetValue(); }
  ///Return the logarithmic unnormalized weights of particles
  arma::vec GetParticleLogWeight(void) { return pParticles.GetLogWeight(); }
  ///Return the unnormalized weights of particls
  arma::vec GetParticleWeight(void) { return pParticles.GetWeight(); }  
  ///Return the unnormalized weights of particls
  double GetParticleWeightN(int n) { return pParticles.GetWeightN(n); }  
  ///Returns the current evolution time of the system.
  long GetTime(void) const {return T;}
  ///Initialise the sampler and its constituent particles.
  void Initialise(void);
  ///Integrate the supplied function with respect to the current particle set.
  double Integrate(double(*pIntegrand)(const Space &,void*), void* pAuxiliary);
  ///Integrate the supplied function over the path path using the supplied width function.
  double IntegratePathSampling(double (*pIntegrand)(long,const particle<Space>&,void*), double (*pWidth)(long,void*), void* pAuxiliary);
  ///Perform one iteration of the simulation algorithm.
  void Iterate(void);
  ///Cancel one iteration of the simulation algorithm.
  void IterateBack(void);
  ///Perform one iteration of the simulation algorithm and return the resulting ess
  double IterateEss(void);
  ///Perform iterations until the specified evolution time is reached
  void IterateUntil(long lTerminate);
  ///Move the particle set by proposing an applying an appropriate move to each particle.
  void MoveParticles(void);
  ///Resample the particle set using the specified resmpling scheme.
  void Resample(ResampleType lMode);
  ///Sets the entire moveset to the one which is supplied
  void SetMoveSet(moveset<Space>& pNewMoveset) {Moves = pNewMoveset;}
  ///Set Resampling Parameters
  void SetResampleParams(ResampleType rtMode, double dThreshold);
  
private:
  ///Duplication of smc::sampler is not currently permitted.
  sampler(const sampler<Space> & sFrom);
  ///Duplication of smc::sampler is not currently permitted.
  sampler<Space> & operator=(const sampler<Space> & sFrom);
  
protected:
  ///Returns the crude normalising constant ratio estimate implied by the weights.
  double CalcLogNC(void) const;
};


/// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
/// store the particle set and to initialise a random number generator.
///
/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
/// \param htHM The history mode to use: set this to SMC_HISTORY_RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
/// \tparam Space The class used to represent a point in the sample space.
template <class Space>
sampler<Space>::sampler(long lSize, HistoryType htHM)
{pRng = new rng();
  N = lSize;
  
   //particle<Space>* pParticles; 
  //pParticles = new particle<Space>;
  
  // //Allocate some storage for internal workspaces
  // dRSWeights = new double[N];
  // ///Structure used internally for resampling.
  // uRSCount  = new unsigned[N];
  // ///Structure used internally for resampling.
  // uRSIndices = new unsigned[N];
  
  //Some workable defaults.	
  htHistoryMode = htHM;
  rtResampleMode = SMC_RESAMPLE_STRATIFIED;
  dResampleThreshold = 0.5 * N;
}

template <class Space>
sampler<Space>::~sampler()
{
  delete pRng;
  
  //if (pParticles)
    //delete pParticles;
  //if (dRSWeights)
  //  delete [] dRSWeights;
  //if (uRSCount)
  //  delete [] uRSCount;
  //if (uRSIndices)
  //  delete [] uRSIndices;
}


template <class Space>
double sampler<Space>::GetESS(void) const
{
	double sum = arma::sum(exp(pParticles.GetLogWeight()));
	double sumsq = arma::sum(exp(2.0*pParticles.GetLogWeight()));
  
  return expl(-log(sumsq) + 2*log(sum));
}

template <class Space>
double sampler<Space>::CalcLogNC(void) const
{
	double dMaxWeight = arma::max(pParticles.GetLogWeight());
	double sum = arma::sum(exp(pParticles.GetLogWeight() - dMaxWeight*arma::ones(N)));
  
  return (dMaxWeight + log(sum));
}


/// At present this function resets the system evolution time to 0 and calls the moveset initialisor to assign each
/// particle in the ensemble.
///
/// Note that the initialisation function must be specified before calling this function.
template <class Space>
void sampler<Space>::Initialise(void)
{
  T = 0;
  dlogNCIt = 0;
  dlogNCPath = 0;
  
  
  // particle<Space> myParticles; 
  // myParticles = Moves.DoInit(pRng);
  // pParticles = myParticles;
  pParticles = Moves.DoInit(pRng);

  //Scaling weights by 1/N (mostly for evidence computation)
	
	pParticles.SetLogWeight(pParticles.GetLogWeight() - log(static_cast<double>(N))*arma::ones(N));
  
   if(htHistoryMode != SMC_HISTORY_NONE) {
    while(History.Pop()!=NULL);
    nResampled = 0;
    History.Push(N, pParticles, 0, historyflags(nResampled));
  }    
  
  //Estimate the normalising constant
  dlogNCIt = CalcLogNC();
  dlogNCPath += dlogNCIt;
  
  //Normalise the weights to sensible values....
  pParticles.SetLogWeight(pParticles.GetLogWeight() - dlogNCIt*arma::ones(N));
  
  //Check if the ESS is below some reasonable threshold and resample if necessary.
  double ESS = GetESS();
  if(ESS < dResampleThreshold) {
  nResampled = 1;
  Resample(rtResampleMode);
  }
  else {
  nResampled = 0;
  }
  //A possible MCMC step could be included here.
    if(Moves.DoMCMC(0,pParticles, pRng))
      nAccepted++;
  
  return;
}

/// This function returns the result of integrating the supplied function under the empirical measure associated with the
/// particle set at the present time. The final argument of the integrand function is a pointer which will be supplied
/// with pAuxiliary to allow for arbitrary additional information to be passed to the function being integrated.
///
/// \param pIntegrand The function to integrate with respect to the particle set 
/// \param pAuxiliary A pointer to any auxiliary data which should be passed to the function

template <class Space>
double sampler<Space>::Integrate(double(*pIntegrand)(const Space&,void*), void * pAuxiliary)
{
  long double rValue = 0;
  long double wSum = 0;
  for(int i =0; i < N; i++)
  {
    rValue += expl(pParticles.GetLogWeightN(i)) * pIntegrand(pParticles.GetValueN(i), pAuxiliary);
    wSum  += expl(pParticles.GetLogWeightN(i));
  }
  
  rValue /= wSum;
  return (double)rValue;
}

/// This function is intended to be used to estimate integrals of the sort which must be evaluated to determine the
/// normalising constant of a distribution obtain using a sequence of potential functions proportional to densities with respect
/// to the initial distribution to define a sequence of distributions leading up to the terminal, interesting distribution.
///
/// In this context, the particle set at each time is used to make an estimate of the path sampling integrand, and a
/// trapezoidal integration is then performed to obtain an estimate of the path sampling integral which is the natural logarithm
/// of the ratio of normalising densities.
///
/// \param pIntegrand  The quantity which we wish to integrate at each time
/// \param pWidth      A pointer to a function which specifies the width of each 

template <class Space>
double sampler<Space>::IntegratePathSampling(double (*pIntegrand)(long,const particle<Space> &,void*), double (*pWidth)(long,void*), void* pAuxiliary)
{
  if(htHistoryMode == SMC_HISTORY_NONE)
    throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "The path sampling integral cannot be computed as the history of the system was not stored.");
  
  History.Push(N, pParticles, nAccepted, historyflags(nResampled));
  double dRes = History.IntegratePathSampling(pIntegrand, pWidth, pAuxiliary);
  History.Pop();
  return dRes;
}

/// The iterate function:
///         -# appends the current particle set to the history if desired
///          -# moves the current particle set
///         -# checks the effective sample size and resamples if necessary
///         -# performs a mcmc step if required
///         -# increments the current evolution time
template <class Space>
void sampler<Space>::Iterate(void)
{
  IterateEss();
  return;
}

template <class Space>
void sampler<Space>::IterateBack(void)
{
  if(htHistoryMode == SMC_HISTORY_NONE)
    throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "An attempt to undo an iteration was made; unforunately, the system history has not been stored.");
  
  //History.Pop(&N, pParticles, &nAccepted, NULL);
  History.Pop(&N, &pParticles, &nAccepted, NULL);
  T--;
  return;
}

template <class Space>
double sampler<Space>::IterateEss(void)
{
  //Initially, the current particle set should be appended to the historical process.
    
  nAccepted = 0;
  
  //Move the particle set.
  MoveParticles();
  
  //Estimate the normalising constant
  dlogNCIt = CalcLogNC();
  dlogNCPath += dlogNCIt;
		
  //Normalise the weights
  pParticles.SetLogWeight(pParticles.GetLogWeight()  - dlogNCIt*arma::ones(N));
  
  //Check if the ESS is below some reasonable threshold and resample if necessary.
  //A mechanism for setting this threshold is required.
  double ESS = GetESS();
  if(ESS < dResampleThreshold) {
    nResampled = 1;
    Resample(rtResampleMode);
  }
  else
    nResampled = 0;
  //A possible MCMC step could be included here.
    if(Moves.DoMCMC(T+1,pParticles, pRng))
      nAccepted++;
  // Increment the evolution time.
  T++;
  
  return ESS;
}

template <class Space>
void sampler<Space>::IterateUntil(long lTerminate)
{
  while(T < lTerminate)
    Iterate();
}

template <class Space>
void sampler<Space>::MoveParticles(void)
{
    Moves.DoMove(T+1,pParticles, pRng);
    //  pParticles[i].Set(pNew.value, pNew.logweight);
}

template <class Space>
void sampler<Space>::Resample(ResampleType lMode)
{
  //Resampling is done in place.
  unsigned uMultinomialCount;
  
  //First obtain a count of the number of children each particle has.
  switch(lMode) {
  case SMC_RESAMPLE_MULTINOMIAL:
    //Sample from a suitable multinomial vector
    dRSWeights = pParticles.GetWeight();
	//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
    //pRng->Multinomial(N,N,dRSWeights,uRSCount.memptr());
    uRSCount = pRng->Multinomial(N,N,dRSWeights);
    break;
    
	
  case SMC_RESAMPLE_RESIDUAL:
    dRSWeights = exp(log((double)N)*arma::ones(N) + pParticles.GetLogWeight() - CalcLogNC()*arma::ones(N));
	
	uRSIndices = arma::zeros<arma::Col<unsigned int> >((int)N);
	//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
	
	for(int i = 0; i < N; ++i)
		uRSIndices(i) = unsigned(floor(dRSWeights(i))); //Reuse temporary storage.
	dRSWeights = dRSWeights - uRSIndices;
	uMultinomialCount = N - arma::sum(uRSIndices);
	
	//pRng->Multinomial(uMultinomialCount,N,dRSWeights,uRSCount.memptr());
	uRSCount = pRng->Multinomial(uMultinomialCount,N,dRSWeights);
	uRSCount += uRSIndices;
	
	break;
	
  case SMC_RESAMPLE_STRATIFIED:
  default:
  {
    // Procedure for stratified sampling
    //Generate a random number between 0 and 1/N times the sum of the weights
    double dRand = pRng->Uniform(0,1.0 / ((double)N));
    
    arma::vec dWeightCumulative = arma::cumsum(exp(pParticles.GetLogWeight() - CalcLogNC()));
	
    int k = 0;
	int j = 0;
	uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
    while(j < N) {
      while((dWeightCumulative(k) - dRand) > ((double)j)/((double)N) && j < N) {
        uRSCount(k)++;
        j++;
        dRand = pRng->Uniform(0,1.0 / ((double)N));
      }
      k++;
    }
    break;
  }
  
  case SMC_RESAMPLE_SYSTEMATIC:
  {
    // Procedure for stratified sampling but with a common RV for each stratum
    //Generate a random number between 0 and 1/N times the sum of the weights
    double dRand = pRng->Uniform(0,1.0 / ((double)N));
    
    int j = 0, k = 0;
	uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
    arma::vec dWeightCumulative = arma::cumsum(exp(pParticles.GetLogWeight() - CalcLogNC()*arma::ones(N)));
    while(k < N) {
      while((dWeightCumulative(k) - dRand) > ((double)j)/((double)N) && j < N) {
        uRSCount(k)++;
        j++;
      }
      k++;
    }
    break;
   }
   
  }
  
  uRSIndices = arma::zeros<arma::Col<unsigned int> >((int)N);
  //Map count to indices to allow in-place resampling
  for (int i=0, j=0; i<N; ++i) {
    if (uRSCount(i)>0) {
      uRSIndices(i) = i;
      while (uRSCount(i)>1) {
        while (uRSCount(j)>0) ++j; // find next free spot
        uRSIndices(j++) = i; // assign index
        --uRSCount(i); // decrement number of remaining offsprings
      }
    }
  }
  Space Current;
  //Perform the replication of the chosen.
  for(int i = 0; i < N ; ++i) {
	if(uRSIndices(i) != static_cast<unsigned int>(i)){
	Current = pParticles.GetValueN((int)uRSIndices(i));
    pParticles.SetValueN(Current,i);
	}
  }
  //Set equal normalised weights
  pParticles.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
} 				
	  
/// This function configures the resampling parameters, allowing the specification of both the resampling
/// mode and the threshold at which resampling is used.
///
/// \param rtMode The resampling mode to be used.
/// \param dThreshold The threshold at which resampling is deemed necesary.
///
/// The rtMode parameter should be set to one of the following:
/// -# SMC_RESAMPLE_MULTINOMIAL to use multinomial resampling  
/// -# SMC_RESAMPLE_RESIDUAL to use residual resampling
/// -# SMC_RESAMPLE_STRATIFIED to use stratified resampling
/// -# SMC_RESAMPLE_SYSTEMATIC to use systematic resampling
///
/// The dThreshold parameter can be set to a value in the range [0,1) corresponding to a fraction of the size of
/// the particle set or it may be set to an integer corresponding to an actual effective sample size.

template <class Space>
void sampler<Space>::SetResampleParams(ResampleType rtMode, double dThreshold)
{
  rtResampleMode = rtMode;
  if(dThreshold < 1)
    dResampleThreshold = dThreshold * N;
  else
    dResampleThreshold = dThreshold;
}

}
#endif
