// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
//
// sampler.h: Rcpp integration of SMC library -- sampler object 
//
// Copyright (C) 2008 - 2009  Adam Johansen
// Copyright (C) 2017		  Adam Johansen, Dirk Eddelbuettel and Leah South
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

#include "history.h"
#include "moveset.h"
#include "algParam.h"
#include "population.h"
#include "smc-exception.h"
#include "staticModelAdapt.h"

///Specifiers for various resampling algorithms:
namespace ResampleType
{
	enum Enum {	MULTINOMIAL = 0, 
				RESIDUAL, 
				STRATIFIED, 
				SYSTEMATIC };
}

///Specifiers for various path sampling methods:
namespace PathSamplingType
{
	enum Enum {	TRAPEZOID2 = 0, 
				TRAPEZOID1, 
				RECTANGLE};
}


///Storage types for the history of the particle system. 
namespace HistoryType
{
	enum Enum {	NONE = 0, 
				RAM};
}

///Class for additional algorithm parameters
class nullParams{
	
};

namespace smc {

	/// A template class for an interacting particle system suitable for SMC sampling
	template <class Space, class Params = nullParams> 
	class sampler
	{
	protected:
		///Number of particles in the system.
		long N;
		///The current evolution time of the system.
		long T;
		
		///The resampling mode which is to be employed.
		ResampleType::Enum rtResampleMode;
		///The effective sample size at which resampling should be used.
		double dResampleThreshold;

		///Structure used internally for resampling.
		arma::vec dRSWeights;
		///Structure used internally for resampling.
		arma::Col<unsigned int> uRSCount;
		///Structure used internally for resampling.
		arma::Col<unsigned int> uRSIndices;

		///The particles within the system.
		population<Space> pPopulation;
		///The set of moves available.
		moveset<Space> Moves;
		///An object for tracking and adapting additional algorithm parameters.
		algParam<Space,Params>* pAlgParams;
		///A flag to track whether the adaptation object needs to be included in this destructor.
		bool pAlgBelong;

		///A mode flag which indicates whether historical information is stored.
		HistoryType::Enum htHistoryMode;
		///The historical process associated with the particle system.
		std::vector<historyelement< population<Space> > > History;
		
		///The number of MCMC moves which have been accepted during this iteration.
		int nAccepted;
		///A flag which tracks whether the ensemble was resampled during this iteration.
		int nResampled;
		///An estimate of the log normalising constant ratio over the entire path.
		double dlogNCPath;
		///An estimate of the log normalising constant ratio over the last step.
		double dlogNCIt;

	public:
		///Create a particle system containing lSize uninitialised particles with the specified history mode.
		sampler(long lSize, HistoryType::Enum htHistoryMode);
		///Create a particle system containing lSize uninitialised particles with the specified history mode and adaptation object.
		sampler(long lSize, HistoryType::Enum htHistoryMode, algParam<Space,Params>* adaptSet);
		///Dispose of a sampler.
		~sampler();
		///Calculates and returns the effective sample size.
		double GetESS(void) const;
		///Returns the number of accepted values from the most recent MCMC iteration.
		int GetAccepted(void) const {return nAccepted;}
		///Returns a flag for whether the ensemble was resampled during the most recent iteration.
		int GetResampled(void) const {return nResampled;}
		/// Returns the effective sample size of the specified particle generation.
		double GetESS(long lGeneration) const;
		///Returns a pointer to the History of the particle system
		const std::vector<historyelement< population<Space> > > & GetHistory(void) const { return History; }
		///Returns a pointer to the algorithm parameters
		algParam<Space,Params> * GetAlgs(void) const { return pAlgParams; }																											   
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
		///Returns the number of evolution times stored in the history.
		long GetHistoryLength(void) const {return History.size();}
		///Returns the values of the particles.
		const std::vector<Space> &  GetPopulationValue(void) const { return pPopulation.GetValue(); }
		///Returns the current population object.
		const population<Space> &  GetPopulation(void) const { return pPopulation; }
		///Returns the value of particle n.
		const Space &  GetPopulationValueN(long n) const { return pPopulation.GetValueN(n); }
		///Returns the logarithmic unnormalized weights of the particles.
		const arma::vec & GetPopulationLogWeight(void) const { return pPopulation.GetLogWeight(); }
		///Returns the unnormalized weights of the particles.
		arma::vec GetPopulationWeight(void) const { return pPopulation.GetWeight(); }  
		///Returns the unnormalized weight of particle n.
		double GetPopulationWeightN(int n) const { return pPopulation.GetWeightN(n); }  
		///Returns the current evolution time of the system.
		long GetTime(void) const {return T;}
		///Initialise the sampler and its constituent particles.
		void Initialise(void);
		///Integrate the supplied function with respect to the current particle set.
		double Integrate(double(*pIntegrand)(const Space &,void*), void* pAuxiliary) const;
		///Integrate the supplied function over the path using the supplied width function and integration method.
		double IntegratePathSampling(PathSamplingType::Enum, double (*pIntegrand)(long,const population<Space>&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary);
		///Integrate the supplied function over the path using the supplied width function and the default integration method (the corrected trapezoid rule).
		double IntegratePathSampling(double (*pIntegrand)(long,const population<Space>&,long,void*), double (*pWidth)(long,void*), void* pAuxiliary) {return IntegratePathSampling(PathSamplingType::TRAPEZOID2, pIntegrand, pWidth, pAuxiliary);}
		///Perform one iteration of the simulation algorithm.
		void Iterate(void);
		///Cancel one iteration of the simulation algorithm.
		void IterateBack(void);
		///Perform one iteration of the simulation algorithm and return the resulting ESS.
		double IterateEss(void);
		///Perform iterations until the specified evolution time is reached.
		void IterateUntil(long lTerminate);
		///Move the particle set by proposing and applying an appropriate move to each particle.
		void MovePopulations(void);
		///Resample the particle set using the specified resampling scheme.
		void Resample(ResampleType::Enum lMode);
		///Sets the entire moveset to the one which is supplied
		void SetMoveSet(moveset<Space>& pNewMoveset) {Moves = pNewMoveset;}
		///Set Resampling Parameters
		void SetResampParams(ResampleType::Enum rtMode, double dThreshold);
		///Dump a specified particle to the specified output stream in a human readable form
		const std::ostream & StreamParticle(std::ostream & os, long n) const;
		///Dump the entire particle set to the specified output stream in a human readable form
		const std::ostream & StreamParticles(std::ostream & os) const;
		///Output a vector indicating the number of accepted MCMC moves at each time instance
		void OstreamMCMCRecordToStream(std::ostream &os) const;
		///Output a 0-1 value vector indicating the times at which resampling occured to an output stream
		void OstreamResamplingRecordToStream(std::ostream &os) const;
		

	private:
		///Duplication of smc::sampler is not currently permitted.
		sampler(const sampler<Space,Params> & sFrom);
		///Duplication of smc::sampler is not currently permitted.
		sampler<Space,Params> & operator=(const sampler<Space,Params> & sFrom);
		///Generate a multinomial random vector with parameters (n,w[1:k]) and store it in X
		void Multinomial(unsigned n, unsigned k, arma::vec w, unsigned int * X);		
		
	protected:
		///Returns the crude normalising constant ratio estimate implied by the weights.
		double CalcLogNC(void) const {return stableLogSumWeights(pPopulation.GetLogWeight());}
	};



	/// The constructor prepares a sampler for use but does not assign any moves to the moveset, set any method for adaptation, initialise the particles
	/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
	/// store the particle set.
	///
	/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
	/// \param htHM The history mode to use: set this to HistoryType::RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
	/// \tparam Space The class used to represent a point in the sample space.
	/// \tparam Params (optional) The class used for any additional parameters.
	template <class Space, class Params>
	sampler<Space,Params>::sampler(long lSize, HistoryType::Enum htHM)
	{
		N = lSize;
		pAlgParams = new algParam<Space,Params>;
		pAlgBelong = 1;
		uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
		
		//Some workable defaults.
		htHistoryMode = htHM;
		rtResampleMode = ResampleType::STRATIFIED;
		dResampleThreshold = 0.5 * N;

	}
	
	/// The constructor prepares a sampler for use but does not assign any moves to the moveset, initialise the particles
	/// or otherwise perform any sampling related tasks. Its main function is to allocate a region of memory in which to
	/// store the particle set.
	///
	/// \param lSize The number of particles present in the ensemble (at time 0 if this is a variable quantity)
	/// \param htHM The history mode to use: set this to HistoryType::RAM to store the whole history of the system and SMC_HISTORY_NONE to avoid doing so.
	/// \param adaptSet The class derived from algParam for parameter adaptation.
	/// \tparam Space The class used to represent a point in the sample space.
	/// \tparam Params The class used for additional algorithm parameters.
	template <class Space, class Params>
	sampler<Space,Params>::sampler(long lSize, HistoryType::Enum htHM, algParam<Space,Params>* adaptSet)
	{
		N = lSize;
		pAlgParams = adaptSet;
		pAlgBelong = 0;
		uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
		
		//Some workable defaults.
		htHistoryMode = htHM;
		rtResampleMode = ResampleType::STRATIFIED;
		dResampleThreshold = 0.5 * N;

	}

	///The destructor deletes the adaptation object only if it was created inside the sampler.
	template <class Space, class Params>
	sampler<Space,Params>::~sampler()
	{
		if(pAlgBelong)
			delete pAlgParams;
	}


	template <class Space, class Params>
	double sampler<Space,Params>::GetESS(void) const
	{
		return expl(2*stableLogSumWeights(pPopulation.GetLogWeight())-stableLogSumWeights(2.0*pPopulation.GetLogWeight()));
	}


	template <class Space, class Params>
	double  sampler<Space,Params>::GetESS(long lGeneration) const
	{
		typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin();
		std::advance(it,lGeneration);
		return it->GetESS(); 
	}

	/// The initialise function:
	///         -# resets the system evolution time to 0 and calls the moveset initialisor to assign each particle in the ensemble.
	///         -# checks the effective sample size and resamples if necessary
	///         -# performs a mcmc step if required
	///         -# appends the particle set to the history if desired
	///
	/// Note that the initialisation function must be specified before calling this function.

	template <class Space, class Params>
	void sampler<Space,Params>::Initialise(void)
	{
		T = 0;
		dlogNCIt = 0;
		dlogNCPath = 0;
		nAccepted = -1;

		std::vector<Space> InitVal(N);
		arma::vec InitWeights(N);
		pPopulation = population<Space>(InitVal,InitWeights);
		Moves.DoInit(pPopulation,N);
		
		//Scaling weights by 1/N (mostly for evidence computation)
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - log(static_cast<double>(N))*arma::ones(N));

		//Estimate the normalising constant
		dlogNCIt = CalcLogNC();
		dlogNCPath += dlogNCIt;

		//Normalise the weights to sensible values....
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - dlogNCIt*arma::ones(N));
		
		//Check if the ESS is below some reasonable threshold and resample if necessary.
		double ESS = GetESS();
		if(ESS < dResampleThreshold) {
			nResampled = 1;
			pAlgParams->updateForMCMC(pPopulation,nAccepted,nResampled);
			Resample(rtResampleMode);
		}
		else {
			nResampled = 0;
			pAlgParams->updateForMCMC(pPopulation,nAccepted,nResampled);
		}
		//A possible MCMC step could be included here.
		nAccepted = Moves.DoMCMC(0,pPopulation, N);

		//Normalise the weights to sensible values....
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N));

		if(htHistoryMode != HistoryType::NONE) {
			History.clear();
			historyelement<population<Space> > histel;
			histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
			History.push_back(histel);
			//History.emplace_back(historyelement<population<Space> >(N, pPopulation, nAccepted, historyflags(nResampled)));
			
		}  
		
		pAlgParams->updateEnd(pPopulation);

		return;
	}

	/// This function returns the result of integrating the supplied function under the empirical measure associated with the
	/// particle set at the present time. The final argument of the integrand function is a pointer which will be supplied
	/// with pAuxiliary to allow for arbitrary additional information to be passed to the function being integrated.
	///
	/// \param pIntegrand The function to integrate with respect to the particle set 
	/// \param pAuxiliary A pointer to any auxiliary data which should be passed to the function

	template <class Space, class Params>
	double sampler<Space,Params>::Integrate(double(*pIntegrand)(const Space&,void*), void * pAuxiliary) const
	{
		long double rValue = 0;
		for(int i =0; i < N; i++)
		{
			rValue += expl(pPopulation.GetLogWeightN(i)) * pIntegrand(pPopulation.GetValueN(i), pAuxiliary);
		}

		rValue /= expl(CalcLogNC());
		return (double)rValue;
	}


	/// This function is intended to be used to estimate integrals of the sort which must be evaluated to determine the
	/// normalising constant of a distribution obtain using a sequence of potential functions proportional to densities with respect
	/// to the initial distribution to define a sequence of distributions leading up to the terminal, interesting distribution.
	///
	/// In this context, the particle set at each time is used to make an estimate of the path sampling integrand, and
	/// numerical integration is then performed to obtain an estimate of the path sampling integral which is the natural logarithm
	/// of the ratio of normalising densities.
	///
	/// The integrand is integrated at every time point in the population history. The results of this integration are
	/// taken to be point-evaluations of the path sampling integrand which are spaced on a grid of intervals given by the
	/// width function. The path sampling integral is then calculated by performing a suitable numerical integration and
	/// the results of this integration is returned.
	///
	/// pAuxiliary is passed to both of the user specified functions to allow the user to pass additional data to either or
	/// both of these functions in a convenient manner. It is safe to use NULL if no such data is used by either function.
	///
	/// \param PStype  The numerical integration method to use
	/// \param pIntegrand  The function to integrated. The first argument is evolution time, the second the population at which the function is to be evaluated, the third is the particle index and the final argument is always pAuxiliary.
	/// \param pWidth      The function which returns the width of the path sampling grid at the specified evolution time. The final argument is always pAuxiliary
	/// \param pAuxiliary  A pointer to auxiliary data to pass to both of the above functions
	/// \tparam Space The class used to represent a point in the sample space.
	/// \tparam Params (optional) The class used for any additional parameters.
	///
	/// The PStype parameter should be set to one of the following:
	/// -# PathSamplingType::RECTANGLE to use the rectangle rule for integration  
	/// -# PathSamplingType::TRAPEZOID1 to use the trapezoidal rule for integration
	/// -# PathSamplingType::TRAPEZOID2 to use the trapezoidal rule for integration with a second order correction

	template <class Space, class Params>
	double sampler<Space,Params>::IntegratePathSampling(PathSamplingType::Enum PStype, double (*pIntegrand)(long,const population<Space> &,long,void*), double (*pWidth)(long,void*), void* pAuxiliary)
	{
		if(htHistoryMode == HistoryType::NONE)
		throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "The path sampling integral cannot be computed as the history of the system was not stored.");

		// historyelement<population<Space> > histel;
		// histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
		// History.push_back(histel);
		
		long lTime = 1;
		long double rValue = 0.0;
		typename std::vector<historyelement<population<Space> > >::const_iterator it;	
		
		switch(PStype) {
		case PathSamplingType::RECTANGLE:
			{
				for(it = ++History.begin(); it!=History.end(); it++){
					rValue += it->Integrate(lTime, pIntegrand, pAuxiliary) * (long double)pWidth(lTime,pAuxiliary);
					lTime++;
				}
				break;
			}
			
			
		case PathSamplingType::TRAPEZOID1:
			{	

				long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
				long double current_expt;
				for(it = ++History.begin(); it!=History.end(); it++){
					current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
					rValue += (long double)pWidth(lTime,pAuxiliary)/2.0 * (previous_expt + current_expt) ;
					lTime++;
					previous_expt = current_expt;
				}
				
				break;
			}

		case PathSamplingType::TRAPEZOID2:
		default:
			{
				long double previous_expt = History.begin()->Integrate(0,pIntegrand,pAuxiliary);
				long double previous_var = History.begin()->Integrate_Var(0,pIntegrand,previous_expt,pAuxiliary);
				long double current_expt;
				long double current_var;
				long double width = 0.0;
				for(it = ++History.begin(); it!=History.end(); it++){
					current_expt = it->Integrate(lTime, pIntegrand, pAuxiliary);
					current_var = it->Integrate_Var(lTime, pIntegrand, current_expt, pAuxiliary);
					width = (long double)pWidth(lTime,pAuxiliary);
					rValue += width/2.0 * (previous_expt + current_expt) - pow(width,2)/12.0*(current_var - previous_var);
					lTime++;
					previous_expt = current_expt;
					previous_var = current_var;
				}	
				
				break;
			}

		}
		
		// History.pop_back();
		
		return ((double)rValue);
	}

	/// The iterate function:
	///         -# moves the current particle set
	///         -# checks the effective sample size and resamples if necessary
	///         -# performs a mcmc step if required
	///         -# appends the current particle set to the history if desired
	///         -# increments the current evolution time
	template <class Space, class Params>
	void sampler<Space,Params>::Iterate(void)
	{
		IterateEss();
		return;
	}

	template <class Space, class Params>
	void sampler<Space,Params>::IterateBack(void)
	{
		if(htHistoryMode == HistoryType::NONE)
		throw SMC_EXCEPTION(SMCX_MISSING_HISTORY, "An attempt to undo an iteration was made; unforunately, the system history has not been stored.");

		History.pop_back();
		historyelement<population<Space> > recent = History.back();
		pPopulation = recent.GetRefs();
		N =recent.GetNumber();
		nAccepted = recent.AcceptCount();
		nResampled = recent.WasResampled();
		T--;
		return;
	}

	template <class Space, class Params>
	double sampler<Space,Params>::IterateEss(void)
	{		
		pAlgParams->updateForMove(pPopulation);
		//Move the particle set.
		MovePopulations();

		//Estimate the normalising constant
		dlogNCIt = CalcLogNC();
		dlogNCPath += dlogNCIt;
		
		//Normalise the weights
		pPopulation.SetLogWeight(pPopulation.GetLogWeight()  - dlogNCIt*arma::ones(N));

		//Check if the ESS is below some reasonable threshold and resample if necessary.
		//A mechanism for setting this threshold is required.
		double ESS = GetESS();
		if(ESS < dResampleThreshold) {
			nResampled = 1;
			pAlgParams->updateForMCMC(pPopulation,nAccepted,nResampled);
			Resample(rtResampleMode);
		}
		else{
		nResampled = 0;
			pAlgParams->updateForMCMC(pPopulation,nAccepted,nResampled);
		}
		//A possible MCMC step could be included here.
		nAccepted = Moves.DoMCMC(T+1,pPopulation,N);
		
		//Normalise the weights to sensible values....
		pPopulation.SetLogWeight(pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N));
		
		// Increment the evolution time.
		T++;

		//Finally, the current particle set should be appended to the historical process.
		if(htHistoryMode != HistoryType::NONE){
			historyelement<population<Space> > histel;
			histel.Set(N, pPopulation, nAccepted, historyflags(nResampled));
			History.push_back(histel);
			//History.emplace_back(historyelement<population<Space> >(N, pPopulation, nAccepted, historyflags(nResampled)));
		}
		
		pAlgParams->updateEnd(pPopulation);

		return ESS;
	}

	template <class Space, class Params>
	void sampler<Space,Params>::IterateUntil(long lTerminate)
	{
		while(T < lTerminate)
		Iterate();
	}

	template <class Space, class Params>
	void sampler<Space,Params>::MovePopulations(void)
	{
		Moves.DoMove(T+1,pPopulation, N);
	}

	template <class Space, class Params>
	void sampler<Space,Params>::Resample(ResampleType::Enum lMode)
	{
		//Resampling is done in place.
		unsigned uMultinomialCount;

		
		//First obtain a count of the number of children each particle has.
		switch(lMode) {
		case ResampleType::MULTINOMIAL:
			//Sample from a suitable multinomial vector
			dRSWeights = pPopulation.GetWeight();
			//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
			Multinomial(N,N,dRSWeights,uRSCount.memptr());
			//uRSCount = pRng->Multinomial(N,N,dRSWeights);
			break;
			
			
		case ResampleType::RESIDUAL:
			dRSWeights = exp(log((double)N)*arma::ones(N) + pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N));
			
			uRSIndices = arma::zeros<arma::Col<unsigned int> >((int)N);
			
			
			for(int i = 0; i < N; ++i)
			uRSIndices(i) = unsigned(floor(dRSWeights(i)));
			dRSWeights = dRSWeights - uRSIndices;
			uMultinomialCount = N - arma::sum(uRSIndices);
			
			//uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
			Multinomial(uMultinomialCount,N,dRSWeights,uRSCount.memptr());
			uRSCount += uRSIndices;
			
			break;
			
		case ResampleType::STRATIFIED:
		default:
			{
				// Procedure for stratified sampling
				//Generate a random number between 0 and 1/N times the sum of the weights
				double dRand = R::runif(0,1/(double)N);
				
				arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - CalcLogNC()));
				
				int k = 0;
				int j = 0;
				uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
				while(j < N) {
					while((dWeightCumulative(k) - dRand) > ((double)j)/((double)N) && j < N) {
						uRSCount(k)++;
						j++;
						dRand = R::runif(0,1/(double)N);
					}
					k++;
				}
				break;
			}

		case ResampleType::SYSTEMATIC:
			{
				// Procedure for stratified sampling but with a common RV for each stratum
				//Generate a random number between 0 and 1/N times the sum of the weights
				double dRand = R::runif(0,1/(double)N);
				
				int j = 0, k = 0;
				uRSCount = arma::zeros<arma::Col<unsigned int> >((int)N);
				arma::vec dWeightCumulative = arma::cumsum(exp(pPopulation.GetLogWeight() - CalcLogNC()*arma::ones(N)));
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
				Current = pPopulation.GetValueN((int)uRSIndices(i));
				pPopulation.SetValueN(Current,i);
			}
		}
		//Set equal normalised weights
		pPopulation.SetLogWeight(- log(static_cast<double>(N))*arma::ones(N));
	} 				
	
	/// This function configures the resampling parameters, allowing the specification of both the resampling
	/// mode and the threshold at which resampling is used.
	///
	/// \param rtMode The resampling mode to be used.
	/// \param dThreshold The threshold at which resampling is deemed necesary.
	///
	/// The rtMode parameter should be set to one of the following:
	/// -# ResampleType::MULTINOMIAL to use multinomial resampling  
	/// -# ResampleType::RESIDUAL to use residual resampling
	/// -# ResampleType::STRATIFIED to use stratified resampling
	/// -# ResampleType::SYSTEMATIC to use systematic resampling
	///
	/// The dThreshold parameter can be set to a value in the range [0,1) corresponding to a fraction of the size of
	/// the particle set or it may be set to an integer corresponding to an actual effective sample size.
	
	template <class Space, class Params>
	void sampler<Space,Params>::SetResampParams(ResampleType::Enum rtMode, double dThreshold)
	{
		rtResampleMode = rtMode;
		if(dThreshold < 1)
		dResampleThreshold = dThreshold * N;
		else
		dResampleThreshold = dThreshold;
	}

	// template <class Space, class Params>
	// void sampler<Space,Params>::SetAdaptSet(algParam<Space,Params>* inParams)
	// {
		// pAlgParams = inParams; //later want to add lSize, history type to this
	// }

	template <class Space, class Params>
	const std::ostream & sampler<Space,Params>::StreamParticle(std::ostream & os, long n) const
	{
		Space val = pPopulation.GetValueN(n);
		double weight = pPopulation.GetWeightN(n);
		os << val << "," << weight;
		return os;
	}

	template <class Space, class Params>
	const std::ostream & sampler<Space,Params>::StreamParticles(std::ostream & os) const
	{
		Space value;
		double weight;
		for(int i = 0; i < N - 1; i++){
			value = pPopulation.GetValueN(i);
			weight = pPopulation.GetWeightN(i);
			os << value << "," << weight << std::endl;
		}
		
		return os;
	}

	/// This function records the MCMC acceptance history to the specified output stream as a list of
	/// the number of moves accepted at each time instant.
	///
	/// \param os The output stream to send the data to.
	template <class Space, class Params>
	void sampler<Space,Params>:: OstreamMCMCRecordToStream(std::ostream &os) const
	{
		os << "Accepted MCMC proposals history:" << std::endl;
		os << "======================" << std::endl;
		for(typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin(); it!=History.end(); it++){
			os << it->AcceptCount() << std::endl;
		}
	}
	/// This function records the resampling history to the specified output stream as a 0-1 valued list which takes
	/// the value 1 for those time instances when resampling occured and 0 otherwise.
	///
	/// \param os The output stream to send the data to.
	template <class Space, class Params>
	void sampler<Space,Params>:: OstreamResamplingRecordToStream(std::ostream &os) const
	{
		os << "Resampling history:" << std::endl;
		os << "======================" << std::endl;
		os << "Flag\t" << "ESS\t" << std::endl;
		for(typename std::vector<historyelement<population<Space> > >::const_iterator it = History.begin(); it!=History.end(); it++){ 
			if(it->WasResampled())
			os << "1\t";
			else
			os << "0\t";

			os << it->GetESS() << std::endl;
		}
	}

	template <class Space, class Params>
	void sampler<Space,Params>::Multinomial(unsigned n, unsigned k, arma::vec w, unsigned int * X) {		
		Rcpp::IntegerVector v(k);
		w = w/arma::sum(w);
		
		double * w_mem = w.memptr();
		
		// R sources:  rmultinom(int n, double* prob, int K, int* rN);
		rmultinom(static_cast<int>(n), const_cast<double*>(w_mem), static_cast<int>(k), &(v[0]));
		
		for (unsigned int i=0; i<k; i++) {
			X[i] = static_cast<unsigned int>(v[i]);
		}	
	}

	
	inline double stableLogSumWeights(const arma::vec & logw){
		long N = logw.n_rows;
		double dMaxWeight = arma::max(logw);
		double sum = arma::sum(exp(logw - dMaxWeight*arma::ones(N)));

		return (dMaxWeight + log(sum));
	}
}


namespace std {
	/// Produce a human-readable display of the state of an smc::sampler class using the stream operator.

	/// \param os The output stream to which the display should be made.
	/// \param s  The sampler which is to be displayed.
	template <class Space, class Params>
	std::ostream & operator<< (std::ostream & os, smc::sampler<Space,Params> & s)
	{
		os << "Sampler Configuration:" << std::endl;
		os << "======================" << std::endl;
		os << "Evolution Time:   " << s.GetTime() << std::endl;
		os << "Particle Set Size: " << s.GetNumber() << std::endl;
		os << "Effective Sample Size: " << s.GetESS() << std::endl;
		os << std::endl;
		os << "Particle Set: " << std::endl;
		s.StreamParticles(os);
		os << std::endl;
		return os;
	}
}


#endif
